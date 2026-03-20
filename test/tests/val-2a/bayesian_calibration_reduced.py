#!/usr/bin/env python3
"""
REDUCED Bayesian Calibration for Val-2a Model (2 Parameters Only)
Uses emcee for MCMC sampling to calibrate only the most influential parameters:
  - Kr_left_fraction (ST = 1.218)
  - gaussian_factor (ST = 0.123)

This provides:
  - ~8-16x faster calibration than 8-parameter version
  - Better parameter identifiability
  - Clearer interpretation of results
  - Lower computational cost for surrogate

Requires: emcee, corner, numpy, pandas, scipy, scikit-learn
Install: pip install emcee corner numpy pandas scipy scikit-learn pyDOE3

Usage:
    # Step 1: Generate ensemble data for surrogate
    python bayesian_calibration_reduced.py --mode generate-ensemble --n-samples 200 --parallel 8

    # Step 2: Build and validate surrogate model
    python bayesian_calibration_reduced.py --mode build-surrogate

    # Step 3: Run Bayesian calibration
    python bayesian_calibration_reduced.py --mode calibrate --walkers 32 --steps 5000

    # Step 4: Analyze results
    python bayesian_calibration_reduced.py --mode analyze
"""

import numpy as np
import pandas as pd
import emcee
import corner
import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.interpolate import interp1d
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel, Matern
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import pickle

# Set plotting style
sns.set_style("whitegrid")
sns.set_context("talk")

# REDUCED parameter set - only 2 most influential
PARAMETERS = {
    'Kr_left_fraction': {
        'bounds': [0.95, 0.9999],
        'log_scale': False,
        'units': 'dimensionless',
        'nominal': 0.99,
        'description': 'Surface cleanup fraction'
    },
    'gaussian_factor': {
        'bounds': [1.0, 2.5],
        'log_scale': False,
        'units': 'dimensionless',
        'nominal': 1.5,
        'description': 'Implantation profile scaling'
    }
}

# Experimental data
EXP_DATA_FILE = 'gold/experiment_data_paper.csv'


def load_experimental_data(exp_file=EXP_DATA_FILE):
    """Load experimental permeation flux data"""
    df = pd.read_csv(exp_file)
    exp_time = df['time (s)'].values
    exp_flux = df['permeation flux (atom/m^2/s)'].values
    return exp_time, exp_flux


def generate_lhs_samples(n_samples=200):
    """Generate Latin Hypercube samples for training data

    Args:
        n_samples: Number of samples

    Returns:
        DataFrame with parameter samples
    """
    from pyDOE3 import lhs

    param_names = list(PARAMETERS.keys())
    n_params = len(param_names)

    # Generate LHS samples in [0, 1]^d
    lhs_samples = lhs(n_params, samples=n_samples, criterion='maximin')

    # Transform to parameter bounds
    samples_dict = {}
    for i, param_name in enumerate(param_names):
        bounds = PARAMETERS[param_name]['bounds']
        if PARAMETERS[param_name]['log_scale']:
            # Log-uniform sampling
            log_bounds = [np.log10(bounds[0]), np.log10(bounds[1])]
            samples_dict[param_name] = 10 ** (log_bounds[0] + lhs_samples[:, i] * (log_bounds[1] - log_bounds[0]))
        else:
            # Uniform sampling
            samples_dict[param_name] = bounds[0] + lhs_samples[:, i] * (bounds[1] - bounds[0])

    return pd.DataFrame(samples_dict)


def _run_single_tmap8(args_tuple):
    """Helper function for parallel TMAP8 execution (must be at module level for pickling)

    Args:
        args_tuple: (idx, params, output_dir, tmap8_exe, input_file_abs, exp_time)

    Returns:
        Dictionary with run results or None if failed
    """
    idx, params, output_dir, tmap8_exe, input_file_abs, exp_time = args_tuple

    run_dir = output_dir / f"run_{idx:04d}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Build command
    args = [str(tmap8_exe), '-i', str(input_file_abs)]
    args.extend([f"{k}={v}" for k, v in params.items()])
    args.extend(['--no-color'])

    try:
        result = subprocess.run(args, cwd=run_dir, capture_output=True, timeout=1800, text=True)

        if result.returncode != 0:
            print(f"Run {idx} failed with return code {result.returncode}")
            print(f"  Error: {result.stderr[:200]}")
            return None

        # Load output
        csv_file = run_dir / 'val-2a_param_reduced_out.csv'
        if not csv_file.exists():
            print(f"Run {idx}: output file not found")
            return None

        df = pd.read_csv(csv_file)

        # Interpolate to experimental times
        f = interp1d(df['time'], df['scaled_recombination_flux_right'],
                    kind='linear', fill_value='extrapolate')
        flux_at_exp_times = f(exp_time)

        return {'run_id': idx, 'flux': flux_at_exp_times, 'params': params}

    except Exception as e:
        print(f"Run {idx} failed: {e}")
        return None


def run_tmap8_ensemble(samples_df, output_dir, tmap8_exe, input_file, n_parallel=8):
    """Run TMAP8 ensemble and extract permeation flux time series

    Args:
        samples_df: DataFrame with parameter samples
        output_dir: Directory for outputs
        tmap8_exe: Path to TMAP8 executable
        input_file: Input file path
        n_parallel: Number of parallel processes

    Returns:
        List of result dictionaries
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    exp_time, _ = load_experimental_data()
    input_file_abs = Path(input_file).resolve()

    # Prepare arguments for parallel execution
    args_list = [
        (idx, row.to_dict(), output_dir, tmap8_exe, input_file_abs, exp_time)
        for idx, row in samples_df.iterrows()
    ]

    # Run in parallel
    results = []
    with ProcessPoolExecutor(max_workers=n_parallel) as executor:
        futures = [executor.submit(_run_single_tmap8, args) for args in args_list]

        for future in tqdm(as_completed(futures), total=len(futures), desc="Running ensemble"):
            result = future.result()
            if result is not None:
                results.append(result)

    return results


class SurrogateModel:
    """Gaussian Process surrogate model for TMAP8 (reduced 2-parameter version)"""

    def __init__(self):
        self.param_names = list(PARAMETERS.keys())
        self.scaler_X = StandardScaler()
        self.gp_models = []  # One GP per time point
        self.exp_times = None
        self.trained = False

    def prepare_training_data(self, ensemble_results):
        """Prepare X, y matrices from ensemble results"""
        # X: parameter values
        X = []
        for result in ensemble_results:
            params = [result['params'][name] for name in self.param_names]
            X.append(params)
        X = np.array(X)

        # y: flux at each experimental time point
        y = np.array([result['flux'] for result in ensemble_results])

        return X, y

    def train(self, X, y, subsample_time=5):
        """Train GP surrogate - separate GP for each time point

        For 2 parameters, we train independent GPs for each time point.
        More robust than single multi-output GP.

        Args:
            X: Parameter values (n_samples, 2)
            y: Flux values (n_samples, n_times)
            subsample_time: Use every Nth time point
        """
        print(f"\nTraining GP surrogate model...")
        print(f"  Training samples: {X.shape[0]}")
        print(f"  Time points: {y.shape[1]} -> {y.shape[1]//subsample_time} (subsampled)")

        # Standardize inputs once
        X_scaled = self.scaler_X.fit_transform(X)

        # Subsample time points
        time_indices = np.arange(0, y.shape[1], subsample_time)
        y_subsampled = y[:, time_indices]

        # Train separate GP for each time point
        print(f"  Training {len(time_indices)} GPs...")
        self.gp_models = []

        kernel_template = ConstantKernel(1.0, (1e-2, 1e2)) * RBF(
            length_scale=[1.0, 1.0],
            length_scale_bounds=(1e-1, 10.0)
        ) + WhiteKernel(noise_level=0.1, noise_level_bounds=(1e-3, 1.0))

        for t_idx in time_indices:
            y_t = y_subsampled[:, list(time_indices).index(t_idx)]

            # Skip if all values are too similar (no variance to learn)
            if np.std(y_t) < 1e10:
                self.gp_models.append(None)
                continue

            gp = GaussianProcessRegressor(
                kernel=kernel_template,
                n_restarts_optimizer=5,
                alpha=1e-8,
                normalize_y=True
            )
            gp.fit(X_scaled, y_t)
            self.gp_models.append(gp)

        self.time_indices = time_indices
        self.trained = True
        print("  Training complete!")

    def predict(self, X, return_std=False):
        """Predict flux time series for given parameters

        Args:
            X: Parameter values (n_samples, 2)
            return_std: If True, return standard deviation

        Returns:
            y_pred: Flux predictions (n_samples, n_times_subsampled)
            y_std (optional): Prediction uncertainties
        """
        if not self.trained:
            raise RuntimeError("Model not trained yet")

        X_scaled = self.scaler_X.transform(X)
        n_samples = X_scaled.shape[0]
        n_times = len(self.gp_models)

        y_pred = np.zeros((n_samples, n_times))
        if return_std:
            y_std = np.zeros((n_samples, n_times))

        # Predict from each GP
        for i, gp in enumerate(self.gp_models):
            if gp is None:
                # Use mean value if GP wasn't trained
                y_pred[:, i] = 0
                if return_std:
                    y_std[:, i] = 1e15
            else:
                if return_std:
                    y_pred[:, i], y_std[:, i] = gp.predict(X_scaled, return_std=True)
                else:
                    y_pred[:, i] = gp.predict(X_scaled)

        if return_std:
            return y_pred, y_std
        else:
            return y_pred

    def validate(self, X_test, y_test):
        """Validate surrogate on test data"""
        y_pred, y_std = self.predict(X_test, return_std=True)

        # Subsample y_test to match training
        y_test_sub = y_test[:, self.time_indices]

        # Calculate metrics
        mse = np.mean((y_pred - y_test_sub)**2)
        rmse = np.sqrt(mse)
        mae = np.mean(np.abs(y_pred - y_test_sub))
        r2 = 1 - np.sum((y_pred - y_test_sub)**2) / np.sum((y_test_sub - np.mean(y_test_sub))**2)

        print(f"\nSurrogate validation metrics:")
        print(f"  RMSE: {rmse:.2e}")
        print(f"  MAE:  {mae:.2e}")
        print(f"  R²:   {r2:.4f}")

        return {'rmse': rmse, 'mae': mae, 'r2': r2}


def log_prior(theta):
    """Log prior probability for parameters

    Args:
        theta: Parameter values [Kr_left_fraction, gaussian_factor]

    Returns:
        log_prior: Log prior probability
    """
    Kr_left_fraction, gaussian_factor = theta

    # Bounds checking
    if not (0.95 <= Kr_left_fraction <= 0.9999):
        return -np.inf
    if not (1.0 <= gaussian_factor <= 2.5):
        return -np.inf

    # Beta prior for Kr_left_fraction: Beta(50, 2) scaled to [0.95, 1.0]
    # Transform to [0, 1]
    x_norm = (Kr_left_fraction - 0.95) / (0.9999 - 0.95)
    log_p_Kr = stats.beta.logpdf(x_norm, a=50, b=2)

    # Uniform prior for gaussian_factor
    log_p_gauss = np.log(1.0 / (2.5 - 1.0))

    return log_p_Kr + log_p_gauss


def log_likelihood(theta, surrogate, exp_time, exp_flux, noise_std=1e16):
    """Log likelihood function using surrogate model

    Args:
        theta: Parameter values
        surrogate: Trained surrogate model
        exp_time: Experimental time points
        exp_flux: Experimental flux values
        noise_std: Measurement noise standard deviation

    Returns:
        log_likelihood: Log likelihood value
    """
    # Predict flux using surrogate
    X = np.array([theta])
    try:
        y_pred = surrogate.predict(X)[0]  # Shape: (n_times_subsampled,)

        # Subsample experimental data to match
        exp_flux_sub = exp_flux[surrogate.time_indices]

        # Gaussian likelihood
        residuals = exp_flux_sub - y_pred
        log_like = -0.5 * np.sum((residuals / noise_std)**2) - len(residuals) * np.log(noise_std * np.sqrt(2 * np.pi))

        return log_like

    except Exception as e:
        print(f"Likelihood evaluation failed: {e}")
        return -np.inf


def log_probability(theta, surrogate, exp_time, exp_flux, noise_std):
    """Log posterior probability"""
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf

    ll = log_likelihood(theta, surrogate, exp_time, exp_flux, noise_std)
    return lp + ll


def run_mcmc(surrogate, exp_time, exp_flux, n_walkers=32, n_steps=5000, noise_std=1e16):
    """Run MCMC calibration using emcee

    Args:
        surrogate: Trained surrogate model
        exp_time: Experimental time points
        exp_flux: Experimental flux values
        n_walkers: Number of MCMC walkers
        n_steps: Number of MCMC steps
        noise_std: Measurement noise standard deviation

    Returns:
        sampler: emcee sampler object
    """
    ndim = 2  # 2 parameters

    # Initialize walkers near nominal values with small perturbations
    nominal = np.array([PARAMETERS['Kr_left_fraction']['nominal'],
                       PARAMETERS['gaussian_factor']['nominal']])
    pos = nominal + 1e-3 * np.random.randn(n_walkers, ndim)

    # Ensure within bounds
    pos[:, 0] = np.clip(pos[:, 0], 0.95, 0.9999)
    pos[:, 1] = np.clip(pos[:, 1], 1.0, 2.5)

    # Setup sampler
    sampler = emcee.EnsembleSampler(
        n_walkers, ndim, log_probability,
        args=(surrogate, exp_time, exp_flux, noise_std)
    )

    # Run MCMC
    print(f"\nRunning MCMC with {n_walkers} walkers for {n_steps} steps...")
    sampler.run_mcmc(pos, n_steps, progress=True)

    return sampler


def analyze_mcmc_results(sampler, param_names, burn_in=1000):
    """Analyze MCMC results

    Args:
        sampler: emcee sampler
        param_names: List of parameter names
        burn_in: Number of burn-in samples to discard

    Returns:
        results: Dictionary with summary statistics
    """
    # Get samples (discard burn-in)
    samples = sampler.get_chain(discard=burn_in, flat=True)

    # Summary statistics
    means = np.mean(samples, axis=0)
    stds = np.std(samples, axis=0)
    medians = np.median(samples, axis=0)
    q16, q84 = np.percentile(samples, [16, 84], axis=0)

    print("\n" + "="*80)
    print("CALIBRATION RESULTS (Reduced 2-Parameter Model)")
    print("="*80)
    print(f"\nBurn-in: {burn_in} steps")
    print(f"Post-burn-in samples: {len(samples)}")
    print()

    results = {}
    for i, name in enumerate(param_names):
        print(f"{name}:")
        print(f"  Mean:   {means[i]:.6f} ± {stds[i]:.6f}")
        print(f"  Median: {medians[i]:.6f}")
        print(f"  95% CI: [{q16[i]:.6f}, {q84[i]:.6f}]")
        print(f"  Nominal: {PARAMETERS[name]['nominal']:.6f}")
        print()

        results[name] = {
            'mean': float(means[i]),
            'std': float(stds[i]),
            'median': float(medians[i]),
            'q16': float(q16[i]),
            'q84': float(q84[i]),
            'nominal': PARAMETERS[name]['nominal']
        }

    # Acceptance fraction
    print(f"Mean acceptance fraction: {np.mean(sampler.acceptance_fraction):.3f}")

    # Autocorrelation time (if converged)
    try:
        tau = sampler.get_autocorr_time()
        print(f"Autocorrelation time: {tau}")
    except Exception:
        print("Autocorrelation time: Not converged")

    print("="*80)

    return results, samples


def plot_mcmc_results(samples, param_names, output_dir):
    """Generate diagnostic plots for MCMC results"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Corner plot
    fig = corner.corner(
        samples,
        labels=param_names,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 12}
    )
    plt.savefig(output_dir / 'corner_plot.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'corner_plot.png'}")

    # Trace plots
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    for i, name in enumerate(param_names):
        axes[i].plot(samples[:, i], alpha=0.3, color='steelblue')
        axes[i].set_ylabel(name)
        axes[i].axhline(PARAMETERS[name]['nominal'], color='red', linestyle='--', label='Nominal')
        axes[i].legend()
    axes[-1].set_xlabel('Sample')
    plt.tight_layout()
    plt.savefig(output_dir / 'trace_plots.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_dir / 'trace_plots.png'}")


def main():
    parser = argparse.ArgumentParser(description='Reduced Bayesian calibration (2 parameters)')
    parser.add_argument('--mode', type=str, required=True,
                       choices=['generate-ensemble', 'build-surrogate', 'calibrate', 'analyze'],
                       help='Execution mode')
    parser.add_argument('--n-samples', type=int, default=200,
                       help='Number of LHS samples for ensemble')
    parser.add_argument('--parallel', type=int, default=8,
                       help='Number of parallel processes')
    parser.add_argument('--tmap8-exe', type=str,
                       default=str(Path.home() / 'projects/TMAP8/tmap8-opt'),
                       help='Path to TMAP8 executable')
    parser.add_argument('--input-file', type=str, default='val-2a_param_reduced.i',
                       help='TMAP8 input file')
    parser.add_argument('--ensemble-dir', type=str, default='calibration_reduced_ensemble',
                       help='Directory for ensemble runs')
    parser.add_argument('--output-dir', type=str, default='calibration_reduced_results',
                       help='Directory for calibration outputs')
    parser.add_argument('--walkers', type=int, default=32,
                       help='Number of MCMC walkers')
    parser.add_argument('--steps', type=int, default=5000,
                       help='Number of MCMC steps')
    parser.add_argument('--burn-in', type=int, default=1000,
                       help='Burn-in steps')
    parser.add_argument('--noise-std', type=float, default=1e16,
                       help='Measurement noise standard deviation')

    args = parser.parse_args()

    if args.mode == 'generate-ensemble':
        print("="*80)
        print("MODE: Generate Ensemble (Reduced 2-Parameter)")
        print("="*80)

        # Generate samples
        samples_df = generate_lhs_samples(args.n_samples)
        print(f"\nGenerated {len(samples_df)} LHS samples")

        # Save samples
        ensemble_dir = Path(args.ensemble_dir)
        ensemble_dir.mkdir(parents=True, exist_ok=True)
        samples_df.to_csv(ensemble_dir / 'lhs_samples.csv', index=False)
        print(f"Saved: {ensemble_dir / 'lhs_samples.csv'}")

        # Run ensemble
        results = run_tmap8_ensemble(
            samples_df, ensemble_dir, args.tmap8_exe, args.input_file, args.parallel
        )

        print(f"\nCompleted {len(results)}/{len(samples_df)} simulations successfully")

        # Save results
        with open(ensemble_dir / 'ensemble_results.pkl', 'wb') as f:
            pickle.dump(results, f)
        print(f"Saved: {ensemble_dir / 'ensemble_results.pkl'}")

    elif args.mode == 'build-surrogate':
        print("="*80)
        print("MODE: Build Surrogate Model")
        print("="*80)

        # Load ensemble results
        ensemble_dir = Path(args.ensemble_dir)
        with open(ensemble_dir / 'ensemble_results.pkl', 'rb') as f:
            results = pickle.load(f)

        print(f"Loaded {len(results)} ensemble results")

        # Create and train surrogate
        surrogate = SurrogateModel()
        X, y = surrogate.prepare_training_data(results)

        # Split train/test
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

        surrogate.train(X_train, y_train, subsample_time=5)

        # Validate
        metrics = surrogate.validate(X_test, y_test)

        # Save surrogate
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        with open(output_dir / 'surrogate_model.pkl', 'wb') as f:
            pickle.dump(surrogate, f)
        print(f"\nSaved: {output_dir / 'surrogate_model.pkl'}")

    elif args.mode == 'calibrate':
        print("="*80)
        print("MODE: Bayesian Calibration")
        print("="*80)

        # Load surrogate
        output_dir = Path(args.output_dir)
        with open(output_dir / 'surrogate_model.pkl', 'rb') as f:
            surrogate = pickle.load(f)
        print("Loaded surrogate model")

        # Load experimental data
        exp_time, exp_flux = load_experimental_data()
        print(f"Loaded experimental data: {len(exp_time)} points")

        # Run MCMC
        sampler = run_mcmc(surrogate, exp_time, exp_flux, args.walkers, args.steps, args.noise_std)

        # Save sampler
        with open(output_dir / 'mcmc_sampler.pkl', 'wb') as f:
            pickle.dump(sampler, f)
        print(f"\nSaved: {output_dir / 'mcmc_sampler.pkl'}")

        # Analyze
        param_names = list(PARAMETERS.keys())
        results, samples = analyze_mcmc_results(sampler, param_names, args.burn_in)

        # Save results
        with open(output_dir / 'calibration_results.json', 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Saved: {output_dir / 'calibration_results.json'}")

        np.save(output_dir / 'mcmc_samples.npy', samples)
        print(f"Saved: {output_dir / 'mcmc_samples.npy'}")

        # Plot
        plot_mcmc_results(samples, param_names, output_dir)

    elif args.mode == 'analyze':
        print("="*80)
        print("MODE: Analyze Results")
        print("="*80)

        output_dir = Path(args.output_dir)

        # Load results
        with open(output_dir / 'calibration_results.json', 'r') as f:
            results = json.load(f)

        samples = np.load(output_dir / 'mcmc_samples.npy')

        print(f"Loaded {len(samples)} samples")
        print("\nCalibrated parameter values:")
        for name, vals in results.items():
            print(f"\n{name}:")
            print(f"  Posterior mean: {vals['mean']:.6f} ± {vals['std']:.6f}")
            print(f"  95% CI: [{vals['q16']:.6f}, {vals['q84']:.6f}]")
            print(f"  Nominal value: {vals['nominal']:.6f}")

        print("\n" + "="*80)


if __name__ == '__main__':
    main()
