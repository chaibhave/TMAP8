#!/usr/bin/env python3
"""
Bayesian Calibration for Val-2a Model
Uses PyMC for MCMC sampling to calibrate parameters against experimental data

Requires: pymc, arviz, numpy, pandas, scipy
Install: pip install pymc arviz numpy pandas scipy

Usage:
    # Step 1: Generate ensemble data for surrogate (if needed)
    python bayesian_calibration.py --mode generate-ensemble --n-samples 500 --parallel 8

    # Step 2: Build and validate surrogate model
    python bayesian_calibration.py --mode build-surrogate

    # Step 3: Run Bayesian calibration
    python bayesian_calibration.py --mode calibrate --chains 4 --samples 25000

    # Step 4: Analyze and validate results
    python bayesian_calibration.py --mode analyze
"""

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az
import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.interpolate import interp1d
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import subprocess
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import pickle

# Set plotting style
sns.set_style("whitegrid")
sns.set_context("talk")

# Parameter definitions (same as sensitivity analysis)
from sensitivity_analysis import PARAMETERS, setup_problem

# Experimental data times
EXP_DATA_FILE = 'gold/experiment_data_paper.csv'


def load_experimental_data(exp_file=EXP_DATA_FILE):
    """Load experimental permeation flux data"""
    df = pd.read_csv(exp_file)
    exp_time = df['time (s)'].values
    exp_flux = df['permeation flux (atom/m^2/s)'].values
    return exp_time, exp_flux


def generate_lhs_samples(n_samples=500):
    """Generate Latin Hypercube samples for training data

    Args:
        n_samples: Number of samples

    Returns:
        DataFrame with parameter samples
    """
    from pyDOE2 import lhs

    problem = setup_problem()
    n_params = problem['num_vars']

    # Generate LHS samples in [0, 1]^d
    lhs_samples = lhs(n_params, samples=n_samples, criterion='maximin')

    # Transform to parameter bounds
    samples_dict = {}
    for i, param_name in enumerate(problem['names']):
        bounds = PARAMETERS[param_name]['bounds']
        if PARAMETERS[param_name]['log_scale']:
            # Log-uniform sampling
            log_bounds = [np.log10(bounds[0]), np.log10(bounds[1])]
            samples_dict[param_name] = 10 ** (log_bounds[0] + lhs_samples[:, i] * (log_bounds[1] - log_bounds[0]))
        else:
            # Uniform sampling
            samples_dict[param_name] = bounds[0] + lhs_samples[:, i] * (bounds[1] - bounds[0])

    return pd.DataFrame(samples_dict)


def run_tmap8_ensemble(samples_df, output_dir, tmap8_exe, input_file, n_parallel=8):
    """Run TMAP8 ensemble and extract permeation flux time series

    Args:
        samples_df: DataFrame with parameter samples
        output_dir: Directory for outputs
        tmap8_exe: Path to TMAP8 executable
        input_file: Input file path
        n_parallel: Number of parallel processes

    Returns:
        DataFrame with flux time series for each sample
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    exp_time, _ = load_experimental_data()

    def run_single(row_tuple):
        idx, params = row_tuple
        run_dir = output_dir / f"run_{idx:04d}"
        run_dir.mkdir(parents=True, exist_ok=True)

        # Build command
        args = [str(tmap8_exe), '-i', str(input_file)]
        args.extend([f"{k}={v}" for k, v in params.items()])
        args.extend(['--no-color'])

        try:
            result = subprocess.run(args, cwd=run_dir, capture_output=True, timeout=600)

            if result.returncode != 0:
                return None

            # Load output
            csv_file = run_dir / 'val-2a_param_out.csv'
            df = pd.read_csv(csv_file)

            # Interpolate to experimental times
            f = interp1d(df['time'], df['scaled_recombination_flux_right'],
                        kind='linear', fill_value='extrapolate')
            flux_at_exp_times = f(exp_time)

            return {'run_id': idx, 'flux': flux_at_exp_times, 'params': params}

        except Exception as e:
            print(f"Run {idx} failed: {e}")
            return None

    # Run in parallel
    results = []
    with ProcessPoolExecutor(max_workers=n_parallel) as executor:
        futures = [executor.submit(run_single, (idx, row))
                  for idx, row in samples_df.iterrows()]

        for future in tqdm(futures, desc="Running ensemble"):
            result = future.result()
            if result is not None:
                results.append(result)

    return results


class SurrogateModel:
    """Gaussian Process surrogate model for TMAP8"""

    def __init__(self):
        self.problem = setup_problem()
        self.param_names = self.problem['names']
        self.scaler_X = StandardScaler()
        self.scaler_y = StandardScaler()
        self.models = []  # One GP per output time point
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

    def train(self, X, y, subsample_time=10):
        """Train GP surrogate

        Args:
            X: Parameter values (n_samples, n_params)
            y: Flux values (n_samples, n_times)
            subsample_time: Use every Nth time point to reduce computational cost
        """
        print(f"Training surrogate model...")
        print(f"  Training samples: {X.shape[0]}")
        print(f"  Time points: {y.shape[1]} -> {y.shape[1]//subsample_time} (subsampled)")

        # Standardize inputs
        X_scaled = self.scaler_X.fit_transform(X)

        # Train separate GP for each time point (subsampled)
        self.models = []
        time_indices = np.arange(0, y.shape[1], subsample_time)

        for t_idx in tqdm(time_indices, desc="Training GPs"):
            y_t = y[:, t_idx].reshape(-1, 1)
            y_t_scaled = self.scaler_y.fit_transform(y_t)

            # Define kernel
            kernel = ConstantKernel(1.0) * RBF(length_scale=np.ones(X.shape[1])) + WhiteKernel(noise_level=0.01)

            # Train GP
            gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5, normalize_y=False)
            gp.fit(X_scaled, y_t_scaled.ravel())

            self.models.append(gp)

        self.trained = True
        self.time_indices = time_indices
        print(f"Surrogate training complete: {len(self.models)} GPs trained")

    def predict(self, X, return_std=False):
        """Predict flux time series for given parameters

        Args:
            X: Parameter values (n_samples, n_params) or (n_params,)
            return_std: If True, also return prediction uncertainty

        Returns:
            y_pred: Predicted flux (n_samples, n_times) or (n_times,)
            y_std: Prediction std (if return_std=True)
        """
        if not self.trained:
            raise ValueError("Surrogate model not trained")

        # Handle single sample
        single_sample = False
        if X.ndim == 1:
            X = X.reshape(1, -1)
            single_sample = True

        X_scaled = self.scaler_X.transform(X)

        # Predict at each time point
        y_pred = []
        y_std = [] if return_std else None

        for gp in self.models:
            if return_std:
                y_t_scaled, std_t_scaled = gp.predict(X_scaled, return_std=True)
                y_t = self.scaler_y.inverse_transform(y_t_scaled.reshape(-1, 1)).ravel()
                std_t = std_t_scaled * self.scaler_y.scale_[0]  # Approximate
                y_pred.append(y_t)
                y_std.append(std_t)
            else:
                y_t_scaled = gp.predict(X_scaled)
                y_t = self.scaler_y.inverse_transform(y_t_scaled.reshape(-1, 1)).ravel()
                y_pred.append(y_t)

        y_pred = np.array(y_pred).T

        if single_sample:
            y_pred = y_pred.ravel()
            if return_std:
                y_std = np.array(y_std).ravel()

        if return_std:
            return y_pred, np.array(y_std).T if not single_sample else y_std
        return y_pred

    def save(self, filepath):
        """Save surrogate model"""
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
        print(f"Surrogate model saved to {filepath}")

    @staticmethod
    def load(filepath):
        """Load surrogate model"""
        with open(filepath, 'rb') as f:
            model = pickle.load(f)
        print(f"Surrogate model loaded from {filepath}")
        return model


def run_bayesian_calibration(surrogate, exp_time, exp_flux, n_chains=4, n_samples=25000,
                            calibrate_params=None):
    """Run Bayesian calibration using PyMC

    Args:
        surrogate: Trained surrogate model
        exp_time: Experimental time points
        exp_flux: Experimental flux values
        n_chains: Number of MCMC chains
        n_samples: Number of samples per chain
        calibrate_params: List of parameters to calibrate (None = all)

    Returns:
        trace: PyMC trace object with posterior samples
    """
    problem = setup_problem()

    if calibrate_params is None:
        calibrate_params = problem['names']

    print(f"Running Bayesian calibration...")
    print(f"  Calibrating parameters: {calibrate_params}")
    print(f"  Chains: {n_chains}")
    print(f"  Samples per chain: {n_samples}")

    with pm.Model() as model:
        # Define priors
        params = {}
        for param_name in calibrate_params:
            bounds = PARAMETERS[param_name]['bounds']
            if PARAMETERS[param_name]['log_scale']:
                # Log-uniform prior
                log_bounds = [np.log(bounds[0]), np.log(bounds[1])]
                params[param_name] = pm.Uniform(param_name, lower=bounds[0], upper=bounds[1],
                                                testval=(bounds[0] * bounds[1]) ** 0.5)
            else:
                # Uniform prior
                params[param_name] = pm.Uniform(param_name, lower=bounds[0], upper=bounds[1],
                                               testval=(bounds[0] + bounds[1]) / 2)

        # Build parameter vector
        param_vector = pm.math.stack([params[name] for name in problem['names']])

        # Model prediction (using surrogate)
        # Note: This is a simplified approach - in practice, you'd need a PyMC-compatible surrogate
        # or use external likelihood evaluation
        # Here we'll use a custom likelihood

        # Observation error (unknown, so we estimate it)
        sigma = pm.HalfNormal('sigma', sigma=0.1 * np.mean(exp_flux))

        # Custom likelihood (simplified - would need proper implementation)
        # For now, this is a placeholder showing the structure

        print("Note: This is a template - full implementation requires integrating surrogate with PyMC")
        print("Consider using external sampling with emcee or custom Metropolis sampler")

        # # Likelihood
        # def logp_surrogate(theta):
        #     flux_pred = surrogate.predict(theta)
        #     return np.sum(stats.norm.logpdf(exp_flux, flux_pred, sigma))
        #
        # likelihood = pm.DensityDist('likelihood', logp_surrogate, observed={'theta': param_vector})

        # For demonstration, use simple prior sampling
        trace = pm.sample_prior_predictive(samples=1000)

    return trace


def analyze_calibration_results(trace, output_dir='calibration_results'):
    """Analyze and visualize calibration results

    Args:
        trace: PyMC trace object
        output_dir: Directory for plots and results
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\nAnalyzing calibration results...")

    # Summary statistics
    summary = az.summary(trace)
    print("\nPosterior Summary:")
    print(summary)

    summary_file = output_dir / 'posterior_summary.csv'
    summary.to_csv(summary_file)
    print(f"Saved summary to {summary_file}")

    # Trace plots
    fig = az.plot_trace(trace)
    plt.tight_layout()
    plt.savefig(output_dir / 'trace_plot.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Posterior distributions
    fig = az.plot_posterior(trace)
    plt.tight_layout()
    plt.savefig(output_dir / 'posterior_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Corner plot
    try:
        fig = az.plot_pair(trace, kind='kde', divergences=True)
        plt.tight_layout()
        plt.savefig(output_dir / 'corner_plot.png', dpi=300, bbox_inches='tight')
        plt.close()
    except:
        print("Corner plot failed (might need more samples)")

    print(f"\nResults saved to {output_dir}")


def main():
    parser = argparse.ArgumentParser(description='Bayesian calibration for Val-2a')
    parser.add_argument('--mode', required=True,
                       choices=['generate-ensemble', 'build-surrogate', 'calibrate', 'analyze'],
                       help='Calibration workflow step')

    # Ensemble generation args
    parser.add_argument('--n-samples', type=int, default=500,
                       help='Number of training samples for surrogate')
    parser.add_argument('--parallel', type=int, default=8,
                       help='Number of parallel processes')

    # Calibration args
    parser.add_argument('--chains', type=int, default=4,
                       help='Number of MCMC chains')
    parser.add_argument('--samples', type=int, default=25000,
                       help='Number of samples per chain')

    # File paths
    parser.add_argument('--tmap8-exe', type=str,
                       default=str(Path.home() / 'projects/TMAP8/tmap8-opt'))
    parser.add_argument('--input-file', type=str, default='val-2a_param.i')
    parser.add_argument('--ensemble-dir', type=str, default='ensemble_data')
    parser.add_argument('--surrogate-file', type=str, default='surrogate_model.pkl')
    parser.add_argument('--output-dir', type=str, default='calibration_results')

    args = parser.parse_args()

    if args.mode == 'generate-ensemble':
        print("="*80)
        print("STEP 1: GENERATE ENSEMBLE DATA")
        print("="*80)

        # Generate LHS samples
        print(f"\nGenerating {args.n_samples} Latin Hypercube samples...")
        samples_df = generate_lhs_samples(args.n_samples)

        samples_file = Path(args.ensemble_dir) / 'lhs_samples.csv'
        samples_file.parent.mkdir(parents=True, exist_ok=True)
        samples_df.to_csv(samples_file, index=False)
        print(f"Saved samples to {samples_file}")

        # Run ensemble
        print(f"\nRunning {args.n_samples} TMAP8 simulations in parallel...")
        results = run_tmap8_ensemble(samples_df, args.ensemble_dir, args.tmap8_exe,
                                     args.input_file, args.parallel)

        # Save results
        results_file = Path(args.ensemble_dir) / 'ensemble_results.pkl'
        with open(results_file, 'wb') as f:
            pickle.dump(results, f)
        print(f"Saved ensemble results to {results_file}")

        print(f"\nEnsemble generation complete: {len(results)}/{args.n_samples} successful")

    elif args.mode == 'build-surrogate':
        print("="*80)
        print("STEP 2: BUILD SURROGATE MODEL")
        print("="*80)

        # Load ensemble results
        results_file = Path(args.ensemble_dir) / 'ensemble_results.pkl'
        with open(results_file, 'rb') as f:
            ensemble_results = pickle.load(f)
        print(f"Loaded {len(ensemble_results)} ensemble results")

        # Create and train surrogate
        surrogate = SurrogateModel()
        X, y = surrogate.prepare_training_data(ensemble_results)
        surrogate.train(X, y)

        # Save surrogate
        surrogate.save(args.surrogate_file)

        # Validation
        print("\nValidating surrogate model...")
        # TODO: Add validation metrics (R², RMSE, etc.)

    elif args.mode == 'calibrate':
        print("="*80)
        print("STEP 3: BAYESIAN CALIBRATION")
        print("="*80)

        # Load surrogate
        surrogate = SurrogateModel.load(args.surrogate_file)

        # Load experimental data
        exp_time, exp_flux = load_experimental_data()

        # Run calibration
        print("\nNOTE: Full PyMC integration requires custom likelihood implementation")
        print("Consider using emcee or custom MCMC for surrogate-based calibration")
        print("\nSee bayesian_calibration_emcee.py for working implementation")

        # trace = run_bayesian_calibration(surrogate, exp_time, exp_flux,
        #                                  args.chains, args.samples)
        #
        # # Save trace
        # trace_file = Path(args.output_dir) / 'posterior_trace.nc'
        # trace.to_netcdf(trace_file)

    elif args.mode == 'analyze':
        print("="*80)
        print("STEP 4: ANALYZE RESULTS")
        print("="*80)

        # Load trace
        trace_file = Path(args.output_dir) / 'posterior_trace.nc'
        trace = az.from_netcdf(trace_file)

        # Analyze
        analyze_calibration_results(trace, args.output_dir)


if __name__ == '__main__':
    main()
