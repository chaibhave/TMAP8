#!/usr/bin/env python3
"""
Uncertainty Propagation for Calibrated Model
Runs TMAP8 with posterior samples to generate prediction uncertainty bands

Usage:
    python uncertainty_propagation.py --n-samples 100 --parallel 8
"""

import numpy as np
import pandas as pd
import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

# Set style
sns.set_style("whitegrid")
sns.set_context("talk")


def run_tmap8_with_params(args_tuple):
    """Run TMAP8 with given parameters

    Args:
        args_tuple: (sample_id, Kr_left_fraction, gaussian_factor, tmap8_exe, input_file, output_dir)

    Returns:
        Dictionary with sample_id and flux time series
    """
    sample_id, kr_frac, gauss_fac, tmap8_exe, input_file, output_dir = args_tuple

    run_dir = output_dir / f"sample_{sample_id:04d}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Build command
    input_file_abs = Path(input_file).resolve()
    cmd = [
        str(tmap8_exe), '-i', str(input_file_abs),
        f'Kr_left_fraction={kr_frac}',
        f'gaussian_factor={gauss_fac}',
        '--no-color'
    ]

    try:
        result = subprocess.run(cmd, cwd=run_dir, capture_output=True, timeout=600, text=True)

        if result.returncode != 0:
            print(f"Sample {sample_id} failed")
            return None

        # Load output
        csv_file = run_dir / 'val-2a_param_reduced_out.csv'
        if not csv_file.exists():
            return None

        df = pd.read_csv(csv_file)

        return {
            'sample_id': sample_id,
            'time': df['time'].values,
            'flux': df['scaled_recombination_flux_right'].values,
            'kr_frac': kr_frac,
            'gauss_fac': gauss_fac
        }

    except Exception as e:
        print(f"Sample {sample_id} error: {e}")
        return None


def run_ensemble(posterior_samples, n_samples, tmap8_exe, input_file, output_dir, n_parallel=8):
    """Run TMAP8 ensemble with posterior samples

    Args:
        posterior_samples: Array of shape (n_total, 2) with [Kr_left_fraction, gaussian_factor]
        n_samples: Number of samples to run
        tmap8_exe: Path to TMAP8 executable
        input_file: TMAP8 input file
        output_dir: Directory for outputs
        n_parallel: Number of parallel processes

    Returns:
        List of result dictionaries
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Randomly select samples from posterior
    n_total = len(posterior_samples)
    if n_samples > n_total:
        n_samples = n_total
        print(f"Warning: Requested {n_samples} samples but only {n_total} available")

    indices = np.random.choice(n_total, size=n_samples, replace=False)
    selected_samples = posterior_samples[indices]

    print(f"Running {n_samples} TMAP8 simulations with posterior samples...")
    print(f"  Kr_left_fraction range: [{selected_samples[:, 0].min():.6f}, {selected_samples[:, 0].max():.6f}]")
    print(f"  gaussian_factor range:  [{selected_samples[:, 1].min():.4f}, {selected_samples[:, 1].max():.4f}]")

    # Prepare arguments
    args_list = [
        (i, params[0], params[1], tmap8_exe, input_file, output_dir)
        for i, params in enumerate(selected_samples)
    ]

    # Run in parallel
    results = []
    with ProcessPoolExecutor(max_workers=n_parallel) as executor:
        futures = [executor.submit(run_tmap8_with_params, args) for args in args_list]

        for future in tqdm(as_completed(futures), total=len(futures), desc="Running ensemble"):
            result = future.result()
            if result is not None:
                results.append(result)

    print(f"\nCompleted {len(results)}/{n_samples} simulations successfully")
    return results


def calculate_uncertainty_bands(results, percentiles=[5, 50, 95]):
    """Calculate uncertainty bands from ensemble results

    Args:
        results: List of result dictionaries
        percentiles: Percentiles to calculate

    Returns:
        Dictionary with time and percentile arrays
    """
    # Use first result's time array as reference
    ref_time = results[0]['time']

    # Interpolate all flux arrays to reference time
    from scipy.interpolate import interp1d

    flux_ensemble = []
    for result in results:
        if len(result['time']) != len(ref_time) or not np.allclose(result['time'], ref_time):
            # Interpolate to reference time
            f = interp1d(result['time'], result['flux'], kind='linear',
                        bounds_error=False, fill_value='extrapolate')
            flux_interp = f(ref_time)
        else:
            flux_interp = result['flux']
        flux_ensemble.append(flux_interp)

    flux_ensemble = np.array(flux_ensemble)  # Shape: (n_samples, n_times)

    # Calculate percentiles
    percentile_data = {}
    for p in percentiles:
        percentile_data[f'p{p}'] = np.percentile(flux_ensemble, p, axis=0)

    return {
        'time': ref_time,
        'percentiles': percentile_data,
        'mean': np.mean(flux_ensemble, axis=0),
        'std': np.std(flux_ensemble, axis=0)
    }


def plot_with_uncertainty(uncertainty_data, exp_data_file, output_file='val-2a_with_uncertainty.png'):
    """Plot model predictions with uncertainty bands

    Args:
        uncertainty_data: Dictionary from calculate_uncertainty_bands
        exp_data_file: Path to experimental data CSV
        output_file: Output plot filename
    """
    # Load experimental data
    exp_data = pd.read_csv(exp_data_file)
    exp_time = exp_data['time (s)'].values
    exp_flux = exp_data['permeation flux (atom/m^2/s)'].values

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    time_hrs = uncertainty_data['time'] / 3600

    # Plot uncertainty bands (95% credible interval)
    ax.fill_between(
        time_hrs,
        uncertainty_data['percentiles']['p5'],
        uncertainty_data['percentiles']['p95'],
        alpha=0.3,
        color='tab:green',
        label='95% Credible Interval'
    )

    # Plot median
    ax.plot(
        time_hrs,
        uncertainty_data['percentiles']['p50'],
        color='tab:green',
        linewidth=2.5,
        label='Calibrated Model (Median)'
    )

    # Plot experimental data
    ax.plot(
        exp_time / 3600,
        exp_flux,
        'k--',
        linewidth=2,
        label='Experimental Data'
    )

    # Formatting
    ax.set_xlabel('Time (hr)', fontsize=14)
    ax.set_ylabel('Deuterium Flux (atom/m²/s)', fontsize=14)
    ax.set_xlim(left=-0.1, right=20000 / 3600)
    ax.set_ylim(bottom=0)
    ax.legend(loc='upper right', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(15, 15))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def plot_comparison_with_uncertainty(uncertainty_data, exp_data_file,
                                     gold_file, claude_file,
                                     output_file='val-2a_comparison_with_uncertainty.png'):
    """Plot comparison including uncertainty bands

    Args:
        uncertainty_data: Dictionary from calculate_uncertainty_bands
        exp_data_file: Path to experimental data
        gold_file: Path to gold simulation CSV
        claude_file: Path to claude-code simulation CSV
        output_file: Output filename
    """
    # Load data
    exp_data = pd.read_csv(exp_data_file)
    exp_time = exp_data['time (s)'].values / 3600
    exp_flux = exp_data['permeation flux (atom/m^2/s)'].values

    gold_data = pd.read_csv(gold_file)
    gold_time = gold_data['time'].values / 3600
    gold_flux = gold_data['scaled_recombination_flux_right'].values

    claude_data = pd.read_csv(claude_file)
    claude_time = claude_data['time'].values / 3600
    claude_flux = claude_data['scaled_recombination_flux_right'].values

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    time_hrs = uncertainty_data['time'] / 3600

    # Plot uncertainty bands FIRST (so they're in background)
    ax.fill_between(
        time_hrs,
        uncertainty_data['percentiles']['p5'],
        uncertainty_data['percentiles']['p95'],
        alpha=0.25,
        color='tab:green',
        label='Calibrated 95% CI',
        zorder=1
    )

    # Plot simulations
    ax.plot(gold_time, gold_flux, '-', color='tab:gray', linewidth=2,
            label='TMAP8 (Gold)', zorder=2)
    ax.plot(claude_time, claude_flux, '-', color='tab:blue', linewidth=2,
            label='TMAP8 (claude-code)', zorder=2)
    ax.plot(time_hrs, uncertainty_data['percentiles']['p50'], '-',
            color='tab:green', linewidth=2.5, label='TMAP8 (Calibrated - Median)', zorder=3)

    # Plot experimental data
    ax.plot(exp_time, exp_flux, 'k--', linewidth=2, label='Experiment', zorder=4)

    # Formatting
    ax.set_xlabel('Time (hr)', fontsize=14)
    ax.set_ylabel('Deuterium Flux (atom/m²/s)', fontsize=14)
    ax.set_xlim(left=-0.1, right=20000 / 3600)
    ax.set_ylim(bottom=0)
    ax.legend(loc='upper right', fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(15, 15))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Uncertainty propagation for calibrated model')
    parser.add_argument('--n-samples', type=int, default=100,
                       help='Number of posterior samples to run')
    parser.add_argument('--parallel', type=int, default=8,
                       help='Number of parallel processes')
    parser.add_argument('--posterior-file', type=str,
                       default='calibration_reduced_results/mcmc_samples.npy',
                       help='Path to posterior samples file')
    parser.add_argument('--tmap8-exe', type=str,
                       default=str(Path.home() / 'projects/TMAP8/tmap8-opt'),
                       help='Path to TMAP8 executable')
    parser.add_argument('--input-file', type=str, default='val-2a_param_reduced.i',
                       help='TMAP8 input file')
    parser.add_argument('--output-dir', type=str, default='uncertainty_propagation',
                       help='Directory for ensemble outputs')

    args = parser.parse_args()

    print("="*80)
    print("UNCERTAINTY PROPAGATION FOR CALIBRATED MODEL")
    print("="*80)
    print()

    # Load posterior samples
    print("Loading posterior samples...")
    posterior_samples = np.load(args.posterior_file)
    print(f"  Loaded {len(posterior_samples)} samples")
    print(f"  Parameters: Kr_left_fraction, gaussian_factor")
    print()

    # Run ensemble
    results = run_ensemble(
        posterior_samples, args.n_samples,
        args.tmap8_exe, args.input_file, args.output_dir, args.parallel
    )

    if len(results) == 0:
        print("ERROR: No successful simulations!")
        return

    # Calculate uncertainty bands
    print("\nCalculating uncertainty bands...")
    uncertainty_data = calculate_uncertainty_bands(results)

    # Save results
    output_dir = Path(args.output_dir)
    np.savez(
        output_dir / 'uncertainty_bands.npz',
        time=uncertainty_data['time'],
        mean=uncertainty_data['mean'],
        std=uncertainty_data['std'],
        **uncertainty_data['percentiles']
    )
    print(f"Saved: {output_dir / 'uncertainty_bands.npz'}")

    # Generate plots
    print("\nGenerating plots...")
    plot_with_uncertainty(uncertainty_data, 'gold/experiment_data_paper.csv')

    # If comparison files exist, make comparison plot
    if Path('gold/val-2a_out.csv').exists() and Path('val-2a_out.csv').exists():
        plot_comparison_with_uncertainty(
            uncertainty_data,
            'gold/experiment_data_paper.csv',
            'gold/val-2a_out.csv',
            'val-2a_out.csv'
        )

    print("\n" + "="*80)
    print("UNCERTAINTY PROPAGATION COMPLETE!")
    print("="*80)
    print(f"Generated {len(results)} model predictions from posterior")
    print(f"Plots saved showing 95% credible interval")
    print("="*80)


if __name__ == '__main__':
    main()
