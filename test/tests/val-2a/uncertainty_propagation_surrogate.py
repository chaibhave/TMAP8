#!/usr/bin/env python3
"""
Fast Uncertainty Propagation Using Surrogate Model
Uses the trained GP surrogate to quickly generate prediction uncertainty bands

Usage:
    python uncertainty_propagation_surrogate.py --n-samples 1000
"""

import numpy as np
import pandas as pd
import pickle
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
from bayesian_calibration_reduced import SurrogateModel

sns.set_style("whitegrid")
sns.set_context("talk")


def main():
    parser = argparse.ArgumentParser(description='Fast uncertainty propagation via surrogate')
    parser.add_argument('--n-samples', type=int, default=1000,
                       help='Number of posterior samples to propagate')
    parser.add_argument('--posterior-file', type=str,
                       default='calibration_reduced_results/mcmc_samples.npy')
    parser.add_argument('--surrogate-file', type=str,
                       default='calibration_reduced_results/surrogate_model.pkl')

    args = parser.parse_args()

    print("="*80)
    print("FAST UNCERTAINTY PROPAGATION (Surrogate-based)")
    print("="*80)
    print()

    # Load surrogate
    print("Loading surrogate model...")
    with open(args.surrogate_file, 'rb') as f:
        surrogate = pickle.load(f)
    print("✓ Surrogate loaded")

    # Load posterior samples
    print("Loading posterior samples...")
    posterior = np.load(args.posterior_file)
    print(f"✓ Loaded {len(posterior)} samples")
    print(f"  Kr_left_fraction: {posterior[:, 0].mean():.6f} ± {posterior[:, 0].std():.6f}")
    print(f"  gaussian_factor:  {posterior[:, 1].mean():.6f} ± {posterior[:, 1].std():.6f}")
    print()

    # Randomly select samples
    n_samples = min(args.n_samples, len(posterior))
    indices = np.random.choice(len(posterior), size=n_samples, replace=False)
    selected_samples = posterior[indices]

    # Predict using surrogate
    print(f"Propagating {n_samples} samples through surrogate...")
    predictions = surrogate.predict(selected_samples)  # Shape: (n_samples, n_times_subsampled)
    print("✓ Predictions complete")

    # Load experimental data for reference time
    exp_data = pd.read_csv('gold/experiment_data_paper.csv')
    exp_time = exp_data['time (s)'].values
    exp_flux = exp_data['permeation flux (atom/m^2/s)'].values

    # Interpolate predictions to full experimental time grid
    ref_time = surrogate.time_indices
    predictions_full = []

    print("Interpolating to full time grid...")
    for pred in predictions:
        # Interpolate from subsampled times to experimental times
        f = interp1d(exp_time[ref_time], pred, kind='cubic',
                    bounds_error=False, fill_value='extrapolate')
        pred_full = f(exp_time)
        predictions_full.append(pred_full)

    predictions_full = np.array(predictions_full)

    # Calculate percentiles
    print("Calculating uncertainty bands...")
    percentiles = {
        'p5': np.percentile(predictions_full, 5, axis=0),
        'p50': np.percentile(predictions_full, 50, axis=0),
        'p95': np.percentile(predictions_full, 95, axis=0),
        'mean': np.mean(predictions_full, axis=0),
        'std': np.std(predictions_full, axis=0)
    }

    # Save results
    output_dir = Path('uncertainty_propagation_surrogate')
    output_dir.mkdir(exist_ok=True)
    np.savez(
        output_dir / 'uncertainty_bands.npz',
        time=exp_time,
        **percentiles
    )
    print(f"✓ Saved: {output_dir / 'uncertainty_bands.npz'}")

    # Generate plots
    print("\nGenerating plots...")

    # Plot 1: Just calibrated model with uncertainty
    fig, ax = plt.subplots(figsize=(10, 6))

    time_hrs = exp_time / 3600

    # 95% credible interval
    ax.fill_between(
        time_hrs,
        percentiles['p5'],
        percentiles['p95'],
        alpha=0.3,
        color='tab:green',
        label='95% Credible Interval'
    )

    # Median
    ax.plot(
        time_hrs,
        percentiles['p50'],
        color='tab:green',
        linewidth=2.5,
        label='Calibrated Model (Median)'
    )

    # Experimental data
    ax.plot(
        time_hrs,
        exp_flux,
        'k--',
        linewidth=2,
        label='Experimental Data'
    )

    ax.set_xlabel('Time (hr)', fontsize=14)
    ax.set_ylabel('Deuterium Flux (atom/m²/s)', fontsize=14)
    ax.set_xlim(left=-0.1, right=20000 / 3600)
    ax.set_ylim(bottom=0)
    ax.legend(loc='upper right', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(15, 15))

    plt.tight_layout()
    plt.savefig('val-2a_with_uncertainty.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: val-2a_with_uncertainty.png")
    plt.close()

    # Plot 2: Comparison with all models
    if Path('gold/val-2a_out.csv').exists() and Path('val-2a_out.csv').exists():
        gold_data = pd.read_csv('gold/val-2a_out.csv')
        claude_data = pd.read_csv('val-2a_out.csv')

        fig, ax = plt.subplots(figsize=(10, 6))

        # Uncertainty band
        ax.fill_between(
            time_hrs,
            percentiles['p5'],
            percentiles['p95'],
            alpha=0.25,
            color='tab:green',
            label='Calibrated 95% CI',
            zorder=1
        )

        # Simulations
        ax.plot(gold_data['time'] / 3600, gold_data['scaled_recombination_flux_right'],
                '-', color='tab:gray', linewidth=2, label='TMAP8 (Gold)', zorder=2)
        ax.plot(claude_data['time'] / 3600, claude_data['scaled_recombination_flux_right'],
                '-', color='tab:blue', linewidth=2, label='TMAP8 (claude-code)', zorder=2)
        ax.plot(time_hrs, percentiles['p50'], '-', color='tab:green', linewidth=2.5,
                label='TMAP8 (Calibrated - Median)', zorder=3)

        # Experimental
        ax.plot(time_hrs, exp_flux, 'k--', linewidth=2, label='Experiment', zorder=4)

        ax.set_xlabel('Time (hr)', fontsize=14)
        ax.set_ylabel('Deuterium Flux (atom/m²/s)', fontsize=14)
        ax.set_xlim(left=-0.1, right=20000 / 3600)
        ax.set_ylim(bottom=0)
        ax.legend(loc='upper right', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(15, 15))

        plt.tight_layout()
        plt.savefig('val-2a_comparison_with_uncertainty.png', dpi=300, bbox_inches='tight')
        print(f"✓ Saved: val-2a_comparison_with_uncertainty.png")
        plt.close()

    print("\n" + "="*80)
    print("UNCERTAINTY PROPAGATION COMPLETE!")
    print("="*80)
    print(f"Propagated {n_samples} posterior samples through surrogate")
    print(f"Generated uncertainty bands in seconds (vs hours for full simulations)")
    print("="*80)


if __name__ == '__main__':
    main()
