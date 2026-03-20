#!/usr/bin/env python3
"""
Plot Sensitivity Analysis Results for Val-2a Calibration

Generates comprehensive visualizations of Sobol sensitivity analysis results

Usage:
    python plot_sensitivity_results.py --input-dir sensitivity_runs
    python plot_sensitivity_results.py --input-dir test_sensitivity --output-dir plots
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
from pathlib import Path

# Set style
sns.set_style("whitegrid")
sns.set_context("talk")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

# Parameter display names (prettier for plots)
PARAM_DISPLAY_NAMES = {
    'diffusivity': 'Diffusivity',
    'Kr_left_max': r'$K_r^{left,max}$',
    'Kr_right': r'$K_r^{right}$',
    'gaussian_factor': 'Gaussian Factor',
    'Kr_left_time_constant': r'$\tau_{cleanup}$',
    'Kr_left_fraction': 'Cleanup Fraction',
    'implantation_depth': 'Implant Depth',
    'implantation_sigma': 'Implant Straggling'
}


def load_sensitivity_data(input_dir):
    """Load sensitivity analysis results"""
    input_dir = Path(input_dir)

    # Load Sobol indices
    si_file = input_dir / 'sensitivity_peak_flux.json'
    with open(si_file, 'r') as f:
        Si = json.load(f)

    # Load simulation results
    results = pd.read_csv(input_dir / 'results.csv')

    # Load parameter samples
    samples = pd.read_csv(input_dir / 'samples.csv')

    return Si, results, samples


def plot_sobol_indices(Si, output_dir):
    """Plot Sobol sensitivity indices (S1 and ST)"""
    param_names = Si['parameter_names']
    display_names = [PARAM_DISPLAY_NAMES.get(p, p) for p in param_names]

    S1 = np.array(Si['S1'])
    S1_conf = np.array(Si['S1_conf'])
    ST = np.array(Si['ST'])
    ST_conf = np.array(Si['ST_conf'])

    # Sort by total-order indices
    sort_idx = np.argsort(ST)[::-1]
    S1 = S1[sort_idx]
    S1_conf = S1_conf[sort_idx]
    ST = ST[sort_idx]
    ST_conf = ST_conf[sort_idx]
    display_names_sorted = [display_names[i] for i in sort_idx]

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    x = np.arange(len(display_names_sorted))
    width = 0.35

    # Plot S1 (first-order)
    bars1 = ax1.barh(x, S1, xerr=S1_conf, height=width,
                     color='steelblue', alpha=0.8, label='First-order (S1)')
    ax1.set_yticks(x)
    ax1.set_yticklabels(display_names_sorted)
    ax1.set_xlabel('Sensitivity Index', fontsize=14, fontweight='bold')
    ax1.set_title('First-Order Sensitivity (Direct Effect)', fontsize=16, fontweight='bold')
    ax1.axvline(0, color='black', linestyle='-', linewidth=0.8)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=12)

    # Add value labels
    for i, (val, conf) in enumerate(zip(S1, S1_conf)):
        if abs(val) > 0.01:
            ax1.text(val + conf + 0.02, i, f'{val:.3f}',
                    va='center', fontsize=10, fontweight='bold')

    # Plot ST (total-order)
    bars2 = ax2.barh(x, ST, xerr=ST_conf, height=width,
                     color='coral', alpha=0.8, label='Total-order (ST)')
    ax2.set_yticks(x)
    ax2.set_yticklabels(display_names_sorted)
    ax2.set_xlabel('Sensitivity Index', fontsize=14, fontweight='bold')
    ax2.set_title('Total-Order Sensitivity (Total Effect)', fontsize=16, fontweight='bold')
    ax2.axvline(0, color='black', linestyle='-', linewidth=0.8)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=12)

    # Add value labels
    for i, (val, conf) in enumerate(zip(ST, ST_conf)):
        if abs(val) > 0.01:
            ax2.text(val + conf + 0.02, i, f'{val:.3f}',
                    va='center', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_dir / 'sobol_indices.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir / 'sobol_indices.png'}")
    plt.close()


def plot_sobol_comparison(Si, output_dir):
    """Plot S1 vs ST comparison"""
    param_names = Si['parameter_names']
    display_names = [PARAM_DISPLAY_NAMES.get(p, p) for p in param_names]

    S1 = np.array(Si['S1'])
    ST = np.array(Si['ST'])

    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot S1 vs ST
    ax.scatter(S1, ST, s=200, alpha=0.7, c=range(len(S1)), cmap='viridis', edgecolors='black', linewidth=2)

    # Add parameter labels
    for i, name in enumerate(display_names):
        ax.annotate(name, (S1[i], ST[i]), xytext=(5, 5), textcoords='offset points',
                   fontsize=11, fontweight='bold', alpha=0.8)

    # Add diagonal line (S1 = ST means no interactions)
    max_val = max(max(abs(S1)), max(abs(ST)))
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=2, label='S1 = ST (no interaction)')

    ax.set_xlabel('First-Order Index (S1)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Total-Order Index (ST)', fontsize=14, fontweight='bold')
    ax.set_title('Sobol Indices Comparison\n(Distance from diagonal = interaction strength)',
                fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=12)
    ax.axhline(0, color='black', linestyle='-', linewidth=0.8)
    ax.axvline(0, color='black', linestyle='-', linewidth=0.8)

    plt.tight_layout()
    plt.savefig(output_dir / 'sobol_comparison.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir / 'sobol_comparison.png'}")
    plt.close()


def plot_parameter_distributions(samples, output_dir):
    """Plot distributions of sampled parameters"""
    n_params = len(samples.columns)
    n_cols = 4
    n_rows = int(np.ceil(n_params / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, n_rows * 4))
    axes = axes.flatten()

    for i, col in enumerate(samples.columns):
        ax = axes[i]
        display_name = PARAM_DISPLAY_NAMES.get(col, col)

        # Histogram
        ax.hist(samples[col], bins=30, alpha=0.7, color='steelblue', edgecolor='black')
        ax.set_xlabel(display_name, fontsize=12, fontweight='bold')
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title(f'{display_name}\n[{samples[col].min():.2e}, {samples[col].max():.2e}]',
                    fontsize=11)
        ax.grid(True, alpha=0.3)

        # Add mean line
        mean_val = samples[col].mean()
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.2e}')
        ax.legend(fontsize=9)

    # Hide extra subplots
    for i in range(n_params, len(axes)):
        axes[i].axis('off')

    plt.tight_layout()
    plt.savefig(output_dir / 'parameter_distributions.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir / 'parameter_distributions.png'}")
    plt.close()


def plot_output_metrics(results, output_dir):
    """Plot distribution of output metrics"""
    successful = results[results['success'] == True]

    metrics = {
        'Peak Flux': ('peak_flux', 'atoms/m²/s'),
        'Total Permeation': ('total_permeation', 'atoms/m²'),
        'Final Inventory': ('final_inventory', 'atoms'),
        'Peak Time': ('peak_time', 's')
    }

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    for i, (title, (col, units)) in enumerate(metrics.items()):
        ax = axes[i]

        data = successful[col]

        # Histogram
        n, bins, patches = ax.hist(data, bins=30, alpha=0.7, color='coral', edgecolor='black')

        # Add statistics
        mean_val = data.mean()
        std_val = data.std()

        ax.axvline(mean_val, color='red', linestyle='--', linewidth=3, label=f'Mean: {mean_val:.2e}')
        ax.axvline(mean_val - std_val, color='orange', linestyle=':', linewidth=2, alpha=0.7)
        ax.axvline(mean_val + std_val, color='orange', linestyle=':', linewidth=2, alpha=0.7,
                  label=f'±1σ: {std_val:.2e}')

        ax.set_xlabel(f'{title} ({units})', fontsize=13, fontweight='bold')
        ax.set_ylabel('Count', fontsize=13, fontweight='bold')
        ax.set_title(f'{title} Distribution\nRange: [{data.min():.2e}, {data.max():.2e}]',
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=11)

    plt.tight_layout()
    plt.savefig(output_dir / 'output_distributions.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir / 'output_distributions.png'}")
    plt.close()


def plot_parameter_output_correlations(samples, results, output_dir):
    """Plot correlations between parameters and peak flux"""
    successful = results[results['success'] == True].reset_index(drop=True)

    # Merge samples with results
    combined = samples.copy()
    combined['peak_flux'] = successful['peak_flux'].values

    n_params = len(samples.columns)
    n_cols = 4
    n_rows = int(np.ceil(n_params / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, n_rows * 4))
    axes = axes.flatten()

    for i, col in enumerate(samples.columns):
        ax = axes[i]
        display_name = PARAM_DISPLAY_NAMES.get(col, col)

        # Scatter plot
        ax.scatter(combined[col], combined['peak_flux'], alpha=0.6, s=50,
                  c=combined['peak_flux'], cmap='viridis', edgecolors='black', linewidth=0.5)

        # Calculate correlation
        corr = np.corrcoef(combined[col], combined['peak_flux'])[0, 1]

        ax.set_xlabel(display_name, fontsize=12, fontweight='bold')
        ax.set_ylabel('Peak Flux (atoms/m²/s)', fontsize=12)
        ax.set_title(f'{display_name}\nCorrelation: {corr:.3f}', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Add trend line if correlation is significant
        if abs(corr) > 0.3:
            z = np.polyfit(combined[col], combined['peak_flux'], 1)
            p = np.poly1d(z)
            x_trend = np.linspace(combined[col].min(), combined[col].max(), 100)
            ax.plot(x_trend, p(x_trend), 'r--', linewidth=2, alpha=0.8)

    # Hide extra subplots
    for i in range(n_params, len(axes)):
        axes[i].axis('off')

    plt.tight_layout()
    plt.savefig(output_dir / 'parameter_output_correlations.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir / 'parameter_output_correlations.png'}")
    plt.close()


def plot_variance_decomposition(Si, output_dir):
    """Plot variance decomposition pie chart"""
    param_names = Si['parameter_names']
    display_names = [PARAM_DISPLAY_NAMES.get(p, p) for p in param_names]
    ST = np.array(Si['ST'])

    # Normalize and handle negative values
    ST_positive = np.maximum(ST, 0)
    total = ST_positive.sum()

    if total > 0:
        ST_normalized = ST_positive / total * 100
    else:
        ST_normalized = np.ones(len(ST)) / len(ST) * 100

    # Sort by contribution
    sort_idx = np.argsort(ST_normalized)[::-1]
    ST_sorted = ST_normalized[sort_idx]
    names_sorted = [display_names[i] for i in sort_idx]

    # Group small contributions
    threshold = 2.0  # percent
    major_mask = ST_sorted >= threshold
    major_values = ST_sorted[major_mask]
    major_names = [names_sorted[i] for i in range(len(names_sorted)) if major_mask[i]]

    if (~major_mask).any():
        minor_sum = ST_sorted[~major_mask].sum()
        major_values = np.append(major_values, minor_sum)
        major_names.append(f'Others (<{threshold}%)')

    # Create pie chart
    fig, ax = plt.subplots(figsize=(12, 12))

    colors = sns.color_palette("husl", len(major_values))
    wedges, texts, autotexts = ax.pie(major_values, labels=major_names, autopct='%1.1f%%',
                                       colors=colors, startangle=90, textprops={'fontsize': 12})

    # Bold percentage text
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(13)

    # Bold labels
    for text in texts:
        text.set_fontweight('bold')
        text.set_fontsize(12)

    ax.set_title('Output Variance Decomposition\n(Based on Total-Order Sobol Indices)',
                fontsize=16, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(output_dir / 'variance_decomposition.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir / 'variance_decomposition.png'}")
    plt.close()


def create_summary_report(Si, results, output_dir):
    """Create text summary report"""
    report_file = output_dir / 'sensitivity_summary.txt'

    param_names = Si['parameter_names']
    S1 = np.array(Si['S1'])
    ST = np.array(Si['ST'])
    ST_conf = np.array(Si['ST_conf'])

    successful = results[results['success'] == True]

    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("SENSITIVITY ANALYSIS SUMMARY REPORT\n")
        f.write("="*80 + "\n\n")

        f.write(f"Simulations: {len(successful)}/{len(results)} successful ({100*len(successful)/len(results):.1f}%)\n\n")

        f.write("Output Metrics:\n")
        f.write("-"*80 + "\n")
        f.write(f"Peak Flux:          {successful['peak_flux'].mean():.3e} ± {successful['peak_flux'].std():.3e} atoms/m²/s\n")
        f.write(f"Total Permeation:   {successful['total_permeation'].mean():.3e} ± {successful['total_permeation'].std():.3e} atoms/m²\n")
        f.write(f"Final Inventory:    {successful['final_inventory'].mean():.3e} ± {successful['final_inventory'].std():.3e} atoms\n\n")

        f.write("Parameter Sensitivity (Total-Order Indices):\n")
        f.write("-"*80 + "\n")

        # Sort by ST
        sort_idx = np.argsort(ST)[::-1]

        for idx in sort_idx:
            param = param_names[idx]
            st = ST[idx]
            conf = ST_conf[idx]

            if st > 0.1:
                importance = "HIGH"
            elif st > 0.05:
                importance = "MODERATE"
            else:
                importance = "LOW"

            f.write(f"{param:25s}  ST = {st:7.4f} ± {conf:7.4f}  [{importance}]\n")

        f.write("\n" + "="*80 + "\n")
        f.write("KEY FINDINGS\n")
        f.write("="*80 + "\n\n")

        f.write("Most Influential Parameters (ST > 0.05):\n")
        for idx in sort_idx:
            if ST[idx] > 0.05:
                f.write(f"  • {param_names[idx]}: ST = {ST[idx]:.4f}\n")

        f.write("\nRecommendation for Calibration:\n")
        f.write("  Focus on parameters with ST > 0.05 for efficient calibration.\n")
        f.write("  Consider fixing parameters with ST < 0.01 to their nominal values.\n")

        f.write("\n" + "="*80 + "\n")

    print(f"Saved: {report_file}")


def main():
    parser = argparse.ArgumentParser(description='Plot sensitivity analysis results')
    parser.add_argument('--input-dir', type=str, default='sensitivity_runs',
                       help='Directory with sensitivity analysis results')
    parser.add_argument('--output-dir', type=str, default='sensitivity_plots',
                       help='Directory for output plots')

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*80)
    print("PLOTTING SENSITIVITY ANALYSIS RESULTS")
    print("="*80)
    print(f"\nInput directory:  {input_dir}")
    print(f"Output directory: {output_dir}\n")

    # Load data
    print("Loading data...")
    Si, results, samples = load_sensitivity_data(input_dir)
    print(f"  Loaded {len(results)} simulation results")
    print(f"  Loaded {len(samples)} parameter samples")
    print(f"  Parameters: {len(Si['parameter_names'])}\n")

    # Generate plots
    print("Generating plots...\n")

    print("1. Sobol indices...")
    plot_sobol_indices(Si, output_dir)

    print("2. Sobol comparison (S1 vs ST)...")
    plot_sobol_comparison(Si, output_dir)

    print("3. Parameter distributions...")
    plot_parameter_distributions(samples, output_dir)

    print("4. Output metric distributions...")
    plot_output_metrics(results, output_dir)

    print("5. Parameter-output correlations...")
    plot_parameter_output_correlations(samples, results, output_dir)

    print("6. Variance decomposition...")
    plot_variance_decomposition(Si, output_dir)

    print("\nCreating summary report...")
    create_summary_report(Si, results, output_dir)

    print("\n" + "="*80)
    print("PLOTTING COMPLETE")
    print("="*80)
    print(f"\nAll plots saved to: {output_dir}")
    print("\nGenerated files:")
    print("  • sobol_indices.png - Bar charts of S1 and ST")
    print("  • sobol_comparison.png - S1 vs ST scatter plot")
    print("  • parameter_distributions.png - Histograms of sampled parameters")
    print("  • output_distributions.png - Histograms of output metrics")
    print("  • parameter_output_correlations.png - Parameter vs peak flux")
    print("  • variance_decomposition.png - Pie chart of variance contributions")
    print("  • sensitivity_summary.txt - Text summary report")
    print("\n" + "="*80)


if __name__ == '__main__':
    main()
