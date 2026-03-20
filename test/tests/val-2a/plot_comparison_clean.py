#!/usr/bin/env python3
"""
Production-ready comparison plot for val-2a validation test
Follows publication standards with proper styling
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns

# Set up publication-quality style
def setup_publication_style():
    """Configure matplotlib for publication-quality figures."""
    plt.rcParams.update({
        # Font settings - use standard fonts for compatibility
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,

        # Figure settings
        "figure.figsize": (8, 5),
        "figure.dpi": 100,
        "savefig.dpi": 300,
        "savefig.format": "png",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.1,

        # Line and marker settings
        "lines.linewidth": 2.0,
        "lines.markersize": 4,

        # Axes settings
        "axes.linewidth": 1.0,
        "axes.grid": True,
        "axes.axisbelow": True,

        # Grid settings
        "grid.linewidth": 0.5,
        "grid.alpha": 0.3,
        "grid.color": "gray",
        "grid.linestyle": "--",

        # Tick settings
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.direction": "in",
        "ytick.direction": "in",

        # Legend settings
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "legend.fancybox": False,
        "legend.edgecolor": "gray",
    })

setup_publication_style()

# Load data
print("Loading simulation and experimental data...")
gold_data = pd.read_csv('gold/val-2a_out.csv')
claude_data = pd.read_csv('val-2a_out.csv')
calib_data = pd.read_csv('val-2a_param_reduced_out.csv')
exp_data = pd.read_csv('gold/experiment_data_paper.csv')

# Extract data
gold_time = gold_data['time'].values / 3600  # Convert to hours
gold_flux = gold_data['scaled_recombination_flux_right'].values

claude_time = claude_data['time'].values / 3600
claude_flux = claude_data['scaled_recombination_flux_right'].values

calib_time = calib_data['time'].values / 3600
calib_flux = calib_data['scaled_recombination_flux_right'].values

exp_time = exp_data['time (s)'].values / 3600
exp_flux = exp_data['permeation flux (atom/m^2/s)'].values

# Create figure
fig, ax = plt.subplots(figsize=(8, 5.5))

# Plot lines with specified styles
ax.plot(gold_time, gold_flux,
        linestyle='-', color='tab:gray', linewidth=2,
        label='TMAP8 (Gold)', zorder=2)

ax.plot(claude_time, claude_flux,
        linestyle='-', color='tab:blue', linewidth=2,
        label='TMAP8 (claude-code)', zorder=2)

ax.plot(calib_time, calib_flux,
        linestyle='-', color='tab:green', linewidth=2.5,
        label='TMAP8 (Calibrated)', zorder=3)

ax.plot(exp_time, exp_flux,
        linestyle='--', color='black', linewidth=2,
        label='Experimental Data', zorder=4)

# Set labels and limits
ax.set_xlabel('Time (hr)', fontsize=11, fontweight='bold')
ax.set_ylabel('Deuterium Flux (atom/m²/s)', fontsize=11, fontweight='bold')
ax.set_xlim(left=-0.1, right=20000 / 3600)
ax.set_ylim(bottom=0)

# Configure legend - place in upper right, away from data
ax.legend(loc='upper right', fontsize=10, framealpha=0.95)

# Calculate RMSPE values
def numerical_solution_on_experiment_input(experiment_input, tmap_input, tmap_output):
    """Interpolate numerical solution to experimental time points"""
    new_tmap_output = np.zeros(len(experiment_input))
    for i in range(len(experiment_input)):
        left_limit = np.argwhere((np.diff(tmap_input < experiment_input[i])))[0][0]
        right_limit = left_limit + 1
        new_tmap_output[i] = (experiment_input[i] - tmap_input[left_limit]) / (
            tmap_input[right_limit] - tmap_input[left_limit]
        ) * (tmap_output[right_limit] - tmap_output[left_limit]) + tmap_output[
            left_limit
        ]
    return new_tmap_output

# Calculate RMSPE for Gold (flux values are already in correct units)
gold_flux_interp = numerical_solution_on_experiment_input(
    exp_data['time (s)'].values, gold_data['time'].values,
    gold_data['scaled_recombination_flux_right'].values
)
RMSE_gold = np.sqrt(np.mean((gold_flux_interp - exp_flux) ** 2))
RMSPE_gold = RMSE_gold * 100 / np.mean(exp_flux)

# Calculate RMSPE for claude-code
claude_flux_interp = numerical_solution_on_experiment_input(
    exp_data['time (s)'].values, claude_data['time'].values,
    claude_data['scaled_recombination_flux_right'].values
)
RMSE_claude = np.sqrt(np.mean((claude_flux_interp - exp_flux) ** 2))
RMSPE_claude = RMSE_claude * 100 / np.mean(exp_flux)

# Calculate RMSPE for Calibrated
calib_flux_interp = numerical_solution_on_experiment_input(
    exp_data['time (s)'].values, calib_data['time'].values,
    calib_data['scaled_recombination_flux_right'].values
)
RMSE_calib = np.sqrt(np.mean((calib_flux_interp - exp_flux) ** 2))
RMSPE_calib = RMSE_calib * 100 / np.mean(exp_flux)

# Add RMSPE text annotations in a properly aligned box
# Position in upper-middle area to avoid data
text_x_pos = 0.42  # Relative position (0-1)
text_y_start = 0.65  # Start from this relative y position

# Create formatted RMSPE text
rmspe_text = (f'RMSPE (Gold) = {RMSPE_gold:.2f}%\n'
              f'RMSPE (claude-code) = {RMSPE_claude:.2f}%\n'
              f'RMSPE (Calibrated) = {RMSPE_calib:.2f}%')

ax.text(text_x_pos, text_y_start, rmspe_text,
        transform=ax.transAxes,
        fontsize=9,
        fontweight='bold',
        verticalalignment='top',
        horizontalalignment='left',
        bbox=dict(boxstyle='round,pad=0.5',
                  facecolor='white',
                  alpha=0.9,
                  edgecolor='gray',
                  linewidth=1.5))

# Configure grid
ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)

# Use scientific notation for y-axis
ax.ticklabel_format(axis='y', style='sci', scilimits=(15, 15))

# Minor ticks
ax.minorticks_on()

# Tight layout
plt.tight_layout()

# Save figure
output_file = 'val-2a_comparison_clean.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file}")

# Also save as PDF for vector graphics
output_pdf = 'val-2a_comparison_clean.pdf'
plt.savefig(output_pdf, format='pdf', bbox_inches='tight')
print(f"✓ Saved: {output_pdf}")

plt.close()

print("\nPlot generation complete!")
