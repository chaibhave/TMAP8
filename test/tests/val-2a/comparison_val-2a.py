import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import pandas as pd
from scipy import special
import scipy.stats as stats
import os

# Changes working directory to script directory (for consistent MooseDocs usage)
script_folder = os.path.dirname(__file__)
os.chdir(script_folder)


def numerical_solution_on_experiment_input(experiment_input, tmap_input, tmap_output):
    """Get new numerical solution based on the experimental input data points

    Args:
        experiment_input (float, ndarray): experimental input data points
        tmap_input (float, ndarray): numerical input data points
        tmap_output (float, ndarray): numerical output data points

    Returns:
        float, ndarray: updated tmap_output based on the data points in experiment_input
    """
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


# Read gold simulation data
if "/tmap8/doc/" in script_folder.lower():  # if in documentation folder
    csv_folder = "../../../../test/tests/val-2a/gold/val-2a_out.csv"
else:  # if in test folder
    csv_folder = "./gold/val-2a_out.csv"
simulation_TMAP4_data = pd.read_csv(csv_folder)
simulation_time_TMAP4 = simulation_TMAP4_data["time"]
simulation_recom_flux_left_TMAP4 = simulation_TMAP4_data[
    "scaled_recombination_flux_left"
]
simulation_recom_flux_right_TMAP4 = simulation_TMAP4_data[
    "scaled_recombination_flux_right"
]

# Read new simulation data
if "/tmap8/doc/" in script_folder.lower():  # if in documentation folder
    csv_folder_new = "../../../../test/tests/val-2a/val-2a_out.csv"
else:  # if in test folder
    csv_folder_new = "./val-2a_out.csv"
try:
    simulation_new_data = pd.read_csv(csv_folder_new)
    simulation_time_new = simulation_new_data["time"]
    simulation_recom_flux_left_new = simulation_new_data[
        "scaled_recombination_flux_left"
    ]
    simulation_recom_flux_right_new = simulation_new_data[
        "scaled_recombination_flux_right"
    ]
    has_new_data = True
except FileNotFoundError:
    has_new_data = False
    print("New simulation data not found, plotting gold results only")

# Read calibrated simulation data
if "/tmap8/doc/" in script_folder.lower():  # if in documentation folder
    csv_folder_calib = "../../../../test/tests/val-2a/val-2a_param_reduced_out.csv"
else:  # if in test folder
    csv_folder_calib = "./val-2a_param_reduced_out.csv"
try:
    simulation_calib_data = pd.read_csv(csv_folder_calib)
    simulation_time_calib = simulation_calib_data["time"]
    simulation_recom_flux_left_calib = simulation_calib_data[
        "scaled_recombination_flux_left"
    ]
    simulation_recom_flux_right_calib = simulation_calib_data[
        "scaled_recombination_flux_right"
    ]
    has_calib_data = True
except FileNotFoundError:
    has_calib_data = False
    print("Calibrated simulation data not found")

# Read experiment data
if "/tmap8/doc/" in script_folder.lower():  # if in documentation folder
    csv_folder = "../../../../test/tests/val-2a/gold/experiment_data_paper.csv"
else:  # if in test folder
    csv_folder = "./gold/experiment_data_paper.csv"
experiment_TMAP4_data = pd.read_csv(csv_folder)
experiment_time_TMAP4 = experiment_TMAP4_data["time (s)"]
experiment_flux_TMAP4 = experiment_TMAP4_data["permeation flux (atom/m^2/s)"]

TMAP4_file_base = "val-2a_comparison"
############################ recombination flux - atom/m$^2$/s ############################
fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

ax.plot(
    simulation_time_TMAP4 / 3600,
    simulation_recom_flux_right_TMAP4,
    linestyle="-",
    label=r"TMAP8 (Gold)",
    c="tab:gray",
    linewidth=2,
)

# Plot new simulation results if available
if has_new_data:
    ax.plot(
        simulation_time_new / 3600,
        simulation_recom_flux_right_new,
        linestyle="-",
        label=r"TMAP8 (claude-code)",
        c="tab:blue",
        linewidth=2,
    )

# Plot calibrated simulation results if available
# if has_calib_data:
#     ax.plot(
#         simulation_time_calib / 3600,
#         simulation_recom_flux_right_calib,
#         linestyle="-",
#         label=r"TMAP8 (Calibrated)",
#         c="tab:green",
#         linewidth=2,
#     )

ax.plot(
    experiment_time_TMAP4 / 3600,
    experiment_flux_TMAP4,
    linestyle="--",
    label=r"Experiment",
    c="k",
    linewidth=1.5,
)

ax.set_xlabel("Time (hr)")
ax.set_ylabel("Deuterium flux (atom/m$^2$/s)")
ax.legend(loc="best")
ax.set_ylim(bottom=0)
ax.set_xlim(left=-0.1, right=2e4 / 3600)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)

# Calculate RMSPE for gold results
tmap_flux_for_rmspe = numerical_solution_on_experiment_input(
    experiment_time_TMAP4, simulation_time_TMAP4, simulation_recom_flux_right_TMAP4
)
RMSE = np.sqrt(np.mean((tmap_flux_for_rmspe - experiment_flux_TMAP4) ** 2))
RMSPE = RMSE * 100 / np.mean(experiment_flux_TMAP4)

# Calculate RMSPE for new results if available
if has_new_data:
    tmap_flux_new_for_rmspe = numerical_solution_on_experiment_input(
        experiment_time_TMAP4, simulation_time_new, simulation_recom_flux_right_new
    )
    RMSE_new = np.sqrt(np.mean((tmap_flux_new_for_rmspe - experiment_flux_TMAP4) ** 2))
    RMSPE_new = RMSE_new * 100 / np.mean(experiment_flux_TMAP4)

# Calculate RMSPE for calibrated results if available
if has_calib_data:
    tmap_flux_calib_for_rmspe = numerical_solution_on_experiment_input(
        experiment_time_TMAP4, simulation_time_calib, simulation_recom_flux_right_calib
    )
    RMSE_calib = np.sqrt(np.mean((tmap_flux_calib_for_rmspe - experiment_flux_TMAP4) ** 2))
    RMSPE_calib = RMSE_calib * 100 / np.mean(experiment_flux_TMAP4)

# Display RMSPE values
text_y_start = 4.5e16
text_y_spacing = 1.5e16
text_x = 1e4 / 3600.0

ax.text(text_x, text_y_start, "RMSPE (Gold) = %.2f%%" % RMSPE, fontweight="bold")
if has_new_data:
    ax.text(text_x, text_y_start - text_y_spacing, "RMSPE (claude-code) = %.2f%%" % RMSPE_new, fontweight="bold", color="tab:blue")
# if has_calib_data:
#     text_offset = 2 * text_y_spacing if has_new_data else text_y_spacing
#     ax.text(text_x, text_y_start - text_offset, "RMSPE (Calibrated) = %.2f%%" % RMSPE_calib, fontweight="bold", color="tab:green")
ax.minorticks_on()
ax.ticklabel_format(axis="y", style="sci", scilimits=(15, 15))
plt.savefig(f"{TMAP4_file_base}.png", bbox_inches="tight", dpi=300)
plt.close(fig)


############################ implantation - atom/m$^2$/s ############################
sigma = 2.4e-9  # m
mu = 14e-9  # m
flux = 4.9e19 * 0.75  # atom/m$^2$/s
coordinate_x = np.arange(0, 25e-9, 0.1e-9)

normal_distribution = flux * 1.5 * stats.norm.pdf(coordinate_x, mu, sigma)
piecewise = np.zeros(len(coordinate_x))
piecewise[80:120] = flux * 0.25 / 4e-9
piecewise[120:160] = flux * 1.00 / 4e-9
piecewise[160:200] = flux * 0.25 / 4e-9

fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

ax.plot(
    coordinate_x,
    normal_distribution,
    linestyle="-",
    label=r"normal distribution",
    c="k",
)
ax.plot(coordinate_x, piecewise, linestyle="-", label=r"piecewise", c="gray")

ax.set_xlabel("Depth (m)")
ax.set_ylabel("Source rate (atom/m$^3$/s)")
ax.legend(loc="best")
ax.set_ylim(bottom=0)
ax.set_xlim(left=0, right=3e-8)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
plt.savefig(f"val-2a_comparison_normal_distribution.png", bbox_inches="tight", dpi=300)
plt.close(fig)
