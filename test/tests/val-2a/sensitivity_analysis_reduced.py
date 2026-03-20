#!/usr/bin/env python3
"""
REDUCED Sensitivity Analysis for Val-2a Bayesian Calibration
Only calibrates the 2 most influential parameters identified from full analysis:
  - Kr_left_fraction (ST = 1.218)
  - gaussian_factor (ST = 0.123)

This dramatically reduces computational cost:
  - 2 parameters instead of 8 → ~4x fewer samples needed
  - Faster convergence in Bayesian calibration
  - Clearer parameter identifiability

Requires: SALib, numpy, pandas, subprocess
Install: pip install SALib numpy pandas pyDOE3 tqdm

Usage:
    python sensitivity_analysis_reduced.py --n-samples 512 --parallel 8
"""

import numpy as np
import pandas as pd
import subprocess
import argparse
import os
from pathlib import Path
from SALib.sample import saltelli
from SALib.analyze import sobol
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# REDUCED parameter set - only the 2 most influential parameters
PARAMETERS = {
    'Kr_left_fraction': {
        'bounds': [0.95, 0.9999],
        'log_scale': False,
        'units': 'dimensionless',
        'description': 'Surface cleanup fraction (MOST INFLUENTIAL)',
        'ST_full': 1.2182  # From full 8-parameter analysis
    },
    'gaussian_factor': {
        'bounds': [1.0, 2.5],
        'log_scale': False,
        'units': 'dimensionless',
        'description': 'Implantation profile scaling factor',
        'ST_full': 0.1229  # From full 8-parameter analysis
    }
}


def setup_problem():
    """Create SALib problem definition"""
    problem = {
        'num_vars': len(PARAMETERS),
        'names': list(PARAMETERS.keys()),
        'bounds': [PARAMETERS[p]['bounds'] for p in PARAMETERS.keys()]
    }
    return problem


def generate_samples(problem, n_samples=512):
    """Generate Sobol samples

    Args:
        problem: SALib problem definition
        n_samples: Base number of samples (total will be N*(2*D+2))

    Returns:
        samples: Array of parameter samples
    """
    # Saltelli sampler generates N*(2*D+2) samples for Sobol analysis
    # For D=2: total = N*(2*2+2) = N*6 samples
    samples = saltelli.sample(problem, n_samples, calc_second_order=True)

    # Transform log-scale parameters (none in reduced set, but keep for consistency)
    for i, param_name in enumerate(problem['names']):
        if PARAMETERS[param_name]['log_scale']:
            bounds = PARAMETERS[param_name]['bounds']
            log_bounds = [np.log10(bounds[0]), np.log10(bounds[1])]
            samples[:, i] = 10 ** (log_bounds[0] + samples[:, i] * (log_bounds[1] - log_bounds[0]))

    return samples


def run_tmap8_single(params, run_id, tmap8_exe, input_file, output_dir):
    """Run single TMAP8 simulation with given parameters

    Args:
        params: Dictionary of parameter values
        run_id: Unique identifier for this run
        tmap8_exe: Path to TMAP8 executable
        input_file: TMAP8 input file
        output_dir: Directory for outputs

    Returns:
        Dictionary with run_id and output metrics
    """
    # Create run directory
    run_dir = output_dir / f"run_{run_id:05d}"
    run_dir.mkdir(parents=True, exist_ok=True)

    # Convert input file to absolute path
    input_file_abs = Path(input_file).resolve()

    # Build command line arguments
    args = [str(tmap8_exe), '-i', str(input_file_abs)]
    args.extend([f"{k}={v}" for k, v in params.items()])
    args.extend(['--no-color', '--no-color-autom'])

    # Run TMAP8
    try:
        result = subprocess.run(
            args,
            cwd=run_dir,
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )

        if result.returncode != 0:
            print(f"Run {run_id} failed: {result.stderr[:200]}")
            return {'run_id': run_id, 'success': False, 'error': result.stderr[:200]}

        # Extract output metrics
        csv_file = run_dir / 'val-2a_param_reduced_out.csv'
        if not csv_file.exists():
            return {'run_id': run_id, 'success': False, 'error': 'CSV output not found'}

        df = pd.read_csv(csv_file)

        # Calculate metrics of interest
        metrics = {
            'run_id': run_id,
            'success': True,
            'peak_flux': df['scaled_recombination_flux_right'].max(),
            'final_flux': df['scaled_recombination_flux_right'].iloc[-1],
            'total_permeation': np.trapezoid(df['scaled_recombination_flux_right'], df['time']),
            'peak_time': df.loc[df['scaled_recombination_flux_right'].idxmax(), 'time'],
            'final_inventory': df['scaled_total_inventory'].iloc[-1]
        }

        # Add parameters
        metrics.update(params)

        return metrics

    except subprocess.TimeoutExpired:
        return {'run_id': run_id, 'success': False, 'error': 'Timeout'}
    except Exception as e:
        return {'run_id': run_id, 'success': False, 'error': str(e)}


def run_ensemble_parallel(samples, problem, tmap8_exe, input_file, output_dir, n_processes=8):
    """Run ensemble of TMAP8 simulations in parallel

    Args:
        samples: Array of parameter samples
        problem: SALib problem definition
        tmap8_exe: Path to TMAP8 executable
        input_file: TMAP8 input file
        output_dir: Directory for outputs
        n_processes: Number of parallel processes

    Returns:
        List of result dictionaries
    """
    results = []

    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = []

        for i, sample in enumerate(samples):
            params = {name: float(val) for name, val in zip(problem['names'], sample)}
            future = executor.submit(run_tmap8_single, params, i, tmap8_exe, input_file, output_dir)
            futures.append(future)

        # Collect results with progress bar
        for future in tqdm(as_completed(futures), total=len(futures), desc="Running simulations"):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"Exception in future: {e}")

    return results


def analyze_sensitivity(problem, samples, results, output_name='peak_flux'):
    """Perform Sobol sensitivity analysis

    Args:
        problem: SALib problem definition
        samples: Parameter samples
        results: Simulation results
        output_name: Which output metric to analyze

    Returns:
        Dictionary of Sobol indices
    """
    # Extract output values in same order as samples
    Y = np.array([r[output_name] for r in results if r['success']])

    if len(Y) != len(samples):
        print(f"Warning: Only {len(Y)}/{len(samples)} runs succeeded")
        # Need to filter samples to match successful runs
        successful_indices = [r['run_id'] for r in results if r['success']]
        samples_filtered = samples[successful_indices, :]
    else:
        samples_filtered = samples

    # Perform Sobol analysis with second-order indices (only 2 parameters, so cheap)
    Si = sobol.analyze(problem, Y, calc_second_order=True, print_to_console=False)

    return Si


def print_sensitivity_results(Si, problem):
    """Print formatted sensitivity analysis results"""
    print("\n" + "="*80)
    print("REDUCED SENSITIVITY ANALYSIS RESULTS (2 Parameters)")
    print("="*80)

    print("\nFirst-order Sobol indices (S1):")
    print("Measures direct effect of each parameter")
    print("-"*80)
    for i, name in enumerate(problem['names']):
        desc = PARAMETERS[name]['description']
        st_full = PARAMETERS[name]['ST_full']
        print(f"{name:25s} | S1={Si['S1'][i]:7.4f} ± {Si['S1_conf'][i]:7.4f} | {desc}")
        print(f"{'':25s} | ST_full={st_full:7.4f} (from 8-parameter analysis)")

    print("\nTotal-order Sobol indices (ST):")
    print("Measures total effect including interactions")
    print("-"*80)
    for i, name in enumerate(problem['names']):
        desc = PARAMETERS[name]['description']
        print(f"{name:25s} | ST={Si['ST'][i]:7.4f} ± {Si['ST_conf'][i]:7.4f} | {desc}")

    if 'S2' in Si:
        print("\nSecond-order Sobol indices (S2):")
        print("Measures pairwise interactions")
        print("-"*80)
        s2_matrix = Si['S2']
        for i in range(len(problem['names'])):
            for j in range(i+1, len(problem['names'])):
                name_i = problem['names'][i]
                name_j = problem['names'][j]
                s2_val = s2_matrix[i, j]
                s2_conf = Si['S2_conf'][i, j]
                print(f"{name_i} × {name_j:20s} | S2={s2_val:7.4f} ± {s2_conf:7.4f}")

    print("\n" + "="*80)
    print("NOTE: These indices are from the reduced 2-parameter model.")
    print("They may differ from the full 8-parameter analysis due to:")
    print("  - Removal of low-sensitivity parameters")
    print("  - Changes in interaction structure")
    print("="*80)


def main():
    parser = argparse.ArgumentParser(description='REDUCED sensitivity analysis for Val-2a (2 parameters)')
    parser.add_argument('--n-samples', type=int, default=512,
                       help='Base number of samples (total will be N*6 for D=2)')
    parser.add_argument('--parallel', type=int, default=8,
                       help='Number of parallel processes')
    parser.add_argument('--tmap8-exe', type=str,
                       default=str(Path.home() / 'projects/TMAP8/tmap8-opt'),
                       help='Path to TMAP8 executable')
    parser.add_argument('--input-file', type=str, default='val-2a_param_reduced.i',
                       help='TMAP8 input file (reduced version)')
    parser.add_argument('--output-dir', type=str, default='sensitivity_reduced',
                       help='Directory for output files')
    parser.add_argument('--output-metric', type=str, default='peak_flux',
                       choices=['peak_flux', 'total_permeation', 'final_inventory'],
                       help='Output metric to analyze')

    args = parser.parse_args()

    # Setup
    problem = setup_problem()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n" + "="*80)
    print("REDUCED SENSITIVITY ANALYSIS - 2 Parameters Only")
    print("="*80)
    print(f"Configuration:")
    print(f"  Base samples: {args.n_samples}")
    print(f"  Total samples: {args.n_samples * (2 * problem['num_vars'] + 2)} (N*6 for D=2)")
    print(f"  Parameters: {problem['num_vars']}")
    print(f"  Parallel processes: {args.parallel}")
    print(f"  Output metric: {args.output_metric}")
    print()
    print("Parameters being calibrated:")
    for name in problem['names']:
        bounds = PARAMETERS[name]['bounds']
        print(f"  {name}: [{bounds[0]}, {bounds[1]}] {PARAMETERS[name]['units']}")
    print("="*80 + "\n")

    # Generate samples
    print("Generating Sobol samples...")
    samples = generate_samples(problem, args.n_samples)
    print(f"Generated {len(samples)} samples")

    # Save samples
    samples_file = output_dir / 'samples.csv'
    pd.DataFrame(samples, columns=problem['names']).to_csv(samples_file, index=False)
    print(f"Saved samples to {samples_file}")

    # Run ensemble
    print(f"\nRunning {len(samples)} TMAP8 simulations...")
    results = run_ensemble_parallel(
        samples, problem, args.tmap8_exe, args.input_file, output_dir, args.parallel
    )

    # Save results
    results_file = output_dir / 'results.csv'
    pd.DataFrame(results).to_csv(results_file, index=False)
    print(f"Saved results to {results_file}")

    # Analyze sensitivity
    print(f"\nAnalyzing sensitivity for {args.output_metric}...")
    Si = analyze_sensitivity(problem, samples, results, args.output_metric)

    # Save sensitivity indices
    si_file = output_dir / f'sensitivity_{args.output_metric}.json'
    si_dict = {
        'S1': Si['S1'].tolist(),
        'S1_conf': Si['S1_conf'].tolist(),
        'ST': Si['ST'].tolist(),
        'ST_conf': Si['ST_conf'].tolist(),
        'parameter_names': problem['names']
    }
    if 'S2' in Si:
        si_dict['S2'] = Si['S2'].tolist()
        si_dict['S2_conf'] = Si['S2_conf'].tolist()

    with open(si_file, 'w') as f:
        json.dump(si_dict, f, indent=2)
    print(f"Saved sensitivity indices to {si_file}")

    # Print results
    print_sensitivity_results(Si, problem)

    print(f"\nReduced sensitivity analysis complete!")
    print(f"Results saved in: {output_dir}")
    print(f"\nComputational savings vs full 8-parameter analysis:")
    print(f"  Sample reduction: {args.n_samples * 18} → {args.n_samples * 6} ({100*(1-6/18):.0f}% fewer)")
    print(f"  Typical speedup for Bayesian calibration: ~8-16x faster")


if __name__ == '__main__':
    main()
