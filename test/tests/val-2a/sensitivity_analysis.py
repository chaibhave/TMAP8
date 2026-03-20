#!/usr/bin/env python3
"""
Sensitivity Analysis for Val-2a Bayesian Calibration
Uses Sobol indices to identify most influential parameters

Requires: SALib, numpy, pandas, subprocess
Install: pip install SALib numpy pandas

Usage:
    python sensitivity_analysis.py --n-samples 1024 --parallel 8
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

# Parameter definitions with prior ranges
PARAMETERS = {
    'diffusivity': {
        'bounds': [1.0e-10, 1.0e-9],
        'log_scale': True,
        'units': 'm^2/s',
        'description': 'Deuterium diffusivity in PCA'
    },
    'Kr_left_max': {
        'bounds': [1.0e-28, 1.0e-26],
        'log_scale': True,
        'units': 'm^4/atom/s',
        'description': 'Maximum upstream recombination coefficient'
    },
    'Kr_right': {
        'bounds': [1.0e-32, 1.0e-30],
        'log_scale': True,
        'units': 'm^4/atom/s',
        'description': 'Downstream recombination coefficient'
    },
    'gaussian_factor': {
        'bounds': [1.0, 2.5],
        'log_scale': False,
        'units': 'dimensionless',
        'description': 'Implantation profile scaling factor'
    },
    'Kr_left_time_constant': {
        'bounds': [1.0e-5, 2.0e-4],
        'log_scale': True,
        'units': '1/s',
        'description': 'Surface cleanup time constant'
    },
    'Kr_left_fraction': {
        'bounds': [0.95, 0.9999],
        'log_scale': False,
        'units': 'dimensionless',
        'description': 'Surface cleanup fraction'
    },
    'implantation_depth': {
        'bounds': [10e-9, 18e-9],
        'log_scale': False,
        'units': 'm',
        'description': 'Average implantation depth'
    },
    'implantation_sigma': {
        'bounds': [1.5e-9, 3.5e-9],
        'log_scale': False,
        'units': 'm',
        'description': 'Implantation depth straggling'
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


def generate_samples(problem, n_samples=1024):
    """Generate Sobol samples

    Args:
        problem: SALib problem definition
        n_samples: Base number of samples (total will be N*(2*D+2))

    Returns:
        samples: Array of parameter samples
    """
    # Saltelli sampler generates N*(2*D+2) samples for Sobol analysis
    samples = saltelli.sample(problem, n_samples, calc_second_order=False)

    # Transform log-scale parameters
    for i, param_name in enumerate(problem['names']):
        if PARAMETERS[param_name]['log_scale']:
            # Convert from linear to log scale
            bounds = PARAMETERS[param_name]['bounds']
            log_bounds = [np.log10(bounds[0]), np.log10(bounds[1])]
            # Samples are in [0,1], scale to log bounds
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
        csv_file = run_dir / 'val-2a_param_out.csv'
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

    # Perform Sobol analysis
    Si = sobol.analyze(problem, Y, calc_second_order=False, print_to_console=False)

    return Si


def print_sensitivity_results(Si, problem):
    """Print formatted sensitivity analysis results"""
    print("\n" + "="*80)
    print("SENSITIVITY ANALYSIS RESULTS")
    print("="*80)

    print("\nFirst-order Sobol indices (S1):")
    print("Measures direct effect of each parameter")
    print("-"*80)
    for i, name in enumerate(problem['names']):
        desc = PARAMETERS[name]['description']
        print(f"{name:25s} | S1={Si['S1'][i]:7.4f} ± {Si['S1_conf'][i]:7.4f} | {desc}")

    print("\nTotal-order Sobol indices (ST):")
    print("Measures total effect including interactions")
    print("-"*80)
    for i, name in enumerate(problem['names']):
        desc = PARAMETERS[name]['description']
        print(f"{name:25s} | ST={Si['ST'][i]:7.4f} ± {Si['ST_conf'][i]:7.4f} | {desc}")

    # Identify most important parameters
    threshold = 0.05
    important_params = [problem['names'][i] for i, st in enumerate(Si['ST']) if st > threshold]

    print(f"\nMost influential parameters (ST > {threshold}):")
    for param in important_params:
        i = problem['names'].index(param)
        print(f"  - {param}: ST={Si['ST'][i]:.4f}")

    print("\n" + "="*80)


def main():
    parser = argparse.ArgumentParser(description='Sensitivity analysis for Val-2a calibration')
    parser.add_argument('--n-samples', type=int, default=512,
                       help='Base number of samples (total will be N*(2*D+2))')
    parser.add_argument('--parallel', type=int, default=8,
                       help='Number of parallel processes')
    parser.add_argument('--tmap8-exe', type=str,
                       default=str(Path.home() / 'projects/TMAP8/tmap8-opt'),
                       help='Path to TMAP8 executable')
    parser.add_argument('--input-file', type=str, default='val-2a_param.i',
                       help='TMAP8 input file')
    parser.add_argument('--output-dir', type=str, default='sensitivity_runs',
                       help='Directory for output files')
    parser.add_argument('--output-metric', type=str, default='peak_flux',
                       choices=['peak_flux', 'total_permeation', 'final_inventory'],
                       help='Output metric to analyze')

    args = parser.parse_args()

    # Setup
    problem = setup_problem()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Sensitivity Analysis Configuration:")
    print(f"  Base samples: {args.n_samples}")
    print(f"  Total samples: {args.n_samples * (2 * problem['num_vars'] + 2)}")
    print(f"  Parameters: {problem['num_vars']}")
    print(f"  Parallel processes: {args.parallel}")
    print(f"  Output metric: {args.output_metric}")
    print()

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
    with open(si_file, 'w') as f:
        json.dump({
            'S1': Si['S1'].tolist(),
            'S1_conf': Si['S1_conf'].tolist(),
            'ST': Si['ST'].tolist(),
            'ST_conf': Si['ST_conf'].tolist(),
            'parameter_names': problem['names']
        }, f, indent=2)
    print(f"Saved sensitivity indices to {si_file}")

    # Print results
    print_sensitivity_results(Si, problem)

    print(f"\nSensitivity analysis complete!")
    print(f"Results saved in: {output_dir}")


if __name__ == '__main__':
    main()
