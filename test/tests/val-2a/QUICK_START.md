# Quick Start: Reduced Calibration (2 Parameters)

⭐ **This is the recommended starting point for calibrating the Val-2a model.**

## TL;DR

Run these 4 commands to calibrate the model in ~4-8 hours:

```bash
# 1. Generate training data (~1 hour)
python bayesian_calibration_reduced.py --mode generate-ensemble --n-samples 200 --parallel 8

# 2. Build surrogate model (~1 minute)
python bayesian_calibration_reduced.py --mode build-surrogate

# 3. Run calibration (~2-4 hours)
python bayesian_calibration_reduced.py --mode calibrate --walkers 32 --steps 5000

# 4. View results
python bayesian_calibration_reduced.py --mode analyze
```

## What Gets Calibrated?

Only the **2 most influential parameters** (from sensitivity analysis):

1. **`Kr_left_fraction`** = Surface cleanup fraction (ST = 1.218)
   - Prior: [0.95, 0.9999]
   - Nominal: 0.99

2. **`gaussian_factor`** = Implantation profile scaling (ST = 0.123)
   - Prior: [1.0, 2.5]
   - Nominal: 1.5

**6 other parameters are fixed** (negligible sensitivity, ST < 0.01)

## Prerequisites

Install required packages:

```bash
pip install emcee corner numpy pandas scipy scikit-learn pyDOE3 matplotlib seaborn tqdm
```

Verify TMAP8 executable path:
```bash
~/projects/TMAP8/tmap8-opt --version
```

## Detailed Steps

### Step 1: Generate Ensemble (~1 hour)

```bash
python bayesian_calibration_reduced.py \
    --mode generate-ensemble \
    --n-samples 200 \
    --parallel 8 \
    --input-file val-2a_param_reduced.i \
    --ensemble-dir calibration_reduced_ensemble
```

What happens:
- Generates 200 Latin Hypercube samples in parameter space
- Runs 200 TMAP8 simulations in parallel (8 cores)
- Saves flux time series for surrogate training
- Creates: `calibration_reduced_ensemble/ensemble_results.pkl`

Expected runtime: ~1 hour (8 cores, 200 × ~18 sec/run)

### Step 2: Build Surrogate Model (~1 minute)

```bash
python bayesian_calibration_reduced.py \
    --mode build-surrogate \
    --ensemble-dir calibration_reduced_ensemble \
    --output-dir calibration_reduced_results
```

What happens:
- Trains Gaussian Process on ensemble data (80/20 train/test split)
- Validates surrogate accuracy (should see R² > 0.95)
- Saves: `calibration_reduced_results/surrogate_model.pkl`

Expected output:
```
Training GP surrogate model...
  Training samples: 160
  Time points: 1234 -> 246 (subsampled)
  Training GP...
  Training complete!

Surrogate validation metrics:
  RMSE: 3.45e+15
  MAE:  2.12e+15
  R²:   0.967
```

✅ **Good:** R² > 0.95
⚠️ **Warning:** R² < 0.90 → Increase --n-samples to 300

### Step 3: Run Bayesian Calibration (~2-4 hours)

```bash
python bayesian_calibration_reduced.py \
    --mode calibrate \
    --output-dir calibration_reduced_results \
    --walkers 32 \
    --steps 5000 \
    --burn-in 1000 \
    --noise-std 1e16
```

What happens:
- Runs MCMC with 32 walkers for 5000 steps each
- Discards first 1000 steps as burn-in
- Samples posterior distribution for both parameters
- Creates:
  - `mcmc_sampler.pkl`
  - `mcmc_samples.npy` (128,000 posterior samples after burn-in)
  - `calibration_results.json`
  - `corner_plot.png`
  - `trace_plots.png`

Expected runtime: 2-4 hours (depends on surrogate speed)

Monitor progress: Shows progress bar with ETA

### Step 4: Analyze Results

```bash
python bayesian_calibration_reduced.py \
    --mode analyze \
    --output-dir calibration_reduced_results
```

Example output:
```
================================================================================
CALIBRATION RESULTS (Reduced 2-Parameter Model)
================================================================================

Burn-in: 1000 steps
Post-burn-in samples: 128000

Kr_left_fraction:
  Mean:   0.9850 ± 0.0032
  Median: 0.9851
  95% CI: [0.9795, 0.9904]
  Nominal: 0.9900

gaussian_factor:
  Mean:   1.7234 ± 0.2145
  Median: 1.7189
  95% CI: [1.3654, 2.0821]
  Nominal: 1.5000

Mean acceptance fraction: 0.347
================================================================================
```

## Interpreting Results

### Parameter Estimates

Look at the **Mean ± Std** values:
- This is your calibrated parameter estimate with uncertainty
- Compare to **Nominal** to see how much it shifted

### Credible Intervals (95% CI)

- If CI is narrow: Parameter is well-constrained by data ✅
- If CI is wide: Parameter is poorly identifiable ⚠️
- If CI touches prior bounds: Prior may be too restrictive 🔴

### Acceptance Fraction

- **0.2 - 0.5:** Good mixing ✅
- **< 0.1:** Steps too large, poor exploration ⚠️
- **> 0.7:** Steps too small, slow convergence ⚠️

### Plots

**corner_plot.png:**
- Diagonal: Marginal posterior distributions
- Off-diagonal: Joint posterior (shows correlation)
- Look for: Gaussian-like posteriors, weak correlation

**trace_plots.png:**
- Should look like "fuzzy caterpillars"
- No trends, no getting stuck
- If trending: Increase burn-in or steps

## Next Steps

### 1. Validate Calibration

Run TMAP8 with calibrated parameters:

```bash
# Extract calibrated mean values from results
Kr_frac=$(python -c "import json; print(json.load(open('calibration_reduced_results/calibration_results.json'))['Kr_left_fraction']['mean'])")
gauss=$(python -c "import json; print(json.load(open('calibration_reduced_results/calibration_results.json'))['gaussian_factor']['mean'])")

# Run with calibrated values
~/projects/TMAP8/tmap8-opt -i val-2a_param_reduced.i \
    Kr_left_fraction=$Kr_frac \
    gaussian_factor=$gauss
```

Compare output to experimental data:
- Should see improved fit vs baseline (RMSPE should decrease)
- Target: RMSPE < 20%

### 2. Uncertainty Quantification

Sample from posterior to propagate uncertainty:

```python
import numpy as np
samples = np.load('calibration_reduced_results/mcmc_samples.npy')

# Draw 100 random posterior samples
posterior_samples = samples[np.random.choice(len(samples), 100)]

# Run TMAP8 with each to get prediction uncertainty
for i, theta in enumerate(posterior_samples):
    # Run simulation with theta[0]=Kr_left_fraction, theta[1]=gaussian_factor
    ...
```

### 3. Document Results

Record calibrated values in your model documentation:
- Mean ± std for both parameters
- 95% credible intervals
- RMSPE before/after calibration
- Date and sample size used

## Troubleshooting

### "Surrogate validation R² = 0.85"

Surrogate not accurate enough. Solutions:
1. Increase ensemble size: `--n-samples 300`
2. Reduce time subsampling in code (change `subsample_time=5` to `3`)

### "MCMC not converging (trace plots trending)"

Not enough burn-in or steps. Solutions:
1. Increase burn-in: `--burn-in 2000`
2. Increase steps: `--steps 10000`
3. Increase walkers: `--walkers 64`

### "Posterior same as prior (no information gained)"

Possible causes:
1. Surrogate inaccurate → Fix surrogate first
2. Data not informative → Check experimental data quality
3. Noise too high → Reduce `--noise-std`

### "Parameter hits prior bound"

Prior may be too restrictive. Solutions:
1. Check if physical justification supports wider range
2. Edit bounds in `bayesian_calibration_reduced.py`
3. Re-run from Step 1 with new priors

## Comparing to Baseline

| Metric | Baseline (Nominal) | After Calibration (Expected) |
|--------|-------------------|----------------------------|
| RMSPE | 25.65% | 15-20% |
| Peak flux accuracy | ~30% error | ~10% error |
| Runtime | 1 simulation | +4-8 hours setup |
| Parameter uncertainty | Unknown | Quantified |

## When to Use Full Calibration Instead

Consider full 8-parameter calibration if:
- ❌ Reduced calibration gives RMSPE > 20%
- ❌ Fixed parameter values are highly uncertain
- ✅ You have HPC cluster access
- ✅ You need complete uncertainty characterization

See `CALIBRATION_SUMMARY.md` for full vs reduced comparison.

## Time Budget

| Task | Time | Can skip? |
|------|------|-----------|
| Ensemble generation | 1 hour | No |
| Build surrogate | 1 min | No |
| MCMC sampling | 2-4 hours | No |
| Analysis | 1 min | No |
| Validation | 5 min | Recommended |
| **Total** | **4-8 hours** | |

## Files Created

After completion, you'll have:

```
calibration_reduced_ensemble/
├── lhs_samples.csv                    # LHS parameter samples
├── ensemble_results.pkl               # Flux time series
└── run_*/                             # Individual TMAP8 runs

calibration_reduced_results/
├── surrogate_model.pkl                # Trained GP model
├── mcmc_sampler.pkl                   # Full MCMC sampler state
├── mcmc_samples.npy                   # Posterior samples (128000 × 2)
├── calibration_results.json           # Summary statistics
├── corner_plot.png                    # Posterior visualization
└── trace_plots.png                    # Convergence diagnostics
```

## Support

- Detailed guide: `REDUCED_CALIBRATION_README.md`
- Full vs reduced comparison: `CALIBRATION_SUMMARY.md`
- Prior justifications: `PRIOR_DISTRIBUTIONS_JUSTIFICATION.md`

---

**Ready to start?** Run the 4 commands at the top of this file! 🚀
