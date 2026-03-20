# Bayesian Calibration Plan for Val-2a Model

## Overview
This document outlines a plan for performing Bayesian calibration and uncertainty quantification (UQ) on the val-2a ion implantation validation case. The goal is to calibrate uncertain model parameters against experimental permeation flux data while treating experimental conditions (beam schedule, sample geometry, temperature) as fixed.

## 1. Parameters for Calibration

### 1.1 Primary Parameters (High Uncertainty)
These parameters have the largest uncertainty and strongest influence on permeation flux:

| Parameter | Current Value | Physical Range | Prior Distribution | Rationale |
|-----------|---------------|----------------|-------------------|-----------|
| `diffusivity` | 3.0×10⁻¹⁰ m²/s | 1×10⁻¹⁰ to 1×10⁻⁹ m²/s | Log-uniform | Diffusivity in PCA steel varies with microstructure, grain boundaries |
| `Kr_left_max` | 1.0×10⁻²⁷ m⁴/atom/s | 1×10⁻²⁸ to 1×10⁻²⁶ m⁴/atom/s | Log-uniform | Surface recombination highly sensitive to surface state |
| `Kr_right` | 2.0×10⁻³¹ m⁴/atom/s | 1×10⁻³² to 1×10⁻³⁰ m⁴/atom/s | Log-uniform | Downstream surface recombination |
| `gaussian_factor` | 1.5 | 1.0 to 2.5 | Uniform | Empirical scaling factor for implantation profile |

### 1.2 Secondary Parameters (Moderate Uncertainty)
These parameters are somewhat constrained by physics but still uncertain:

| Parameter | Current Value | Physical Range | Prior Distribution | Rationale |
|-----------|---------------|----------------|-------------------|-----------|
| `Kr_left_time_constant` | 6.0×10⁻⁵ s⁻¹ | 1×10⁻⁵ to 2×10⁻⁴ s⁻¹ | Log-uniform | Surface cleanup rate - depends on sputtering dynamics |
| `Kr_left_fraction` | 0.9999 | 0.95 to 0.9999 | Beta(α=50, β=2) | Fraction of surface cleanup - near 1 but uncertain |
| `implantation_depth` | 14×10⁻⁹ m | 10×10⁻⁹ to 18×10⁻⁹ m | Normal(μ=14, σ=2) nm | Ion implantation depth from SRIM calculations |
| `implantation_sigma` | 2.4×10⁻⁹ m | 1.5×10⁻⁹ to 3.5×10⁻⁹ m | Normal(μ=2.4, σ=0.5) nm | Width of Gaussian implantation profile |

### 1.3 Fixed Parameters
These are well-known from experimental conditions:

- `thickness` = 5×10⁻⁴ m (measured)
- `implantation_flux` = 3.675×10¹⁹ atom/m²/s (measured, 75% retention)
- `end_time` = 20,000 s (experimental duration)
- Beam schedule times: beam_on_1_end, beam_off_1_end, etc. (recorded)

## 2. Observable Data

### 2.1 Primary Observable
- **Permeation flux** (`scaled_recombination_flux_right`): Time series of deuterium flux at downstream surface
- **Experimental data**: 271 data points from `gold/experiment_data_paper.csv`
- **Key features**:
  - Peak flux during each beam-on period (~3×10¹⁷ atoms/m²/s)
  - Decay rates during beam-off periods
  - Overall time evolution

### 2.2 Secondary Observables (for validation)
- **Total inventory** (`scaled_total_inventory`): Total deuterium in sample
- **Upstream flux** (`scaled_recombination_flux_left`): Re-emission from implantation surface

## 3. Bayesian Calibration Approach

### 3.1 Mathematical Framework

**Forward model:**
```
y = f(θ, x) + ε
```
where:
- `y` = observed permeation flux (experimental data)
- `f(θ, x)` = TMAP8 simulation output
- `θ` = parameter vector to calibrate
- `x` = experimental conditions (fixed)
- `ε` = measurement error ~ N(0, σ²)

**Bayes' theorem:**
```
p(θ|y) ∝ p(y|θ) × p(θ)
```
where:
- `p(θ|y)` = posterior distribution (what we want)
- `p(y|θ)` = likelihood (how well model matches data)
- `p(θ)` = prior distribution (initial beliefs about parameters)

**Likelihood function:**
```
p(y|θ) = ∏ᵢ N(yᵢ | f(θ, xᵢ), σ²)
```
Assumes Gaussian measurement error with unknown variance σ².

### 3.2 Sampling Strategy

**MCMC (Markov Chain Monte Carlo):**
- Algorithm: Adaptive Metropolis or DREAM(ZS)
- Chain length: 50,000 - 100,000 samples
- Burn-in: 20% of samples
- Thinning: Every 10th sample to reduce autocorrelation
- Parallel chains: 4-8 chains for convergence diagnostics

**Surrogate modeling (for efficiency):**
- Build Gaussian Process (GP) surrogate of TMAP8 model
- Training samples: 200-500 Latin Hypercube samples
- Use surrogate for rapid likelihood evaluations
- Periodically validate surrogate accuracy

## 4. Implementation in MOOSE/TMAP8

### 4.1 MOOSE Stochastic Tools Module

MOOSE provides several approaches for Bayesian calibration:

#### Option 1: Parameter Study + External Analysis
1. Use `Sampler` objects to generate parameter samples
2. Run ensemble of TMAP8 simulations
3. Export results to CSV
4. Perform Bayesian analysis in Python (PyMC3, emcee, or Stan)

**Pros:** Most flexible, access to mature Bayesian libraries
**Cons:** Requires external tooling, less integrated

#### Option 2: MOOSE-integrated Surrogate + Optimization
1. Use `Trainers` to build surrogate models (Polynomial Chaos, Gaussian Process)
2. Use `AdaptiveImportanceSampler` or `PMCMCSSampler` for Bayesian inference
3. Analyze posteriors with built-in statistics

**Pros:** Fully integrated, efficient with surrogates
**Cons:** Less mature than external libraries, limited MCMC diagnostics

#### Recommended: Hybrid Approach
1. Use MOOSE for parameter sampling and surrogate building
2. Export surrogate or run ensemble in MOOSE
3. Use Python for MCMC sampling and posterior analysis
4. Validate with select full TMAP8 runs

### 4.2 Required Input File Modifications

Create new input file: `val-2a_calibration.i`

```moose
# Add parameter definitions that can be controlled by sampler
diffusivity_param = 3.0e-10  # Will be overridden by sampler

[Distributions]
  [diffusivity_dist]
    type = Uniform
    lower_bound = 1.0e-10
    upper_bound = 1.0e-9
  []

  [Kr_left_max_dist]
    type = Uniform
    lower_bound = 1.0e-28
    upper_bound = 1.0e-26
  []

  # ... (other parameter distributions)
[]

[Samplers]
  [parameter_sampler]
    type = LatinHypercube
    num_rows = 500  # Number of samples for surrogate training
    distributions = 'diffusivity_dist Kr_left_max_dist ...'
    execute_on = INITIAL
  []
[]

[MultiApps]
  [runner]
    type = SamplerFullSolveMultiApp
    sampler = parameter_sampler
    input_files = 'val-2a.i'
  []
[]

[Transfers]
  [parameters_to_sub]
    type = SamplerParameterTransfer
    to_multi_app = runner
    sampler = parameter_sampler
    parameters = 'Materials/diffusivity_material/prop_values ...'
  []

  [results_from_sub]
    type = SamplerPostprocessorTransfer
    from_multi_app = runner
    sampler = parameter_sampler
    to_vector_postprocessor = storage
    from_postprocessor = scaled_recombination_flux_right
  []
[]

[Trainers]
  [gp_trainer]
    type = GaussianProcessTrainer
    sampler = parameter_sampler
    results_vpp = storage
    results_vector = results:scaled_recombination_flux_right
  []
[]

[Outputs]
  csv = true
  execute_on = FINAL
[]
```

### 4.3 Python Analysis Scripts

Create `bayesian_analysis.py`:

```python
import numpy as np
import pandas as pd
import pymc as pm
import arviz as az

# 1. Load MOOSE ensemble results or surrogate
# 2. Define likelihood function
# 3. Set up PyMC model with priors
# 4. Run MCMC sampling
# 5. Analyze posteriors (trace plots, corner plots, convergence)
# 6. Validate calibrated parameters
```

## 5. Implementation Steps

### Phase 1: Setup and Initial Exploration (Week 1)
1. **Create sampling infrastructure**
   - Modify val-2a.i to accept parameter overrides
   - Create val-2a_calibration.i for parameter studies
   - Test with small parameter sweep (10-20 samples)

2. **Sensitivity analysis**
   - Run Sobol sensitivity analysis to identify most influential parameters
   - Use `SobolSampler` with ~1000 samples
   - Focus calibration on high-sensitivity parameters

3. **Define priors**
   - Review literature for parameter ranges
   - Consult with domain experts
   - Document prior choices with justification

### Phase 2: Surrogate Model Development (Week 2)
1. **Generate training data**
   - Latin Hypercube sampling: 500 samples across parameter space
   - Run full TMAP8 ensemble (parallelize with mpirun)
   - Extract time series data at experimental measurement times

2. **Build surrogate model**
   - Train Gaussian Process or Polynomial Chaos Expansion
   - Validate surrogate accuracy (R² > 0.95, RMSE < 5% of data range)
   - Identify regions needing more samples

3. **Refine surrogate**
   - Add samples in high-error regions
   - Re-train and validate
   - Document surrogate accuracy

### Phase 3: Bayesian Calibration (Week 3)
1. **Set up MCMC sampling**
   - Implement likelihood function (compare model to experimental data)
   - Estimate measurement error variance (or treat as unknown parameter)
   - Configure MCMC: 4 chains × 25,000 samples each

2. **Run calibration**
   - Start MCMC sampling
   - Monitor convergence (Gelman-Rubin statistic, effective sample size)
   - Adjust if needed (longer chains, different initial values)

3. **Validate with full model**
   - Sample 50 parameter sets from posterior
   - Run full TMAP8 simulations (not surrogate)
   - Verify predictions match experimental data

### Phase 4: Analysis and Validation (Week 4)
1. **Posterior analysis**
   - Plot marginal distributions for each parameter
   - Create corner plots showing correlations
   - Calculate credible intervals (50%, 90%, 95%)
   - Identify parameter identifiability issues

2. **Model predictions**
   - Generate posterior predictive distributions
   - Create prediction bands for permeation flux
   - Compare uncertainty pre/post calibration

3. **Validation**
   - Cross-validation: hold out subset of data, calibrate on rest
   - Physics checks: ensure calibrated parameters are physically reasonable
   - Sensitivity: verify calibration reduces parameter uncertainty

4. **Documentation**
   - Write calibration report with results
   - Document final parameter distributions
   - Provide recommendations for future experiments

## 6. Computational Requirements

### 6.1 Resources Needed

**For surrogate training (500 samples):**
- CPU time: 500 runs × 55 sec = 7.6 hours (serial)
- With 8 cores: ~1 hour wall time
- Storage: ~500 MB (CSV outputs)

**For MCMC (100,000 surrogate evaluations):**
- CPU time: ~1-2 hours (surrogate is fast)
- Storage: ~100 MB (MCMC chains)

**For validation (50 full runs):**
- CPU time: 50 runs × 55 sec = 45 minutes
- Storage: ~50 MB

**Total:** ~4-5 hours wall time, ~1 GB storage

### 6.2 Parallelization Strategy
- Use `MultiApps` with `SamplerFullSolveMultiApp` for ensemble runs
- Run with `mpirun -n 8` for training data generation
- Use parallel MCMC chains (4-8 chains)

## 7. Expected Outcomes

### 7.1 Calibrated Parameters
- **Posterior distributions** for all calibrated parameters
- **Credible intervals** (90% and 95%)
- **Parameter correlations** revealing physical dependencies
- **Reduced uncertainty** compared to priors (quantified)

### 7.2 Model Improvements
- **Better fit to data**: RMSPE reduced from 25.65% to target <15%
- **Quantified uncertainty**: Prediction bands around model output
- **Identifiability**: Which parameters can/cannot be constrained by data

### 7.3 Physical Insights
- **Dominant mechanisms**: Which parameters control peak flux vs. decay?
- **Surface vs. bulk**: Relative importance of surface recombination vs. diffusion
- **Time-dependent effects**: Evolution of surface state during implantation

### 7.4 Recommendations
- **Parameter refinement**: Updated values with uncertainty for future simulations
- **Experimental design**: Which measurements would most reduce uncertainty?
- **Model structure**: Is current model adequate or are new physics needed?

## 8. Potential Challenges and Mitigation

### 8.1 Computational Cost
**Challenge:** Full TMAP8 runs take ~1 minute, may need thousands for MCMC

**Mitigation:**
- Build accurate surrogate model (GP or PCE)
- Use adaptive sampling to focus on high-likelihood regions
- Parallelize ensemble runs across multiple cores
- Consider reduced-order models for screening

### 8.2 Parameter Identifiability
**Challenge:** Multiple parameter combinations may fit data equally well

**Mitigation:**
- Start with sensitivity analysis to identify important parameters
- Calibrate in stages: fix some parameters, calibrate others
- Use informative priors from literature/physics
- Analyze posterior correlations to identify trade-offs

### 8.3 Model Inadequacy
**Challenge:** Model structure may not capture all physics

**Mitigation:**
- Include model discrepancy term in likelihood
- Validate residuals are random (no systematic errors)
- If needed, refine model physics before calibration
- Use cross-validation to detect overfitting

### 8.4 Measurement Uncertainty
**Challenge:** Experimental data has unknown measurement error

**Mitigation:**
- Treat measurement variance as unknown parameter
- Use replicate measurements if available
- Consider heteroscedastic errors (varying with flux magnitude)
- Validate with independent data if possible

## 9. Software and Tools

### 9.1 Required Software
- **TMAP8** (current installation)
- **MOOSE** with Stochastic Tools module
- **Python 3.8+** with:
  - PyMC (Bayesian inference)
  - ArviZ (MCMC diagnostics)
  - scikit-learn (surrogate models)
  - matplotlib, seaborn (plotting)
  - pandas, numpy (data handling)

### 9.2 Optional Tools
- **emcee**: Alternative MCMC sampler in Python
- **Stan**: More sophisticated MCMC with gradient information
- **UQPy**: UQ library with various sampling methods
- **ChaosPy**: Polynomial Chaos Expansion tools

## 10. References and Resources

### 10.1 MOOSE Documentation
- Stochastic Tools Module: https://mooseframework.inl.gov/modules/stochastic_tools/
- Parameter Studies: https://mooseframework.inl.gov/modules/stochastic_tools/examples/parameter_study.html
- Surrogate Modeling: https://mooseframework.inl.gov/modules/stochastic_tools/surrogates/

### 10.2 Bayesian Calibration Literature
- Kennedy and O'Hagan (2001): "Bayesian calibration of computer models"
- Higdon et al. (2004): "Computer model calibration using high-dimensional output"
- Gelman et al. (2013): "Bayesian Data Analysis" (textbook)

### 10.3 UQ in Nuclear/Materials Science
- Permann et al. (2020): "MOOSE: Enabling massively parallel multiphysics simulation"
- NEAMS uncertainty quantification guidelines
- ASME V&V 20-2009: Standard for verification and validation

## 11. Next Steps

1. **Immediate (Day 1-3)**:
   - Review this plan with team/advisor
   - Set up Python environment (PyMC, ArviZ, etc.)
   - Test parameter override in val-2a.i

2. **Short-term (Week 1)**:
   - Run sensitivity analysis (Sobol indices)
   - Create calibration input file
   - Generate initial training dataset (100 samples)

3. **Medium-term (Week 2-3)**:
   - Build and validate surrogate model
   - Implement Bayesian calibration
   - Run MCMC and analyze posteriors

4. **Long-term (Week 4+)**:
   - Validate calibrated model
   - Write calibration report
   - Apply to other validation cases (val-2b, val-2c)

---

## Appendix A: Example Python Workflow

```python
import pymc as pm
import numpy as np
import pandas as pd

# Load experimental data
exp_data = pd.read_csv('gold/experiment_data_paper.csv')
exp_time = exp_data['time (s)'].values
exp_flux = exp_data['permeation flux (atom/m^2/s)'].values
exp_std = 0.1 * exp_flux  # Assume 10% measurement error

# Load surrogate model (from MOOSE)
surrogate = load_surrogate('gp_model.pkl')

# Define Bayesian model
with pm.Model() as model:
    # Priors
    diffusivity = pm.Uniform('diffusivity', 1e-10, 1e-9)
    Kr_left_max = pm.Uniform('Kr_left_max', 1e-28, 1e-26)
    Kr_right = pm.Uniform('Kr_right', 1e-32, 1e-30)
    gaussian_factor = pm.Uniform('gaussian_factor', 1.0, 2.5)

    # Combine parameters
    theta = pm.math.stack([diffusivity, Kr_left_max, Kr_right, gaussian_factor])

    # Likelihood (using surrogate)
    model_flux = surrogate.predict(theta, exp_time)
    likelihood = pm.Normal('obs', mu=model_flux, sigma=exp_std, observed=exp_flux)

    # Sample posterior
    trace = pm.sample(10000, tune=2000, chains=4, target_accept=0.9)

# Analyze results
summary = az.summary(trace)
az.plot_trace(trace)
az.plot_posterior(trace)
```

## Appendix B: Expected Timeline

| Phase | Duration | Key Deliverables |
|-------|----------|------------------|
| Setup | 3 days | Input files, sensitivity analysis |
| Surrogate | 5 days | Trained surrogate, validation report |
| Calibration | 7 days | Posterior samples, convergence diagnostics |
| Analysis | 5 days | Calibration report, parameter recommendations |
| **Total** | **~3 weeks** | Complete calibrated model with uncertainty |

---

**Document prepared by:** claude-code
**Date:** 2026-02-14
**TMAP8 Version:** val-2a (ion implantation validation case)
**Status:** Planning document - ready for implementation
