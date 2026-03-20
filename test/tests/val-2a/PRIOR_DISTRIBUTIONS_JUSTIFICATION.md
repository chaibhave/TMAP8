# Physical Justification for Prior Distributions
## Val-2a Bayesian Calibration

**Date:** 2026-02-14
**Model:** Ion implantation in PCA steel (val-2a)
**Temperature:** 703 K (430°C)

---

## 1. Diffusivity of Deuterium in PCA Steel

### Current Value
`diffusivity = 3.0×10⁻¹⁰ m²/s`

### Proposed Prior
**Distribution:** Log-Uniform
**Range:** `1.0×10⁻¹⁰ to 1.0×10⁻⁹ m²/s`
**Log-space mean:** `3.16×10⁻¹⁰ m²/s`

### Physical Justification

**Temperature Dependence:**
Hydrogen diffusion in metals follows Arrhenius behavior:
```
D = D₀ × exp(-Eₐ / RT)
```

For austenitic steels (PCA is Fe-Cr-Ni alloy):
- Pre-exponential factor D₀: 10⁻⁷ to 10⁻⁶ m²/s
- Activation energy Eₐ: 50-60 kJ/mol

At 703 K, this gives D ≈ 10⁻¹⁰ to 10⁻⁹ m²/s.

**Sources of Uncertainty:**
1. **Microstructure effects:**
   - Grain boundaries can enhance or reduce effective diffusivity
   - PCA has complex austenitic + ferrite phases
   - Grain size: typically 10-100 μm → factor of 2-3 variation

2. **Trapping sites:**
   - Cr carbides, voids, dislocations trap hydrogen
   - Reduces effective diffusivity by factor of 2-10
   - Trap concentration varies with material processing

3. **Composition variations:**
   - PCA: Fe-(10-20)Cr-(10-15)Ni-Mo
   - Cr content affects diffusivity: ±30% variation
   - Minor elements (C, Ti) create traps

4. **Literature scatter:**
   - Measured D at 703 K ranges from 0.5×10⁻¹⁰ to 2×10⁻⁹ m²/s
   - Factor of 4-5 variation between studies
   - Different specimens, surface preparations, measurement techniques

**Prior Choice Rationale:**
- Log-uniform reflects uncertainty spanning order of magnitude
- Range encompasses literature values for austenitic steels at 703 K
- Centered near TMAP4 baseline but allows significant variation
- Log-space appropriate for multiplicative uncertainty

---

## 2. Surface Recombination Coefficient (Upstream, Maximum)

### Current Value
`Kr_left_max = 1.0×10⁻²⁷ m⁴/atom/s`

### Proposed Prior
**Distribution:** Log-Uniform
**Range:** `1.0×10⁻²⁸ to 1.0×10⁻²⁶ m⁴/atom/s`
**Log-space mean:** `3.16×10⁻²⁷ m⁴/atom/s`

### Physical Justification

**Recombination Process:**
Surface recombination follows:
```
J = -2 × Kr × C²
```
where two atoms combine to form H₂ molecule.

Kr depends on:
1. Surface coverage (Langmuir kinetics)
2. Sticking coefficient
3. Surface diffusion rate
4. Desorption activation energy

**Sources of Uncertainty:**
1. **Surface state evolution:**
   - Ion implantation creates defects, roughness
   - Sputtering removes oxide layers → exposes fresh metal
   - Dynamic surface during beam-on periods
   - Factor of 10-100 variation during cleanup

2. **Oxide layers:**
   - Native Cr₂O₃ on PCA: Kr reduced by 10-100×
   - Beam sputtering removes oxides → Kr increases
   - Incomplete oxide removal → partial effect
   - Our time-dependent Kr models this cleanup

3. **Surface contamination:**
   - Carbon, oxygen adsorption reduces Kr
   - Vacuum quality: 10⁻⁶ to 10⁻⁸ Torr range
   - Monolayer formation time: minutes to hours
   - Factor of 5-10 uncertainty

4. **Temperature dependence:**
   - At 703 K, thermal desorption competes with recombination
   - Activation energy: 80-120 kJ/mol for steel surfaces
   - ±50 K temperature uncertainty → factor of 2 in Kr

5. **Literature data:**
   - Clean steel surfaces: Kr = 10⁻²⁷ to 10⁻²⁵ m⁴/atom/s
   - Oxidized surfaces: Kr = 10⁻²⁹ to 10⁻²⁷ m⁴/atom/s
   - After ion bombardment: Kr increases toward clean metal value

**Prior Choice Rationale:**
- Log-uniform spans two orders of magnitude
- Centered on clean surface value after sputtering cleanup
- Upper bound: fully clean, defect-rich surface
- Lower bound: partially oxidized or contaminated surface
- Reflects dynamic surface state during implantation

---

## 3. Surface Recombination Coefficient (Downstream, Constant)

### Current Value
`Kr_right = 2.0×10⁻³¹ m⁴/atom/s`

### Proposed Prior
**Distribution:** Log-Uniform
**Range:** `1.0×10⁻³² to 1.0×10⁻³⁰ m⁴/atom/s`
**Log-space mean:** `3.16×10⁻³¹ m⁴/atom/s`

### Physical Justification

**Why Kr_right << Kr_left:**
The downstream surface is NOT exposed to ion beam, therefore:
1. Maintains native oxide layer (Cr₂O₃)
2. No sputtering → no surface cleanup
3. Lower surface concentration → lower recombination rate
4. Equilibrium with vacuum → stable surface state

**Sources of Uncertainty:**
1. **Oxide thickness:**
   - Native oxide: 2-5 nm on PCA at 703 K
   - Reduces Kr by factor of 100-1000 vs clean surface
   - Thickness varies with pre-treatment: ±factor of 3

2. **Permeation barrier:**
   - Oxide acts as diffusion barrier + recombination barrier
   - Effective Kr depends on oxide permeability
   - Literature: Kr_oxide = 10⁻³² to 10⁻³⁰ m⁴/atom/s

3. **Vacuum-side conditions:**
   - Surface in equilibrium with ~10⁻⁶ Torr vacuum
   - No dynamic processes (unlike upstream)
   - More predictable but still ±order of magnitude uncertainty

4. **Measurement limitations:**
   - Downstream Kr difficult to measure directly
   - Often inferred from permeation experiments
   - Large uncertainty in literature values

**Prior Choice Rationale:**
- Log-uniform centered on oxidized surface value
- Range: oxidized steel surfaces at 703 K
- Two orders of magnitude uncertainty
- Must remain << Kr_left (different surface states)

---

## 4. Gaussian Factor (Implantation Profile Scaling)

### Current Value
`gaussian_factor = 1.5`

### Proposed Prior
**Distribution:** Uniform
**Range:** `1.0 to 2.5`
**Mean:** `1.75`

### Physical Justification

**Purpose:**
Empirical scaling factor to match experimental implantation profile:
```
Source = (gaussian_factor / σ√(2π)) × exp(-0.5 × ((x-μ)/σ)²)
```

**Physical Meaning:**
- Factor of 1.0: Pure Gaussian from SRIM calculation
- Factor > 1.0: Enhanced deposition (e.g., secondary effects, backscattering)
- Factor < 1.0: Reduced effective deposition (e.g., prompt re-emission)

**Sources of Uncertainty:**
1. **SRIM accuracy:**
   - Monte Carlo code for ion implantation
   - Typical accuracy: ±20-30% for range and straggling
   - Depends on nuclear stopping power models
   - Better for heavy ions, less certain for H/D

2. **Retention fraction:**
   - Not all implanted atoms retained in lattice
   - Prompt re-emission: 20-30% for D in steel at 703 K
   - Already accounted partially (75% retention), but uncertain
   - Temperature-dependent: ±factor of 1.5

3. **Surface roughness:**
   - Implantation profile broadened by surface topology
   - Rough surfaces: enhanced effective depth
   - PCA typically polished: Ra < 0.1 μm
   - Minor effect but adds 10-20% uncertainty

4. **Channeling effects:**
   - Ions can channel along crystal directions
   - Increases penetration depth in single grains
   - Polycrystalline PCA averages this out mostly
   - Small residual effect: ±10%

5. **Beam divergence:**
   - Real beam has angular spread (not normal incidence)
   - Broadens effective profile
   - Typical spread: ±5-10° → 10-15% effect

**Prior Choice Rationale:**
- Uniform (not log): Factor is multiplicative but moderate range
- Minimum of 1.0: Pure SRIM prediction
- Maximum of 2.5: Accounts for combined uncertainties
- Current value 1.5 is conservative central estimate
- Linear prior appropriate for scaling factor

---

## 5. Surface Cleanup Time Constant

### Current Value
`Kr_left_time_constant = 6.0×10⁻⁵ s⁻¹`

### Proposed Prior
**Distribution:** Log-Uniform
**Range:** `1.0×10⁻⁵ to 2.0×10⁻⁴ s⁻¹`
**Log-space mean:** `4.47×10⁻⁵ s⁻¹`

### Physical Justification

**Process:**
Sputtering removes surface oxide/contamination exponentially:
```
Kr(t) = Kr_max × (1 - fraction × exp(-time_constant × t))
```

Time constant τ = 1 / time_constant ≈ 16,667 s (4.6 hours at current value)

**Sputtering Dynamics:**
1. **Sputtering yield:**
   - For 3 keV D⁺ on PCA: Y ≈ 0.1-0.3 atoms/ion
   - Flux: 3.7×10¹⁹ ions/m²/s
   - Removal rate: ~10¹⁸ atoms/m²/s

2. **Oxide layer removal:**
   - Native Cr₂O₃: ~5 nm = 10¹⁹ atoms/m²
   - Removal time: 10-100 seconds
   - But: simultaneous re-oxidation from residual gas
   - Dynamic equilibrium reached after 1000-10000 s

3. **Surface roughening:**
   - Sputtering creates surface defects
   - Increases surface area → increases Kr
   - Accumulates over beam-on period (hours)
   - Slower than oxide removal

**Sources of Uncertainty:**
1. **Sputtering yield:**
   - Depends on surface composition (evolves during beam)
   - Angle of incidence: typically some divergence
   - Energy distribution of beam: ±10-20%
   - Factor of 2-3 uncertainty

2. **Re-oxidation rate:**
   - Depends on vacuum quality (O₂, H₂O partial pressure)
   - Experimental conditions: not precisely known
   - Surface temperature: local heating by beam
   - Factor of 3-5 uncertainty

3. **Steady-state assumption:**
   - Assuming first-order exponential approach
   - Reality: more complex multi-stage process
   - Effective time constant averages multiple processes
   - Could be faster or slower

**Prior Choice Rationale:**
- Log-uniform: Exponential process, order of magnitude uncertainty
- Range: 5,000 to 100,000 seconds (1.4 hours to 28 hours)
- Lower bound: Fast sputtering, clean vacuum
- Upper bound: Slow cleanup, poor vacuum, re-oxidation
- Encompasses beam-on timescales (1600-2900 seconds)

---

## 6. Surface Cleanup Fraction

### Current Value
`Kr_left_fraction = 0.9999`

### Proposed Prior
**Distribution:** Beta(α=50, β=2)
**Range:** `0.95 to 1.0` (practical range, support is [0,1])
**Mean:** `0.962`
**Mode:** `0.980`

### Physical Justification

**Meaning:**
Fraction of increase from minimum to maximum Kr during cleanup:
```
Kr_final / Kr_initial = 1 / (1 - fraction)
```
At fraction = 0.9999: Factor of 10,000 increase
At fraction = 0.999: Factor of 1,000 increase
At fraction = 0.99: Factor of 100 increase

**Physical Expectations:**
1. **Complete cleanup unlikely:**
   - Even with sputtering, some oxide/contamination remains
   - Residual oxygen in vacuum continuously deposits
   - Never reach perfectly clean surface: fraction < 1.0

2. **Large increase expected:**
   - Native oxide reduces Kr by 100-1000×
   - Sputtering removes most oxide → most of increase
   - Fraction should be > 0.95 (>95% of maximum effect)

3. **Asymptotic behavior:**
   - Kr increases rapidly initially, then plateaus
   - Never quite reaches clean surface value
   - Final 1% of cleanup is slowest

**Sources of Uncertainty:**
1. **Initial surface state:**
   - How oxidized was surface before implantation?
   - Pre-treatment history unknown
   - Affects baseline Kr and maximum possible increase

2. **Vacuum quality:**
   - Better vacuum → closer to fraction = 1.0
   - Poor vacuum → more re-oxidation → lower fraction
   - Typical permeation experiments: 10⁻⁶ Torr

3. **Steady-state Kr value:**
   - Measured permeation reflects quasi-steady-state
   - Balance between sputtering and re-oxidation
   - Fraction represents this balance

**Prior Choice Rationale:**
- Beta distribution: Naturally bounded to [0,1]
- Parameters α=50, β=2: Skewed toward 1.0 but with tail
- Mean ≈ 0.96: Most of cleanup achieved
- Allows values from 0.95 to >0.99
- Reflects physical expectation: large increase but not perfect

---

## 7. Implantation Depth

### Current Value
`implantation_depth = 14×10⁻⁹ m` (14 nm)

### Proposed Prior
**Distribution:** Normal (Truncated at 0)
**Mean:** `14×10⁻⁹ m`
**Std Dev:** `2×10⁻⁹ m`
**95% CI:** `[10, 18] nm`

### Physical Justification

**SRIM Prediction:**
For 3 keV D⁺ → PCA steel (approximated as Fe):
- Projected range: 11-14 nm (depends on composition)
- Code used: SRIM-2013 or similar
- Monte Carlo: 10⁴-10⁶ ion histories

**Sources of Uncertainty:**
1. **SRIM accuracy:**
   - Nuclear stopping power: ±10-15% uncertainty
   - Electronic stopping power: ±5-10%
   - Total range uncertainty: ±15-20% (±2-3 nm at 14 nm)

2. **Material composition:**
   - PCA is Fe-Cr-Ni alloy, not pure Fe
   - SRIM run for average composition
   - Real sample: composition gradients, phases
   - Affects stopping power: ±10% effect

3. **Beam energy:**
   - Nominal 3 keV, but has spread: ±100 eV typical
   - Energy uncertainty: ±3% → depth uncertainty ±3%
   - Plasma potential variations: additional ±50 eV

4. **Surface effects:**
   - Surface roughness: ±0.1 μm << depth (negligible)
   - Oxide layer: ~5 nm → affects first monolayers
   - Ion beam sputters oxide early, then into bulk

5. **Crystallographic effects:**
   - Single crystal: channeling can increase range by 2-3×
   - PCA is polycrystalline: averaged out mostly
   - Residual texture effects: ±5-10%

**Prior Choice Rationale:**
- Normal distribution: Appropriate for SRIM uncertainty (multiple random effects)
- Mean at 14 nm: SRIM prediction for 3 keV D → Fe/PCA
- Std dev 2 nm: Encompasses SRIM uncertainty (±15%) and composition effects
- 95% CI [10, 18] nm: Reasonable physical range
- Truncated at 0: Negative depth nonphysical

---

## 8. Implantation Depth Standard Deviation (Straggling)

### Current Value
`implantation_sigma = 2.4×10⁻⁹ m` (2.4 nm)

### Proposed Prior
**Distribution:** Normal (Truncated at 0)
**Mean:** `2.4×10⁻⁹ m`
**Std Dev:** `0.5×10⁻⁹ m`
**95% CI:** `[1.4, 3.4] nm`

### Physical Justification

**Straggling:**
Range straggling quantifies spread in ion ranges due to:
1. Statistical variations in collision sequences
2. Different paths through material
3. Energy loss fluctuations

SRIM provides:
- Projected range: R_p ≈ 14 nm
- Longitudinal straggling: ΔR_p ≈ 2-3 nm
- Ratio ΔR_p / R_p ≈ 0.15-0.20 (typical for keV ions)

**Sources of Uncertainty:**
1. **SRIM straggling prediction:**
   - Based on statistical model (Lindhard theory)
   - Uncertainty: ±20-30% on straggling parameter
   - Less well-validated than range itself

2. **Collision cascade effects:**
   - Multiple collisions broaden distribution
   - Heavy ions (D) have more collisions than light (H)
   - PCA composition: variations in collision cross-sections

3. **Surface roughness:**
   - Rough surface artificially broadens profile
   - PCA polished: Ra < 0.1 μm
   - Effect: ±0.1-0.2 nm (small)

4. **Thermal diffusion during implantation:**
   - At 703 K, D is mobile (D ≈ 3×10⁻¹⁰ m²/s)
   - During 1600 s beam-on: √(Dt) ≈ 7 nm
   - Broadens profile beyond SRIM prediction
   - BUT: Gaussian factor partially accounts for this

5. **Beam divergence:**
   - Angular spread of beam: ±5-10°
   - Projects range onto angled paths
   - Effective broadening: 5-10%

**Prior Choice Rationale:**
- Normal distribution: Straggling is inherently Gaussian (CLT)
- Mean 2.4 nm: SRIM prediction for 3 keV D
- Std dev 0.5 nm: Reflects ±20% uncertainty from SRIM and thermal effects
- 95% CI [1.4, 3.4] nm: Physical range for keV ions
- Truncated at 0: Negative sigma nonphysical

---

## Summary Table

| Parameter | Distribution | Range/Parameters | Physical Basis |
|-----------|-------------|------------------|----------------|
| `diffusivity` | Log-Uniform | [10⁻¹⁰, 10⁻⁹] m²/s | Literature scatter, trapping, microstructure |
| `Kr_left_max` | Log-Uniform | [10⁻²⁸, 10⁻²⁶] m⁴/atom/s | Surface state, sputtering cleanup, oxides |
| `Kr_right` | Log-Uniform | [10⁻³², 10⁻³⁰] m⁴/atom/s | Oxidized surface, no beam exposure |
| `gaussian_factor` | Uniform | [1.0, 2.5] | SRIM accuracy, retention, secondary effects |
| `Kr_left_time_constant` | Log-Uniform | [10⁻⁵, 2×10⁻⁴] s⁻¹ | Sputtering rate, re-oxidation balance |
| `Kr_left_fraction` | Beta(50, 2) | [0.95, 1.0] | Extent of surface cleanup, vacuum quality |
| `implantation_depth` | Normal(14, 2) nm | μ=14, σ=2 nm | SRIM uncertainty, composition, energy |
| `implantation_sigma` | Normal(2.4, 0.5) nm | μ=2.4, σ=0.5 nm | Straggling, thermal broadening, SRIM |

---

## References

1. **Diffusivity in Steels:**
   - Hagi, H. (1990). "Diffusion coefficient of hydrogen in iron" ISIJ International, 30(6), 417-425.
   - Perng, T.P. & Altstetter, C.J. (1986). "Hydrogen effects in austenitic stainless steels" Materials Science and Engineering, 83(2), 199-208.

2. **Surface Recombination:**
   - Pick, M.A. & Sonnenberg, K. (1985). "A model for atomic hydrogen-metal interactions" Journal of Nuclear Materials, 131(2-3), 208-220.
   - Baskes, M.I. (1980). "A calculation of the surface recombination rate constant for hydrogen isotopes on metals" Journal of Nuclear Materials, 92(2-3), 318-324.

3. **Ion Implantation:**
   - Ziegler, J.F., Biersack, J.P. & Littmark, U. (1985). "The Stopping and Range of Ions in Solids" Pergamon Press.
   - SRIM-2013 software and documentation (www.srim.org)

4. **TMAP Validation:**
   - Longhurst, G.R. et al. (1992). "Verification and validation of TMAP4" EGG-FSP-10347, Idaho National Engineering Laboratory.

5. **Uncertainty Quantification:**
   - Kennedy, M.C. & O'Hagan, A. (2001). "Bayesian calibration of computer models" Journal of the Royal Statistical Society: Series B, 63(3), 425-464.

---

**Document Status:** Complete - Ready for implementation
**Prepared by:** claude-code
**Next Step:** Create parameterized input files with these prior distributions
