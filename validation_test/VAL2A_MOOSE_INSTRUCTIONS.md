# MOOSE Simulation Setup Instructions: val-2a Ion Implantation Experiment

## Experiment Overview
This simulation recreates the TMAP4 val-2a validation case: an ion implantation experiment conducted at INEL in 1985. The experiment measures deuterium permeation through a Primary Candidate Alloy (PCA - modified 316 stainless steel) sample during pulsed ion beam exposure.

**Reference**: TMAP4 Verification and Validation Report, Section 2a, page 49-50

## Physical Problem Description

### Geometry
- **1D slab geometry** (through-thickness transport)
- **Sample thickness**: 0.5 mm (5×10⁻⁴ m)
- **Sample diameter**: 2.5 cm (for reference only; 1D model)
- **Domain**: 0 ≤ x ≤ 5×10⁻⁴ m

### Material
- Primary Candidate Alloy (PCA) - modified 316 stainless steel
- Temperature: 703 K (constant, isothermal)

## Physics Requirements

### 1. Diffusion Physics
Implement transient hydrogen diffusion using Fick's law:
```
∂C/∂t = ∇·(D∇C) + S
```

Where:
- C = deuterium concentration (atoms/m³)
- D = diffusion coefficient (m²/s)
- S = volumetric source term (atoms/m³/s)

**Material Properties**:
- Diffusion coefficient: D = 3.0×10⁻¹⁰ m²/s (constant)
- No temperature dependence needed (isothermal case)

### 2. Ion Implantation Source
Implement a **subsurface volumetric source** representing ion implantation:

**Implantation parameters**:
- Average implantation depth: 11 nm (1.1×10⁻⁸ m) from surface
- Implantation depth spread: ±5.4 nm (Gaussian distribution recommended)
- Incident flux: 6.53×10¹⁹ D/m²/s
- Retention fraction: 75% (25% immediately re-emitted)
- Effective source flux: 4.9×10¹⁹ D/m²/s

**Time-dependent beam schedule** (use PiecewiseLinear or similar):
```
Time (s)          Flux (D/m²/s)
0 - 6420          4.9×10¹⁹ (beam on)
6420 - 9420       0.0 (beam off)
9420 - 12480      4.9×10¹⁹ (beam on)
12480 - 14940     0.0 (beam off)
14940 - 18180     4.9×10¹⁹ (beam on)
18180 - end       0.0 (beam off)
```

**Implementation approach**:
- Apply source over first 5-6 mesh elements near x=0
- Source intensity should peak at ~11 nm depth
- Use a Gaussian or box function centered at implantation depth
- Total integrated source must equal the retention-corrected flux

### 3. Boundary Conditions

#### Left Boundary (x = 0, Implantation Surface)
**Recombination-limited flux BC** (mixed BC):

The flux leaving the surface is governed by:
```
J = k_r * C_s² - k_d * P^0.5
```

Where:
- J = flux (D/m²/s)
- k_r = recombination rate coefficient (m⁴/D/s)
- k_d = dissociation rate coefficient (D/m²/s/Pa^0.5)
- C_s = surface concentration (D/m³)
- P = upstream pressure (Pa)

**Time-dependent rate coefficients** (surface cleanup model):
```
k_d(t) = 8.959×10¹⁸ × (1 - 0.9999×exp(-6.0×10⁻⁵×t)) [D/m²/s/Pa^0.5]
k_r(t) = 1.0×10⁻²⁷ × (1 - 0.9999×exp(-6.0×10⁻⁵×t)) [m⁴/D/s]
```

**Upstream pressure schedule** (Pa):
```
Time (s)          Pressure (Pa)
0 - 6420          4.0×10⁻⁵
6420 - 9420       9.0×10⁻⁶
9420 - 12480      4.0×10⁻⁵
12480 - 14940     9.0×10⁻⁶
14940 - 18180     4.0×10⁻⁵
18180 - end       9.0×10⁻⁶
```

**Implementation notes**:
- This BC may require a custom IntegratedBC or NeumannBC
- The recombination term is nonlinear (quadratic in concentration)
- May need to use ParsedFunction for time-dependent coefficients

#### Right Boundary (x = L, Downstream Surface)
**Similar recombination BC but with different parameters**:
```
k_d = 1.7918×10¹⁵ [D/m²/s/Pa^0.5] (constant)
k_r = 2.0×10⁻³¹ [m⁴/D/s] (constant, essentially irreversible)
P_downstream = 2.0×10⁻⁶ Pa (constant)
```

### 4. Initial Conditions
- **Initial concentration**: C(x,0) = 0.0 D/m³ (everywhere)
- **Initial temperature**: T = 703 K (constant)

## Mesh Requirements

### Spatial Discretization
The mesh must be **highly refined near the implantation zone** and can be coarser in the bulk:

**Recommended mesh structure** (from TMAP4 input):
```
Region              Number of Elements    Element Size
0 - 20 nm           5 elements            4 nm each
20 - 30 nm          1 element             10 nm
30 - 130 nm         1 element             100 nm
130 - 1130 nm       1 element             1 μm
1130 - 11130 nm     1 element             10 μm
11130 - 500000 nm   10 elements           ~48.8 μm each
```

**Total nodes**: 21 (20 elements)

**Implementation**:
- Use GeneratedMesh with bias or
- Use MultiAppGeometricRestart with refinement zones or
- Manually specify mesh blocks with different sizing

### Temporal Discretization
- **Total simulation time**: 19,200 s (320 minutes, ~5.3 hours)
- **Initial time step**: 20 s
- **Time stepping**: Adaptive recommended, but fixed Δt = 20 s acceptable
- **Output frequency**: Every 60 s (3 outputs per time step)

## Required MOOSE Objects

### Modules/Actions
```
[GlobalParams]
  temperature = 703
[]

[Mesh]
  # 1D mesh with variable spacing (see above)
[]

[Variables]
  [./concentration]
    initial_condition = 0.0
  [../]
[]
```

### Kernels
1. **TimeDerivative**: ∂C/∂t term
2. **MatDiffusion** or **Diffusion**: D∇²C term
3. **BodyForce** or **CoupledForce**: Volumetric implantation source S(x,t)

### Materials
1. **Diffusivity**: D = 3.0×10⁻¹⁰ m²/s
2. **Temperature**: T = 703 K (constant)

### BCs
1. **Left boundary**: Custom recombination BC (see equations above)
2. **Right boundary**: Custom recombination BC with different parameters

### Functions
1. **Implantation flux schedule**: PiecewiseLinear for beam on/off
2. **Upstream pressure schedule**: PiecewiseLinear
3. **Surface cleanup function**: ParsedFunction for exponential time dependence
4. **Implantation depth profile**: ParsedFunction or similar for Gaussian distribution

### AuxVariables & AuxKernels
Consider adding:
- **mobile_concentration**: concentration in mobile lattice sites
- **flux**: to compute and visualize flux at boundaries
- **implantation_source**: to visualize source distribution

### Postprocessors
**Required outputs for validation**:
1. **PointValue**: Concentration at specific locations (x = 0, x = L)
2. **SideFluxIntegral** or custom: Downstream permeation flux (primary validation metric)
3. **SideFluxIntegral**: Upstream flux (for mass balance checking)
4. **ElementIntegralVariablePostprocessor**: Total D inventory in sample
5. **TimePostprocessor**: For time tracking

**Key metric**: Downstream permeation flux vs. time should match Figure 9 in TMAP4 V&V report

## Validation Criteria

### Qualitative Checks
1. **Breakthrough time**: Should show delayed permeation breakthrough (~few thousand seconds)
2. **Steady-state flux**: When beam is on continuously, flux should approach steady value
3. **Beam-off transients**: Flux should decay when beam turns off
4. **Multiple cycles**: Three beam-on periods should show consistent behavior

### Quantitative Comparison
Compare your MOOSE results against TMAP4 results shown in **Figure 9** (page 1261):
- Permeation flux magnitude during beam-on periods
- Time to reach quasi-steady state
- Decay rate during beam-off periods
- Overall transient shape

**Expected behavior**:
- Initial breakthrough around 1000-2000 s
- Quasi-steady flux: ~2-4×10¹⁸ D/m²/s (order of magnitude, read from figure)
- Rapid rise when beam turns on
- Gradual decay when beam turns off

### Known Modeling Limitations
From TMAP4 report:
1. Surface cleanup modeled as exponential in time (not ion fluence)
2. 1D approximation (ignores edge effects on 2.5 cm diameter sample)
3. No chamber wall pumping effects
4. No surface contamination layers explicitly modeled

## Implementation Steps

### Step 1: Create Basic Input File
- Set up 1D mesh with proper refinement
- Add concentration variable
- Add diffusion kernel with constant D

### Step 2: Add Implantation Source
- Implement subsurface source at ~11 nm depth
- Add time-dependent flux schedule
- Verify source term integrates correctly

### Step 3: Implement Boundary Conditions
- Start with simplified BCs (e.g., fixed concentration) to test
- Implement full recombination BC on downstream side
- Add time-dependent recombination BC on upstream side

### Step 4: Add Postprocessors
- Track downstream flux (primary output)
- Track inventory for mass balance
- Add any diagnostic outputs

### Step 5: Run and Validate
- Run full 19,200 s simulation
- Compare permeation transient to Figure 9
- Check mass balance
- Verify physically reasonable results

### Step 6: Document Results
- Plot permeation flux vs. time
- Compare quantitatively with TMAP4
- Document any differences and potential causes
- Consider sensitivity to mesh refinement, time step size

## Additional Resources

### Required Files
- TMAP4 input file: `VAL-2A.INP` (lines 2136-2219 in TMAP4_V&V.md)
- Validation figure: Figure 9 (page 1261, file: `_page_57_Figure_0.jpeg`)

### Useful MOOSE Resources
- **TMAP8** (MOOSE-based tritium transport code): May have similar BCs already implemented
- **Chemical Reactions module**: For surface chemistry
- **ADKernels**: Consider using automatic differentiation for better convergence

### Physics References
- TRIM code: For ion implantation depth calculations
- Sieverts' law: For hydrogen solubility (C ∝ √P)
- McNabb-Foster model: If adding trapping (not needed here)

## Expected Deliverables

1. **MOOSE input file** (`.i` file) implementing the simulation
2. **Plot** of downstream permeation flux vs. time
3. **Comparison plot** overlaying MOOSE results with TMAP4/experimental data
4. **Brief validation report** documenting agreement and any discrepancies
5. **Convergence studies** (mesh, time step) if time permits

## Contact/Questions

When implementing, consider:
- What MOOSE version are you using?
- Do you have access to TMAP8 or similar hydrogen transport applications?
- Do you need custom BCs or can you use existing MOOSE objects?
- What output format is preferred for comparison plots?

---

**Document Version**: 1.0
**Created**: 2026-02-14
**Based on**: TMAP4 V&V Report, Section 2a (val-2a.inp)
