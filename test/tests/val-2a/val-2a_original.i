# Validation Problem #2a from TMAP4/TMAP7 V&V document
# Ion Implantation Experiment on Primary Candidate Alloy (PCA)
# Deuterium permeation through 0.5 mm steel sample

# Geometry
thickness = 5e-4 # m - sample thickness (0.5 mm)

# Material properties
diffusivity = 3.0e-10 # m^2/s - deuterium diffusivity in PCA

# Implantation parameters
implantation_depth = 14e-9 # m - average implantation depth (11 nm from TMAP4 instructions, 14 nm from doc)
implantation_sigma = 2.4e-9 # m - standard deviation of implantation depth
implantation_flux = ${fparse 4.9e19 * 0.75} # atom/m^2/s - 75% retention, 25% re-emitted
gaussian_factor = 1.5 # factor to match experimental data

# Boundary condition parameters - Left (upstream, implantation side)
Kr_left_max = 1.0e-27 # m^4/atom/s - maximum recombination coefficient
Kr_left_time_constant = 6.0e-5 # 1/s - time constant for surface cleanup
Kr_left_fraction = 0.9999 # fraction of cleanup

# Boundary condition parameters - Right (downstream, permeation side)
Kr_right = 2.0e-31 # m^4/atom/s - downstream recombination coefficient (constant)

# Time parameters
end_time = 20000 # s - total simulation time

# Beam schedule times (from experimental paper via documentation)
beam_on_1_end = 5820
beam_off_1_end = 9056
beam_on_2_end = 12062
beam_off_2_end = 14572
beam_on_3_end = 17678

[Mesh]
  [generated_mesh]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0
    xmax = ${thickness}
    nx = 100
    bias_x = 1.1 # bias towards left (implantation) side
  []
[]

[Variables]
  [concentration]
    initial_condition = 0.0
  []
[]

[AuxVariables]
  [Kr_left_aux]
    initial_condition = 0.0
  []
[]

[Kernels]
  [time_derivative]
    type = TimeDerivative
    variable = concentration
  []
  [diffusion]
    type = ADMatDiffusion
    variable = concentration
    diffusivity = diffusivity_property
  []
  [source]
    type = BodyForce
    variable = concentration
    function = implantation_source_function
  []
[]

[AuxKernels]
  [Kr_left_calc]
    type = FunctionAux
    variable = Kr_left_aux
    function = Kr_left_function
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[BCs]
  [left_recombination]
    type = ADMatNeumannBC
    variable = concentration
    boundary = left
    value = 1
    boundary_material = flux_left
  []
  [right_recombination]
    type = ADMatNeumannBC
    variable = concentration
    boundary = right
    value = 1
    boundary_material = flux_right
  []
[]

[Functions]
  # Beam schedule function - piecewise constant for instant beam on/off
  [beam_flux_function]
    type = PiecewiseConstant
    x = '0                  ${beam_on_1_end}  ${beam_off_1_end}  ${beam_on_2_end}  ${beam_off_2_end}  ${beam_on_3_end}'
    y = '${implantation_flux} 0                  ${implantation_flux}  0                  ${implantation_flux}  0'
    direction = left
  []

  # Gaussian distribution for implantation depth profile
  [implantation_distribution]
    type = ParsedFunction
    expression = '(${gaussian_factor} / (${implantation_sigma} * sqrt(2.0 * 3.14159265359))) * exp(-0.5 * ((x - ${implantation_depth}) / ${implantation_sigma})^2)'
  []

  # Combined source function
  [implantation_source_function]
    type = CompositeFunction
    functions = 'implantation_distribution beam_flux_function'
  []

  # Time-dependent recombination coefficient on left (upstream) surface
  [Kr_left_function]
    type = ParsedFunction
    expression = '${Kr_left_max} * (1.0 - ${Kr_left_fraction} * exp(-${Kr_left_time_constant} * t))'
  []

  # Adaptive max time step - small during beam-on (high flux), larger during beam-off (decay)
  [max_dt_function]
    type = ParsedFunction
    expression = 'if(t<${beam_on_1_end}, 4,
                  if(t<${beam_off_1_end}, 300,
                  if(t<${beam_on_2_end}, 4,
                  if(t<${beam_off_2_end}, 300,
                  if(t<${beam_on_3_end}, 4, 300)))))'
  []
[]

[Materials]
  [diffusivity_material]
    type = ADGenericConstantMaterial
    prop_names = 'diffusivity_property'
    prop_values = '${diffusivity}'
  []

  # Left boundary flux material (time-dependent Kr)
  [recombination_rate_left]
    type = ADGenericFunctionMaterial
    prop_names = 'Kr_left'
    prop_values = 'Kr_left_function'
  []

  [flux_left]
    type = ADDerivativeParsedMaterial
    coupled_variables = 'concentration'
    property_name = 'flux_left'
    material_property_names = 'Kr_left'
    expression = '-2.0 * Kr_left * concentration^2'
  []

  # Right boundary flux material (constant Kr)
  [recombination_rate_right]
    type = ADGenericConstantMaterial
    prop_names = 'Kr_right'
    prop_values = '${Kr_right}'
  []

  [flux_right]
    type = ADDerivativeParsedMaterial
    coupled_variables = 'concentration'
    property_name = 'flux_right'
    material_property_names = 'Kr_right'
    expression = '-2.0 * Kr_right * concentration^2'
  []
[]

[Postprocessors]
  [time]
    type = TimePostprocessor
  []

  # Recombination flux on left (upstream) surface
  [recombination_flux_left]
    type = ADSideAverageMaterialProperty
    boundary = 'left'
    property = flux_left
    outputs = none
  []

  [scaled_recombination_flux_left]
    type = ScalePostprocessor
    scaling_factor = -1.0
    value = recombination_flux_left
  []

  # Recombination flux on right (downstream) surface
  [recombination_flux_right]
    type = ADSideAverageMaterialProperty
    boundary = 'right'
    property = flux_right
    outputs = none
  []

  [scaled_recombination_flux_right]
    type = ScalePostprocessor
    scaling_factor = -1.0
    value = recombination_flux_right
  []

  # Total inventory
  [total_inventory]
    type = ElementIntegralVariablePostprocessor
    variable = concentration
  []

  [scaled_total_inventory]
    type = ScalePostprocessor
    scaling_factor = 1.0
    value = total_inventory
  []

  # Recombination coefficients for diagnostics
  [Kr_left_pp]
    type = ElementAverageValue
    variable = Kr_left_aux
  []

  [Kr_right_pp]
    type = ADElementAverageMaterialProperty
    mat_prop = Kr_right
  []

  # Adaptive maximum time step size - small during beam-on, larger during beam-off
  [max_dt_pp]
    type = FunctionValuePostprocessor
    function = max_dt_function
    execute_on = 'INITIAL TIMESTEP_END'
    outputs = none
  []

  # Beam flux output for verification
  [beam_flux_pp]
    type = FunctionValuePostprocessor
    function = beam_flux_function
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  nl_max_its = 20

  start_time = 0
  end_time = ${end_time}

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 4
    optimal_iterations = 4
    growth_factor = 1.05
    cutback_factor = 0.5
    # timestep_limiting_postprocessor = max_dt_pp
  []

  automatic_scaling = true
[]

[Outputs]
  file_base = 'val-2a_out'
  [csv]
    type = CSV
  []
  [exodus]
    type = Exodus
    output_material_properties = true
  []
[]
