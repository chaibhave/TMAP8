 Based on your experience working on this problem, what information / skills / hooks / .MD files do you think I could have provided
   to
    improve your workflow? Refer to saved conversations in this folder for context, you can also look at the git diff to see what
  changes you made.

⏺ I'll analyze the project to understand the work that was done and provide recommendations for improving the workflow.

⏺ Searched for 3 patterns (ctrl+o to expand)

⏺ Bash(git diff --stat)
  ⎿   moose                                  |   0
      test/tests/val-2a/comparison_val-2a.py |  98 ++++++++-
      test/tests/val-2a/val-2a.i             | 349 +++++++++++++++++++--------------
      3 files changed, 292 insertions(+), 155 deletions(-)

⏺ Bash(git diff test/tests/val-2a/comparison_val-2a.py | head -100)
  ⎿  diff --git a/test/tests/val-2a/comparison_val-2a.py b/test/tests/val-2a/comparison_val-2a.py
     index 7ad49f12..7fa7c4f9 100644
     --- a/test/tests/val-2a/comparison_val-2a.py
     … +36 lines (ctrl+o to expand)

⏺ Bash(git diff test/tests/val-2a/val-2a.i | head -100)
  ⎿  diff --git a/test/tests/val-2a/val-2a.i b/test/tests/val-2a/val-2a.i
     index 85d41569..93b2bf1a 100644
     --- a/test/tests/val-2a/val-2a.i
     … +39 lines (ctrl+o to expand)

⏺ Recalled 1 memory, searched for 1 pattern, read 6 files (ctrl+o to expand)

⏺ Based on my analysis of the work done (extensive Bayesian calibration framework for TMAP8 validation), here are my recommendations for
  improving future workflows:

  1. Project-Specific CLAUDE.md File

  Create /Users/bhavcv/projects/TMAP8/CLAUDE.md with:

  # TMAP8 Project Context

  ## Project Overview
  - TMAP8: Tritium Migration Analysis Program (MOOSE-based)
  - Primary language: C++ (MOOSE framework) + Python (analysis scripts)
  - Input files: `.i` format (MOOSE/libMesh syntax)

  ## Key Conventions
  - Executable location: `~/projects/TMAP8/tmap8-opt`
  - Test structure: `test/tests/<test-name>/<test-name>.i`
  - Gold files: `test/tests/<test-name>/gold/`
  - Units: Always use SI with explicit conversion comments
  - Parameter names: Use physics-based naming (e.g., `Kr_left_fraction`, not `param1`)

  ## Common Workflows
  - Run simulation: `~/projects/TMAP8/tmap8-opt -i <input>.i [param=value ...]`
  - Typical runtime: ~18 seconds for val-2a test
  - Output format: CSV files with time series data

  ## Physics Domain Knowledge
  - Hydrogen isotope transport in materials
  - Ion implantation: typically 10-50 nm depth with Gaussian profile
  - Recombination coefficients: Kr units are m^4/(atom·s)
  - Diffusivity: D units are m^2/s
  - Soret effect: often negligible at room temperature

  ## Documentation Standards
  - Use `moose-docs` skill for MOOSE class documentation
  - Scientific workflows need: README + QUICK_START + detailed docs
  - Always include computational cost estimates in calibration docs

  2. Custom Skills to Create

  a) tmap8-run Skill

  name: tmap8-run
  description: Run TMAP8 simulations with proper error handling and output parsing
  tools: [Bash, Read, Grep]
  This would handle common patterns like:
  - Checking executable exists
  - Running with timeout
  - Parsing CSV outputs
  - Comparing to gold files

  b) bayesian-calibration-setup Skill

  name: bayesian-calibration-setup
  description: Set up Bayesian calibration workflow for MOOSE-based simulations
  This would automate:
  - Creating parameterized input files
  - Setting up sensitivity analysis
  - Creating calibration scripts with proper priors
  - Generating documentation structure

  c) scientific-workflow-docs Skill

  name: scientific-workflow-docs
  description: Generate comprehensive documentation for scientific computing workflows
  Would create:
  - README (overview)
  - QUICK_START (TL;DR commands)
  - Detailed methodology docs
  - Troubleshooting guides
  - File organization diagrams

  3. Git Hooks

  Pre-commit Hook (.git/hooks/pre-commit)

  #!/bin/bash
  # Prevent committing large simulation output directories
  if git diff --cached --name-only | grep -qE "calibration.*ensemble/|sensitivity_results/"; then
      echo "❌ Error: Attempting to commit large simulation output directory"
      echo "Add to .gitignore instead"
      exit 1
  fi

  # Warn if committing .i file without gold file update
  if git diff --cached --name-only | grep -q "\.i$"; then
      echo "⚠️   Warning: MOOSE input file changed. Consider updating gold files."
  fi

  Post-checkout Hook

  #!/bin/bash
  # Remind about submodule updates for MOOSE
  if git diff --name-only HEAD@{1} HEAD | grep -q "^moose$"; then
      echo "📦 MOOSE submodule changed. Run: git submodule update --init --recursive"
  fi

  4. Context Files

  /test/tests/val-2a/CONTEXT.md

  # Val-2a Test Context

  ## Physical System
  - Material: Primary Candidate Alloy (PCA) steel
  - Isotope: Deuterium
  - Sample thickness: 0.5 mm
  - Implantation energy: Creates 11-14 nm depth profile
  - Temperature: Room temperature (~300K)

  ## Experimental Data Source
  - Original paper: [citation]
  - Experiment date: [date]
  - Key measurements: Permeation flux vs time
  - Beam schedule: 3 on/off cycles

  ## Known Issues
  - Baseline model underpredicts peaks (~25% RMSPE)
  - Surface recombination dominates (not bulk diffusion)
  - Soret effect negligible at this temperature

  ## Parameter Sensitivity (from analysis)
  - **HIGH**: Kr_left_fraction (ST=1.22), gaussian_factor (ST=0.12)
  - **LOW**: All other 6 parameters (ST<0.01)

  5. Memory File Enhancement

  Add to /Users/bhavcv/.claude/projects/-Users-bhavcv-projects-TMAP8/memory/MEMORY.md:

  # TMAP8 Project Memory

  ## MOOSE Input File Patterns
  - Always use `[Mesh]`, `[Variables]`, `[Kernels]`, `[BCs]` blocks
  - Units: Define with `fparse` for clarity: `${fparse 3e-10 # m^2/s}`
  - Parameters: Define at top for easy CLI override
  - Postprocessors: Use for flux/inventory tracking

  ## Bayesian Calibration Workflow
  1. Sensitivity analysis first (Sobol with N=512 samples)
  2. Identify important parameters (ST > 0.1)
  3. Fix insensitive parameters at nominal values
  4. LHS sampling for surrogate training
  5. GP surrogate (subsample time by 5x for efficiency)
  6. MCMC with emcee (32 walkers, 5000 steps, 1000 burn-in)

  ## Common Parameter Ranges (Physical)
  - Diffusivity D: 1e-12 to 1e-9 m^2/s (steel, room temp)
  - Recombination Kr: 1e-32 to 1e-26 m^4/(atom·s)
  - Solubility: 1e18 to 1e22 atoms/m^3/Pa^0.5
  - Implantation depth: 5-50 nm (typical ion energies)

  ## User Preferences
  - Prefers thorough documentation with quick-start guides
  - Likes computational cost estimates in docs
  - Values comparison of full vs reduced approaches

  6. Additional Documentation Templates

  COMPUTATIONAL_WORKFLOW_TEMPLATE.md

  A template for documenting scientific computing workflows with sections for:
  - TL;DR commands
  - Prerequisites (dependencies, executables)
  - Step-by-step workflow
  - Expected outputs with examples
  - Computational cost breakdown
  - Troubleshooting decision tree
  - Comparison tables

  Why These Would Help

  1. CLAUDE.md: Would have given me immediate context about MOOSE syntax, units conventions, and typical workflows
  2. tmap8-run skill: Would have standardized simulation execution instead of me figuring out the command pattern each time
  3. Git hooks: Would have prevented potential issues with large files and reminded about gold file updates
  4. CONTEXT.md: Would have provided physics domain knowledge upfront, reducing time spent understanding the problem
  5. Memory file: Would have preserved learnings about parameter ranges and workflow patterns across sessions
  6. Scientific workflow skill: Would have automated the creation of the comprehensive documentation structure I ended up creating

  The most impactful would be #1 (CLAUDE.md) and #4 (test-specific CONTEXT.md) as they provide immediate, always-available context
  without requiring tool invocations.


One issue I noticed with your iterations was that you never questioned that PiecewiseLinear might be generating a different ion
  implantation behavior than PiecewiseConstant. What could I have done to ensure you spotted that quicker / didn't miss that in the
  first place? Make sure the solution is general and not specific to this instance

⏺ Excellent question. This highlights a critical gap: I modified computational model code without verifying behavioral equivalence. Here
  are general strategies to catch this class of error:

  1. Mandatory Verification Protocol in CLAUDE.md

  Add to /Users/bhavcv/projects/TMAP8/CLAUDE.md:

  ## Critical Rule: Computational Model Changes

  ⚠️  **BLOCKING REQUIREMENT**: Before declaring any simulation input file change complete:

  1. **Baseline Capture**: Run original version, save output with `_baseline` suffix
  2. **Modified Run**: Run modified version with same parameters
  3. **Comparison**: Quantitative diff of key outputs (must show < 1% difference unless intentional)
  4. **Document Intent**: State explicitly what SHOULD change vs what should stay same

  ### Example Workflow
  ```bash
  # Before modifying test.i:
  ./tmap8-opt -i test.i          # generates test_out.csv
  cp test_out.csv test_baseline.csv

  # After modifications:
  ./tmap8-opt -i test.i          # generates test_out.csv
  python compare_outputs.py test_baseline.csv test_out.csv

  # Must show: "✓ Outputs match within 0.1%" OR
  #           "⚠ Outputs differ by X% - EXPECTED because [reason]"

  Red Flags That Require Verification

  - Changing function types (PiecewiseLinear ↔ PiecewiseConstant ↔ others)
  - Modifying mesh parameters (nx, bias, geometry)
  - Changing numerical schemes (time integration, solvers)
  - Updating boundary condition types
  - Refactoring without explicit physics changes

  ## **2. Git Pre-Commit Hook with Automated Checks**

  `.git/hooks/pre-commit`:
  ```bash
  #!/bin/bash

  # Detect MOOSE input file changes
  changed_inputs=$(git diff --cached --name-only --diff-filter=M | grep '\.i$')

  if [ -n "$changed_inputs" ]; then
      echo "═══════════════════════════════════════════════════════════"
      echo "⚠️   MOOSE INPUT FILE MODIFIED"
      echo "═══════════════════════════════════════════════════════════"
      echo "Changed files: $changed_inputs"
      echo ""
      echo "Required verification (unless new file):"
      echo "  1. Baseline output exists with '_baseline' or '_gold' suffix?"
      echo "  2. Quantitative comparison performed?"
      echo "  3. Differences documented (if any)?"
      echo ""
      echo "Run verification:"
      echo "  ./scripts/verify_model_change.sh <test-name>"
      echo ""
      read -p "Verification complete? (y/n) " -n 1 -r
      echo
      if [[ ! $REPLY =~ ^[Yy]$ ]]; then
          echo "❌ Commit aborted. Complete verification first."
          exit 1
      fi
  fi

  3. Verification Helper Script

  /scripts/verify_model_change.sh:
  #!/bin/bash
  # Usage: verify_model_change.sh <test-name>

  TEST=$1
  BASELINE="${TEST}_baseline_out.csv"
  CURRENT="${TEST}_out.csv"

  if [ ! -f "$BASELINE" ]; then
      echo "❌ No baseline found. Run: mv ${CURRENT} ${BASELINE}"
      exit 1
  fi

  # Run comparison (example with Python)
  python3 - <<EOF
  import pandas as pd
  import numpy as np

  baseline = pd.read_csv('$BASELINE')
  current = pd.read_csv('$CURRENT')

  # Check column alignment
  if not all(baseline.columns == current.columns):
      print("❌ Column mismatch! Output structure changed.")
      exit(1)

  # Compute relative differences for all numeric columns
  max_diff = 0
  for col in baseline.select_dtypes(include=[np.number]).columns:
      rel_diff = np.max(np.abs(baseline[col] - current[col]) / (np.abs(baseline[col]) + 1e-30))
      max_diff = max(max_diff, rel_diff)
      if rel_diff > 0.01:  # 1% threshold
          print(f"⚠️   {col}: {rel_diff*100:.2f}% difference")

  if max_diff < 0.001:
      print(f"✅ VERIFIED: Outputs identical (<0.1% difference)")
  elif max_diff < 0.01:
      print(f"✅ VERIFIED: Outputs nearly identical ({max_diff*100:.2f}% difference)")
  else:
      print(f"❌ SIGNIFICANT CHANGE: {max_diff*100:.2f}% difference")
      print("   Document reason in commit message or revert changes.")
      exit(1)
  EOF

  4. Task-Based Verification Checklist

  Create a model-verification skill that I must invoke:

  # ~/.claude/skills/model-verification.md
  ---
  name: model-verification
  description: Verification checklist for computational model changes
  ---

  When modifying simulation input files, follow this checklist:

  ## Pre-Modification
  - [ ] Capture baseline output: `./run_simulation.sh > baseline.out`
  - [ ] Document baseline metrics (runtime, key outputs, checksums)
  - [ ] Identify what SHOULD change vs what should stay constant

  ## Post-Modification
  - [ ] Run modified simulation with identical inputs
  - [ ] Quantitative comparison of outputs (< 1% diff unless intended)
  - [ ] Check side effects:
    - [ ] Runtime changed significantly?
    - [ ] Convergence behavior different?
    - [ ] Memory usage different?
    - [ ] Output file structure/columns same?

  ## Red Flags (require extra scrutiny)
  - [ ] Function type changes (PiecewiseLinear vs PiecewiseConstant vs Function)
  - [ ] Time discretization changes (dt, scheme)
  - [ ] Mesh changes (nx, bias, element type)
  - [ ] Boundary condition type changes
  - [ ] Solver/preconditioner changes

  ## Documentation
  - [ ] If outputs differ: Explain WHY in commit message
  - [ ] If outputs same: Note "verified behaviorally equivalent"
  - [ ] Update gold files if intentional change

  ## Sign-off
  I have verified that [describe what you verified].
  Behavioral changes: [NONE | list specific intended changes]

  5. Automated Regression Testing Structure

  /test/tests/val-2a/test_regression.py:
  """
  Regression test: Ensures changes don't break existing behavior
  Run automatically by git hooks or CI
  """
  import subprocess
  import pandas as pd
  import numpy as np

  def test_val2a_output_stable():
      """Verify val-2a produces consistent output"""

      # Run simulation
      result = subprocess.run([
          '~/projects/TMAP8/tmap8-opt',
          '-i', 'val-2a.i'
      ], capture_output=True)

      assert result.returncode == 0, "Simulation failed"

      # Load outputs
      current = pd.read_csv('val-2a_out.csv')
      gold = pd.read_csv('gold/val-2a_out.csv')

      # Check flux values within 1%
      flux_diff = np.max(np.abs(
          current['scaled_recombination_flux_right'] -
          gold['scaled_recombination_flux_right']
      )) / np.max(gold['scaled_recombination_flux_right'])

      assert flux_diff < 0.01, f"Output changed by {flux_diff*100}%"

  6. Memory Reminder

  Add to /Users/bhavcv/.claude/projects/-Users-bhavcv-projects-TMAP8/memory/MEMORY.md:

  ## Critical Workflow Rule: Verify Before/After

  **NEVER modify a simulation input file without verification:**

  1. Run original → save baseline
  2. Make changes
  3. Run modified → compare outputs
  4. Only proceed if: outputs identical OR differences explained

  **Common gotchas that change behavior silently:**
  - Function types (PiecewiseLinear vs PiecewiseConstant vs ...)
  - Interpolation schemes (linear vs constant vs ...)
  - Mesh refinement (changes discretization error)
  - Solver tolerances
  - Output timesteps

  **Template response when modifying .i files:**
  > "Before marking this complete, I need to verify behavioral equivalence:
  > 1. Running baseline...
  > 2. Running modified...
  > 3. Comparing outputs...
  > Result: [identical | differs by X% due to Y]"

  7. AskUserQuestion Integration

  Modify my workflow to automatically ask when I detect model changes:

  # In CLAUDE.md

  ## Automated Verification Questions

  When you detect yourself modifying a computational model file (.i, .cfg, .inp):

  MUST use AskUserQuestion:
  - "I'm about to change X to Y in the model. Should I verify behavioral equivalence first?"
    - Options: "Yes, run baseline comparison" | "No, this is a new feature" | "Already verified manually"

  This forces me to PAUSE and consider whether my change is behavior-preserving.

  Why These Would Have Caught the Issue

  The PiecewiseLinear vs PiecewiseConstant issue would have been caught by:

  1. Pre-commit hook: Would have stopped me and asked "did you verify?"
  2. Verification script: Would have shown that flux output changed when function type changed
  3. Checklist: "Function type changes" would have been flagged as high-risk
  4. Memory reminder: Would have reminded me that interpolation schemes affect behavior
  5. AskUserQuestion: Would have forced me to pause and question the change

  The General Principle

  Treat computational model modifications like code refactoring: prove behavioral equivalence unless the change is explicitly intended to
   alter behavior.

  This is the scientific computing equivalent of unit tests - you need verification that your "refactoring" doesn't change the physics.

  The most powerful single addition would be #1 (CLAUDE.md rule) + #2 (git hook), as they create a forcing function that prevents me from
   proceeding without verification. The hook physically blocks the commit, and the CLAUDE.md rule gives me the protocol to follow.

  In this instance you didn't change the function type, you were writing a new input file from scratch. I wonder if it would be
  helpful for you to have a way to ensure the code behind a given input file block (aka the class this block actually creates an
  object of) does what you conceptually think it does.

⏺ Ah, this is a much deeper and more important issue! You're right - I wasn't changing a function type, I was writing from scratch and
  assumed PiecewiseLinear would create a step function when given step-like data points, without verifying that assumption.

  The core problem: I used a MOOSE class without verifying that its implementation matches my conceptual understanding.

  Here are general strategies to ensure I verify "does this class actually do what I think it does?":

  1. Documentation Lookup Protocol in CLAUDE.md

  ## BLOCKING REQUIREMENT: Class Documentation Verification

  Before using ANY framework class in a simulation input file:

  ### Step 1: State Your Intent
  Write what you THINK the class does:
  > "I'm using PiecewiseLinear to create a step function for the beam schedule.
  > I expect it to hold constant values between time points."

  ### Step 2: Read Documentation
  ```bash
  # For MOOSE classes:
  ./tmap8-opt --dump
  grep -A 20 "class_name"

  # Or use moose-docs skill
  /moose-docs <ClassName>

  Step 3: Verify with Minimal Example

  Create a trivial test case that demonstrates the behavior:
  [Problem]
    type = FEProblem
    solve = false
  []

  [Functions]
    [test_func]
      type = PiecewiseLinear  # Class under test
      x = '0 1 2'
      y = '0 1 1'
    []
  []

  [Postprocessors]
    [eval_at_0.5]
      type = FunctionValuePostprocessor
      function = test_func
      point = '0 0 0'
      execute_on = 'initial'
    []
  []

  [Executioner]
    type = Transient
    num_steps = 1
  []

  [Outputs]
    csv = true
  []

  Expected: value at t=0.5 is ???
  Actual: [run and check]
  Conclusion: Does behavior match expectation?

  Step 4: Document in Code Comments

  # Using PiecewiseLinear - verified it does LINEAR INTERPOLATION between points
  # NOT step function behavior. For steps, use PiecewiseConstant.

  ## **2. Class Behavior Verification Skill**

  Create `~/.claude/skills/verify-class-behavior.md`:

  ```markdown
  ---
  name: verify-class-behavior
  description: Verify that a framework class behaves as conceptually expected
  ---

  # Class Behavior Verification Protocol

  ## Usage
  Use this when: Writing simulation input with a class you haven't used before OR using a class in a new way.

  ## Protocol

  1. **Conceptual Description**
     Write down what you think the class does in plain English.

  2. **Documentation Check**
     - Read official docs
     - Look for: parameters, algorithms, edge cases, interpolation methods
     - Red flags: "linear", "interpolation", "smooth" when you want discrete/step behavior

  3. **Minimal Test Case**
     Create smallest possible input that uses ONLY this class:
     - 3-5 data points maximum
     - Simple, predictable expected output
     - Can be run in seconds

  4. **Verification Questions**
     - What happens between data points? (interpolation? constant? extrapolation?)
     - What happens outside domain? (error? clamp? extrapolate?)
     - What happens at exact data points?
     - Does it match my conceptual understanding?

  5. **Visualization (if applicable)**
     Plot or tabulate the actual function values:
     ```python
     import numpy as np
     import matplotlib.pyplot as plt

     # Load output from minimal test
     # Plot actual vs expected behavior

  Example: PiecewiseLinear

  Conceptual: "I want a step function for beam on/off"

  Documentation: "PiecewiseLinear performs LINEAR INTERPOLATION between points"

  Test:
  x = [0, 1, 1.01, 2]
  y = [1, 1, 0,    0]

  Question: What is f(1.005)?
  - Expected (if step): 0 or 1
  - Actual: 0.5 (linear interpolation!)

  Conclusion: ❌ WRONG CLASS for step function. Use PiecewiseConstant instead.

  ## **3. Automatic Minimal Test Generation**

  Add helper script `/scripts/test_function_behavior.py`:

  ```python
  #!/usr/bin/env python3
  """
  Generate and run minimal MOOSE input to verify function behavior
  Usage: python test_function_behavior.py <function_type> <x_values> <y_values>
  """

  import sys
  import subprocess
  import pandas as pd

  def generate_test_input(func_type, x_vals, y_vals, test_points):
      """Generate minimal MOOSE input to test a function"""

      input_template = f"""
  [Mesh]
    type = GeneratedMesh
    dim = 1
    nx = 1
  []

  [Variables]
    [dummy]
    []
  []

  [Problem]
    solve = false
  []

  [Functions]
    [test_func]
      type = {func_type}
      x = '{" ".join(map(str, x_vals))}'
      y = '{" ".join(map(str, y_vals))}'
    []
  []

  [Postprocessors]
  """

      # Create postprocessors to evaluate at test points
      for i, t in enumerate(test_points):
          input_template += f"""
    [eval_t_{t}]
      type = FunctionValuePostprocessor
      function = test_func
      execute_on = 'initial'
      indirect_value = {t}
    []
  """

      input_template += """
  [Executioner]
    type = Steady
  []

  [Outputs]
    csv = true
    execute_on = 'initial'
  []
  """

      return input_template

  def run_test(func_type, x_vals, y_vals):
      """Run test and show function values at key points"""

      # Test at: data points, midpoints, and boundaries
      test_points = []
      test_points.extend(x_vals)  # At data points
      for i in range(len(x_vals)-1):
          test_points.append((x_vals[i] + x_vals[i+1]) / 2)  # Midpoints

      input_str = generate_test_input(func_type, x_vals, y_vals, test_points)

      # Write and run
      with open('test_function.i', 'w') as f:
          f.write(input_str)

      result = subprocess.run(
          ['~/projects/TMAP8/tmap8-opt', '-i', 'test_function.i'],
          capture_output=True, text=True
      )

      if result.returncode != 0:
          print(f"❌ Test failed:\n{result.stderr}")
          return

      # Parse output
      df = pd.read_csv('test_function_out.csv')

      print(f"\n{'='*60}")
      print(f"Testing: {func_type}")
      print(f"Data points: x={x_vals}, y={y_vals}")
      print(f"{'='*60}")
      print(f"\n{'Time':>10} | {'Value':>12} | {'Note'}")
      print(f"{'-'*10}|{'-'*14}|{'-'*30}")

      for col in df.columns:
          if col.startswith('eval_t_'):
              t = float(col.replace('eval_t_', ''))
              val = df[col].iloc[0]

              # Check if this is at a data point or between
              if t in x_vals:
                  note = "← data point"
              else:
                  note = "← between points"

              print(f"{t:>10.4f} | {val:>12.6f} | {note}")

      print(f"\n{'='*60}")

  if __name__ == '__main__':
      # Example usage
      print("\n🔍 Testing PiecewiseLinear with step-like data:")
      run_test('PiecewiseLinear',
               x_vals=[0, 5820, 5821, 9056],
               y_vals=[1, 1, 0, 0])

      print("\n🔍 Testing PiecewiseConstant with same data:")
      run_test('PiecewiseConstant',
               x_vals=[0, 5820, 9056],
               y_vals=[1, 1, 0])

  4. Visual Function Inspector Tool

  /scripts/plot_function_object.py:

  #!/usr/bin/env python3
  """
  Visualize what a MOOSE Function object actually produces
  Helps catch misunderstandings before using in production
  """

  import matplotlib.pyplot as plt
  import numpy as np
  import subprocess
  import pandas as pd

  def visualize_moose_function(func_type, x_data, y_data, title=""):
      """Plot what MOOSE function actually evaluates to"""

      # Generate dense evaluation points
      x_min, x_max = min(x_data), max(x_data)
      x_eval = np.linspace(x_min, x_max, 1000)

      # Create MOOSE input that evaluates function at these points
      # [implementation details...]

      fig, ax = plt.subplots(figsize=(10, 6))

      # Plot actual behavior
      ax.plot(x_eval, y_eval, 'b-', linewidth=2, label=f'{func_type} (actual)')

      # Plot data points
      ax.plot(x_data, y_data, 'ro', markersize=8, label='Data points')

      # Plot what you might expect (step function)
      y_step = np.zeros_like(x_eval)
      for i in range(len(x_data)-1):
          mask = (x_eval >= x_data[i]) & (x_eval < x_data[i+1])
          y_step[mask] = y_data[i]
      ax.plot(x_eval, y_step, 'g--', alpha=0.5, linewidth=2, label='Expected (step)')

      ax.set_xlabel('x (time)')
      ax.set_ylabel('y (value)')
      ax.set_title(f'Function Behavior: {func_type}\n{title}')
      ax.legend()
      ax.grid(True, alpha=0.3)

      plt.tight_layout()
      plt.savefig(f'function_verification_{func_type}.png', dpi=150)
      print(f"📊 Saved: function_verification_{func_type}.png")

      # Analysis
      if not np.allclose(y_eval, y_step, rtol=0.01):
          print("\n⚠️   WARNING: Actual behavior differs from step function!")
          print(f"   Max deviation: {np.max(np.abs(y_eval - y_step)):.3f}")
          print(f"   This is LINEAR INTERPOLATION, not step/constant.")
          return False
      else:
          print("\n✅ Behavior matches step function")
          return True

  5. Add to CLAUDE.md: Common Class Gotchas

  ## Common MOOSE Class Misconceptions

  ### Functions
  | Class Name | What You Might Think | What It Actually Does |
  |------------|---------------------|---------------------|
  | `PiecewiseLinear` | Step function if points close together | LINEAR interpolation always |
  | `PiecewiseConstant` | Linear interpolation | Step function (right-continuous) |
  | `ParsedFunction` | Simple math | Full muParser, can be slow in loops |
  | `Smooth` functions | Small smoothing | Can significantly alter peaks/transients |

  ### Time Integration
  | Class | What You Might Think | What It Actually Does |
  |-------|---------------------|---------------------|
  | `Implicit Euler` | Just "implicit" | First-order accurate, can be diffusive |
  | `BDF2` | "Better Euler" | Second-order, but can overshoot |

  ### ALWAYS verify with minimal test case first!

  6. Checklist for New Class Usage

  Add to task workflow:

  When using a framework class for the first time:

  - [ ] Read documentation for this specific class
  - [ ] Note: interpolation method, boundary behavior, algorithmic details
  - [ ] Create 5-line test case with 3 data points
  - [ ] Verify behavior at: data points, midpoints, boundaries
  - [ ] If behavior != expectation: find correct class or adjust approach
  - [ ] Add comment in code: "Using X because [behavior] verified"

  7. Integration with Existing Tools

  Extend the moose-docs skill:

  # When moose-docs is called, also:

  1. Show common gotchas for this class type
  2. Generate a minimal test template
  3. Suggest: "Would you like me to verify behavior with a test case first?"

  Why This Would Have Caught the Issue

  If I had followed this protocol:

  1. State intent: "I want a step function for beam on/off schedule"
  2. Check docs: Would have seen "PiecewiseLinear: linear interpolation between points"
  3. Minimal test:
  x = [0, 1, 1.01, 2]
  y = [1, 1, 0, 0]
  3. Evaluate at t=1.005 → would get 0.5, not 1.0
  4. Visualization: Plot would show linear ramp, not step
  5. Conclusion: "This is wrong class! Need PiecewiseConstant or smaller dt transition"

  The General Principle

  Never use a framework class based on name alone. Always verify that the implementation matches your conceptual model through minimal
  testing.

  This is like the programming principle "don't assume API behavior - test it", but applied to configuration-driven frameworks where the
  "API" is a class name in an input file.

  Most impactful single addition: The visualization tool (#4) combined with the protocol requirement (#1), because seeing a plot of "what
   you get vs what you expected" immediately reveals conceptual mismatches.
