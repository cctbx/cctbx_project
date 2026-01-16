"""
Dry Run Manager for PHENIX AI Agent Testing.

Provides simulated program execution using pre-recorded logs and outputs,
allowing full workflow testing without running actual PHENIX programs.
"""

from __future__ import absolute_import, division, print_function
import os
import shutil
import yaml


class DryRunManager:
  """
  Manages dry run test scenarios.

  Provides dummy log files and output files for simulated program execution,
  enabling full agent workflow testing without running real PHENIX programs.

  Usage:
      manager = DryRunManager("xray_basic")
      log_text, error_text, output_files = manager.get_dummy_result(
          "phenix.xtriage data.mtz",
          working_dir="/path/to/workdir"
      )
  """

  def __init__(self, scenario_name, scenarios_dir=None):
    """
    Initialize the dry run manager.

    Args:
        scenario_name: Name of the test scenario (folder name)
        scenarios_dir: Path to scenarios directory (auto-detected if None)
    """
    self.scenario_name = scenario_name
    self.scenarios_dir = scenarios_dir or self._find_scenarios_dir()
    self.scenario_dir = os.path.join(self.scenarios_dir, scenario_name)

    if not os.path.isdir(self.scenario_dir):
      raise ValueError(
        "Scenario '%s' not found in %s" % (scenario_name, self.scenarios_dir)
      )

    self.scenario = self._load_scenario()
    self.program_counts = {}  # Track invocation count per program
    self.step_index = 0       # Current position in steps list

  def _find_scenarios_dir(self):
    """Find the test scenarios directory."""
    # Look relative to this file
    this_dir = os.path.dirname(os.path.abspath(__file__))

    # Try tests/scenarios relative to parent
    parent_dir = os.path.dirname(this_dir)
    scenarios_dir = os.path.join(parent_dir, "tests", "scenarios")

    if os.path.isdir(scenarios_dir):
      return scenarios_dir

    # Try relative to cwd
    if os.path.isdir("tests/scenarios"):
      return os.path.abspath("tests/scenarios")

    raise ValueError(
      "Could not find test scenarios directory. "
      "Expected at: %s" % scenarios_dir
    )

  def _load_scenario(self):
    """Load scenario configuration from YAML."""
    yaml_path = os.path.join(self.scenario_dir, "scenario.yaml")

    if not os.path.isfile(yaml_path):
      raise ValueError(
        "Scenario file not found: %s" % yaml_path
      )

    with open(yaml_path, 'r') as f:
      return yaml.safe_load(f)

  def get_initial_files(self):
    """
    Get list of initial input files for this scenario.

    Returns:
        list: Absolute paths to input files
    """
    files = []
    for rel_path in self.scenario.get('initial_files', []):
      abs_path = os.path.join(self.scenario_dir, rel_path)
      if os.path.exists(abs_path):
        files.append(abs_path)
      else:
        print("Warning: Initial file not found: %s" % abs_path)
    return files

  def get_dummy_result(self, command, working_dir):
    """
    Get dummy log and output files for a command.

    Args:
        command: The command string (e.g., "phenix.xtriage data.mtz")
        working_dir: Directory where output files should be copied

    Returns:
        tuple: (log_text, error_text, output_files)

    Raises:
        ValueError: If no dummy data exists for this program call
    """
    # Extract program name from command
    program = command.split()[0]  # e.g., "phenix.xtriage"

    # Track how many times this program has been called
    count = self.program_counts.get(program, 0) + 1
    self.program_counts[program] = count

    # Find matching step in scenario
    step = self._find_step(program, count)

    if step is None:
      # No more steps for this program - might be at end of workflow
      available = self._get_available_steps()
      raise ValueError(
        "No dummy data for %s (call #%d). "
        "Available steps: %s" % (program, count, available)
      )

    # Read log file
    log_path = os.path.join(self.scenario_dir, step['log'])
    if not os.path.isfile(log_path):
      raise ValueError("Log file not found: %s" % log_path)

    with open(log_path, 'r') as f:
      log_text = f.read()

    # Get error text if specified
    error_text = ""
    if 'error' in step:
      error_path = os.path.join(self.scenario_dir, step['error'])
      if os.path.isfile(error_path):
        with open(error_path, 'r') as f:
          error_text = f.read()

    # Copy output files to working directory
    output_files = []
    for output_rel in step.get('outputs', []):
      src = os.path.join(self.scenario_dir, output_rel)
      if os.path.isfile(src):
        dst = os.path.join(working_dir, os.path.basename(src))
        shutil.copy(src, dst)
        output_files.append(dst)
      else:
        print("Warning: Output file not found: %s" % src)

    self.step_index += 1
    return log_text, error_text, output_files

  def _find_step(self, program, count):
    """
    Find the step definition for a program invocation.

    Args:
        program: Program name (e.g., "phenix.xtriage")
        count: Which invocation of this program (1, 2, 3, ...)

    Returns:
        dict: Step definition or None if not found
    """
    steps = self.scenario.get('steps', [])

    # Count how many times we've seen this program in steps
    seen = 0
    for step in steps:
      if step.get('program') == program:
        seen += 1
        if seen == count:
          return step

    return None

  def _get_available_steps(self):
    """Get summary of available steps for error messages."""
    steps = self.scenario.get('steps', [])
    programs = [s.get('program', 'unknown') for s in steps]
    return programs

  def get_expected_outcome(self):
    """
    Get expected outcome for validation.

    Returns:
        dict: Expected outcome configuration
    """
    return self.scenario.get('expected', {})

  def reset(self):
    """Reset state for a new run."""
    self.program_counts = {}
    self.step_index = 0


def list_scenarios(scenarios_dir=None):
  """
  List available test scenarios.

  Args:
      scenarios_dir: Path to scenarios directory (auto-detected if None)

  Returns:
      list: Names of available scenarios
  """
  if scenarios_dir is None:
    # Find scenarios dir
    this_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(this_dir)
    scenarios_dir = os.path.join(parent_dir, "tests", "scenarios")

  if not os.path.isdir(scenarios_dir):
    return []

  scenarios = []
  for name in os.listdir(scenarios_dir):
    scenario_path = os.path.join(scenarios_dir, name)
    yaml_path = os.path.join(scenario_path, "scenario.yaml")
    if os.path.isdir(scenario_path) and os.path.isfile(yaml_path):
      scenarios.append(name)

  return sorted(scenarios)


def get_scenario_info(scenario_name, scenarios_dir=None):
  """
  Get information about a scenario.

  Args:
      scenario_name: Name of the scenario
      scenarios_dir: Path to scenarios directory

  Returns:
      dict: Scenario metadata
  """
  manager = DryRunManager(scenario_name, scenarios_dir)
  scenario = manager.scenario

  return {
    'name': scenario.get('name', scenario_name),
    'description': scenario.get('description', ''),
    'experiment_type': scenario.get('experiment_type', 'unknown'),
    'num_steps': len(scenario.get('steps', [])),
    'programs': [s.get('program') for s in scenario.get('steps', [])],
    'expected': scenario.get('expected', {}),
  }
