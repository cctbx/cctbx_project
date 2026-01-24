#!/usr/bin/env python
"""
Program Validator - Validates program configuration completeness.

This tool checks that a program is properly configured across all necessary
files in the PHENIX AI Agent codebase.

Usage:
    python agent/program_validator.py phenix.map_symmetry
    python agent/program_validator.py --all
    python agent/program_validator.py --list
"""

import os
import re
import sys
import yaml


def get_knowledge_dir():
  """Find the knowledge directory."""
  # Try relative to this file
  this_dir = os.path.dirname(os.path.abspath(__file__))
  parent_dir = os.path.dirname(this_dir)
  knowledge_dir = os.path.join(parent_dir, "knowledge")
  if os.path.isdir(knowledge_dir):
    return knowledge_dir

  # Try relative to cwd
  if os.path.isdir("knowledge"):
    return "knowledge"

  return None


def get_code_dirs():
  """Find code directories."""
  this_dir = os.path.dirname(os.path.abspath(__file__))
  parent_dir = os.path.dirname(this_dir)

  dirs = {
    "agent": os.path.join(parent_dir, "agent"),
    "phenix_ai": os.path.join(parent_dir, "phenix_ai"),
  }

  return {k: v for k, v in dirs.items() if os.path.isdir(v)}


def load_yaml_file(filepath):
  """Load a YAML file."""
  with open(filepath, 'r') as f:
    return yaml.safe_load(f)


def check_programs_yaml(program_name, knowledge_dir):
  """Check if program is defined in programs.yaml."""
  issues = []
  info = {}

  filepath = os.path.join(knowledge_dir, "programs.yaml")
  if not os.path.exists(filepath):
    return ["programs.yaml not found"], info

  programs = load_yaml_file(filepath)
  if not programs:
    return ["programs.yaml is empty"], info

  if program_name not in programs:
    issues.append(f"Program '{program_name}' not defined in programs.yaml")
    return issues, info

  prog_def = programs[program_name]
  info["definition"] = prog_def

  # Check required fields
  required_fields = ["description", "category", "inputs", "command"]
  for field in required_fields:
    if field not in prog_def:
      issues.append(f"Missing required field '{field}' in programs.yaml")

  # Check experiment_types
  if "experiment_types" not in prog_def:
    issues.append("Missing 'experiment_types' - should be [xray], [cryoem], or both")

  # Check if outputs are defined
  if "outputs" not in prog_def:
    issues.append("No 'outputs' section - should define expected output files/metrics")
  else:
    if "metrics" in prog_def["outputs"]:
      info["metrics"] = prog_def["outputs"]["metrics"]

  # Check log_parsing
  if "log_parsing" not in prog_def:
    issues.append("No 'log_parsing' section - metrics won't be extracted from YAML patterns")

  return issues, info


def check_workflows_yaml(program_name, knowledge_dir):
  """Check if program is included in workflows.yaml."""
  issues = []
  info = {"phases": []}

  filepath = os.path.join(knowledge_dir, "workflows.yaml")
  if not os.path.exists(filepath):
    return ["workflows.yaml not found"], info

  workflows = load_yaml_file(filepath)
  if not workflows:
    return ["workflows.yaml is empty"], info

  # Search for program in all workflows and phases
  found = False
  for workflow_name, workflow in workflows.items():
    if not isinstance(workflow, dict):
      continue
    phases = workflow.get("phases", {})
    for phase_name, phase in phases.items():
      if not isinstance(phase, dict):
        continue
      programs = phase.get("programs", [])
      for prog in programs:
        prog_name = prog if isinstance(prog, str) else prog.get("program", "")
        if prog_name == program_name:
          found = True
          info["phases"].append(f"{workflow_name}/{phase_name}")

  if not found:
    issues.append(f"Program '{program_name}' not found in any workflow phase in workflows.yaml")

  return issues, info


def check_log_parsers(program_name, code_dirs):
  """Check if program has log parsing support."""
  issues = []
  info = {}

  phenix_ai_dir = code_dirs.get("phenix_ai")
  if not phenix_ai_dir:
    return ["phenix_ai directory not found"], info

  filepath = os.path.join(phenix_ai_dir, "log_parsers.py")
  if not os.path.exists(filepath):
    return ["log_parsers.py not found"], info

  with open(filepath, 'r') as f:
    content = f.read()

  # Get short name (e.g., "map_symmetry" from "phenix.map_symmetry")
  short_name = program_name.replace("phenix.", "").replace(".", "_")

  # Check for extraction function
  func_pattern = f"def _extract_{short_name}_metrics"
  if func_pattern not in content:
    issues.append(f"No extraction function '_extract_{short_name}_metrics()' in log_parsers.py")
  else:
    info["has_extraction_function"] = True

  # Check if it's called in extract_all_metrics
  call_pattern = f"_extract_{short_name}_metrics"
  in_extract_all = re.search(
    rf'def extract_all_metrics.*?(?=\ndef [a-z])',
    content,
    re.DOTALL
  )
  if in_extract_all and call_pattern not in in_extract_all.group():
    issues.append(f"Extraction function not called in extract_all_metrics()")

  # Check detect_program
  if short_name not in content.split("def detect_program")[1].split("def ")[0] if "def detect_program" in content else "":
    # More lenient check
    detect_section = re.search(r'def detect_program.*?(?=\ndef )', content, re.DOTALL)
    if detect_section and short_name.replace("_", "") not in detect_section.group().lower().replace("_", ""):
      issues.append(f"Program may not be detected in detect_program()")

  return issues, info


def check_workflow_state(program_name, code_dirs):
  """Check if program has history tracking."""
  issues = []
  info = {}

  agent_dir = code_dirs.get("agent")
  if not agent_dir:
    return ["agent directory not found"], info

  filepath = os.path.join(agent_dir, "workflow_state.py")
  if not os.path.exists(filepath):
    return ["workflow_state.py not found"], info

  with open(filepath, 'r') as f:
    content = f.read()

  short_name = program_name.replace("phenix.", "").replace(".", "_")
  done_flag = f"{short_name}_done"

  # Check for done flag in info dict
  if done_flag not in content:
    issues.append(f"No '{done_flag}' tracking flag in workflow_state.py")
  else:
    info["has_tracking"] = True

  return issues, info


def check_workflow_engine(program_name, code_dirs):
  """Check if program is in workflow engine context (if needed)."""
  issues = []
  info = {}

  agent_dir = code_dirs.get("agent")
  if not agent_dir:
    return [], info  # Not critical

  filepath = os.path.join(agent_dir, "workflow_engine.py")
  if not os.path.exists(filepath):
    return [], info

  with open(filepath, 'r') as f:
    content = f.read()

  short_name = program_name.replace("phenix.", "").replace(".", "_")
  done_flag = f"{short_name}_done"

  # Only warn if tracked in workflow_state but not in engine
  if done_flag in open(os.path.join(agent_dir, "workflow_state.py")).read():
    if done_flag not in content:
      issues.append(f"'{done_flag}' tracked in workflow_state but not in workflow_engine context (may be fine)")
      info["optional_warning"] = True

  return issues, info


def check_session_summary(program_name, code_dirs):
  """Check if program metrics appear in session summary."""
  issues = []
  info = {}

  agent_dir = code_dirs.get("agent")
  if not agent_dir:
    return ["agent directory not found"], info

  filepath = os.path.join(agent_dir, "session.py")
  if not os.path.exists(filepath):
    return ["session.py not found"], info

  with open(filepath, 'r') as f:
    content = f.read()

  short_name = program_name.replace("phenix.", "").replace(".", "_")

  # Check _get_step_metric
  if short_name not in content:
    issues.append(f"Program '{short_name}' not referenced in session.py - metrics may not display in summary")

  return issues, info


def check_directive_extractor(program_name, code_dirs):
  """Check if program has tutorial/directive patterns."""
  issues = []
  info = {}

  agent_dir = code_dirs.get("agent")
  if not agent_dir:
    return [], info

  filepath = os.path.join(agent_dir, "directive_extractor.py")
  if not os.path.exists(filepath):
    return [], info

  with open(filepath, 'r') as f:
    content = f.read()

  # This is optional - only note if missing
  if program_name not in content:
    short_name = program_name.replace("phenix.", "")
    if short_name not in content:
      info["note"] = "No tutorial patterns - users can't trigger this program by name"

  return issues, info


def validate_program(program_name):
  """Validate that a program is fully configured."""
  print(f"\n{'='*60}")
  print(f"Validating: {program_name}")
  print('='*60)

  knowledge_dir = get_knowledge_dir()
  code_dirs = get_code_dirs()

  if not knowledge_dir:
    print("ERROR: Could not find knowledge directory")
    return False

  all_issues = []
  all_info = {}

  # Check each component
  checks = [
    ("programs.yaml", check_programs_yaml, (program_name, knowledge_dir)),
    ("workflows.yaml", check_workflows_yaml, (program_name, knowledge_dir)),
    ("log_parsers.py", check_log_parsers, (program_name, code_dirs)),
    ("workflow_state.py", check_workflow_state, (program_name, code_dirs)),
    ("workflow_engine.py", check_workflow_engine, (program_name, code_dirs)),
    ("session.py", check_session_summary, (program_name, code_dirs)),
    ("directive_extractor.py", check_directive_extractor, (program_name, code_dirs)),
  ]

  for name, check_func, args in checks:
    issues, info = check_func(*args)
    all_issues.extend(issues)
    all_info[name] = info

    status = "✓" if not issues else "✗"
    print(f"\n{status} {name}")

    if issues:
      for issue in issues:
        print(f"  - {issue}")
    elif info:
      for key, value in info.items():
        if key == "note":
          print(f"  Note: {value}")
        elif key == "phases":
          print(f"  Found in: {', '.join(value)}")

  # Summary
  print(f"\n{'-'*60}")
  if all_issues:
    print(f"RESULT: {len(all_issues)} issue(s) found")
    return False
  else:
    print("RESULT: Program is fully configured ✓")
    return True


def list_programs():
  """List all defined programs."""
  knowledge_dir = get_knowledge_dir()
  if not knowledge_dir:
    print("ERROR: Could not find knowledge directory")
    return

  filepath = os.path.join(knowledge_dir, "programs.yaml")
  programs = load_yaml_file(filepath)

  print("\nDefined Programs:")
  print("-" * 40)
  for name, defn in sorted(programs.items()):
    category = defn.get("category", "unknown")
    exp_types = defn.get("experiment_types", [])
    print(f"  {name}")
    print(f"    category: {category}, types: {exp_types}")


def validate_all():
  """Validate all programs."""
  knowledge_dir = get_knowledge_dir()
  if not knowledge_dir:
    print("ERROR: Could not find knowledge directory")
    return

  filepath = os.path.join(knowledge_dir, "programs.yaml")
  programs = load_yaml_file(filepath)

  results = {}
  for program_name in sorted(programs.keys()):
    results[program_name] = validate_program(program_name)

  # Final summary
  print("\n" + "=" * 60)
  print("SUMMARY")
  print("=" * 60)
  passed = sum(1 for v in results.values() if v)
  total = len(results)
  print(f"Passed: {passed}/{total}")

  if passed < total:
    print("\nPrograms with issues:")
    for name, passed in results.items():
      if not passed:
        print(f"  - {name}")


def main():
  if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(1)

  arg = sys.argv[1]

  if arg == "--list":
    list_programs()
  elif arg == "--all":
    validate_all()
  elif arg.startswith("--"):
    print(f"Unknown option: {arg}")
    print(__doc__)
    sys.exit(1)
  else:
    # Assume it's a program name
    program_name = arg
    if not program_name.startswith("phenix."):
      program_name = f"phenix.{program_name}"
    success = validate_program(program_name)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
  main()
