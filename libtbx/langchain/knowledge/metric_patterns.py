"""
Centralized metric pattern loading from programs.yaml.

This module provides a single source of truth for metric extraction patterns,
eliminating duplication between log_parsers.py and session.py.

Usage:
    from libtbx.langchain.knowledge.metric_patterns import (
        extract_metrics_for_program,
        get_all_metric_patterns,
        get_metric_display_config,
    )

    # Extract metrics from a log using YAML-defined patterns
    metrics = extract_metrics_for_program(log_text, "phenix.map_symmetry")

    # Get display configuration for formatting
    config = get_metric_display_config("phenix.map_symmetry", "symmetry_type")
"""

import re
from functools import lru_cache

# Silence unused import warning - lru_cache is used as decorator
assert lru_cache is not None


@lru_cache(maxsize=1)
def get_all_metric_patterns():
  """
  Load all metric patterns from programs.yaml.

  Returns:
      dict: Program name -> metric name -> pattern configuration

  Pattern configuration includes:
      - pattern: Compiled regex pattern
      - type: 'float', 'int', 'string', or 'bool'/'boolean'
      - display_name: Human-readable name
      - summary_format: Format string for summary display
      - no_match_pattern: Pattern to detect "no result" cases
      - no_match_value: Value to use when no_match_pattern matches
      - fallback_pattern: Alternative pattern to try
      - extract: 'first' or 'last' for multiple matches
  """
  from libtbx.langchain.knowledge.yaml_loader import load_programs

  programs = load_programs()
  if not programs:
    return {}

  all_patterns = {}

  for prog_name, prog_def in programs.items():
    if not isinstance(prog_def, dict):
      continue

    log_parsing = prog_def.get('log_parsing', {})
    if not log_parsing:
      continue

    prog_patterns = {}

    for metric_name, metric_def in log_parsing.items():
      if not isinstance(metric_def, dict):
        continue

      pattern_str = metric_def.get('pattern')
      if not pattern_str:
        continue

      try:
        # Compile primary pattern
        compiled_pattern = re.compile(pattern_str, re.IGNORECASE)

        config = {
          'pattern': compiled_pattern,
          'pattern_str': pattern_str,  # Keep original for debugging
          'type': metric_def.get('type', 'string'),
          'display_name': metric_def.get('display_name', _format_metric_name(metric_name)),
          'summary_format': metric_def.get('summary_format', '{value}'),
          'extract': metric_def.get('extract', 'first'),
        }

        # Optional fields
        if 'no_match_pattern' in metric_def:
          config['no_match_pattern'] = re.compile(
            metric_def['no_match_pattern'], re.IGNORECASE
          )
          config['no_match_value'] = metric_def.get('no_match_value')

        if 'fallback_pattern' in metric_def:
          config['fallback_pattern'] = re.compile(
            metric_def['fallback_pattern'], re.IGNORECASE
          )

        prog_patterns[metric_name] = config

      except re.error as e:
        # Log error but continue - don't fail on one bad pattern
        print(f"Warning: Invalid regex pattern for {prog_name}.{metric_name}: {e}")
        continue

    if prog_patterns:
      all_patterns[prog_name] = prog_patterns

  return all_patterns


def _format_metric_name(name):
  """Convert metric_name to 'Metric Name'."""
  return name.replace('_', ' ').title()


def _convert_value(value_str, value_type):
  """Convert extracted string value to appropriate type."""
  if value_str is None:
    return None

  value_str = value_str.strip()

  if value_type == 'float':
    try:
      return float(value_str)
    except ValueError:
      return None
  elif value_type == 'int':
    try:
      return int(value_str)
    except ValueError:
      return None
  elif value_type in ('bool', 'boolean'):
    # For boolean patterns, match existing means True
    return True
  else:
    return value_str


def extract_metrics_for_program(log_text, program_name):
  """
  Extract metrics for a specific program using YAML patterns.

  Args:
      log_text: Log file content
      program_name: Full program name (e.g., 'phenix.map_symmetry')

  Returns:
      dict: metric_name -> extracted value
  """
  if not log_text:
    return {}

  all_patterns = get_all_metric_patterns()
  prog_patterns = all_patterns.get(program_name, {})

  if not prog_patterns:
    return {}

  metrics = {}

  for metric_name, config in prog_patterns.items():
    value = None

    # Try primary pattern
    pattern = config['pattern']
    matches = list(pattern.finditer(log_text))

    if matches:
      # Handle extract: first vs last
      if config.get('extract') == 'last':
        match = matches[-1]
      else:
        match = matches[0]

      # Get first non-None group (patterns may have alternation)
      for group in match.groups():
        if group is not None:
          value = _convert_value(group, config['type'])
          break

    # Try fallback pattern if primary didn't match
    if value is None and 'fallback_pattern' in config:
      fallback_match = config['fallback_pattern'].search(log_text)
      if fallback_match:
        for group in fallback_match.groups():
          if group is not None:
            value = _convert_value(group, config['type'])
            break

    # Handle no-match pattern (e.g., "No symmetry found")
    if value is None and 'no_match_pattern' in config:
      if config['no_match_pattern'].search(log_text):
        value = config.get('no_match_value')

    # For boolean type, presence of match means True
    if config['type'] in ('bool', 'boolean'):
      if matches:
        value = True
      # If no_match_pattern matched, value was already set

    if value is not None:
      metrics[metric_name] = value

  return metrics


def get_metric_display_config(program_name, metric_name):
  """
  Get display configuration for a metric.

  Args:
      program_name: Full program name
      metric_name: Metric name

  Returns:
      dict with display_name, summary_format, etc. or empty dict
  """
  all_patterns = get_all_metric_patterns()
  prog_patterns = all_patterns.get(program_name, {})
  return prog_patterns.get(metric_name, {})


def get_program_metrics(program_name):
  """
  Get list of metrics defined for a program.

  Args:
      program_name: Full program name

  Returns:
      list of metric names
  """
  all_patterns = get_all_metric_patterns()
  prog_patterns = all_patterns.get(program_name, {})
  return list(prog_patterns.keys())


def format_metric_value(program_name, metric_name, value):
  """
  Format a metric value using YAML-defined format.

  Args:
      program_name: Full program name
      metric_name: Metric name
      value: The metric value

  Returns:
      Formatted string
  """
  config = get_metric_display_config(program_name, metric_name)
  if not config:
    return str(value)

  summary_format = config.get('summary_format', '{value}')
  try:
    return summary_format.format(value=value)
  except (ValueError, KeyError):
    return str(value)


def extract_metrics_generic(log_text):
  """
  Extract metrics from log text without knowing the program.

  Tries all program patterns and returns combined results.
  Useful when program is unknown or for generic metric detection.

  Args:
      log_text: Log file content

  Returns:
      dict: metric_name -> extracted value (from first matching program)
  """
  all_patterns = get_all_metric_patterns()

  # Try to detect program first
  program = None
  log_lower = log_text[:2000].lower() if log_text else ""

  # Common program signatures
  program_signatures = [
    ('phenix.xtriage', ['xtriage', 'relative wilson']),
    ('phenix.mtriage', ['mtriage', 'fsc_model']),
    ('phenix.refine', ['phenix.refine', 'r_work', 'r_free']),
    ('phenix.real_space_refine', ['real_space_refine', 'cc_mask']),
    ('phenix.map_symmetry', ['map_symmetry', 'ncs type']),
    ('phenix.molprobity', ['molprobity', 'clashscore']),
    ('phenix.phaser', ['phaser', 'tfz=', 'llg=']),
  ]

  for prog_name, signatures in program_signatures:
    if any(sig in log_lower for sig in signatures):
      program = prog_name
      break

  if program:
    return extract_metrics_for_program(log_text, program)

  # Fallback: try all patterns and return combined results
  combined = {}
  for prog_name in all_patterns:
    metrics = extract_metrics_for_program(log_text, prog_name)
    # Don't overwrite existing metrics
    for k, v in metrics.items():
      if k not in combined:
        combined[k] = v

  return combined


# =============================================================================
# Testing / CLI
# =============================================================================

def _test_patterns():
  """Test pattern extraction."""
  print("Testing metric pattern extraction...")

  # Test map_symmetry
  map_symmetry_log = """
  Best NCS type is:
    SCORE    CC   OPERATORS     SYMMETRY
     3.40   0.93    14          D7 (a)  Best NCS type

  Final symmetry obtained:
  NCS type: D7 (a)
  Correlation of symmetry-related regions: 0.93   Copies: 14
  """

  metrics = extract_metrics_for_program(map_symmetry_log, "phenix.map_symmetry")
  print(f"map_symmetry metrics: {metrics}")
  assert 'symmetry_type' in metrics, "Should extract symmetry_type"
  assert metrics.get('symmetry_type') == 'D7 (a)', f"Expected 'D7 (a)', got {metrics.get('symmetry_type')}"

  # Test no symmetry case
  no_sym_log = "No suitable symmetry found"
  metrics = extract_metrics_for_program(no_sym_log, "phenix.map_symmetry")
  print(f"no symmetry metrics: {metrics}")
  assert metrics.get('symmetry_type') == 'None', "Should extract 'None' for no symmetry"

  # Test refine
  refine_log = """
  Final R-work = 0.2150
  Final R-free = 0.2580
  """
  metrics = extract_metrics_for_program(refine_log, "phenix.refine")
  print(f"refine metrics: {metrics}")
  assert 'r_free' in metrics, "Should extract r_free"
  assert abs(metrics.get('r_free', 0) - 0.258) < 0.001, "R-free should be 0.258"

  print("All pattern tests passed!")


if __name__ == "__main__":
  import sys

  if len(sys.argv) > 1 and sys.argv[1] == "--test":
    _test_patterns()
  else:
    # Show all patterns
    patterns = get_all_metric_patterns()
    print(f"Loaded patterns for {len(patterns)} programs:")
    for prog_name, prog_patterns in sorted(patterns.items()):
      print(f"\n{prog_name}:")
      for metric_name, config in prog_patterns.items():
        print(f"  {metric_name}: {config['display_name']} ({config['type']})")
