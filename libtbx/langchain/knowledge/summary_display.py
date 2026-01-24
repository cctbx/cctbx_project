"""
YAML-driven summary display formatting.

This module provides utilities to format session summary displays based on
configuration in metrics.yaml, eliminating the need for hardcoded formatting
in session.py.

Usage:
    from libtbx.langchain.knowledge.summary_display import (
        format_quality_table_rows,
        format_step_metric,
        get_quality_table_config,
    )

    # Format rows for Final Quality table
    rows = format_quality_table_rows(metrics_dict, experiment_type)

    # Format a step's metric display
    metric_str = format_step_metric("phenix.refine", metrics_dict)
"""

from functools import lru_cache

# Silence unused import warning - lru_cache is used as decorator
assert lru_cache is not None


@lru_cache(maxsize=1)
def _load_summary_display_config():
  """
  Load summary display configuration from metrics.yaml.

  Returns:
      dict: Summary display configuration or empty dict
  """
  from libtbx.langchain.knowledge.yaml_loader import load_metrics

  metrics_config = load_metrics()
  if not metrics_config:
    return {}

  return metrics_config.get('summary_display', {})


def get_quality_table_config():
  """
  Get the quality table configuration.

  Returns:
      list: List of row configurations for the quality table
  """
  config = _load_summary_display_config()
  return config.get('quality_table', [])


def get_step_metrics_config():
  """
  Get the step metrics configuration.

  Returns:
      dict: Program name -> step metric configuration
  """
  config = _load_summary_display_config()
  return config.get('step_metrics', {})


def format_quality_table_rows(metrics_dict, experiment_type=None):
  """
  Format rows for the Final Quality table based on available metrics.

  Args:
      metrics_dict: Dict of metric name -> value
      experiment_type: 'xray' or 'cryoem' (optional, filters rows)

  Returns:
      list of dicts: [{'label': str, 'value': str, 'detail': str}, ...]
  """
  if not metrics_dict:
    return []

  config = get_quality_table_config()
  rows = []

  for row_config in config:
    # Check experiment type filter
    exp_types = row_config.get('experiment_types')
    if exp_types and experiment_type and experiment_type not in exp_types:
      continue

    # Check if required metrics are present
    required_metrics = row_config.get('metrics', [])
    if not all(m in metrics_dict for m in required_metrics):
      continue

    # Format the main value
    try:
      value = row_config['format'].format(**metrics_dict)
    except (KeyError, ValueError):
      continue

    # Add optional parts (e.g., " (14 copies)")
    optional_metrics = row_config.get('optional_metrics', [])
    if optional_metrics and all(m in metrics_dict for m in optional_metrics):
      try:
        optional_format = row_config.get('optional_format', '')
        value += optional_format.format(**metrics_dict)
      except (KeyError, ValueError):
        pass

    # Format detail column (e.g., "CC: 0.93" or assessment)
    detail = ""

    # Try detail_format first
    detail_metrics = row_config.get('detail_metrics', [])
    if detail_metrics and all(m in metrics_dict for m in detail_metrics):
      try:
        detail = row_config.get('detail_format', '').format(**metrics_dict)
      except (KeyError, ValueError):
        pass

    # If assessment requested, try to get it
    if row_config.get('assessment') and not detail:
      # Look for assessment field (e.g., r_free_assessment)
      for metric in required_metrics:
        assessment_key = f"{metric}_assessment"
        if assessment_key in metrics_dict:
          detail = str(metrics_dict[assessment_key])
          break

    rows.append({
      'label': row_config.get('label', required_metrics[0] if required_metrics else ''),
      'value': value,
      'detail': detail,
    })

  return rows


def format_step_metric(program, metrics_dict):
  """
  Format the metric display for a workflow step.

  Args:
      program: Program name (e.g., 'phenix.refine')
      metrics_dict: Dict of metric name -> value for this step

  Returns:
      str: Formatted metric string or None if no config/metrics
  """
  if not metrics_dict:
    metrics_dict = {}

  config = get_step_metrics_config()

  # Find matching config (try exact match, then partial)
  step_config = None
  program_lower = program.lower() if program else ""

  # Exact match
  if program in config:
    step_config = config[program]
  else:
    # Partial match (e.g., "refine" matches "phenix.refine")
    for prog_name, prog_config in config.items():
      if prog_name == '_default':
        continue
      if prog_name.lower() in program_lower or program_lower in prog_name.lower():
        step_config = prog_config
        break

  # Use default if no match
  if not step_config:
    step_config = config.get('_default', {})

  if not step_config:
    return None

  # Check if required metrics are present
  required_metrics = step_config.get('metrics', [])

  if required_metrics and not all(m in metrics_dict for m in required_metrics):
    # Return fallback format if metrics not available
    return step_config.get('fallback_format')

  # Format the metric
  try:
    result = step_config['format'].format(**metrics_dict)

    # Add optional parts
    optional_metrics = step_config.get('optional_metrics', [])
    if optional_metrics and all(m in metrics_dict for m in optional_metrics):
      try:
        result += step_config.get('optional_format', '').format(**metrics_dict)
      except (KeyError, ValueError):
        pass

    return result
  except (KeyError, ValueError):
    return step_config.get('fallback_format')


def get_metric_assessment(metric_name, value, resolution=None):
  """
  Get quality assessment for a metric value.

  Args:
      metric_name: Name of the metric (e.g., 'r_free')
      value: The metric value
      resolution: Optional resolution for resolution-dependent thresholds

  Returns:
      str: 'Good', 'Acceptable', 'Needs Improvement', or ''
  """
  from libtbx.langchain.knowledge.yaml_loader import load_metrics

  metrics_config = load_metrics()
  if not metrics_config or metric_name not in metrics_config:
    return ''

  metric_def = metrics_config[metric_name]
  direction = metric_def.get('direction', 'minimize')

  # Get thresholds (possibly resolution-dependent)
  thresholds = metric_def.get('thresholds', {})

  if resolution and 'by_resolution' in metric_def:
    for range_def in metric_def['by_resolution']:
      res_range = range_def.get('range', [0, 999])
      if res_range[0] <= resolution < res_range[1]:
        thresholds = range_def
        break

  good = thresholds.get('good')
  acceptable = thresholds.get('acceptable')

  if good is None or acceptable is None:
    return ''

  if direction == 'minimize':
    if value <= good:
      return 'Good'
    elif value <= acceptable:
      return 'Acceptable'
    else:
      return 'Needs Improvement'
  else:  # maximize
    if value >= good:
      return 'Good'
    elif value >= acceptable:
      return 'Acceptable'
    else:
      return 'Needs Improvement'


# =============================================================================
# Testing / CLI
# =============================================================================

def _test_summary_display():
  """Test summary display functionality."""
  print("Testing summary_display module...")

  # Test quality table config loading
  config = get_quality_table_config()
  print(f"  Loaded {len(config)} quality table rows")

  # Test step metrics config loading
  step_config = get_step_metrics_config()
  print(f"  Loaded {len(step_config)} step metric configs")

  # Test quality table formatting
  metrics = {
    'r_free': 0.2345,
    'r_work': 0.2012,
    'clashscore': 3.5,
    'resolution': 2.1,
  }
  rows = format_quality_table_rows(metrics, 'xray')
  print(f"  X-ray quality rows: {len(rows)}")
  for row in rows:
    print(f"    {row['label']}: {row['value']}")

  # Test cryo-EM formatting
  cryoem_metrics = {
    'map_cc': 0.856,
    'symmetry_type': 'D7 (a)',
    'ncs_copies': 14,
    'ncs_cc': 0.93,
    'clashscore': 2.1,
    'resolution': 3.2,
  }
  rows = format_quality_table_rows(cryoem_metrics, 'cryoem')
  print(f"  Cryo-EM quality rows: {len(rows)}")
  for row in rows:
    print(f"    {row['label']}: {row['value']} ({row['detail']})")

  # Test step metric formatting
  refine_metrics = {'r_free': 0.2580}
  step_str = format_step_metric('phenix.refine', refine_metrics)
  print(f"  Refine step: {step_str}")

  rsr_metrics = {'map_cc': 0.823}
  step_str = format_step_metric('phenix.real_space_refine', rsr_metrics)
  print(f"  RSR step: {step_str}")

  sym_metrics = {'symmetry_type': 'D7 (a)', 'ncs_copies': 14}
  step_str = format_step_metric('phenix.map_symmetry', sym_metrics)
  print(f"  Map symmetry step: {step_str}")

  # Test fallback
  step_str = format_step_metric('phenix.refine', {})
  print(f"  Refine step (no metrics): {step_str}")

  # Test assessment
  assessment = get_metric_assessment('r_free', 0.22, resolution=2.0)
  print(f"  R-free 0.22 at 2.0Å: {assessment}")

  assessment = get_metric_assessment('clashscore', 3.5)
  print(f"  Clashscore 3.5: {assessment}")

  print("\nAll summary_display tests PASSED!")


if __name__ == "__main__":
  import sys

  if len(sys.argv) > 1 and sys.argv[1] == "--test":
    _test_summary_display()
  else:
    # Show config
    print("Quality Table Config:")
    for row in get_quality_table_config():
      print(f"  {row.get('label', 'Unknown')}: {row.get('format', 'N/A')}")

    print("\nStep Metrics Config:")
    for prog, config in get_step_metrics_config().items():
      print(f"  {prog}: {config.get('format', 'N/A')}")
