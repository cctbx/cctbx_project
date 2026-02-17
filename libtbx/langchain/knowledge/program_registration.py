"""
Auto-registration of programs from YAML configuration.

This module automatically generates tracking flags and detection patterns
from programs.yaml, eliminating the need to manually add done flags when
adding new programs.

Usage:
    from libtbx.langchain.knowledge.program_registration import (
        get_trackable_programs,
        get_initial_history_flags,
        detect_programs_in_history,
    )

    # Get initial flags dict with all done flags set to False
    flags = get_initial_history_flags()

    # Detect which programs have been run
    updated_flags = detect_programs_in_history(history, flags)
"""

from functools import lru_cache

# Silence unused import warning - lru_cache is used as decorator
assert lru_cache is not None


@lru_cache(maxsize=1)
def get_trackable_programs():
  """
  Get all programs with run_once: true (from done_tracking block).

  Returns:
      dict: program_name -> {
          'short_name': str,       # e.g., 'map_symmetry'
          'done_flag': str,        # e.g., 'map_symmetry_done'
          'detection_patterns': list,  # Patterns to detect in history
      }
  """
  try:
    from libtbx.langchain.knowledge.yaml_loader import load_programs
  except ImportError:
    from knowledge.yaml_loader import load_programs

  programs = load_programs()
  if not programs:
    return {}

  trackable = {}

  for prog_name, prog_def in programs.items():
    if not isinstance(prog_def, dict):
      continue

    # Read from done_tracking block (single source of truth)
    tracking = prog_def.get('done_tracking', {})
    if not tracking.get('run_once', False):
      continue

    # Use explicit flag name from YAML
    done_flag = tracking.get('flag')
    if not done_flag:
      continue

    # Generate short name for detection patterns
    short_name = prog_name.replace('phenix.', '').replace('.', '_')

    trackable[prog_name] = {
      'short_name': short_name,
      'done_flag': done_flag,
      'detection_patterns': [
        short_name,                          # e.g., 'map_symmetry'
        short_name.replace('_', ''),         # e.g., 'mapsymmetry'
        prog_name,                           # e.g., 'phenix.map_symmetry'
        prog_name.replace('.', ''),          # e.g., 'phenixmapsymmetry'
      ],
    }

  return trackable


def get_initial_history_flags():
  """
  Get initial history info dict with all auto-generated done flags.

  This returns ONLY the auto-generated flags. The caller should merge
  these with any manually-defined flags.

  Returns:
      dict: flag_name -> False for all trackable programs
  """
  trackable = get_trackable_programs()
  return {info['done_flag']: False for info in trackable.values()}


def detect_programs_in_history(history, existing_flags=None):
  """
  Detect which trackable programs have been run SUCCESSFULLY in history.

  Only marks programs as done if they completed without errors.
  Failed runs don't count - the program can be retried.

  Args:
      history: List of history entries (dicts with 'program', 'command', 'result', etc.)
      existing_flags: Optional existing flags dict to update

  Returns:
      dict: Updated flags with detected programs marked as True
  """
  if existing_flags is None:
    flags = get_initial_history_flags()
  else:
    flags = dict(existing_flags)

  if not history:
    return flags

  trackable = get_trackable_programs()

  for entry in history:
    # Extract text to search in
    if isinstance(entry, str):
      combined = entry.lower()
      result = ""
    elif isinstance(entry, dict):
      prog = entry.get('program', '')
      cmd = entry.get('command', '')
      combined = f"{prog} {cmd}".lower()
      result = str(entry.get('result', '')).upper()
    else:
      continue

    # Check if this run FAILED - if so, don't mark as done
    # This allows the program to be retried
    # Look for failure PATTERNS, not just the word (avoid matching "No ERROR detected")
    failure_patterns = [
      'FAILED',           # Common failure indicator
      'SORRY:',           # Phenix error prefix
      'SORRY ',           # Phenix error prefix with space
      'ERROR:',           # Error with colon
      'ERROR ',           # Error as prefix
      ': ERROR',          # Error after colon
      'TRACEBACK',        # Python exception
      'EXCEPTION',        # Exception indicator
    ]
    is_failed = any(pattern in result for pattern in failure_patterns)
    if is_failed:
      continue  # Skip failed runs - program can be retried

    # Remove special characters for matching
    combined_clean = combined.replace('_', '').replace('.', '').replace('-', '')

    # Check each trackable program
    for prog_name, prog_info in trackable.items():
      if flags.get(prog_info['done_flag']):
        continue  # Already detected

      for pattern in prog_info['detection_patterns']:
        pattern_clean = pattern.lower().replace('_', '').replace('.', '').replace('-', '')
        if pattern_clean in combined_clean:
          flags[prog_info['done_flag']] = True
          break

  return flags


def get_all_done_flags():
  """
  Get list of all auto-generated done flag names.

  Returns:
      list: List of flag names (e.g., ['xtriage_done', 'mtriage_done', ...])
  """
  trackable = get_trackable_programs()
  return [info['done_flag'] for info in trackable.values()]


def is_program_trackable(program_name):
  """
  Check if a program is auto-trackable (has done_tracking.run_once: true).

  Args:
      program_name: Full program name (e.g., 'phenix.map_symmetry')

  Returns:
      bool: True if program is trackable
  """
  trackable = get_trackable_programs()
  return program_name in trackable


def get_done_flag_for_program(program_name):
  """
  Get the done flag name for a program.

  Args:
      program_name: Full program name (e.g., 'phenix.map_symmetry')

  Returns:
      str or None: Flag name (e.g., 'map_symmetry_done') or None if not trackable
  """
  trackable = get_trackable_programs()
  prog_info = trackable.get(program_name)
  return prog_info['done_flag'] if prog_info else None


def get_program_done_flag_map():
  """
  Get a mapping of ALL program names to their workflow done flags.

  Reads from done_tracking blocks in programs.yaml — this is the
  SINGLE SOURCE OF TRUTH for program→done_flag mapping.

  Programs without done_tracking (e.g., phenix.map_correlations) are not
  included — they don't have workflow done flags.

  Returns:
      dict: {program_name: done_flag_name}
  """
  try:
    from libtbx.langchain.knowledge.yaml_loader import load_programs
  except ImportError:
    from knowledge.yaml_loader import load_programs

  programs = load_programs()
  if not programs:
    return {}

  result = {}
  for prog_name, prog_def in programs.items():
    if not isinstance(prog_def, dict):
      continue
    tracking = prog_def.get('done_tracking')
    if tracking and tracking.get('flag'):
      result[prog_name] = tracking['flag']

  return result


# =============================================================================
# Testing / CLI
# =============================================================================

def _test_registration():
  """Test auto-registration functionality."""
  print("Testing program auto-registration...")

  # Test trackable programs
  trackable = get_trackable_programs()
  print(f"Found {len(trackable)} trackable programs:")
  for prog_name, info in sorted(trackable.items()):
    print(f"  {prog_name}: {info['done_flag']}")

  # Test initial flags
  flags = get_initial_history_flags()
  print(f"\nInitial flags ({len(flags)}):")
  for flag, value in sorted(flags.items()):
    print(f"  {flag}: {value}")

  # Test detection
  test_history = [
    {"program": "phenix.xtriage", "command": "phenix.xtriage data.mtz"},
    {"program": "phenix.mtriage", "command": "phenix.mtriage map.mrc"},
    {"program": "phenix.map_symmetry", "command": "phenix.map_symmetry map.mrc"},
  ]

  detected = detect_programs_in_history(test_history, flags)
  print(f"\nAfter detecting in history:")
  for flag, value in sorted(detected.items()):
    print(f"  {flag}: {value}")

  # Verify
  assert detected.get('xtriage_done') == True, "Should detect xtriage"
  assert detected.get('mtriage_done') == True, "Should detect mtriage"
  assert detected.get('map_symmetry_done') == True, "Should detect map_symmetry"

  print("\nAll program_registration tests PASSED!")


if __name__ == "__main__":
  import sys

  if len(sys.argv) > 1 and sys.argv[1] == "--test":
    _test_registration()
  else:
    # Show trackable programs
    trackable = get_trackable_programs()
    print(f"Trackable programs ({len(trackable)}):")
    for prog_name, info in sorted(trackable.items()):
      print(f"  {prog_name}")
      print(f"    flag: {info['done_flag']}")
      print(f"    patterns: {info['detection_patterns']}")
