from __future__ import absolute_import, division, print_function
import os
import json

TEMPLATE_PATH = os.path.join(os.path.dirname(__file__), "command_templates.json")


class TemplateBuilder(object):
  """
  Deterministic command constructor.
  Converts abstract IntentJSON into concrete CLI strings.
  """

  def __init__(self, template_path=TEMPLATE_PATH):
    if not os.path.exists(template_path):
      raise RuntimeError("Command templates not found at: %s" % template_path)

    with open(template_path, 'r') as f:
      self.templates = json.load(f)

  def get_programs(self):
    """Return list of supported programs."""
    return list(self.templates.keys())

  def build_command(self, program, files, strategy=None, log=None):
    """
    Builds a command from an explicit Intent.

    Args:
      program (str): e.g., "phenix.refine"
      files (dict): e.g., {"model": "file.pdb", "data": "file.mtz"}
                    For multiple files: {"half_map": ["map_1.ccp4", "map_2.ccp4"]}
      strategy (dict): e.g., {"use_simulated_annealing": True, "nproc": 8}
      log: Optional logging function

    Returns:
      str: The full shell command.
    """
    if log is None:
      log = lambda x: None

    if program not in self.templates:
      raise ValueError("Unknown program: %s" % program)

    template = self.templates[program]
    cmd_parts = [program]

    # 1. Handle File Inputs
    for slot_name, slot_def in template.get("file_slots", {}).items():
      file_path = files.get(slot_name)

      if not file_path and slot_def.get("required"):
        raise ValueError("Missing required file slot '%s' for %s" % (slot_name, program))

      if file_path:
        flag = slot_def.get("flag", "")
        is_multiple = slot_def.get("multiple", False)

        # Handle multiple files (e.g., half_map=file1.ccp4 half_map=file2.ccp4)
        if is_multiple and isinstance(file_path, list):
          for fp in file_path:
            if flag and not flag.endswith("="):
              cmd_parts.append("%s %s" % (flag, fp))
            else:
              cmd_parts.append("%s%s" % (flag, fp))
        else:
          # Single file
          if flag and not flag.endswith("="):
            # Case: "-model file.pdb" (space between flag and value)
            cmd_parts.append("%s %s" % (flag, file_path))
          else:
            # Case: "data=file.mtz" OR just "file.pdb" (positional)
            cmd_parts.append("%s%s" % (flag, file_path))

    # 2. Handle Defaults (can be overridden by strategy)
    defaults = dict(template.get("defaults", {}))

    # 3. Handle Strategy Flags
    if strategy:
      strat_defs = template.get("strategy_flags", {})
      for strat_key, strat_val in strategy.items():
        if strat_key in strat_defs:
          mapping = strat_defs[strat_key]

          if isinstance(strat_val, bool):
            key_str = "true" if strat_val else "false"
            if key_str in mapping:
              arg = mapping[key_str]
              if arg:
                cmd_parts.append(arg)
                log("DEBUG: Applied boolean strategy %s=%s -> %s" % (strat_key, strat_val, arg))

          elif "format" in mapping and strat_val is not None:
            # Handle numeric/string values with format
            formatted = mapping["format"].format(strat_val)
            cmd_parts.append(formatted)
            log("DEBUG: Applied format strategy %s=%s -> %s" % (strat_key, strat_val, formatted))

            # If this overrides a default, remove the default
            # e.g., nproc=8 overrides control.nproc=4
            for default_key in list(defaults.keys()):
              if strat_key in default_key.lower() or default_key.lower().endswith(strat_key):
                del defaults[default_key]
                log("DEBUG: Strategy %s=%s overrides default %s" % (strat_key, strat_val, default_key))
        else:
          # Unknown strategy key - DO NOT pass through raw
          # This prevents malformed parameters like "nproc=4" instead of "control.nproc=4"
          log("WARNING: Unknown strategy flag '%s' for %s - IGNORED (not in template)" % (strat_key, program))

    # 4. Add remaining defaults
    for key, val in defaults.items():
      cmd_parts.append("%s=%s" % (key, val))

    command = " ".join(cmd_parts)
    log("DEBUG: Built command: %s" % command)
    return command

  def build_command_for_program(self, program, available_files, categorized_files=None, log=None):
    """
    THE FALLBACK: Auto-selects files based on extensions and patterns.

    Args:
      program (str): The program to build a command for
      available_files (list): List of available file paths
      categorized_files (dict): Optional pre-categorized files from workflow_state
                                (e.g., {'full_map': [...], 'half_map': [...], ...})
      log: Optional logging function

    Returns:
      str: Command string, or None if required files not found
    """
    if log is None:
      log = lambda x: None

    if program not in self.templates:
      log("ERROR: Unknown program %s" % program)
      return None

    template = self.templates[program]
    selected_files = {}
    used_files = set()  # Track used files to avoid duplicates

    for slot_name, slot_def in template.get("file_slots", {}).items():
      valid_exts = slot_def.get("extensions", [])
      exclude_patterns = slot_def.get("exclude_patterns", [])
      prefer_patterns = slot_def.get("prefer_patterns", [])
      priority_patterns = slot_def.get("priority_patterns", [])
      is_multiple = slot_def.get("multiple", False)

      # If we have pre-categorized files and this slot has a matching category
      if categorized_files and slot_name in categorized_files:
        candidates = list(categorized_files[slot_name])
        log("DEBUG: Using pre-categorized files for slot '%s': %s" % (slot_name, candidates))
        
        # If category is empty, skip this slot (don't fall back to extension matching)
        if not candidates:
          if slot_def.get("required"):
            log("ERROR: No candidates for required slot '%s'" % slot_name)
            return None
          continue
      else:
        # No categorized files for this slot - filter by extension
        candidates = [f for f in available_files
                      if any(f.lower().endswith(ext) for ext in valid_exts)]

      # Remove already used files
      candidates = [f for f in candidates if f not in used_files]

      # Apply exclude patterns
      for pattern in exclude_patterns:
        candidates = [f for f in candidates if pattern.lower() not in f.lower()]

      if not candidates:
        if slot_def.get("required"):
          log("ERROR: No candidates for required slot '%s'" % slot_name)
          return None
        continue

      # Score and sort candidates
      def score_file(f):
        score = 0
        basename = os.path.basename(f).lower()

        # Priority patterns (highest)
        for i, pattern in enumerate(priority_patterns):
          if pattern.lower() in basename:
            score += 1000 - i * 10  # Earlier patterns = higher score

        # Prefer patterns
        for pattern in prefer_patterns:
          if pattern.lower() in basename:
            score += 100

        # General heuristics
        if "with_ligand" in basename:
          score += 50
        if "_refine_" in basename:
          score += 40
        if "phaser" in basename.lower():
          score += 30

        # Prefer higher cycle numbers for refined models (most recent)
        # e.g., refine_003.pdb scores higher than refine_001.pdb
        import re
        cycle_match = re.search(r'refine[_.]?(\d+)', basename)
        if cycle_match:
          cycle_num = int(cycle_match.group(1))
          score += cycle_num * 5  # Each cycle adds 5 points

        # Prefer longer names (usually more specific)
        score += len(basename) * 0.1

        return score

      candidates.sort(key=score_file, reverse=True)

      if is_multiple:
        # For multiple file slots (e.g., half_map), use all candidates
        selected_files[slot_name] = candidates
        for f in candidates:
          used_files.add(f)
        log("DEBUG: Auto-selected %d files for slot '%s': %s" % (len(candidates), slot_name, candidates))
      else:
        # Single file slot - pick the best one
        best_file = candidates[0]
        selected_files[slot_name] = best_file
        used_files.add(best_file)
        log("DEBUG: Auto-selected '%s' for slot '%s'" % (best_file, slot_name))

    try:
      return self.build_command(program, selected_files, strategy={}, log=log)
    except ValueError as e:
      log("ERROR: Failed to build command: %s" % e)
      return None

  def validate_intent(self, program, files, available_files):
    """
    Validate that an intent's files exist in the available files list.

    Args:
      program (str): The program name
      files (dict): File selections from intent
      available_files (list): Files available on client

    Returns:
      tuple: (is_valid, error_message)
    """
    if program not in self.templates:
      return False, "Unknown program: %s" % program

    available_basenames = set(os.path.basename(f) for f in available_files)
    available_set = set(available_files) | available_basenames

    missing = []
    for slot_name, file_path in files.items():
      basename = os.path.basename(file_path)
      if file_path not in available_set and basename not in available_set:
        missing.append("%s=%s" % (slot_name, file_path))

    if missing:
      return False, "Files not found: %s" % ", ".join(missing)

    return True, None

  def get_required_slots(self, program):
    """
    Get the required file slots for a program.

    Args:
      program (str): The program name

    Returns:
      list: Names of required file slots
    """
    if program not in self.templates:
      return []

    template = self.templates[program]
    required = []
    for slot_name, slot_def in template.get("file_slots", {}).items():
      if slot_def.get("required"):
        required.append(slot_name)
    return required

  def get_strategy_options(self, program):
    """
    Get available strategy options for a program.

    Args:
      program (str): The program name

    Returns:
      dict: Strategy flag definitions
    """
    if program not in self.templates:
      return {}

    return self.templates[program].get("strategy_flags", {})
