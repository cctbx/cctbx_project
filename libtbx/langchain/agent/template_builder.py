"""
Template Builder for PHENIX AI Agent.

Builds command strings from program, files, and strategy.

This module now supports two backends:
1. YAML-driven ProgramRegistry (preferred)
2. Legacy JSON command_templates.json (fallback)

Usage:
    builder = TemplateBuilder()
    cmd = builder.build_command("phenix.refine", {"model": "x.pdb", "mtz": "x.mtz"})
"""

from __future__ import absolute_import, division, print_function
import os
import json
import re

# Import ProgramRegistry
from libtbx.langchain.agent.program_registry import ProgramRegistry

TEMPLATE_PATH = os.path.join(os.path.dirname(__file__), "command_templates.json")


class TemplateBuilder(object):
    """
    Deterministic command constructor.
    Converts abstract IntentJSON into concrete CLI strings.

    Uses YAML-driven ProgramRegistry.
    """

    def __init__(self, template_path=TEMPLATE_PATH, use_yaml=True):
        """
        Initialize the builder.

        Args:
            template_path: Path to JSON templates (legacy fallback)
            use_yaml: If True, use YAML-driven ProgramRegistry
        """
        self.use_yaml = use_yaml
        self._registry = None
        self._json_templates = None
        self._template_path = template_path

        if self.use_yaml:
            self._registry = ProgramRegistry(use_yaml=True)
        else:
            self._load_json_templates()

    def _load_json_templates(self):
        """Load legacy JSON templates."""
        if not os.path.exists(self._template_path):
            raise RuntimeError("Command templates not found at: %s" % self._template_path)

        with open(self._template_path, 'r') as f:
            self._json_templates = json.load(f)

    @property
    def templates(self):
        """For backward compatibility - return JSON templates."""
        if self._json_templates is None:
            self._load_json_templates()
        return self._json_templates

    def get_programs(self):
        """Return list of supported programs."""
        if self.use_yaml:
            return self._registry.get_all_programs()
        return list(self.templates.keys())

    def build_command(self, program, files, strategy=None, log=None):
        """
        Builds a command from an explicit Intent.

        Args:
            program (str): e.g., "phenix.refine"
            files (dict): e.g., {"model": "file.pdb", "mtz": "file.mtz"}
                          For multiple files: {"half_map": ["map_1.ccp4", "map_2.ccp4"]}
            strategy (dict): e.g., {"simulated_annealing": True, "nproc": 8}
            log: Optional logging function

        Returns:
            str: The full shell command.
        """
        if log is None:
            log = lambda x: None

        if self.use_yaml:
            return self._build_command_yaml(program, files, strategy, log)
        else:
            return self._build_command_json(program, files, strategy, log)

    def _build_command_yaml(self, program, files, strategy, log):
        """Build command using YAML-driven registry."""
        try:
            cmd = self._registry.build_command(program, files, strategy, log)
            log("DEBUG: Built command (YAML): %s" % cmd)
            return cmd
        except ValueError as e:
            raise ValueError(str(e))

    def _build_command_json(self, program, files, strategy, log):
        """Build command using JSON templates (legacy)."""
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
                        cmd_parts.append("%s %s" % (flag, file_path))
                    else:
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
                        formatted = mapping["format"].format(strat_val)
                        cmd_parts.append(formatted)
                        log("DEBUG: Applied format strategy %s=%s -> %s" % (strat_key, strat_val, formatted))

                        # Override defaults
                        for default_key in list(defaults.keys()):
                            if strat_key in default_key.lower() or default_key.lower().endswith(strat_key):
                                del defaults[default_key]
                else:
                    log("WARNING: Unknown strategy flag '%s' for %s - IGNORED" % (strat_key, program))

        # 4. Add remaining defaults
        for key, val in defaults.items():
            cmd_parts.append("%s=%s" % (key, val))

        command = " ".join(cmd_parts)
        log("DEBUG: Built command (JSON): %s" % command)
        return command

    def build_command_for_program(self, program, available_files, categorized_files=None,
                                   log=None, context=None):
        """
        THE FALLBACK: Auto-selects files based on extensions and patterns.

        Args:
            program (str): The program to build a command for
            available_files (list): List of available file paths
            categorized_files (dict): Optional pre-categorized files from workflow_state
            log: Optional logging function
            context (dict): Optional context for invariant validation (resolution, etc.)

        Returns:
            str: Command string, or None if required files not found
        """
        if log is None:
            log = lambda x: None

        # Get template (from YAML or JSON)
        if self.use_yaml:
            prog_def = self._registry.get_program(program)
            if not prog_def:
                log("ERROR: Unknown program %s" % program)
                return None
            # Convert YAML format to match JSON structure for file selection
            file_slots = self._yaml_inputs_to_slots(prog_def)
        else:
            if program not in self.templates:
                log("ERROR: Unknown program %s" % program)
                return None
            file_slots = self.templates[program].get("file_slots", {})

        selected_files = {}
        used_files = set()

        for slot_name, slot_def in file_slots.items():
            valid_exts = slot_def.get("extensions", [])
            exclude_patterns = slot_def.get("exclude_patterns", [])
            prefer_patterns = slot_def.get("prefer_patterns", [])
            priority_patterns = slot_def.get("priority_patterns", [])
            is_multiple = slot_def.get("multiple", False)
            is_required = slot_def.get("required", False)

            # If we have pre-categorized files and this slot has a matching category
            if categorized_files and slot_name in categorized_files:
                candidates = list(categorized_files[slot_name])
                log("DEBUG: Using pre-categorized files for slot '%s': %s" % (slot_name, candidates))

                if not candidates:
                    if is_required:
                        log("ERROR: No candidates for required slot '%s'" % slot_name)
                        return None
                    continue
            else:
                # Filter by extension
                candidates = [f for f in available_files
                              if any(f.lower().endswith(ext) for ext in valid_exts)]

            # Remove already used files
            candidates = [f for f in candidates if f not in used_files]

            # Apply exclude patterns
            for pattern in exclude_patterns:
                candidates = [f for f in candidates if pattern.lower() not in f.lower()]

            if not candidates:
                if is_required:
                    log("ERROR: No candidates for required slot '%s'" % slot_name)
                    return None
                continue

            # Score and sort candidates
            candidates.sort(key=lambda f: self._score_file(f, priority_patterns, prefer_patterns), reverse=True)

            if is_multiple:
                selected_files[slot_name] = candidates
                for f in candidates:
                    used_files.add(f)
                log("DEBUG: Auto-selected %d files for slot '%s'" % (len(candidates), slot_name))
            else:
                best_file = candidates[0]
                selected_files[slot_name] = best_file
                used_files.add(best_file)
                log("DEBUG: Auto-selected '%s' for slot '%s'" % (best_file, slot_name))

        try:
            # Apply invariant validation and fixes (e.g., auto-fill resolution)
            strategy = {}
            if context:
                selected_files, strategy, warnings = self.validate_and_fix(
                    program, selected_files, strategy,
                    log=log, context=context
                )
                for warning in warnings:
                    log("FALLBACK: %s" % warning)

            return self.build_command(program, selected_files, strategy, log=log)
        except ValueError as e:
            log("ERROR: Failed to build command: %s" % e)
            return None

    def _yaml_inputs_to_slots(self, prog_def):
        """Convert YAML inputs format to JSON file_slots format."""
        slots = {}
        inputs = prog_def.get("inputs", {})

        for section, is_required in [("required", True), ("optional", False)]:
            for slot_name, slot_def in inputs.get(section, {}).items():
                slots[slot_name] = {
                    "extensions": slot_def.get("extensions", []),
                    "flag": slot_def.get("flag", ""),
                    "required": is_required,
                    "multiple": slot_def.get("multiple", False),
                    "exclude_patterns": slot_def.get("exclude_patterns", []),
                    "prefer_patterns": slot_def.get("prefer_patterns", []),
                    "priority_patterns": slot_def.get("priority_patterns", []),
                }

        return slots

    def _score_file(self, f, priority_patterns, prefer_patterns):
        """Score a file for selection."""
        score = 0
        basename = os.path.basename(f).lower()

        # Priority patterns (highest)
        for i, pattern in enumerate(priority_patterns):
            if pattern.lower() in basename:
                score += 1000 - i * 10

        # Prefer patterns
        for pattern in prefer_patterns:
            if pattern.lower() in basename:
                score += 100

        # General heuristics
        if "with_ligand" in basename:
            score += 50
        if "_refine_" in basename:
            score += 40
        if "phaser" in basename:
            score += 30

        # Prefer higher cycle numbers
        cycle_match = re.search(r'refine[_.]?(\d+)', basename)
        if cycle_match:
            score += int(cycle_match.group(1)) * 5

        # Prefer longer names
        score += len(basename) * 0.1

        return score

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
        if self.use_yaml:
            prog = self._registry.get_program(program)
            if not prog:
                return False, "Unknown program: %s" % program
        else:
            if program not in self.templates:
                return False, "Unknown program: %s" % program

        available_basenames = set(os.path.basename(f) for f in available_files)
        available_set = set(available_files) | available_basenames

        missing = []
        for slot_name, file_path in files.items():
            if isinstance(file_path, list):
                for fp in file_path:
                    basename = os.path.basename(fp)
                    if fp not in available_set and basename not in available_set:
                        missing.append("%s=%s" % (slot_name, fp))
            else:
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
        if self.use_yaml:
            return self._registry.get_required_inputs(program)

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
        if self.use_yaml:
            return self._registry.get_strategy_flags(program)

        if program not in self.templates:
            return {}

        return self.templates[program].get("strategy_flags", {})

    def get_program_description(self, program):
        """
        Get human-readable description of a program.

        Args:
            program (str): The program name

        Returns:
            str: Description
        """
        if self.use_yaml:
            return self._registry.get_description(program)

        if program not in self.templates:
            return ""

        return self.templates[program].get("description", "")

    def get_program_hints(self, program):
        """
        Get usage hints for a program.

        Args:
            program (str): The program name

        Returns:
            list: Hint strings
        """
        if self.use_yaml:
            return self._registry.get_hints(program)

        if program not in self.templates:
            return []

        return self.templates[program].get("hints", [])

    def validate_and_fix(self, program, files, strategy, log=None, context=None):
        """
        Validate program invariants and apply fixes if needed.

        This is the single place where program-specific constraints are enforced.
        Invariants are defined in programs.yaml and checked here before
        building the command.

        Args:
            program (str): The program name
            files (dict): Selected files
            strategy (dict): Strategy options
            log: Optional logging function
            context (dict): Optional context for auto-fills (resolution, etc.)

        Returns:
            tuple: (files, strategy, warnings)
                - files: Possibly modified files dict
                - strategy: Possibly modified strategy dict
                - warnings: List of warning messages about fixes applied
        """
        if log is None:
            log = lambda x: None
        if context is None:
            context = {}

        warnings = []

        if not self.use_yaml:
            return files, strategy, warnings

        # Get program definition
        prog_def = self._registry.get_program(program)
        if not prog_def:
            return files, strategy, warnings

        # Get invariants
        invariants = prog_def.get("invariants", [])
        if not invariants:
            return files, strategy, warnings

        # Check each invariant
        for invariant in invariants:
            name = invariant.get("name", "unnamed")
            check = invariant.get("check", {})

            if not self._check_invariant(check, files, strategy):
                # Invariant violated - apply fix
                fix = invariant.get("fix", {})
                message = invariant.get("message", "Invariant '%s' violated, fix applied" % name)

                files, strategy, fix_msg = self._apply_fix(fix, files, strategy, context)
                if fix_msg:
                    message = fix_msg  # Use more specific message from fix
                warnings.append(message)
                log("INVARIANT: %s - %s" % (name, message))

        return files, strategy, warnings

    def _check_invariant(self, check, files, strategy):
        """
        Check if an invariant condition is satisfied.

        Supports:
            - has_file: [list of file slots] - True if any file slot is filled
            - strategy_equals: {key: value} - True if strategy[key] == value
            - any_of: [list of conditions] - True if any condition is met
            - all_of: [list of conditions] - True if all conditions are met

        Args:
            check: Check definition from YAML
            files: Current files dict
            strategy: Current strategy dict

        Returns:
            bool: True if check passes, False if violated
        """
        if not check:
            return True

        # Handle any_of (OR logic)
        if "any_of" in check:
            for sub_check in check["any_of"]:
                if self._check_invariant(sub_check, files, strategy):
                    return True
            return False

        # Handle all_of (AND logic)
        if "all_of" in check:
            for sub_check in check["all_of"]:
                if not self._check_invariant(sub_check, files, strategy):
                    return False
            return True

        # Handle has_file
        if "has_file" in check:
            file_slots = check["has_file"]
            if isinstance(file_slots, str):
                file_slots = [file_slots]
            for slot in file_slots:
                if files.get(slot):
                    return True
            return False

        # Handle has_strategy - check if strategy key exists and has a value
        if "has_strategy" in check:
            strategy_keys = check["has_strategy"]
            if isinstance(strategy_keys, str):
                strategy_keys = [strategy_keys]
            for key in strategy_keys:
                val = strategy.get(key)
                if val is not None and val != "":
                    return True
            return False

        # Handle strategy_equals
        if "strategy_equals" in check:
            for key, expected_value in check["strategy_equals"].items():
                if strategy.get(key) == expected_value:
                    return True
            return False

        return True

    def _apply_fix(self, fix, files, strategy, context=None):
        """
        Apply a fix to files and/or strategy.

        Supports:
            - set_strategy: {key: value} - Set strategy values
            - set_file: {slot: value} - Set file values
            - remove_file: [list] - Remove file slots
            - auto_fill_resolution: true - Auto-fill resolution from context

        Args:
            fix: Fix definition from YAML
            files: Current files dict
            strategy: Current strategy dict
            context: Optional context dict with resolution, etc.

        Returns:
            tuple: (files, strategy, message) - Modified dicts and optional message
        """
        # Make copies to avoid mutating originals
        files = dict(files)
        strategy = dict(strategy)
        message = None
        context = context or {}

        if "set_strategy" in fix:
            for key, value in fix["set_strategy"].items():
                strategy[key] = value

        if "set_file" in fix:
            for slot, value in fix["set_file"].items():
                files[slot] = value

        if "remove_file" in fix:
            for slot in fix["remove_file"]:
                files.pop(slot, None)

        # Auto-fill resolution from context
        if fix.get("auto_fill_resolution"):
            resolution = None
            source = None

            # Priority order for resolution sources
            if context.get("session_resolution"):
                resolution = context["session_resolution"]
                source = "session"
            elif context.get("resolution"):
                resolution = context["resolution"]
                source = "workflow"
            elif context.get("workflow_state", {}).get("resolution"):
                resolution = context["workflow_state"]["resolution"]
                source = "workflow_state"

            if resolution:
                strategy["resolution"] = round(float(resolution), 1)
                message = "Auto-filled resolution=%.1fÅ from %s" % (resolution, source)

        return files, strategy, message


# =============================================================================
# TESTING
# =============================================================================

if __name__ == "__main__":
    print("Testing TemplateBuilder...")
    print()

    builder = TemplateBuilder()
    print("Using YAML:", builder.use_yaml)
    print()

    # Test get_programs
    programs = builder.get_programs()
    print("Programs available: %d" % len(programs))
    print()

    # Test build_command
    print("Testing build_command for phenix.refine:")
    files = {"model": "test.pdb", "mtz": "test.mtz"}
    try:
        cmd = builder.build_command("phenix.refine", files)
        print("  Command:", cmd)
    except Exception as e:
        print("  Error:", e)
    print()

    # Test with strategy
    print("Testing with strategy:")
    files = {"model": "test.pdb", "mtz": "test.mtz"}
    strategy = {"generate_rfree_flags": True}
    try:
        cmd = builder.build_command("phenix.refine", files, strategy)
        print("  Command:", cmd)
    except Exception as e:
        print("  Error:", e)
    print()

    # Test get_required_slots
    print("Required slots for phenix.refine:")
    print("  ", builder.get_required_slots("phenix.refine"))
    print()

    # Test description
    print("Description for phenix.refine:")
    print("  ", builder.get_program_description("phenix.refine"))
