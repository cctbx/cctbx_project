"""
Program Registry for PHENIX AI Agent.

This module provides a unified interface to program definitions,
reading from programs.yaml (YAML-driven) or command_templates.json (legacy).

The registry provides:
- Program descriptions and hints
- Input/output requirements
- Command building
- Log parsing patterns

Usage:
    from libtbx.langchain.agent.program_registry import ProgramRegistry

    registry = ProgramRegistry()

    # Get program info
    info = registry.get_program("phenix.refine")

    # Build command
    cmd = registry.build_command("phenix.refine", files, strategy)

    # Get all programs for experiment type
    programs = registry.get_programs_for_experiment("xray")
"""

from __future__ import absolute_import, division, print_function

import os
import re

# Import YAML loader
from libtbx.langchain.knowledge.yaml_loader import (
    get_program,
    get_all_programs,
    get_programs_for_experiment,
    get_program_description,
    get_program_inputs,
    get_program_outputs,
    get_program_log_patterns,
    get_program_hints,
)


class ProgramRegistry:
    """
    Registry providing program information from YAML or JSON templates.

    This class abstracts the source of program definitions, allowing
    gradual migration from JSON to YAML.
    """

    def __init__(self, use_yaml=True, template_path=None):
        """
        Initialize the registry.

        Args:
            use_yaml: If True, prefer YAML definitions
            template_path: Path to command_templates.json (fallback)
        """
        self.use_yaml = use_yaml
        self._json_templates = None
        self._template_path = template_path

        if not self.use_yaml:
            self._load_json_templates()

    def _load_json_templates(self):
        """Load legacy JSON templates."""
        import json

        if self._template_path is None:
            self._template_path = os.path.join(
                os.path.dirname(__file__), "command_templates.json"
            )

        if os.path.exists(self._template_path):
            with open(self._template_path, 'r') as f:
                self._json_templates = json.load(f)
        else:
            self._json_templates = {}

    # =========================================================================
    # PROGRAM INFORMATION
    # =========================================================================

    def get_program(self, program_name):
        """
        Get full definition for a program.

        Args:
            program_name: Name of program (e.g., "phenix.refine")

        Returns:
            dict: Program definition or None
        """
        if self.use_yaml:
            return get_program(program_name)
        else:
            return self._json_templates.get(program_name)

    def get_all_programs(self):
        """
        Get list of all known program names.

        Returns:
            list: Program names
        """
        if self.use_yaml:
            return get_all_programs()
        else:
            return list(self._json_templates.keys())

    def get_description(self, program_name):
        """
        Get human-readable description.

        Args:
            program_name: Name of program

        Returns:
            str: Description
        """
        if self.use_yaml:
            return get_program_description(program_name)
        else:
            prog = self._json_templates.get(program_name, {})
            return prog.get("description", "")

    def get_hints(self, program_name):
        """
        Get usage hints for a program.

        Args:
            program_name: Name of program

        Returns:
            list: Hint strings
        """
        if self.use_yaml:
            return get_program_hints(program_name)
        else:
            prog = self._json_templates.get(program_name, {})
            return prog.get("hints", [])

    def get_programs_for_experiment(self, experiment_type):
        """
        Get programs valid for an experiment type.

        Args:
            experiment_type: "xray" or "cryoem"

        Returns:
            list: Program names
        """
        if self.use_yaml:
            return get_programs_for_experiment(experiment_type)
        else:
            # In JSON, all programs are potentially valid
            # Filter based on known categorization
            xray_only = ["phenix.xtriage", "phenix.phaser", "phenix.refine",
                        "phenix.autobuild", "phenix.autosol"]
            cryoem_only = ["phenix.mtriage", "phenix.real_space_refine",
                         "phenix.dock_in_map", "phenix.resolve_cryo_em",
                         "phenix.map_sharpening"]
            both = ["phenix.predict_and_build", "phenix.molprobity",
                   "phenix.ligandfit", "phenix.pdbtools"]

            if experiment_type == "xray":
                return xray_only + both
            elif experiment_type == "cryoem":
                return cryoem_only + both
            else:
                return self.get_all_programs()

    # =========================================================================
    # INPUT/OUTPUT REQUIREMENTS
    # =========================================================================

    def get_required_inputs(self, program_name):
        """
        Get required input file types.

        Args:
            program_name: Name of program

        Returns:
            list: Required input slot names
        """
        if self.use_yaml:
            inputs = get_program_inputs(program_name)
            required = inputs.get("required", {})
            return list(required.keys())
        else:
            prog = self._json_templates.get(program_name, {})
            required = []
            for slot_name, slot_def in prog.get("file_slots", {}).items():
                if slot_def.get("required"):
                    required.append(slot_name)
            return required

    def get_input_priorities(self, program_name, input_name):
        """
        Get file category priorities for an input slot.

        Returns the priority order of file categories for selecting
        files for this input, plus any categories to exclude.

        Args:
            program_name: Name of program
            input_name: Name of input slot (e.g., "model", "mtz")

        Returns:
            dict: {
                "categories": [list of categories in priority order],
                "prefer_subcategories": [list of subcategories to prefer within parent],
                "fallback_categories": [list of fallback categories if primary not found],
                "exclude_categories": [list of categories to never use]
            }
        """
        if self.use_yaml:
            prog_def = get_program(program_name)
            if not prog_def:
                return {"categories": [], "prefer_subcategories": [],
                        "fallback_categories": [], "exclude_categories": []}

            input_priorities = prog_def.get("input_priorities", {})
            slot_priorities = input_priorities.get(input_name, {})

            categories = slot_priorities.get("categories", [])
            prefer_subcategories = slot_priorities.get("prefer_subcategories", [])
            fallback_categories = slot_priorities.get("fallback_categories", [])
            exclude_categories = list(slot_priorities.get("exclude_categories", []))

            # If program requires_full_map and this is a map input,
            # automatically exclude half_map
            if prog_def.get("requires_full_map") and input_name in ("map", "full_map"):
                if "half_map" not in exclude_categories:
                    exclude_categories.append("half_map")

            return {
                "categories": categories,
                "prefer_subcategories": prefer_subcategories,
                "fallback_categories": fallback_categories,
                "exclude_categories": exclude_categories
            }
        else:
            # Legacy JSON doesn't have this feature
            return {"categories": [], "prefer_subcategories": [],
                    "fallback_categories": [], "exclude_categories": []}

    def requires_full_map(self, program_name):
        """
        Check if a program requires a full map (cannot use half-maps).

        Args:
            program_name: Name of program

        Returns:
            bool: True if program requires full map
        """
        if self.use_yaml:
            prog_def = get_program(program_name)
            if prog_def:
                return prog_def.get("requires_full_map", False)
        return False

    def get_user_advice_keywords(self, program_name):
        """
        Get keywords that indicate user wants this program.

        These are keywords in user advice that should cause
        this program to be preferred when it's valid.

        Args:
            program_name: Name of program

        Returns:
            list: Keywords (all lowercase)
        """
        if self.use_yaml:
            prog_def = get_program(program_name)
            if not prog_def:
                return []
            return prog_def.get("user_advice_keywords", [])
        else:
            return []

    def get_optional_inputs(self, program_name):
        """
        Get optional input file types.

        Args:
            program_name: Name of program

        Returns:
            list: Optional input slot names
        """
        if self.use_yaml:
            inputs = get_program_inputs(program_name)
            optional = inputs.get("optional", {})
            return list(optional.keys())
        else:
            prog = self._json_templates.get(program_name, {})
            optional = []
            for slot_name, slot_def in prog.get("file_slots", {}).items():
                if not slot_def.get("required"):
                    optional.append(slot_name)
            return optional

    def get_input_extensions(self, program_name, slot_name):
        """
        Get valid file extensions for an input slot.

        Args:
            program_name: Name of program
            slot_name: Name of input slot (e.g., "mtz", "model")

        Returns:
            list: Valid extensions (e.g., [".mtz", ".sca"])
        """
        if self.use_yaml:
            inputs = get_program_inputs(program_name)
            for section in ["required", "optional"]:
                if slot_name in inputs.get(section, {}):
                    return inputs[section][slot_name].get("extensions", [])
            return []
        else:
            prog = self._json_templates.get(program_name, {})
            slot_def = prog.get("file_slots", {}).get(slot_name, {})
            return slot_def.get("extensions", [])

    def get_output_metrics(self, program_name):
        """
        Get metrics this program produces.

        Args:
            program_name: Name of program

        Returns:
            list: Metric names
        """
        if self.use_yaml:
            outputs = get_program_outputs(program_name)
            return outputs.get("metrics", [])
        else:
            # JSON doesn't track output metrics
            return []

    # =========================================================================
    # COMMAND BUILDING
    # =========================================================================

    def get_command_template(self, program_name):
        """
        Get command template string.

        Args:
            program_name: Name of program

        Returns:
            str: Template like "phenix.refine {model} {mtz}"
        """
        if self.use_yaml:
            prog = get_program(program_name)
            if prog:
                return prog.get("command", program_name)
            return program_name
        else:
            # JSON doesn't have explicit templates
            return program_name

    def get_strategy_flags(self, program_name):
        """
        Get available strategy flags for a program.

        Args:
            program_name: Name of program

        Returns:
            dict: Strategy flag definitions
        """
        if self.use_yaml:
            prog = get_program(program_name)
            if prog:
                return prog.get("strategy_flags", {})
            return {}
        else:
            prog = self._json_templates.get(program_name, {})
            return prog.get("strategy_flags", {})

    def get_defaults(self, program_name):
        """
        Get default parameter values.

        Args:
            program_name: Name of program

        Returns:
            dict: Default parameters
        """
        if self.use_yaml:
            prog = get_program(program_name)
            if prog:
                return prog.get("defaults", {})
            return {}
        else:
            prog = self._json_templates.get(program_name, {})
            return prog.get("defaults", {})

    def get_invariants(self, program_name):
        """
        Get program invariants (validation rules and auto-fills).

        Args:
            program_name: Name of program

        Returns:
            list: List of invariant definitions
        """
        if self.use_yaml:
            prog = get_program(program_name)
            if prog:
                return prog.get("invariants", [])
            return []
        else:
            # JSON templates don't have invariants
            return []

    def build_command(self, program_name, files, strategy=None, log=None,
                      file_sources=None, strategy_sources=None):
        """
        Build a command string.

        Args:
            program_name: Name of program
            files: Dict mapping slot names to file paths
            strategy: Optional dict of strategy flags
            log: Optional logging function
            file_sources: Optional dict mapping slot names to source strings
                         (e.g., {"model": "llm_selected", "data_mtz": "best_files"})
            strategy_sources: Optional dict mapping strategy keys to source strings
                             (e.g., {"resolution": "invariant", "nproc": "llm_strategy"})

        Returns:
            str: Complete command string
        """
        if log is None:
            log = lambda x: None
        if file_sources is None:
            file_sources = {}
        if strategy_sources is None:
            strategy_sources = {}

        prog = self.get_program(program_name)
        if not prog:
            raise ValueError("Unknown program: %s" % program_name)

        # Special handling for Phaser multi-ensemble mode
        if prog.get("multi_ensemble") and program_name == "phenix.phaser":
            return self._build_phaser_command(prog, files, strategy, log)

        # Start with command template from YAML if available
        if self.use_yaml:
            cmd_template = prog.get("command", program_name)

            # Get input definitions to check for flags and multiple
            inputs = prog.get("inputs", {})
            all_inputs = {}
            all_inputs.update(inputs.get("required", {}))
            all_inputs.update(inputs.get("optional", {}))

            # Substitute file placeholders
            cmd = cmd_template
            for slot_name, file_path in files.items():
                placeholder = "{%s}" % slot_name
                if placeholder in cmd:
                    # Get the flag for this input
                    slot_def = all_inputs.get(slot_name, {})
                    flag = slot_def.get("flag", "")
                    is_multiple = slot_def.get("multiple", False)

                    if is_multiple and isinstance(file_path, list):
                        # Multiple files - add flag prefix to each
                        if flag:
                            replacement = " ".join("%s%s" % (flag, fp) for fp in file_path)
                        else:
                            replacement = " ".join(file_path)
                    else:
                        # Single file
                        if isinstance(file_path, list):
                            file_path = file_path[0] if file_path else ""
                        # Apply flag if defined
                        if flag:
                            replacement = "%s%s" % (flag, file_path)
                        else:
                            replacement = str(file_path)

                    cmd = cmd.replace(placeholder, replacement)

            # Append files that weren't matched to any placeholder in the template.
            # This handles optional inputs (e.g., partial_model for autosol MR-SAD)
            # that have a flag defined in YAML but no placeholder in the command template.
            for slot_name, file_path in files.items():
                placeholder = "{%s}" % slot_name
                if placeholder not in cmd_template:
                    # This file had no placeholder - append it using its flag
                    slot_def = all_inputs.get(slot_name, {})
                    flag = slot_def.get("flag", "")
                    if flag and file_path:
                        is_multiple = slot_def.get("multiple", False)
                        if is_multiple and isinstance(file_path, list):
                            for fp in file_path:
                                cmd += " %s%s" % (flag, fp)
                        else:
                            if isinstance(file_path, list):
                                file_path = file_path[0] if file_path else ""
                            cmd += " %s%s" % (flag, file_path)

            # Remove any unfilled placeholders (optional files not provided)
            import re
            cmd = re.sub(r'\{[a-z_]+\}', '', cmd)
            cmd = ' '.join(cmd.split())  # Clean up extra spaces

            cmd_parts = cmd.split()
        else:
            # Legacy JSON path - build from scratch
            cmd_parts = [program_name]

            # Get input definitions
            all_inputs = prog.get("file_slots", {})

            # Add files
            for slot_name, slot_def in all_inputs.items():
                file_path = files.get(slot_name)

                if not file_path:
                    is_required = slot_def.get("required", False)
                    if is_required:
                        raise ValueError("Missing required file '%s' for %s" % (slot_name, program_name))
                    continue

                flag = slot_def.get("flag", "")
                is_multiple = slot_def.get("multiple", False)

                if is_multiple and isinstance(file_path, list):
                    for fp in file_path:
                        if flag and not flag.endswith("="):
                            cmd_parts.append("%s %s" % (flag, fp))
                        else:
                            cmd_parts.append("%s%s" % (flag, fp))
                else:
                    if flag and not flag.endswith("="):
                        cmd_parts.append("%s %s" % (flag, file_path))
                    else:
                        cmd_parts.append("%s%s" % (flag, file_path))

        # Handle defaults (can be overridden by strategy)
        defaults = dict(self.get_defaults(program_name))

        # Handle strategy flags
        if strategy:
            strategy_defs = self.get_strategy_flags(program_name)

            for key, value in strategy.items():
                if key not in strategy_defs:
                    # Check if this looks like a PHENIX parameter (contains dots)
                    # These can be passed through directly for error recovery
                    if '.' in key or '=' in key:
                        # Pass through as key=value
                        cmd_parts.append("%s=%s" % (key, value))
                        log("PASSTHROUGH: Adding %s=%s (not in strategy_defs)" % (key, value))
                    else:
                        log("WARNING: Unknown strategy '%s' for %s" % (key, program_name))
                    continue

                flag_def = strategy_defs[key]

                if self.use_yaml:
                    # YAML format: {"flag": "...", "type": "..."}
                    flag_template = flag_def.get("flag", "")
                    if "{value}" in flag_template:
                        cmd_parts.append(flag_template.replace("{value}", str(value)))
                    elif value:  # Boolean true
                        cmd_parts.append(flag_template)
                else:
                    # JSON format: {"true": "...", "false": "..."} or {"format": "..."}
                    if isinstance(value, bool):
                        key_str = "true" if value else "false"
                        if key_str in flag_def:
                            arg = flag_def[key_str]
                            if arg:
                                cmd_parts.append(arg)
                    elif "format" in flag_def and value is not None:
                        cmd_parts.append(flag_def["format"].format(value))

                # Remove overridden defaults
                for default_key in list(defaults.keys()):
                    if key in default_key.lower():
                        del defaults[default_key]

        # Add remaining defaults
        default_keys_used = []
        for key, val in defaults.items():
            cmd_parts.append("%s=%s" % (key, val))
            default_keys_used.append(key)

        # --- Provenance summary ---
        # Log where each piece of the command came from
        command_str = " ".join(cmd_parts)
        log("BUILD_PROVENANCE: --- Command provenance for %s ---" % program_name)
        log("BUILD_PROVENANCE: Command: %s" % command_str)

        # File provenance
        for slot_name, file_path in files.items():
            if isinstance(file_path, list):
                basename = ", ".join(os.path.basename(f) for f in file_path)
            else:
                basename = os.path.basename(str(file_path))
            source = file_sources.get(slot_name, "unknown")
            # Get the flag used
            inputs_def = prog.get("inputs", {}) if self.use_yaml else {}
            all_inp = {}
            all_inp.update(inputs_def.get("required", {}))
            all_inp.update(inputs_def.get("optional", {}))
            slot_def = all_inp.get(slot_name, prog.get("file_slots", {}).get(slot_name, {}))
            flag = slot_def.get("flag", "(positional)")
            log("BUILD_PROVENANCE:   file %-20s = %-30s [source: %s, flag: %s]" % (
                slot_name, basename, source, flag or "(positional)"))

        # Strategy provenance
        if strategy:
            for key, value in strategy.items():
                source = strategy_sources.get(key, "unknown")
                log("BUILD_PROVENANCE:   strat %-19s = %-30s [source: %s]" % (key, value, source))

        # Default provenance
        for key in default_keys_used:
            log("BUILD_PROVENANCE:   default %-17s = %-30s [source: yaml_default]" % (
                key, defaults.get(key, "?")))

        log("BUILD_PROVENANCE: --- end ---")

        return command_str

    def _build_phaser_command(self, prog, files, strategy, log):
        """
        Build command for phenix.phaser.

        Args:
            prog: Program definition from YAML
            files: Dict with mtz, model, sequence (optional)
            strategy: Optional strategy dict
            log: Logging function

        Returns:
            str: Complete command string
        """
        import os

        cmd_parts = ["phenix.phaser"]

        # Get the reflection data file
        data_mtz = files.get("data_mtz") or files.get("mtz")  # Backward compat
        if not data_mtz:
            raise ValueError("Missing required data_mtz file for phenix.phaser")
        cmd_parts.append(str(data_mtz))

        # Get model
        model = files.get("model")
        if isinstance(model, list):
            model = model[0] if model else None
        if not model:
            raise ValueError("Missing required model file for phenix.phaser")
        cmd_parts.append(model)

        # Get sequence - try to match to model name
        sequences = files.get("sequence", [])
        if isinstance(sequences, str):
            sequences = [sequences]

        if sequences:
            # Try to find a sequence that matches the model name
            model_base = os.path.splitext(os.path.basename(model))[0].lower()
            matched_seq = None

            for seq in sequences:
                seq_base = os.path.splitext(os.path.basename(seq))[0].lower()
                if seq_base == model_base or model_base in seq_base or seq_base in model_base:
                    matched_seq = seq
                    break

            if matched_seq:
                cmd_parts.append(matched_seq)
                log("DEBUG: Using matching sequence %s for model %s" % (matched_seq, model))
            else:
                # No match found - use first sequence
                cmd_parts.append(sequences[0])
                log("DEBUG: No matching sequence found for %s, using %s" % (model, sequences[0]))

        # Add defaults
        defaults = dict(prog.get("defaults", {}))

        # Handle strategy flags
        if strategy:
            strategy_defs = prog.get("strategy_flags", {})
            for key, value in strategy.items():
                if key in strategy_defs:
                    flag_def = strategy_defs[key]
                    flag_template = flag_def.get("flag", "")
                    if "{value}" in flag_template:
                        cmd_parts.append(flag_template.replace("{value}", str(value)))
                    elif value:
                        cmd_parts.append(flag_template)

                    # Remove overridden defaults
                    for default_key in list(defaults.keys()):
                        if key in default_key.lower():
                            del defaults[default_key]

        # Add remaining defaults
        for key, val in defaults.items():
            cmd_parts.append("%s=%s" % (key, val))

        return " ".join(cmd_parts)

    # =========================================================================
    # LOG PARSING
    # =========================================================================

    def get_log_patterns(self, program_name):
        """
        Get log parsing patterns for a program.

        Args:
            program_name: Name of program

        Returns:
            dict: metric_name -> pattern definition
        """
        if self.use_yaml:
            return get_program_log_patterns(program_name)
        else:
            # JSON doesn't have log patterns
            return {}

    def extract_metrics_from_log(self, log_text, program_name=None):
        """
        Extract metrics from a log file.

        Args:
            log_text: Log file content
            program_name: Optional program name (auto-detected if not provided)

        Returns:
            dict: Extracted metrics
        """
        if not log_text:
            return {}

        # Auto-detect program if not provided
        if not program_name:
            program_name = self._detect_program(log_text)

        metrics = {"program": program_name}

        # Get patterns for this program
        patterns = self.get_log_patterns(program_name)

        for metric_name, pattern_def in patterns.items():
            if isinstance(pattern_def, dict):
                pattern = pattern_def.get("pattern", "")
                value_type = pattern_def.get("type", "float")
                extract_method = pattern_def.get("extract", "first")
            else:
                pattern = pattern_def
                value_type = "float"
                extract_method = "first"

            if not pattern:
                continue

            matches = re.findall(pattern, log_text, re.IGNORECASE)
            if not matches:
                continue

            # Get appropriate match
            raw_value = matches[-1] if extract_method == "last" else matches[0]

            # Convert type
            if value_type == "float":
                try:
                    metrics[metric_name] = float(raw_value)
                except (ValueError, TypeError):
                    pass
            elif value_type == "boolean":
                metrics[metric_name] = True
            else:
                metrics[metric_name] = str(raw_value).strip()

        return metrics

    def _detect_program(self, log_text):
        """
        Detect which program generated a log.

        Delegates to phenix_ai.log_parsers.detect_program() for consistency.
        """
        try:
            from phenix_ai.log_parsers import detect_program
        except ImportError:
            # Fallback for standalone testing
            try:
                from libtbx.phenix_ai.log_parsers import detect_program
            except ImportError:
                # Last resort - define inline
                def detect_program(log_text):
                    if not log_text:
                        return None
                    log_lower = log_text.lower()
                    if "phenix.refine" in log_lower:
                        return "phenix.refine"
                    elif "phenix.phaser" in log_lower:
                        return "phenix.phaser"
                    return None
        return detect_program(log_text)

    # =========================================================================
    # PROMPT GENERATION
    # =========================================================================

    def format_programs_for_prompt(self, program_names):
        """
        Format program descriptions for LLM prompt.

        Args:
            program_names: List of program names

        Returns:
            str: Formatted text
        """
        lines = []
        for name in program_names:
            desc = self.get_description(name)
            if desc:
                lines.append("- %s: %s" % (name, desc))
            else:
                lines.append("- %s" % name)
        return "\n".join(lines)

    def format_program_hints(self, program_name):
        """
        Format hints for a specific program.

        Args:
            program_name: Name of program

        Returns:
            str: Formatted hints
        """
        hints = self.get_hints(program_name)
        if not hints:
            return ""

        lines = ["Hints for %s:" % program_name]
        for hint in hints:
            lines.append("  - %s" % hint)
        return "\n".join(lines)


# =============================================================================
# SINGLETON INSTANCE
# =============================================================================

_registry = None

def get_registry():
    """Get the global ProgramRegistry instance."""
    global _registry
    if _registry is None:
        _registry = ProgramRegistry()
    return _registry


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def get_program_info(program_name):
    """Get program definition."""
    return get_registry().get_program(program_name)

def get_all_program_names():
    """Get all program names."""
    return get_registry().get_all_programs()

def get_program_desc(program_name):
    """Get program description."""
    return get_registry().get_description(program_name)

def build_program_command(program_name, files, strategy=None):
    """Build command for a program."""
    return get_registry().build_command(program_name, files, strategy)


# =============================================================================
# MAIN (for testing)
# =============================================================================

if __name__ == "__main__":
    print("Testing ProgramRegistry...")
    print()

    registry = ProgramRegistry()

    # Test basic info
    print("All programs (%d):" % len(registry.get_all_programs()))
    for name in sorted(registry.get_all_programs())[:5]:
        print("  - %s: %s" % (name, registry.get_description(name)[:50]))
    print()

    # Test experiment filtering
    print("X-ray programs:")
    for name in registry.get_programs_for_experiment("xray")[:5]:
        print("  - %s" % name)
    print()

    # Test input requirements
    print("phenix.refine inputs:")
    print("  Required:", registry.get_required_inputs("phenix.refine"))
    print("  Optional:", registry.get_optional_inputs("phenix.refine"))
    print()

    # Test command building
    print("Command building test:")
    files = {"model": "test.pdb", "data_mtz": "test.mtz"}
    try:
        cmd = registry.build_command("phenix.refine", files)
        print("  Command:", cmd)
    except Exception as e:
        print("  Error:", e)
    print()

    # Test hints
    print("phenix.refine hints:")
    for hint in registry.get_hints("phenix.refine"):
        print("  -", hint)
