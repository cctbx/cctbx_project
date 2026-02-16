"""
Command Builder - Unified Command Generation for PHENIX AI Agent.

This module consolidates all command generation logic into a single entry point.

Usage:
    from libtbx.langchain.agent.command_builder import CommandBuilder, CommandContext

    builder = CommandBuilder()
    context = CommandContext.from_state(state)
    command = builder.build("phenix.refine", available_files, context)

Pipeline:
    1. _select_files()     - Choose input files based on priorities
    2. _build_strategy()   - Build strategy flags from context/LLM hints
    3. _apply_invariants() - Apply program-specific rules and auto-fills
    4. _assemble_command() - Build final command string
"""

from __future__ import absolute_import, division, print_function

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple

# Centralized pattern utilities
from libtbx.langchain.agent.pattern_manager import extract_cycle_number, extract_all_numbers

# Silence unused import warning - Tuple is used in type hints
assert Tuple is not None

# Import registry for YAML access
from libtbx.langchain.agent.program_registry import ProgramRegistry


# =============================================================================
# COMMAND CONTEXT
# =============================================================================

@dataclass
class CommandContext:
    """
    All context needed for command generation.

    This consolidates all the scattered context that was previously passed
    as separate parameters (best_files, rfree_mtz, context dict, etc.)
    """

    # Cycle info
    cycle_number: int = 1

    # Experiment info
    experiment_type: str = ""  # 'xray' or 'cryoem'
    resolution: Optional[float] = None

    # File tracking
    best_files: Dict[str, str] = field(default_factory=dict)
    rfree_mtz: Optional[str] = None  # Locked R-free MTZ (X-ray only)
    rfree_resolution: Optional[float] = None  # Resolution at which R-free flags were generated (if limited)
    categorized_files: Dict[str, List[str]] = field(default_factory=dict)

    # Workflow info
    workflow_state: str = ""
    history: List[Dict] = field(default_factory=list)

    # LLM suggestions (optional - used when LLM provides hints)
    llm_files: Optional[Dict[str, Any]] = None  # {slot: path} from LLM
    llm_strategy: Optional[Dict[str, Any]] = None  # {flag: value} from LLM

    # Error recovery strategies (file-keyed)
    # Format: {file_path: {flags: {...}, program_scope: [...], reason: "...",
    #                       selected_label: "...", selected_label_pair: "..."}}
    recovery_strategies: Dict[str, Dict] = field(default_factory=dict)

    # User directives (program_settings, constraints, etc.)
    directives: Dict[str, Any] = field(default_factory=dict)

    # Logging callback
    log: Any = None  # Callable for logging, or None

    @classmethod
    def from_state(cls, state: dict) -> 'CommandContext':
        """
        Create context from graph state dict.

        Args:
            state: Graph state dictionary containing session_info, workflow_state, etc.

        Returns:
            CommandContext instance
        """
        session_info = state.get("session_info", {})
        workflow_state = state.get("workflow_state", {})

        # Get resolution from multiple sources
        resolution = (
            state.get("session_resolution") or
            state.get("resolution") or
            workflow_state.get("resolution")
        )

        return cls(
            cycle_number=state.get("cycle_number", 1),
            experiment_type=session_info.get("experiment_type", ""),
            resolution=resolution,
            best_files=session_info.get("best_files", {}),
            rfree_mtz=session_info.get("rfree_mtz"),
            rfree_resolution=session_info.get("rfree_resolution"),
            categorized_files=workflow_state.get("categorized_files", {}),
            workflow_state=workflow_state.get("state", ""),
            history=state.get("history", []),
            llm_files=state.get("corrected_files"),  # From LLM after path correction
            llm_strategy=state.get("strategy"),  # From LLM
            recovery_strategies=session_info.get("recovery_strategies", {}),
            directives=state.get("directives", {}),
        )

    def _log(self, msg: str):
        """Log a message if logging is enabled."""
        if self.log:
            self.log(msg)


# =============================================================================
# COMMAND BUILDER
# =============================================================================

class CommandBuilder:
    """
    Single entry point for all command generation.

    This replaces the fragmented system of:
    - template_builder.build_command()
    - template_builder.build_command_for_program()
    - program_registry.build_command()
    - graph_nodes file selection logic

    Usage:
        builder = CommandBuilder()
        context = CommandContext.from_state(state)
        command = builder.build("phenix.refine", available_files, context)

        # Get file selection details (for event logging)
        selections = builder.get_selection_details()
    """

    # Map input slot names to best_files categories
    SLOT_TO_BEST_CATEGORY = {
        "model": "model",
        "pdb_file": "model",
        "map": "map",
        "full_map": "map",
        "data_mtz": "data_mtz",
        "map_coeffs_mtz": "map_coeffs_mtz",
        "hkl_file": "data_mtz",
        "data": "data_mtz",
        "sequence": "sequence",
        "seq_file": "sequence",
        "ligand_cif": "ligand_cif",
        "ligand": "ligand_cif",
    }

    # Slot name aliases - maps common LLM names to canonical program input names
    # When the LLM uses "data", we should match it to the "data_mtz" input slot
    SLOT_ALIASES = {
        "data": "data_mtz",
        "hkl_file": "data_mtz",
        "mtz": "data_mtz",          # Old name maps to new
        "pdb": "model",
        "pdb_file": "model",
        "seq_file": "sequence",
        "ligand": "ligand_cif",
        "ligand_file": "ligand",
        "map_file": "map",
        "full_map": "map",
    }

    # Parameter name each program uses for specifying data labels from an MTZ file.
    # When a recovery strategy identifies the correct label (e.g., "FTOXD3"),
    # we need to know which parameter to set for the current program.
    #
    # All PHENIX programs follow the pattern: <scope>.input.xray_data.obs_labels
    # The parameter always takes just the main label name (e.g., "FTOXD3"),
    # not the F,SIGF pair. PHENIX finds the sigma column automatically.
    #
    # Canonical reference: knowledge/recoverable_errors.yaml data_label_parameters
    DATA_LABEL_PARAMETERS = {
        "phenix.xtriage": {
            "parameter": "scaling.input.xray_data.obs_labels",
        },
        "phenix.refine": {
            "parameter": "miller_array.labels.name",
        },
        "phenix.phaser": {
            "parameter": "phaser.hklin.labin",
        },
        # NOTE: phenix.autosol intentionally excluded - it handles labels internally
        "phenix.autobuild": {
            "parameter": "autobuild.input.xray_data.obs_labels",
        },
        "phenix.maps": {
            "parameter": "maps.input.xray_data.obs_labels",
        },
    }
    DATA_LABEL_DEFAULT_PARAMETER = "obs_labels"

    def __init__(self):
        """Initialize builder with program registry."""
        self._registry = ProgramRegistry()
        # Track file selection details for event logging
        self._selection_details = {}
        # Track prerequisite program needed (set when build blocked by missing input)
        self._prerequisite_program = None
        self._prerequisite_reason = None

    def get_prerequisite(self):
        """
        Get prerequisite program needed before the requested program can run.

        Returns:
            tuple: (program_name, reason) or (None, None) if no prerequisite needed
        """
        return self._prerequisite_program, self._prerequisite_reason

    def get_selection_details(self):
        """
        Get details about file selections made during last build().

        Returns:
            Dict mapping input slot names to selection details:
            {
                "model": {"selected": "/path/to/file.pdb", "reason": "best_files", ...},
                "data_mtz": {"selected": "/path/to/data.mtz", "reason": "rfree_locked", ...},
            }
        """
        return self._selection_details

    def _record_selection(self, slot, filepath, reason, **kwargs):
        """Record a file selection for later retrieval."""
        self._selection_details[slot] = {
            "selected": os.path.basename(filepath) if filepath else None,
            "path": filepath,
            "reason": reason,
            **kwargs
        }

    def build(self, program: str, available_files: List[str],
              context: CommandContext) -> Optional[str]:
        """
        Build a command for the given program.

        This is the SINGLE ENTRY POINT for all command generation.

        Args:
            program: Program name (e.g., "phenix.refine")
            available_files: List of available file paths
            context: CommandContext with all needed information

        Returns:
            Command string, or None if command cannot be built
        """
        # Clear previous selection details and prerequisite
        self._selection_details = {}
        self._prerequisite_program = None
        self._prerequisite_reason = None

        self._log(context, "BUILD: Starting command generation for %s" % program)

        # 1. Select files
        files = self._select_files(program, available_files, context)
        if files is None:
            self._log(context, "BUILD: Failed to select required files")
            return None

        # 2. Build strategy (and track sources)
        strategy = self._build_strategy(program, context)
        strategy_sources = {}
        if context.llm_strategy:
            for k in context.llm_strategy:
                if k in strategy:
                    strategy_sources[k] = "llm_strategy"
        # Track auto-generated keys
        for k in strategy:
            if k not in strategy_sources:
                strategy_sources[k] = "auto_generated"

        # 2.5. Apply file-specific recovery strategies
        # This adds flags from error recovery (e.g., obs_labels for ambiguous data)
        pre_recovery_keys = set(strategy.keys())
        strategy = self._apply_recovery_strategies(program, files, strategy, context)
        for k in strategy:
            if k not in pre_recovery_keys:
                strategy_sources[k] = "recovery"

        # 3. Apply invariants (may modify files and strategy)
        pre_invariant_keys = set(strategy.keys())
        files, strategy = self._apply_invariants(program, files, strategy, context)
        for k in strategy:
            if k not in pre_invariant_keys:
                strategy_sources[k] = "invariant"

        # 3.5 Check required strategies (e.g., resolution for map_to_model)
        # If an invariant declares has_strategy and auto-fill couldn't provide it,
        # the command cannot succeed — block it now rather than let it fail at runtime
        # Exception: invariants with prerequisite=X signal that program X should run first
        invariants = self._registry.get_invariants(program)
        for inv in invariants:
            required_key = inv.get("check", {}).get("has_strategy")
            if required_key and required_key not in strategy:
                msg = inv.get("message", "%s required" % required_key)
                prerequisite = inv.get("prerequisite")
                if prerequisite:
                    # This invariant can be satisfied by running a prerequisite program first
                    self._log(context, "BUILD: BLOCKED - %s (prerequisite: %s)" % (msg, prerequisite))
                    self._prerequisite_program = prerequisite
                    self._prerequisite_reason = msg
                    return None
                else:
                    self._log(context, "BUILD: BLOCKED - %s (not in strategy or context)" % msg)
                    return None

        # Build file provenance dict from _selection_details
        file_sources = {}
        for slot in (files or {}):
            detail = self._selection_details.get(slot, {})
            file_sources[slot] = detail.get("reason", "unknown")

        # Log final selections
        if files:
            files_summary = ", ".join("%s=%s" % (k, os.path.basename(str(v)))
                                      for k, v in files.items())
            self._log(context, "BUILD: Final files: {%s}" % files_summary)

        # 4. Assemble final command (provenance is logged inside registry.build_command)
        command = self._assemble_command(program, files, strategy,
                                         context=context,
                                         file_sources=file_sources,
                                         strategy_sources=strategy_sources)

        # Post-assembly validation: reject commands with no input files
        # (e.g., "phenix.mtriage" with no map is useless and may hang)
        if command:
            has_any_file = bool(files) and any(v for v in files.values() if v)
            if not has_any_file:
                self._log(context, "BUILD: Rejected command with no input files: %s" % command[:80])
                return None

            # Also reject commands that are just the program name with no arguments
            # (can happen when optional file placeholders are all stripped)
            cmd_parts = command.strip().split()
            if len(cmd_parts) < 2:
                self._log(context, "BUILD: Rejected bare command '%s' (no arguments)" % command)
                return None

        if command:
            self._log(context, "BUILD: Command = %s" % command[:100])

        # 5. Log any LLM file requests that were truly rejected
        if context.llm_files:
            files_dict = files or {}
            rejected = []
            for s, p in context.llm_files.items():
                if not p:
                    continue
                # Check original slot name AND its canonical alias
                canonical = self.SLOT_ALIASES.get(s, s)
                if s not in files_dict and canonical not in files_dict:
                    # Also check fuzzy matches (e.g., "map" matched to "full_map")
                    fuzzy_matched = any(s in k or k in s for k in files_dict)
                    if not fuzzy_matched:
                        rejected.append("%s=%s" % (s, os.path.basename(str(p))))
            if rejected:
                self._log(context, "BUILD: REJECTED LLM file requests: %s" % ", ".join(rejected))

        return command

    def _log(self, context: CommandContext, msg: str):
        """Log a message through context."""
        context._log(msg)

    # =========================================================================
    # STEP 1: FILE SELECTION
    # =========================================================================

    def _select_files(self, program: str, available_files: List[str],
                      context: CommandContext) -> Optional[Dict[str, Any]]:
        """
        Select input files for the program.

        Priority order:
        1. LLM hints (if provided and valid)
        2. Locked rfree_mtz (for MTZ slots in X-ray refinement)
        3. Best files (from BestFilesTracker)
        4. Category-based selection (from YAML input_priorities)
        5. Extension-based fallback

        Args:
            program: Program name
            available_files: List of available file paths
            context: CommandContext

        Returns:
            Dict of {slot_name: file_path}, or None if required files missing
        """
        # Build lookup tables
        basename_to_path = {os.path.basename(f): f for f in available_files}

        # Get program definition
        prog_def = self._registry.get_program(program)
        if not prog_def:
            self._log(context, "BUILD: Unknown program %s" % program)
            return None

        inputs = prog_def.get("inputs", {})
        required_inputs = list(inputs.get("required", {}).keys())
        optional_inputs = list(inputs.get("optional", {}).keys())
        all_inputs = required_inputs + optional_inputs

        selected_files = {}

        # For refinement programs, best_files should take priority over LLM suggestions
        # This prevents the LLM from selecting the original input model instead of
        # the latest refined model
        is_refinement_program = program in ("phenix.refine", "phenix.real_space_refine")

        # Pre-populate with best_files for refinement (can be overridden by explicit LLM choices)
        if is_refinement_program and context.best_files:
            # Model slot
            if context.best_files.get("model") and os.path.exists(context.best_files["model"]):
                for slot in ["model", "pdb", "pdb_file"]:
                    if slot in all_inputs:
                        selected_files[slot] = context.best_files["model"]
                        self._record_selection(slot, context.best_files["model"], "best_files")
                        self._log(context, "BUILD: Using best_model for %s: %s" % (
                            slot, os.path.basename(context.best_files["model"])))
                        break

        # Process LLM hints (if any)
        if context.llm_files:
            llm_summary = ", ".join("%s=%s" % (k, os.path.basename(str(v)) if v else "None")
                                    for k, v in context.llm_files.items())
            self._log(context, "BUILD: LLM requested files: {%s}" % llm_summary)

            for llm_slot, filepath in context.llm_files.items():
                if not filepath:
                    continue

                # Map LLM slot name to canonical input slot name
                # e.g., "data" -> "mtz", "pdb" -> "model"
                canonical_slot = self.SLOT_ALIASES.get(llm_slot, llm_slot)

                # Check if canonical slot exists in program inputs
                if canonical_slot not in all_inputs:
                    # Try the original slot name as fallback
                    if llm_slot in all_inputs:
                        canonical_slot = llm_slot
                    else:
                        # Fuzzy match: try input names containing or contained by the LLM slot
                        # e.g., LLM says "map" -> match "full_map"; LLM says "model" -> match "partial_model"
                        fuzzy_match = None
                        for inp_name in all_inputs:
                            # Match "map" to "full_map", "half_map", etc.
                            if (llm_slot in inp_name or inp_name in llm_slot) and llm_slot != inp_name:
                                # Skip slots that explicitly disable auto-fill
                                inp_def = (inputs.get("required", {}).get(inp_name) or
                                          inputs.get("optional", {}).get(inp_name) or {})
                                if inp_def.get("auto_fill") is False:
                                    self._log(context, "BUILD: Fuzzy match '%s'->'%s' skipped (auto_fill=false)" % (
                                        llm_slot, inp_name))
                                    continue
                                fuzzy_match = inp_name
                                break
                        if fuzzy_match:
                            self._log(context, "BUILD: Fuzzy-matched LLM slot '%s' to '%s'" % (
                                llm_slot, fuzzy_match))
                            canonical_slot = fuzzy_match
                        else:
                            self._log(context, "BUILD: LLM slot '%s' not found in program inputs (tried alias '%s')" % (
                                llm_slot, canonical_slot))
                            continue

                # Validate and correct path
                corrected = self._correct_file_path(filepath, basename_to_path, available_files)
                if corrected:
                    # Handle list case (multiple files)
                    corrected_str = corrected[0] if isinstance(corrected, list) else corrected

                    # Validate extension matches slot requirements
                    slot_def = inputs.get("required", {}).get(canonical_slot) or inputs.get("optional", {}).get(canonical_slot)
                    if slot_def:
                        slot_extensions = slot_def.get("extensions", [])
                        if slot_extensions:
                            if not any(corrected_str.lower().endswith(ext) for ext in slot_extensions):
                                self._log(context, "BUILD: LLM file rejected (wrong extension): %s=%s (expected %s)" % (
                                    canonical_slot, os.path.basename(corrected_str), slot_extensions))
                                continue  # Skip this file

                    # For refinement, check if LLM is trying to use a non-best model
                    if is_refinement_program and canonical_slot in ("model", "pdb", "pdb_file"):
                        best_model = context.best_files.get("model")
                        if best_model and corrected_str != best_model and os.path.exists(best_model):
                            # LLM chose a different model - log and keep best_model
                            self._log(context, "BUILD: LLM suggested %s=%s but using best_model=%s instead" % (
                                canonical_slot, os.path.basename(corrected_str), os.path.basename(best_model)))
                            self._record_selection(canonical_slot, best_model, "best_files_override",
                                llm_suggested=os.path.basename(corrected_str))
                            continue  # Skip LLM's choice, keep best_model

                    # Log if we're using an alias
                    if llm_slot != canonical_slot:
                        self._log(context, "BUILD: Mapped LLM slot '%s' to '%s'" % (llm_slot, canonical_slot))

                    # Validate against input_priorities exclude_categories
                    # This prevents the LLM from choosing e.g. PHASER.1.mtz for autosol
                    # when autosol needs original data (not phased MTZ)
                    priorities = self._registry.get_input_priorities(program, canonical_slot)
                    exclude_cats = priorities.get("exclude_categories", [])
                    if exclude_cats and context.categorized_files:
                        corrected_base = os.path.basename(corrected_str)
                        excluded = False
                        for exc in exclude_cats:
                            cat_files = context.categorized_files.get(exc, [])
                            # Check both full path and basename for robustness
                            cat_basenames = [os.path.basename(f) for f in cat_files]
                            if corrected_str in cat_files or corrected_base in cat_basenames:
                                excluded = True
                                self._log(context,
                                    "BUILD: LLM file rejected (in excluded category '%s'): %s=%s" % (
                                        exc, canonical_slot, corrected_base))
                                break
                        if excluded:
                            continue  # Skip this file, let auto-fill find the right one

                    # Universal safety check: reject any file in intermediate categories
                    # These are NEVER valid inputs regardless of program
                    if context.categorized_files:
                        corrected_base = os.path.basename(corrected_str) if isinstance(corrected_str, str) else corrected_str
                        intermediate_cats = [
                            "intermediate", "map_to_model_intermediate",
                            "intermediate_mr", "autobuild_temp",
                            "carryover_temp", "trace_fragment"
                        ]
                        is_intermediate = False
                        for icat in intermediate_cats:
                            cat_files = context.categorized_files.get(icat, [])
                            cat_basenames = [os.path.basename(f) for f in cat_files]
                            if corrected_str in cat_files or corrected_base in cat_basenames:
                                is_intermediate = True
                                self._log(context,
                                    "BUILD: LLM file rejected (intermediate '%s'): %s=%s" % (
                                        icat, canonical_slot, corrected_base))
                                break
                        if is_intermediate:
                            continue  # Skip intermediate file

                    selected_files[canonical_slot] = corrected
                    self._record_selection(canonical_slot, corrected_str, "llm_selected")
                else:
                    self._log(context, "BUILD: LLM file rejected (not found): %s=%s" % (
                        llm_slot, os.path.basename(str(filepath))))
        else:
            self._log(context, "BUILD: No LLM file hints, will auto-select")

        # Auto-fill missing required inputs
        for input_name in required_inputs:
            if input_name in selected_files:
                continue  # Already selected by LLM

            input_def = inputs["required"][input_name]
            file_found = self._find_file_for_slot(
                program, input_name, input_def,
                available_files, context, basename_to_path
            )

            if file_found:
                selected_files[input_name] = file_found
                # Handle both single files and lists of files
                if isinstance(file_found, list):
                    basenames = [os.path.basename(f) for f in file_found]
                    self._log(context, "BUILD: Auto-filled %s=%s" % (
                        input_name, ", ".join(basenames)))
                else:
                    self._log(context, "BUILD: Auto-filled %s=%s" % (
                        input_name, os.path.basename(file_found)))

        # Check all required inputs are filled
        missing = [inp for inp in required_inputs if inp not in selected_files]
        if missing:
            self._log(context, "BUILD: Missing required inputs: %s" % ", ".join(missing))
            return None

        # Auto-fill optional inputs
        for input_name in optional_inputs:
            if input_name in selected_files:
                continue

            input_def = inputs["optional"][input_name]

            # Skip auto-fill if input explicitly disables it (e.g., partial_model)
            if input_def.get("auto_fill") is False:
                continue

            file_found = self._find_file_for_slot(
                program, input_name, input_def,
                available_files, context, basename_to_path
            )

            if file_found:
                selected_files[input_name] = file_found

        # Record any selections not already tracked
        # (Don't overwrite provenance from LLM/best_files/rfree_locked/etc.)
        for slot, filepath in selected_files.items():
            if slot not in self._selection_details:
                if isinstance(filepath, list):
                    self._record_selection(slot, filepath[0] if filepath else None, "auto_selected")
                else:
                    self._record_selection(slot, filepath, "auto_selected")

        # POST-SELECTION VALIDATION: Remove redundant half_map if full_map is present
        # Programs like map_to_model and resolve_cryo_em only need half-maps when
        # no full map is available. Having both causes confusion.
        # EXCEPTION: If the "full_map" is actually a categorized half_map (mis-selected
        # via generic 'map' category fallback), remove it instead and keep the half_maps.
        if "full_map" in selected_files and "half_map" in selected_files:
            full_map_path = selected_files["full_map"]
            if isinstance(full_map_path, list):
                full_map_path = full_map_path[0] if full_map_path else ""
            # Check if this "full_map" is actually a half map
            half_map_files = context.categorized_files.get("half_map", [])
            full_map_files = context.categorized_files.get("full_map", [])
            if full_map_path in half_map_files and full_map_path not in full_map_files:
                # The "full_map" is a mis-selected half map — remove it, keep half_maps
                self._log(context, "BUILD: Removing mis-selected full_map=%s (it's actually a half-map)" %
                          os.path.basename(str(full_map_path)))
                del selected_files["full_map"]
            else:
                # Genuine full_map — remove redundant half_maps
                self._log(context, "BUILD: Removing redundant half_map (full_map already selected)")
                del selected_files["half_map"]

        # Also check for programs with requires_full_map that shouldn't get any half-maps
        prog_def = self._registry.get_program(program)
        if prog_def and prog_def.get("requires_full_map"):
            # Check if any selected file is actually a half-map
            for slot in list(selected_files.keys()):
                filepath = selected_files[slot]
                if isinstance(filepath, str):
                    check_files = [filepath]
                elif isinstance(filepath, list):
                    check_files = filepath
                else:
                    continue
                for f in check_files:
                    if f in context.categorized_files.get("half_map", []):
                        self._log(context, "BUILD: Removing %s=%s (half-map in requires_full_map program)" %
                                  (slot, os.path.basename(f)))
                        del selected_files[slot]
                        break

        # Final check: if program has no required inputs, at least one optional
        # input must be filled (prevents useless commands like "phenix.mtriage" alone)
        if not required_inputs and not selected_files:
            self._log(context, "BUILD: No files selected for %s (all inputs optional but none available)" % program)
            return None

        return selected_files

    def _correct_file_path(self, filepath: Any, basename_to_path: Dict[str, str],
                           available_files: List[str]) -> Optional[str]:
        """Correct/validate a file path from LLM."""
        if isinstance(filepath, list):
            # Handle multiple files
            corrected = []
            for fp in filepath:
                c = self._correct_single_path(fp, basename_to_path, available_files)
                if c:
                    corrected.append(c)
            return corrected if corrected else None
        else:
            return self._correct_single_path(filepath, basename_to_path, available_files)

    def _correct_single_path(self, filepath: str, basename_to_path: Dict[str, str],
                             available_files: List[str]) -> Optional[str]:
        """Correct a single file path."""
        if not filepath:
            return None
        basename = os.path.basename(filepath)
        if basename in basename_to_path:
            return basename_to_path[basename]
        if filepath in available_files:
            return filepath
        return None

    def _find_file_for_slot(self, program: str, input_name: str, input_def: dict,
                            available_files: List[str], context: CommandContext,
                            basename_to_path: Dict[str, str]) -> Optional[str]:
        """
        Find the best file for an input slot.

        Priority:
        1. Locked rfree_mtz (for MTZ in X-ray refinement)
        2. Best files (from BestFilesTracker)
        3. Category-based selection (input_priorities)
        4. Extension/pattern fallback
        """
        extensions = input_def.get("extensions", [])
        is_multiple = input_def.get("multiple", False)

        # PRIORITY 1: Locked R-free data_mtz (X-ray refinement only)
        # Some programs (e.g., autosol) need original data, not rfree-locked MTZ
        priorities = self._registry.get_input_priorities(program, input_name)
        skip_rfree = priorities.get("skip_rfree_lock", False)
        if not skip_rfree and input_name in ("data_mtz", "hkl_file", "data") and context.rfree_mtz:
            if os.path.exists(context.rfree_mtz):
                self._log(context, "BUILD: Using LOCKED R-free data_mtz for %s" % input_name)
                self._record_selection(input_name, context.rfree_mtz, "rfree_locked")
                return context.rfree_mtz

        # Check if input_priorities specifies a specific subcategory
        # If so, skip generic best_files and use category-based selection instead
        priority_categories = priorities.get("categories", [])

        # Subcategories that are more specific than generic best_files
        # These indicate the program needs a specific type of file, not just any "best" file
        specific_subcategories = {"original_data_mtz", "phased_data_mtz", "half_map", "full_map",
                                  "optimized_full_map", "intermediate_map",
                                  "ligand_fit_output", "with_ligand", "processed_predicted",
                                  "refine_map_coeffs", "denmod_map_coeffs", "predict_build_map_coeffs"}
        uses_specific_subcategory = any(cat in specific_subcategories for cat in priority_categories)

        # PRIORITY 2: Best files (but respect exclude_categories and skip if specific subcategory needed)
        if not uses_specific_subcategory:
            best_category = self.SLOT_TO_BEST_CATEGORY.get(input_name, input_name)
            best_path = context.best_files.get(best_category)
            if best_path and os.path.exists(best_path):
                # Verify extension matches
                if any(best_path.lower().endswith(ext) for ext in extensions):
                    # Check if best file is in an excluded category
                    exclude_categories = priorities.get("exclude_categories", [])
                    excluded = any(best_path in context.categorized_files.get(exc, [])
                                   for exc in exclude_categories)
                    if not excluded:
                        self._log(context, "BUILD: Using best_%s for %s" % (best_category, input_name))
                        self._record_selection(input_name, best_path, "best_files")
                        return best_path
                    else:
                        self._log(context, "BUILD: Skipping best_%s (in excluded category)" % best_category)
        else:
            self._log(context, "BUILD: Skipping best_files for %s (program needs specific subcategory: %s)" %
                     (input_name, priority_categories))

        # PRIORITY 2.5: Files with recovery strategies
        # When a recoverable error identified the correct labels for a file,
        # prefer that file over others (prevents auto-fill from picking a
        # different MTZ during forced retry).
        if context.recovery_strategies and extensions:
            for recovery_file in context.recovery_strategies:
                recovery_base = os.path.basename(recovery_file)
                # Check if this recovery file is available and has the right extension
                for avail in available_files:
                    if (os.path.basename(avail) == recovery_base and
                            any(avail.lower().endswith(ext) for ext in extensions)):
                        self._log(context, "BUILD: Using file with recovery strategy for %s: %s" %
                                 (input_name, recovery_base))
                        self._record_selection(input_name, avail, "recovery_preferred")
                        return avail

        # PRIORITY 3: Category-based selection (from input_priorities)
        #
        # This supports semantic parent categories with subcategory preferences:
        #   categories: [model]  - Request by parent category
        #   prefer_subcategories: [refined, phaser_output]  - Preference order within parent
        #   fallback_categories: [search_model]  - Try these if primary is empty
        #   exclude_categories: [ligand]  - Never use these
        prefer_subcategories = priorities.get("prefer_subcategories", [])
        fallback_categories = priorities.get("fallback_categories", [])
        exclude_categories = priorities.get("exclude_categories", [])

        def is_excluded(f):
            """Check if file is in any excluded category."""
            return any(f in context.categorized_files.get(exc, [])
                      for exc in exclude_categories)

        def find_in_categories(categories, prefer_subs=None):
            """
            Find best file from given categories.
            If prefer_subs is provided, check those subcategories first (in order),
            then fall back to the full category.
            """
            # If subcategory preferences are specified, check those first
            if prefer_subs:
                for subcat in prefer_subs:
                    subcat_files = context.categorized_files.get(subcat, [])
                    valid_files = [f for f in subcat_files
                                  if not is_excluded(f)
                                  and any(f.lower().endswith(ext) for ext in extensions)]
                    if valid_files:
                        self._log(context, "BUILD: Found file in preferred subcategory '%s'" % subcat)
                        return self._get_most_recent_file(valid_files)

            # Then check the parent categories
            for cat in categories:
                cat_files = context.categorized_files.get(cat, [])
                valid_files = [f for f in cat_files
                              if not is_excluded(f)
                              and any(f.lower().endswith(ext) for ext in extensions)]
                if valid_files:
                    self._log(context, "BUILD: Found file in category '%s'" % cat)
                    return self._get_most_recent_file(valid_files)

            return None

        if priority_categories:
            # Try primary categories with subcategory preferences
            result = find_in_categories(priority_categories, prefer_subcategories)
            if result:
                return result

            # Try fallback categories if specified
            if fallback_categories:
                self._log(context, "BUILD: No files in primary categories, trying fallback")
                result = find_in_categories(fallback_categories)
                if result:
                    return result

        # PRIORITY 4: Extension-based fallback with recency preference
        # SPECIAL: For half_map slot, only use files explicitly categorized as half-maps
        # Don't fall back to extension matching - this prevents the same file from being
        # used as both full_map and half_map
        if input_name == "half_map":
            half_map_files = context.categorized_files.get("half_map", [])
            if half_map_files:
                if is_multiple:
                    return half_map_files
                return self._get_most_recent_file(half_map_files)
            return None  # No half-maps available

        # SPECIAL: If program requires a specific subcategory (like refined_mtz),
        # don't fall back to generic extension matching.
        # This prevents ligandfit from using input MTZ when it needs refined MTZ with map coefficients.
        if uses_specific_subcategory:
            self._log(context, "BUILD: No files found for required subcategory %s, skipping extension fallback" %
                     priority_categories)
            return None

        candidates = [f for f in available_files
                      if any(f.lower().endswith(ext) for ext in extensions)]

        # Apply exclude_categories from input_priorities (same as PRIORITY 3)
        # This prevents ligand files from being selected as search models, etc.
        if exclude_categories:
            candidates = [f for f in candidates
                         if not any(f in context.categorized_files.get(exc, [])
                                   for exc in exclude_categories)]

        # Apply exclude patterns from input_def
        exclude_patterns = input_def.get("exclude_patterns", [])
        if exclude_patterns:
            candidates = [f for f in candidates
                          if not any(pat.lower() in os.path.basename(f).lower()
                                    for pat in exclude_patterns)]

        # Apply priority_patterns from input_def (prefer files matching these patterns)
        priority_patterns = input_def.get("priority_patterns", [])
        if priority_patterns and candidates:
            prioritized = [f for f in candidates
                          if any(pat.lower() in os.path.basename(f).lower()
                                for pat in priority_patterns)]
            if prioritized:
                self._log(context, "BUILD: Prioritizing files matching patterns: %s" % priority_patterns)
                candidates = prioritized

        # Apply prefer_patterns from input_def (prefer files matching these patterns)
        prefer_patterns = input_def.get("prefer_patterns", [])
        if prefer_patterns and candidates:
            preferred = [f for f in candidates
                        if any(pat.lower() in os.path.basename(f).lower()
                              for pat in prefer_patterns)]
            if preferred:
                self._log(context, "BUILD: Preferring files matching patterns: %s" % prefer_patterns)
                candidates = preferred

        # NOTE: Half-map vs full-map categorization is handled by _categorize_files()
        # in workflow_state.py. If a single map was promoted to full_map, trust that decision.
        # We no longer filter out half-map-named files here.

        if candidates:
            if is_multiple:
                return candidates
            return self._get_most_recent_file(candidates)

        return None

    def _get_most_recent_file(self, file_list: List[str]) -> Optional[str]:
        """
        Get the most recent file from a list based on cycle numbering.

        Uses centralized cycle extraction from pattern_manager for consistency
        with graph_nodes and other parts of the agent.

        For files with numbered suffixes (e.g., rsr_011_real_space_refined_002.pdb),
        first sorts by primary cycle number, then by all numbers to handle
        within-cycle variants (rsr_001_..._002.pdb > rsr_001_..._001.pdb).
        """
        if not file_list:
            return None
        if len(file_list) == 1:
            return file_list[0]

        def sort_key(path):
            # Primary: cycle number (for cross-cycle comparison)
            cycle = extract_cycle_number(path, default=0)
            # Secondary: all numbers as tuple (for within-cycle comparison)
            all_nums = extract_all_numbers(path)
            try:
                mtime = os.path.getmtime(path) if os.path.exists(path) else 0
            except OSError:
                mtime = 0
            return (cycle, all_nums, mtime)

        sorted_files = sorted(file_list, key=sort_key, reverse=True)
        return sorted_files[0]

    # =========================================================================
    # STEP 2: BUILD STRATEGY
    # =========================================================================

    def _build_strategy(self, program: str, context: CommandContext) -> Dict[str, Any]:
        """
        Build strategy flags from context, LLM hints, AND user directives.

        Priority (highest wins):
        1. User directives (program_settings) - explicit user request
        2. LLM suggestions - from planning step
        3. Auto-generated defaults (output_prefix, etc.)
        """
        strategy = {}

        # Start with LLM suggestions if provided
        if context.llm_strategy:
            strategy.update(context.llm_strategy)

        # CRITICAL: Resolution must come from verified sources (mtriage, xtriage, log analysis),
        # NOT from LLM guessing. The LLM often hallucinate resolution values.
        # Strip LLM-provided resolution if there's no verified source to confirm it.
        if "resolution" in strategy and context.llm_strategy and "resolution" in context.llm_strategy:
            # Check if we have a verified resolution from context
            verified_resolution = context.resolution  # From session_resolution or workflow_state
            if verified_resolution is None:
                # No verified source — LLM made this up
                hallucinated = strategy.pop("resolution")
                self._log(context, "BUILD: Stripped LLM-hallucinated resolution=%.1f "
                          "(no verified source — must run mtriage first)" % float(hallucinated))

        # Override with user directive program_settings (higher priority)
        if context.directives:
            prog_settings = context.directives.get("program_settings", {})

            # Track which normalized keys we've already set (prevents duplicates
            # when both "cycles" and "macro_cycles" map to the same key)
            directive_keys_set = set()

            # Check for program-specific settings first, then default
            for settings_key in [program, "default"]:
                settings = prog_settings.get(settings_key, {})
                if isinstance(settings, dict):
                    for key, value in settings.items():
                        # Map common aliases
                        # "cycles" and "macro_cycles" both mean the same thing
                        normalized_key = key
                        if key == "cycles":
                            # For real_space_refine, "cycles" means "macro_cycles"
                            if program == "phenix.real_space_refine":
                                normalized_key = "macro_cycles"
                            # For phenix.refine, "cycles" maps to strategy_flags "cycles"

                        # Skip if we've already set this normalized key
                        if normalized_key in directive_keys_set:
                            continue
                        directive_keys_set.add(normalized_key)

                        strategy[normalized_key] = value
                        self._log(context, "BUILD: Directive override: %s=%s (from %s)" %
                                  (normalized_key, value, settings_key))

        # Auto-fill output_prefix for refinement programs
        if program in ("phenix.refine", "phenix.real_space_refine"):
            if "output_prefix" not in strategy:
                strategy["output_prefix"] = self._generate_output_prefix(program, context)

        # Sanitize autosol atom_type: LLM often puts multiple atoms in one field
        # e.g. "Se, S" or "Se S" instead of atom_type="Se" + additional_atom_types="S"
        # The space/comma causes command-line splitting: "autosol.atom_type=Se," "S" → crash
        if program == "phenix.autosol" and "atom_type" in strategy:
            raw = str(strategy["atom_type"]).strip()
            # Split on comma, space, semicolon, slash, "and", "or"
            import re as _re
            parts = [p.strip().rstrip(",") for p in _re.split(r'[,;\s/]+|(?:\band\b)|(?:\bor\b)', raw) if p.strip()]
            # Filter to valid atom symbols (1-2 chars, capitalized)
            valid_atoms = []
            for p in parts:
                cleaned = p.strip().rstrip(",").strip()
                if cleaned and 1 <= len(cleaned) <= 2 and cleaned[0].isupper():
                    valid_atoms.append(cleaned)
            if valid_atoms:
                strategy["atom_type"] = valid_atoms[0]
                if len(valid_atoms) > 1:
                    # Move extras to additional_atom_types
                    extras = ",".join(valid_atoms[1:])
                    if "additional_atom_types" not in strategy:
                        strategy["additional_atom_types"] = extras
                        self._log(context, "BUILD: Split atom_type '%s' → atom_type=%s, additional_atom_types=%s" %
                                  (raw, valid_atoms[0], extras))
                    else:
                        self._log(context, "BUILD: Trimmed atom_type '%s' → atom_type=%s (additional_atom_types already set)" %
                                  (raw, valid_atoms[0]))
                elif raw != valid_atoms[0]:
                    self._log(context, "BUILD: Cleaned atom_type '%s' → '%s'" % (raw, valid_atoms[0]))

        return strategy

    def _generate_output_prefix(self, program: str, context: CommandContext) -> str:
        """Generate unique output prefix based on history."""
        if program == "phenix.refine":
            prefix_base = "refine"
            matching_runs = [h for h in context.history
                           if isinstance(h, dict) and
                           h.get("program") == "phenix.refine" and
                           str(h.get("result", "")).startswith("SUCCESS")]
        else:  # phenix.real_space_refine
            prefix_base = "rsr"
            matching_runs = [h for h in context.history
                           if isinstance(h, dict) and
                           h.get("program") == "phenix.real_space_refine" and
                           str(h.get("result", "")).startswith("SUCCESS")]

        run_count = len(matching_runs)
        next_number = run_count + 1

        if matching_runs:
            cycle_nums = [str(h.get("cycle_number", "?")) for h in matching_runs]
            self._log(context, "BUILD: Found %d previous %s runs in cycles: %s" % (
                run_count, program, ", ".join(cycle_nums)))

        prefix = "%s_%03d" % (prefix_base, next_number)
        self._log(context, "BUILD: output_prefix=%s" % prefix)

        return prefix

    # =========================================================================
    # STEP 2.5: APPLY RECOVERY STRATEGIES
    # =========================================================================

    def _apply_recovery_strategies(self, program: str, files: Dict[str, Any],
                                   strategy: Dict[str, Any],
                                   context: CommandContext) -> Dict[str, Any]:
        """
        Merge recovery strategies for files being used.

        When a recoverable error occurred (e.g., ambiguous data labels),
        the ErrorAnalyzer stored a recovery strategy keyed by file path.
        This method checks if any of our selected files have recovery
        strategies and merges their flags into the command strategy.

        For data label recovery, the parameter name is translated per-program
        using DATA_LABEL_PARAMETERS (e.g., xtriage uses
        scaling.input.xray_data.obs_labels, refine uses xray_data.labels).

        Args:
            program: Program name being built
            files: Dict of selected files {slot: path}
            strategy: Current strategy dict
            context: CommandContext with recovery_strategies

        Returns:
            Updated strategy dict with recovery flags merged in
        """
        if not context.recovery_strategies:
            return strategy

        # Debug: log what we're working with
        self._log(context,
                  f"BUILD: Recovery strategies for: {list(context.recovery_strategies.keys())}")
        self._log(context,
                  f"BUILD: Selected files: {files}")

        # Build a robust set of identifiers for selected files
        # Include raw paths, absolute paths, and basenames
        selected_identifiers = set()
        for slot, value in files.items():
            if isinstance(value, str):
                selected_identifiers.add(value)                       # Raw path
                selected_identifiers.add(os.path.abspath(value))      # Absolute
                selected_identifiers.add(os.path.basename(value))     # Basename
            elif isinstance(value, list):
                for v in value:
                    if isinstance(v, str):
                        selected_identifiers.add(v)
                        selected_identifiers.add(os.path.abspath(v))
                        selected_identifiers.add(os.path.basename(v))

        self._log(context, f"BUILD: File identifiers: {selected_identifiers}")

        # Check each recovery strategy
        for recovery_file, recovery in context.recovery_strategies.items():
            # Robust matching: Check if ANY identifier matches
            recovery_abs = os.path.abspath(recovery_file)
            recovery_base = os.path.basename(recovery_file)

            file_match = (
                recovery_file in selected_identifiers or
                recovery_abs in selected_identifiers or
                recovery_base in selected_identifiers
            )

            self._log(context,
                      f"BUILD: Checking {recovery_base}: match={file_match}")

            if file_match:
                # Check program scope (if specified and non-empty)
                scope = recovery.get("program_scope", [])
                if scope and program not in scope:
                    # This recovery is scoped to different programs
                    self._log(context,
                              f"BUILD: Skipping recovery for {recovery_base} "
                              f"(scoped to {scope}, current is {program})")
                    continue

                # If recovery has a selected_label, translate the parameter
                # name for the current program instead of using stored flags
                selected_label = recovery.get("selected_label", "")

                if selected_label:
                    # Use program-specific parameter name
                    # All PHENIX programs just need the main label (e.g., "FTOXD3")
                    # PHENIX finds the sigma column automatically
                    label_config = self.DATA_LABEL_PARAMETERS.get(program)
                    if label_config:
                        param_name = label_config["parameter"]
                    else:
                        # Program not in DATA_LABEL_PARAMETERS — it handles
                        # labels internally (e.g., autosol). Skip label recovery.
                        self._log(context,
                                  f"BUILD: Skipping label recovery for {program} "
                                  f"(not in DATA_LABEL_PARAMETERS, handles labels internally)")
                        continue

                    if param_name in strategy:
                        self._log(context,
                                  f"BUILD: Recovery label {param_name}={selected_label} "
                                  f"overrides existing {strategy[param_name]}")
                    strategy[param_name] = selected_label

                    self._log(context,
                              f"BUILD: Applied data label for {program}: "
                              f"{param_name}={selected_label}")
                else:
                    # No selected_label - fall back to raw flags
                    flags = recovery.get("flags", {})
                    if flags:
                        for key, value in flags.items():
                            if key in strategy:
                                self._log(context,
                                          f"BUILD: Recovery flag {key}={value} "
                                          f"overrides existing {strategy[key]}")
                            strategy[key] = value

                        self._log(context,
                                  f"BUILD: Applied recovery flags for "
                                  f"{recovery_base}: {flags}")

                # Log the reason for transparency
                reason = recovery.get("reason", "")
                if reason:
                    self._log(context, f"BUILD: Recovery reason: {reason}")
            else:
                self._log(context,
                          f"BUILD: No match for {recovery_base} in selected files")

        return strategy

    # =========================================================================
    # STEP 3: APPLY INVARIANTS
    # =========================================================================

    def _apply_invariants(self, program: str, files: Dict[str, Any],
                          strategy: Dict[str, Any],
                          context: CommandContext) -> Tuple[Dict, Dict]:
        """
        Apply program-specific invariants (auto-fills, validations).

        This handles things like:
        - Auto-fill resolution from context
        - Ensure R-free flag generation on first refinement
        - Apply twin law if detected

        Args:
            program: Program name
            files: Selected files
            strategy: Strategy flags
            context: CommandContext

        Returns:
            Tuple of (modified_files, modified_strategy)
        """
        invariants = self._registry.get_invariants(program)

        for inv in invariants:
            name = inv.get("name", "unnamed")
            check = inv.get("check", {})
            fix = inv.get("fix", {})

            # Check file_matches condition (e.g., for denmod MTZ labels)
            file_matches = check.get("file_matches", {})
            if file_matches:
                match_found = False
                for slot, patterns in file_matches.items():
                    selected_file = files.get(slot, "")
                    if selected_file:
                        basename = os.path.basename(selected_file).lower()
                        for pattern in patterns:
                            # Convert glob pattern to simple substring match
                            pattern_clean = pattern.lower().replace("*", "")
                            if pattern_clean in basename:
                                match_found = True
                                break
                    if match_found:
                        break

                if match_found:
                    # Apply the fix
                    set_strategy = fix.get("set_strategy", {})
                    for key, value in set_strategy.items():
                        strategy[key] = value
                        self._log(context, "BUILD: Invariant '%s' - set %s=%s" % (name, key, value))

            # Auto-fill resolution
            if fix.get("auto_fill_resolution") and "resolution" not in strategy:
                if context.resolution:
                    # Round to 1 decimal place
                    strategy["resolution"] = round(context.resolution, 1)
                    self._log(context, "BUILD: Auto-filled resolution=%.1f from context" %
                              context.resolution)
                else:
                    self._log(context, "BUILD: WARNING - resolution required for %s "
                              "but not available in session" % program)

            # Auto-fill output_prefix (if not already set)
            if fix.get("auto_fill_output_prefix") and "output_prefix" not in strategy:
                prefix = fix.get("auto_fill_output_prefix")
                if isinstance(prefix, str):
                    strategy["output_prefix"] = self._generate_output_prefix(program, context)

        # Check for programs that need R-free resolution matching
        # If R-free flags were generated at a limited resolution, certain programs
        # (like polder) MUST use the same resolution limit or they will crash
        program_def = self._registry.get_program(program)
        if program_def:
            fixes = program_def.get("fixes", {})
            if fixes.get("auto_fill_rfree_resolution") and context.rfree_resolution:
                if "high_resolution" not in strategy:
                    strategy["high_resolution"] = context.rfree_resolution
                    self._log(context,
                        "BUILD: R-free flags limited to %.1fÅ - setting high_resolution=%.1f for %s" %
                        (context.rfree_resolution, context.rfree_resolution, program))

            # Auto-fill output_name with format string (e.g., for pdbtools)
            # This allows specifying output filenames like "{protein_base}_with_ligand.pdb"
            if "output_name" in fixes and "output_name" not in strategy:
                output_fix = fixes["output_name"]
                if isinstance(output_fix, dict) and output_fix.get("source") == "auto":
                    format_str = output_fix.get("format", "")
                    if format_str:
                        # Build substitution dict from files
                        subs = {}
                        for slot, filepath in files.items():
                            if filepath:
                                basename = os.path.basename(filepath)
                                name_no_ext = os.path.splitext(basename)[0]
                                subs[slot] = basename
                                subs[slot + "_base"] = name_no_ext
                        try:
                            output_name = format_str.format(**subs)
                            strategy["output_name"] = output_name
                            self._log(context, "BUILD: Auto-set output_name=%s" % output_name)
                        except KeyError as e:
                            self._log(context, "BUILD: Warning - could not format output_name: %s" % e)

        # X-ray specific: R-free flags on first refinement
        if program == "phenix.refine" and context.experiment_type == "xray":
            refine_count = sum(1 for h in context.history
                              if isinstance(h, dict) and
                              h.get("program") == "phenix.refine" and
                              str(h.get("result", "")).startswith("SUCCESS"))
            if refine_count == 0:
                if "generate_rfree_flags" not in strategy:
                    strategy["generate_rfree_flags"] = True
                    self._log(context, "BUILD: First refinement - will generate R-free flags")

        # predict_and_build: resolution is required for building (stop_after_predict=False)
        # If no resolution available, must use stop_after_predict=True (prediction only)
        if program == "phenix.predict_and_build":
            # Check if we're trying to run the full workflow (building, not just prediction)
            wants_full_workflow = strategy.get("stop_after_predict") is False or "stop_after_predict" not in strategy

            if wants_full_workflow and not context.resolution:
                # No resolution available - must use prediction-only mode
                strategy["stop_after_predict"] = True
                self._log(context, "BUILD: No resolution available - setting stop_after_predict=True (prediction only). Run xtriage/mtriage first to get resolution for full workflow.")

        return files, strategy

    # =========================================================================
    # STEP 4: ASSEMBLE COMMAND
    # =========================================================================

    def _assemble_command(self, program: str, files: Dict[str, Any],
                          strategy: Dict[str, Any],
                          context: Optional[CommandContext] = None,
                          file_sources: Optional[Dict[str, str]] = None,
                          strategy_sources: Optional[Dict[str, str]] = None) -> Optional[str]:
        """
        Assemble the final command string.

        Uses the existing registry.build_command() for now, which handles:
        - Template substitution
        - Strategy flag formatting
        - Default values

        Args:
            program: Program name
            files: Selected files
            strategy: Strategy flags
            context: CommandContext (for logging)
            file_sources: Dict mapping slot names to provenance strings
            strategy_sources: Dict mapping strategy keys to provenance strings

        Returns:
            Command string
        """
        log_fn = context._log if context else lambda x: None
        try:
            return self._registry.build_command(
                program, files, strategy,
                log=log_fn,
                file_sources=file_sources,
                strategy_sources=strategy_sources)
        except Exception as e:
            return None


# =============================================================================
# MODULE-LEVEL UTILITIES
# =============================================================================

_builder_instance = None

def get_command_builder() -> CommandBuilder:
    """Get singleton CommandBuilder instance."""
    global _builder_instance
    if _builder_instance is None:
        _builder_instance = CommandBuilder()
    return _builder_instance


# =============================================================================
# TESTS
# =============================================================================

if __name__ == "__main__":
    print("Testing CommandBuilder...")

    # Test CommandContext creation
    state = {
        "cycle_number": 3,
        "session_info": {
            "experiment_type": "cryoem",
            "best_files": {"model": "/path/to/best.pdb"},
            "rfree_mtz": None,
        },
        "workflow_state": {
            "resolution": 2.8,
            "categorized_files": {
                "pdb": ["/path/to/model.pdb"],
                "full_map": ["/path/to/map.ccp4"],
            },
            "state": "cryoem_refined",
        },
        "history": [
            {"cycle_number": 1, "program": "phenix.mtriage", "result": "SUCCESS"},
            {"cycle_number": 2, "program": "phenix.real_space_refine", "result": "SUCCESS: done"},
        ],
        "corrected_files": {"model": "/path/to/model.pdb"},
        "strategy": {"resolution": 2.8},
    }

    context = CommandContext.from_state(state)
    print("  Context created:")
    print("    cycle_number:", context.cycle_number)
    print("    experiment_type:", context.experiment_type)
    print("    resolution:", context.resolution)
    print("    workflow_state:", context.workflow_state)
    print("    llm_files:", context.llm_files)
    print("    llm_strategy:", context.llm_strategy)

    # Test builder instantiation
    builder = get_command_builder()
    print("  Builder created:", type(builder).__name__)

    print("\nPhase 1 implementation complete!")
    print("Next: Phase 2 - Add compatibility layer")
