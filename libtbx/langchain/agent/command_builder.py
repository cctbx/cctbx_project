"""
Command Builder - Unified Command Generation for PHENIX AI Agent.

This module consolidates all command generation logic into a single entry point.

Usage:
    from agent.command_builder import CommandBuilder, CommandContext

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
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple

# Silence unused import warning - Tuple is used in type hints
assert Tuple is not None

# Import registry for YAML access
try:
    from agent.program_registry import ProgramRegistry
except ImportError:
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
    categorized_files: Dict[str, List[str]] = field(default_factory=dict)

    # Workflow info
    workflow_state: str = ""
    history: List[Dict] = field(default_factory=list)

    # LLM suggestions (optional - used when LLM provides hints)
    llm_files: Optional[Dict[str, Any]] = None  # {slot: path} from LLM
    llm_strategy: Optional[Dict[str, Any]] = None  # {flag: value} from LLM

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
            categorized_files=workflow_state.get("categorized_files", {}),
            workflow_state=workflow_state.get("state", ""),
            history=state.get("history", []),
            llm_files=state.get("corrected_files"),  # From LLM after path correction
            llm_strategy=state.get("strategy"),  # From LLM
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
    """

    # Map input slot names to best_files categories
    SLOT_TO_BEST_CATEGORY = {
        "model": "model",
        "pdb_file": "model",
        "map": "map",
        "full_map": "map",
        "mtz": "mtz",
        "hkl_file": "mtz",
        "data": "mtz",
        "map_coefficients": "map_coefficients",
        "sequence": "sequence",
        "seq_file": "sequence",
        "ligand_cif": "ligand_cif",
        "ligand": "ligand_cif",
    }

    def __init__(self):
        """Initialize builder with program registry."""
        self._registry = ProgramRegistry()

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
        self._log(context, "BUILD: Starting command generation for %s" % program)

        # 1. Select files
        files = self._select_files(program, available_files, context)
        if files is None:
            self._log(context, "BUILD: Failed to select required files")
            return None

        # 2. Build strategy
        strategy = self._build_strategy(program, context)

        # 3. Apply invariants (may modify files and strategy)
        files, strategy = self._apply_invariants(program, files, strategy, context)

        # Log final selections
        if files:
            files_summary = ", ".join("%s=%s" % (k, os.path.basename(str(v)))
                                      for k, v in files.items())
            self._log(context, "BUILD: Final files: {%s}" % files_summary)

        # 4. Assemble final command
        command = self._assemble_command(program, files, strategy)

        if command:
            self._log(context, "BUILD: Command = %s" % command[:100])

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

        # Process LLM hints first (if any)
        if context.llm_files:
            llm_summary = ", ".join("%s=%s" % (k, os.path.basename(str(v)) if v else "None")
                                    for k, v in context.llm_files.items())
            self._log(context, "BUILD: LLM requested files: {%s}" % llm_summary)

            for slot, filepath in context.llm_files.items():
                if filepath and slot in all_inputs:
                    # Validate and correct path
                    corrected = self._correct_file_path(filepath, basename_to_path, available_files)
                    if corrected:
                        selected_files[slot] = corrected
                    else:
                        self._log(context, "BUILD: LLM file rejected (not found): %s=%s" % (
                            slot, os.path.basename(str(filepath))))
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
            file_found = self._find_file_for_slot(
                program, input_name, input_def,
                available_files, context, basename_to_path
            )

            if file_found:
                selected_files[input_name] = file_found

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

        # PRIORITY 1: Locked R-free MTZ (X-ray refinement only)
        if input_name in ("mtz", "hkl_file", "data") and context.rfree_mtz:
            if os.path.exists(context.rfree_mtz):
                self._log(context, "BUILD: Using LOCKED R-free MTZ for %s" % input_name)
                return context.rfree_mtz

        # PRIORITY 2: Best files
        best_category = self.SLOT_TO_BEST_CATEGORY.get(input_name, input_name)
        best_path = context.best_files.get(best_category)
        if best_path and os.path.exists(best_path):
            # Verify extension matches
            if any(best_path.lower().endswith(ext) for ext in extensions):
                self._log(context, "BUILD: Using best_%s for %s" % (best_category, input_name))
                return best_path

        # PRIORITY 3: Category-based selection (from input_priorities)
        priorities = self._registry.get_input_priorities(program, input_name)
        priority_categories = priorities.get("categories", [])
        exclude_categories = priorities.get("exclude_categories", [])

        if priority_categories:
            for cat in priority_categories:
                cat_files = context.categorized_files.get(cat, [])
                for f in cat_files:
                    # Check exclusions
                    excluded = any(f in context.categorized_files.get(exc, [])
                                   for exc in exclude_categories)
                    if not excluded and any(f.lower().endswith(ext) for ext in extensions):
                        return self._get_most_recent_file([f] + [
                            cf for cf in cat_files
                            if cf != f and any(cf.lower().endswith(ext) for ext in extensions)
                        ])

        # PRIORITY 4: Extension-based fallback with recency preference
        candidates = [f for f in available_files
                      if any(f.lower().endswith(ext) for ext in extensions)]

        # Apply exclude patterns from input_def
        exclude_patterns = input_def.get("exclude_patterns", [])
        if exclude_patterns:
            candidates = [f for f in candidates
                          if not any(pat.lower() in os.path.basename(f).lower()
                                    for pat in exclude_patterns)]

        if candidates:
            if is_multiple:
                return candidates
            return self._get_most_recent_file(candidates)

        return None

    def _get_most_recent_file(self, file_list: List[str]) -> Optional[str]:
        """
        Get the most recent file from a list based on filename numbering.

        For files with numbered suffixes (e.g., rsr_001_real_space_refined_002.pdb),
        prefers higher numbers.
        """
        if not file_list:
            return None
        if len(file_list) == 1:
            return file_list[0]

        def get_file_number(path):
            basename = os.path.basename(path)
            numbers = re.findall(r'(\d+)', basename)
            if numbers:
                return int(numbers[-1])
            return 0

        def sort_key(path):
            try:
                mtime = os.path.getmtime(path) if os.path.exists(path) else 0
            except OSError:
                mtime = 0
            return (get_file_number(path), mtime)

        sorted_files = sorted(file_list, key=sort_key, reverse=True)
        return sorted_files[0]

    # =========================================================================
    # STEP 2: BUILD STRATEGY
    # =========================================================================

    def _build_strategy(self, program: str, context: CommandContext) -> Dict[str, Any]:
        """
        Build strategy flags from context and LLM hints.

        Args:
            program: Program name
            context: CommandContext

        Returns:
            Dict of strategy flags
        """
        strategy = {}

        # Start with LLM suggestions if provided
        if context.llm_strategy:
            strategy.update(context.llm_strategy)

        # Auto-fill output_prefix for refinement programs
        if program in ("phenix.refine", "phenix.real_space_refine"):
            if "output_prefix" not in strategy:
                strategy["output_prefix"] = self._generate_output_prefix(program, context)

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
            fix = inv.get("fix", {})

            # Auto-fill resolution
            if fix.get("auto_fill_resolution") and "resolution" not in strategy:
                if context.resolution:
                    # Round to 1 decimal place
                    strategy["resolution"] = round(context.resolution, 1)
                    self._log(context, "BUILD: Auto-filled resolution=%.1f from context" %
                              context.resolution)

            # Auto-fill output_prefix (if not already set)
            if fix.get("auto_fill_output_prefix") and "output_prefix" not in strategy:
                prefix = fix.get("auto_fill_output_prefix")
                if isinstance(prefix, str):
                    strategy["output_prefix"] = self._generate_output_prefix(program, context)

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

        return files, strategy

    # =========================================================================
    # STEP 4: ASSEMBLE COMMAND
    # =========================================================================

    def _assemble_command(self, program: str, files: Dict[str, Any],
                          strategy: Dict[str, Any]) -> Optional[str]:
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

        Returns:
            Command string
        """
        try:
            return self._registry.build_command(program, files, strategy)
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
