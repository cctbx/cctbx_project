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

try:
  from libtbx.langchain.agent.file_utils import matches_exclude_pattern
except ImportError:
  from agent.file_utils import matches_exclude_pattern

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

  # Pre-extracted HETATM residues from best model (for polder selection
  # on server where client PDB files are not on disk).
  # Format: [[chain, resseq, resname], ...]
  model_hetatm_residues: Optional[List] = None

  # Whether files are on the local filesystem (False = server mode).
  # When False, skip all disk I/O for client files (content guards,
  # file reading, etc.) — they will always fail.
  files_local: bool = True

  # Logging callback
  log: Any = None  # Callable for logging, or None

  def file_preferences(self) -> Dict[str, Any]:
    """Return directives.file_preferences dict (empty dict if not set)."""
    return self.directives.get("file_preferences", {}) if self.directives else {}

  def excluded_files(self) -> List[str]:
    """Return list of basenames/paths explicitly excluded by the user."""
    return self.file_preferences().get("exclude", [])

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
      model_hetatm_residues=session_info.get("model_hetatm_residues"),
      files_local=state.get("files_local", True),
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

  @staticmethod
  def _best_path(value):
    """
        Safely extract a single path string from a best_files value.

        best_files values are normally strings, but a client may store
        multi-file entries (e.g. half_map) as a list.  All single-file
        slots should use only the first element in that case; multi-file
        slots (half_map) bypass best_files lookup entirely via the
        specific_subcategories guard, so they never reach this helper.

        Returns None if value is falsy or an empty list.
        """
    if not value:
      return None
    if isinstance(value, list):
      return value[0] if value else None
    return value

  def _file_is_available(self, path):
    """
        Check if a file path is known to be available.

        Replaces os.path.exists() for file-availability guards throughout
        the command builder.  This is critical for server execution where
        client file paths are not on the local filesystem, causing
        os.path.exists() to always return False and silently skipping
        best_files, rfree_mtz, and other validated file selections.

        The check uses basename membership in available_files (set at the
        start of build()) plus a fallback to os.path.exists() for paths
        not in the set (e.g., companion files discovered on the client).

        Args:
            path: File path to check

        Returns:
            bool: True if the file is known to be available
        """
    if not path:
      return False
    bn = os.path.basename(str(path))
    if bn in self._available_basenames:
      return True
    # Fallback: disk check (works on client, no-op on server)
    return os.path.exists(str(path))

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
  # Loaded from: knowledge/recoverable_errors.yaml data_label_parameters
  _data_label_cache = None

  @classmethod
  def _get_data_label_parameters(cls):
    """Load data label parameters from YAML (cached)."""
    if cls._data_label_cache is None:
      try:
        try:
          from libtbx.langchain.knowledge.yaml_loader import load_recoverable_errors
        except ImportError:
          from knowledge.yaml_loader import load_recoverable_errors
        config = load_recoverable_errors()
        dlp = config.get("data_label_parameters", {})
        # Extract default parameter without mutating the cached YAML
        cls._data_label_default = dlp.get("default", {}).get("parameter", "obs_labels")
        cls._data_label_cache = {k: v for k, v in dlp.items() if k != "default"}
      except Exception:
        import warnings
        warnings.warn(
          "Could not load data_label_parameters from YAML. "
          "Using hardcoded fallback — this data may be stale.",
          DeprecationWarning, stacklevel=2
        )
        cls._data_label_cache = {
          "phenix.xtriage": {"parameter": "scaling.input.xray_data.obs_labels"},
          "phenix.refine": {"parameter": "miller_array.labels.name"},
          "phenix.phaser": {"parameter": "phaser.hklin.labin"},
          # NOTE: phenix.autobuild does NOT accept obs_labels — removed
          "phenix.maps": {"parameter": "maps.input.xray_data.obs_labels"},
        }
        cls._data_label_default = "obs_labels"
    return cls._data_label_cache

  def __init__(self):
    """Initialize builder with program registry."""
    self._registry = ProgramRegistry()
    # Track file selection details for event logging
    self._selection_details = {}
    # Track prerequisite program needed (set when build blocked by missing input)
    self._prerequisite_program = None
    self._prerequisite_reason = None
    # Basename lookup for _file_is_available() — populated at start of build()
    self._available_basenames = set()

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
    # filepath may be a list for multiple-file slots (e.g. half_map)
    display = None
    if filepath:
      first = filepath[0] if isinstance(filepath, list) else filepath
      display = os.path.basename(first) if first else None
    self._selection_details[slot] = {
      "selected": display,
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

    # Build basename lookup for _file_is_available() — enables server
    # execution where os.path.exists() would always return False.
    self._available_basenames = {os.path.basename(f) for f in available_files if f}
    # Also include best_files and rfree_mtz — these are client-validated
    # paths that may not appear verbatim in available_files.
    if context.best_files:
      for v in context.best_files.values():
        p = self._best_path(v)
        if p:
          self._available_basenames.add(os.path.basename(str(p)))
    if context.rfree_mtz:
      self._available_basenames.add(os.path.basename(str(context.rfree_mtz)))

    self._log(context, "BUILD: Starting command generation for %s" % program)

    # 0.5 Reference model exclusion (Fix C, Tier 1).
    # When the strategy contains reference_model.file=FILENAME, that
    # file is a RESTRAINT source (not a model to refine).  Exclude it
    # from primary model selection so the command builder doesn't pick
    # it as the positional model argument alongside the real model.
    # This prevents "Wrong number of models" from phenix.refine.
    if context.llm_strategy:
      _ref_file = (context.llm_strategy.get("reference_model.file")
                   or context.llm_strategy.get("reference_model_file"))
      if _ref_file:
        _ref_basename = os.path.basename(str(_ref_file))
        # Inject into the exclude list (context is per-build, safe to modify)
        if context.directives is None:
          context.directives = {}
        _fp = context.directives.setdefault("file_preferences", {})
        _excl = _fp.setdefault("exclude", [])
        if _ref_basename not in _excl:
          _excl.append(_ref_basename)
          self._log(context,
            "BUILD: Excluding reference model '%s' from primary model selection"
            % _ref_basename)

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
      # Respect only_for_experiment_type: skip blocking check for invariants
      # that don't apply to the current session type.
      only_for = inv.get("check", {}).get("only_for_experiment_type")
      if only_for:
        session_exp_type = getattr(context, "experiment_type", None)
        if session_exp_type and session_exp_type != only_for:
          continue
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

    # 3.7 Resolve file paths in strategy values.
    # Strategy entries like reference_model.file=4pf4.pdb or
    # secondary_structure.input.file_name=curated_ss.param contain
    # relative filenames that need absolute paths for PHENIX.
    # Detect values that look like filenames (by extension) and resolve
    # via the same basename → full-path lookup used for file selection.
    _FILE_EXTENSIONS = frozenset({
      '.pdb', '.cif', '.mtz', '.params', '.eff',
      '.dat', '.fa', '.fasta', '.seq', '.phil',
      '.ncs_spec', '.param',
    })
    if strategy:
      _basename_to_path = {
        os.path.basename(f): f for f in available_files if f}
      for key in list(strategy.keys()):
        val = str(strategy[key])
        val_lower = val.lower()
        if any(val_lower.endswith(ext) for ext in _FILE_EXTENSIONS):
          # Looks like a filename — try to resolve
          resolved = self._correct_single_path(
            val, _basename_to_path, available_files)
          if resolved:
            strategy[key] = resolved
            self._log(context,
              "BUILD: Resolved strategy path %s=%s"
              % (key, os.path.basename(resolved)))
          elif os.path.isfile(val):
            strategy[key] = os.path.abspath(val)
          elif os.path.isfile(
              os.path.join(
                os.path.dirname(available_files[0])
                if available_files else ".",
                val)):
            strategy[key] = os.path.abspath(
              os.path.join(
                os.path.dirname(available_files[0]),
                val))

    # 4. Assemble final command (provenance is logged inside registry.build_command)
    command = self._assemble_command(program, files, strategy,
                    context=context,
                    file_sources=file_sources,
                    strategy_sources=strategy_sources)

    # Post-assembly: inject recovery-sourced strategy entries that the
    # template-based build_command silently dropped.
    #
    # build_command only emits strategy entries matching strategy_flags
    # keys in programs.yaml.  Recovery-injected params like
    # scaling.input.xray_data.obs_labels are fully-qualified PHIL paths
    # that don't match short strategy_flags names (unit_cell, space_group).
    # They get silently dropped → the recovery fails.
    #
    # Fix: append any recovery-sourced strategy entry whose key (or leaf)
    # is not already present in the command.
    if command and strategy_sources:
      for key, value in strategy.items():
        if strategy_sources.get(key) != "recovery":
          continue
        # Check if already in command (key or leaf name)
        leaf = key.split(".")[-1] if "." in key else key
        if key in command or leaf in command:
          continue
        # Append as key=value
        token = '%s=%s' % (key, value)
        command = command + ' ' + token
        self._log(context,
             "BUILD: Appended recovery param: %s" % token)

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

  def _should_exclude(self, path, exclude_categories, desired_categories,
            categorized_files):
    """Check if a file should be excluded from selection.

        The exclude_categories mechanism prevents using the wrong type of file
        (e.g., raw Fobs MTZ where map coefficients are needed).  However, a
        file can legitimately belong to BOTH an excluded category and a desired
        category due to categorizer bugs or dual-categorization.

        Rule: a file is excluded only if it IS in an excluded category AND is
        NOT in any desired category.  Membership in a desired category is
        positive evidence that overrides the exclusion.

        Example: refine_001_001.mtz may be in both data_mtz (excluded) and
        map_coeffs_mtz (desired).  Without this rule, ligandfit can never
        find its map coefficients when dual-categorization occurs.

        Args:
            path: File path to check
            exclude_categories: List of category names to exclude
            desired_categories: List of category names we're looking for
            categorized_files: Dict of category -> [file paths]

        Returns:
            bool: True if the file should be excluded
        """
    if not exclude_categories:
      return False
    in_excluded = any(path in categorized_files.get(exc, [])
             for exc in exclude_categories)
    if not in_excluded:
      return False
    # File is in an excluded category.  Check if it's also in a desired
    # category — if so, the positive identification wins.
    if desired_categories:
      in_desired = any(path in categorized_files.get(cat, [])
              for cat in desired_categories)
      if in_desired:
        return False  # Desired category overrides exclusion
    return True

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
      # Model slot - verify best_model isn't in an excluded category
      best_model_path = self._best_path(context.best_files.get("model"))
      if best_model_path and self._file_is_available(best_model_path):
        best_model_excluded = False
        # Check against program's exclude_categories for model slot
        priorities = self._registry.get_input_priorities(program, "model")
        if priorities:
          exclude_cats = priorities.get("exclude_categories", [])
          if exclude_cats and context.categorized_files:
            best_bn = os.path.basename(best_model_path)
            for cat in exclude_cats:
              cat_files = context.categorized_files.get(cat, [])
              if any(os.path.basename(f) == best_bn for f in cat_files):
                self._log(context,
                  "BUILD: best_model %s is in excluded category '%s', skipping" % (
                  best_bn, cat))
                best_model_excluded = True
                break

        if not best_model_excluded:
          for slot in ["model", "pdb", "pdb_file"]:
            if slot in all_inputs:
              selected_files[slot] = best_model_path
              self._record_selection(slot, best_model_path, "best_files")
              self._log(context, "BUILD: Using best_model for %s: %s" % (
                slot, os.path.basename(best_model_path)))
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
        # Special handling for ligand slot: phenix.ligandfit
        # accepts either a file path (atp.pdb, atp.cif) or a
        # 3-letter residue code (ATP, HEM, NAG).  The code
        # path is not a file and won't be in available_files.
        # Note: SLOT_ALIASES maps "ligand" → "ligand_cif",
        # so canonical_slot may be either name.
        if canonical_slot in ("ligand", "ligand_cif"):
          filepath_str = str(filepath).strip()
          # Use the program's actual input name, not the alias
          _actual_slot = ("ligand"
            if "ligand" in all_inputs else canonical_slot)
          # Case 1: 2-4 letter residue code (no path sep, no dots
          # or a dot followed by 3-letter ext only)
          is_residue_code = (
            len(filepath_str) <= 4
            and filepath_str.replace("_", "").isalnum()
            and os.sep not in filepath_str
            and "/" not in filepath_str
            and "." not in filepath_str
          )
          if is_residue_code:
            selected_files[_actual_slot] = filepath_str
            self._log(context,
              "BUILD: Accepted ligand residue code: %s"
              % filepath_str)
            self._record_selection(
              _actual_slot, filepath_str,
              "llm_selected_code")
            continue
          # Case 2: file exists on disk but not in available_files
          if os.path.isfile(filepath_str):
            selected_files[_actual_slot] = filepath_str
            self._log(context,
              "BUILD: Accepted ligand file (on disk): %s"
              % os.path.basename(filepath_str))
            self._record_selection(
              _actual_slot, filepath_str,
              "llm_selected_disk")
            continue

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
            best_model = self._best_path(context.best_files.get("model"))
            if best_model and corrected_str != best_model and self._file_is_available(best_model):
              # If the LLM chose a with_ligand file, trust it — it's the
              # model-with-ligand that should be refined after ligandfit.
              # The tracker may still point to the old ligand-free model
              # because pdbtools-combined files have no R-free metrics yet.
              llm_categories = []
              if context.categorized_files:
                llm_bn = os.path.basename(corrected_str)
                for cat, files in context.categorized_files.items():
                  if any(os.path.basename(f) == llm_bn for f in files):
                    llm_categories.append(cat)
              llm_is_with_ligand = any(
                c in ("with_ligand", "ligand_fit_output")
                for c in llm_categories
              )
              # Also recognize pdbtools output (*_modified.pdb) — this is
              # the default naming when combining model + ligand
              if not llm_is_with_ligand:
                llm_bn_lower = os.path.basename(corrected_str).lower()
                if '_modified' in llm_bn_lower and llm_bn_lower.endswith('.pdb'):
                  llm_is_with_ligand = True
              if llm_is_with_ligand:
                self._log(context, "BUILD: LLM chose with_ligand model=%s — trusting LLM (contains ligand)" %
                     os.path.basename(corrected_str))
                # Fall through: LLM's choice will be used as-is
              else:
                # Check if best_model is in an excluded category
                best_model_ok = True
                priorities = self._registry.get_input_priorities(program, "model")
                if priorities:
                  exclude_cats = priorities.get("exclude_categories", [])
                  if exclude_cats and context.categorized_files:
                    best_bn = os.path.basename(best_model)
                    for cat in exclude_cats:
                      if any(os.path.basename(f) == best_bn
                         for f in context.categorized_files.get(cat, [])):
                        best_model_ok = False
                        break

                if best_model_ok:
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

          # Content-based guard: validate PDB file suitability
          # Skip on server where client files aren't on disk.
          if corrected_str.lower().endswith('.pdb') and context.files_local:
            try:
              try:
                from libtbx.langchain.agent.workflow_state import _pdb_is_small_molecule
              except ImportError:
                from agent.workflow_state import _pdb_is_small_molecule
              is_small = _pdb_is_small_molecule(corrected_str)

              # Reject small-molecule PDB from model/protein slots
              if canonical_slot in ("model", "protein", "pdb_file") and is_small:
                self._log(context,
                  "BUILD: LLM file rejected (small molecule in model slot): %s=%s" % (
                    canonical_slot, os.path.basename(corrected_str)))
                continue

              # Reject protein PDB from ligand slot
              if canonical_slot == "ligand" and not is_small:
                try:
                  try:
                    from libtbx.langchain.agent.workflow_state import _pdb_is_protein_model
                  except ImportError:
                    from agent.workflow_state import _pdb_is_protein_model
                  if _pdb_is_protein_model(corrected_str):
                    self._log(context,
                      "BUILD: LLM file rejected (protein model in ligand slot): %s=%s" % (
                        canonical_slot, os.path.basename(corrected_str)))
                    continue
                except Exception:
                  pass
            except Exception:
              pass

          # Apply exclude_patterns from input definition to LLM selections
          # (same check that _find_file_for_slot applies for auto-fill)
          if canonical_slot in all_inputs:
            # Look up input_def from required or optional
            _idef = inputs.get("required", {}).get(canonical_slot,
                inputs.get("optional", {}).get(canonical_slot, {}))
            _excl = _idef.get("exclude_patterns", [])
            if _excl:
              if matches_exclude_pattern(os.path.basename(corrected_str), _excl):
                self._log(context,
                  "BUILD: LLM file rejected (matches exclude_patterns): %s=%s" % (
                    canonical_slot, os.path.basename(corrected_str)))
                continue

          # For inputs with multiple: true (e.g. half_map for resolve_cryo_em),
          # append to a list rather than overwriting.  This handles the case
          # where the LLM assigns half_map_1=file1, half_map_2=file2 — both
          # fuzzy-match to "half_map" and would otherwise overwrite each other.
          slot_def = (inputs.get("required", {}).get(canonical_slot) or
                      inputs.get("optional", {}).get(canonical_slot) or {})
          if slot_def.get("multiple") and canonical_slot in selected_files:
            existing = selected_files[canonical_slot]
            if isinstance(existing, list):
              existing.append(corrected_str)
            else:
              # existing might be a string or a one-element list from
              # _correct_file_path — normalize to a flat string first
              prev = existing[0] if isinstance(existing, list) else existing
              selected_files[canonical_slot] = [prev, corrected_str]
            self._record_selection(canonical_slot, corrected_str, "llm_selected_append")
          else:
            # For multiple:true slots, store the first value as a string
            # (the append branch above will convert to list when the second
            # value arrives).  Use corrected_str, not corrected, to ensure
            # we always store a plain string — never a one-element list.
            if slot_def.get("multiple"):
              selected_files[canonical_slot] = corrected_str
            else:
              selected_files[canonical_slot] = corrected
            self._record_selection(canonical_slot, corrected_str, "llm_selected")
        else:
          self._log(context, "BUILD: LLM file rejected (not found): %s=%s" % (
            llm_slot, os.path.basename(str(filepath))))
    else:
      self._log(context, "BUILD: No LLM file hints, will auto-select")

    # Apply file_preferences from user directives (D1 fix)
    # Lower priority than LLM hints but higher than auto-fill.
    # Keys: model, sequence, data_mtz, map_coeffs_mtz. 'exclude' handled in _find_file_for_slot.
    file_prefs = context.file_preferences() if callable(getattr(context, 'file_preferences', None)) else {}
    PREF_SLOT_MAP = {
      "model":          ["model", "pdb_file", "pdb"],
      "sequence":       ["sequence", "seq_file"],
      "data_mtz":       ["data_mtz", "hkl_file"],
      "map_coeffs_mtz": ["map_coeffs_mtz"],
    }
    for pref_key, candidate_slots in PREF_SLOT_MAP.items():
      pref_val = file_prefs.get(pref_key)
      if not pref_val:
        continue
      for slot in candidate_slots:
        if slot not in all_inputs or slot in selected_files:
          continue
        corrected = self._correct_file_path(pref_val, basename_to_path, available_files)
        if corrected:
          corrected_str = corrected[0] if isinstance(corrected, list) else corrected
          selected_files[slot] = corrected_str
          self._record_selection(slot, corrected_str, "file_preference")
          self._log(context, "BUILD: file_preference %s=%s" % (
            slot, os.path.basename(corrected_str)))
        else:
          self._log(context, "BUILD: file_preference %s=%s not found in available files" % (
            pref_key, os.path.basename(str(pref_val))))
        break  # Try only the first matching slot for this preference key

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
    # Special case: the 'ligand' slot for phenix.ligandfit can be
    # satisfied by a residue code in the strategy dict (e.g.
    # strategy={ligand: "ATP"}) instead of a file path.
    missing = [inp for inp in required_inputs if inp not in selected_files]
    if missing and "ligand" in missing and context.llm_strategy:
      ligand_val = context.llm_strategy.get("ligand", "")
      if ligand_val:
        ligand_str = str(ligand_val).strip()
        # Accept 2-4 letter residue codes or file paths on disk
        is_code = (
          len(ligand_str) <= 4
          and ligand_str.replace("_", "").isalnum()
          and "." not in ligand_str
        )
        if is_code or os.path.isfile(ligand_str):
          selected_files["ligand"] = ligand_str
          self._log(context,
            "BUILD: Accepted ligand from strategy: %s"
            % ligand_str)
          missing.remove("ligand")
    # Also check if ligand_cif alias was filled but
    # "ligand" is the actual required input name
    if missing and "ligand" in missing:
      if "ligand_cif" in selected_files:
        selected_files["ligand"] = selected_files.pop(
          "ligand_cif")
        missing.remove("ligand")
    if missing:
      self._log(context, "BUILD: Missing required inputs: %s" % ", ".join(missing))
      self._last_missing_slots = missing
      return None
    self._last_missing_slots = None

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

    # SUPPLEMENT: For multiple:true slots, the LLM often assigns only one
    # file (e.g. half_map=file2.ccp4).  The slot is now "filled", so
    # auto-fill above skips it.  But programs like resolve_cryo_em and
    # mtriage need ALL matching files (both half-maps).  Backfill from the
    # category to collect any files the LLM missed.
    all_input_defs = {}
    all_input_defs.update(inputs.get("required", {}))
    all_input_defs.update(inputs.get("optional", {}))
    for input_name, input_def in all_input_defs.items():
      if not input_def.get("multiple"):
        continue
      if input_name not in selected_files:
        continue  # truly empty — auto-fill already handled it
      existing = selected_files[input_name]
      existing_list = existing if isinstance(existing, list) else [existing]
      # Use _find_file_for_slot to collect ALL category matches
      all_found = self._find_file_for_slot(
        program, input_name, input_def,
        available_files, context, basename_to_path)
      if not isinstance(all_found, list):
        continue  # single file or None — nothing to supplement
      if len(all_found) <= len(existing_list):
        continue  # no extra files available
      # Merge: add files not already selected.
      # Use os.path.realpath to handle symlinks robustly.
      existing_real = set(os.path.realpath(f) for f in existing_list)
      added = []
      for f in all_found:
        if os.path.realpath(f) not in existing_real:
          existing_list.append(f)
          existing_real.add(os.path.realpath(f))
          added.append(os.path.basename(f))
      if added:
        selected_files[input_name] = existing_list
        self._log(context,
          "BUILD: Supplemented %s with %d file(s): %s"
          % (input_name, len(added), ", ".join(added)))

    # POST-SELECTION VALIDATION: Remove redundant half_map if full_map is present
    # Programs like map_to_model only need half-maps when no full map is available.
    # Having both causes confusion.
    # EXCEPTION 1: If the "full_map" is actually a categorized half_map (mis-selected
    # via generic 'map' category fallback), remove it instead and keep the half_maps.
    # EXCEPTION 2: Programs with keep_half_maps_with_full_map: true (e.g. map_to_model)
    # genuinely need both simultaneously.
    # EXCEPTION 3: Programs with prefers_half_maps: true (e.g. resolve_cryo_em,
    # mtriage, map_sharpening) need half-maps as source truth — drop the full_map.
    if "full_map" in selected_files and "half_map" in selected_files:
      prog_def = self._registry.get_program(program)
      if prog_def and prog_def.get("keep_half_maps_with_full_map"):
        # Program explicitly needs both — don't remove either
        self._log(context, "BUILD: Keeping half_maps alongside full_map (%s uses both)" % program)
      elif prog_def and prog_def.get("prefers_half_maps"):
        # Program needs half-maps as source truth (e.g. resolve_cryo_em)
        # Drop the full_map, keep the half-maps
        self._log(context, "BUILD: Keeping half_maps, removing full_map (%s prefers half-maps)" % program)
        del selected_files["full_map"]
      else:
        full_map_path = selected_files["full_map"]
        if isinstance(full_map_path, list):
          full_map_path = full_map_path[0] if full_map_path else ""
        # Check if this "full_map" is actually a half map.
        # A genuine full map is in full_map OR optimized_full_map category;
        # a mis-selected half map is in half_map and NOT in either full_map category.
        half_map_files = context.categorized_files.get("half_map", [])
        genuine_full_map_files = (
          context.categorized_files.get("full_map", []) +
          context.categorized_files.get("optimized_full_map", [])
        )
        if full_map_path in half_map_files and full_map_path not in genuine_full_map_files:
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

    # POST-SELECTION: Remove redundant map/map_coeffs overlap.
    # Programs like map_correlations accept both map and map_coeffs inputs.
    # When both are available, use the appropriate one for the experiment type:
    #   X-ray  → prefer map_coeffs (MTZ), drop the map
    #   CryoEM → prefer map, drop the map_coeffs
    if "full_map" in selected_files and "data_mtz" in selected_files:
      exp_type = getattr(context, "experiment_type", "") or ""
      if exp_type == "xray":
        self._log(context,
          "BUILD: Dropping full_map (X-ray: prefer map_coeffs over map)")
        del selected_files["full_map"]
      elif exp_type == "cryoem":
        self._log(context,
          "BUILD: Dropping data_mtz (cryo-EM: prefer map over map_coeffs)")
        del selected_files["data_mtz"]

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
      if self._file_is_available(context.rfree_mtz):
        self._log(context, "BUILD: Using LOCKED R-free data_mtz for %s" % input_name)
        self._record_selection(input_name, context.rfree_mtz, "rfree_locked")
        return context.rfree_mtz

    # Check if input_priorities specifies a specific subcategory
    # If so, skip generic best_files and use category-based selection instead
    priority_categories = priorities.get("categories", [])

    # Subcategories that are more specific than generic best_files
    # These indicate the program needs a specific type of file, not just any "best" file
    specific_subcategories = {"original_data_mtz", "phased_data_mtz", "half_map", "full_map",
                 "optimized_full_map",
                 "ligand_fit_output", "with_ligand", "processed_predicted",
                 "refine_map_coeffs", "denmod_map_coeffs", "predict_build_map_coeffs"}
    uses_specific_subcategory = any(cat in specific_subcategories for cat in priority_categories)

    # PRIORITY 2: Best files (but respect exclude_categories and skip if specific subcategory needed)
    if not uses_specific_subcategory:
      best_category = self.SLOT_TO_BEST_CATEGORY.get(input_name, input_name)
      best_path = self._best_path(context.best_files.get(best_category))
      if best_path and self._file_is_available(best_path):
        # Verify extension matches
        if any(best_path.lower().endswith(ext) for ext in extensions):
          # Check if best file is in an excluded category
          exclude_categories = priorities.get("exclude_categories", [])
          excluded = self._should_exclude(
            best_path, exclude_categories, priority_categories,
            context.categorized_files)
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
      """Check if file is in any excluded category OR is user-excluded via file_preferences."""
      if self._should_exclude(f, exclude_categories, priority_categories,
                  context.categorized_files):
        return True
      # D1: user-specified exclude list (basename comparison)
      user_excl = context.excluded_files() if callable(getattr(context, 'excluded_files', None)) else []
      if user_excl:
        fbn = os.path.basename(f).lower()
        if any(os.path.basename(e).lower() == fbn for e in user_excl):
          self._log(context, "BUILD: file_preference exclude removed (cat): %s" % os.path.basename(f))
          return True
      return False

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
      if is_multiple:
        # For multiple-file slots (e.g. half_map for resolve_cryo_em / mtriage),
        # collect ALL valid files from every priority category — not just the
        # most-recent one.  find_in_categories() always returns a single file so
        # it must not be used for is_multiple slots.
        all_valid = []
        seen = set()
        for cat in priority_categories:
          cat_files = context.categorized_files.get(cat, [])
          for f in cat_files:
            if f not in seen and not is_excluded(f) and                                 any(f.lower().endswith(ext) for ext in extensions):
              all_valid.append(f)
              seen.add(f)
        if all_valid:
          self._log(context, "BUILD: Found %d files for multiple slot '%s' in categories %s" % (
            len(all_valid), input_name, priority_categories))
          self._record_selection(input_name, all_valid, "category_multiple")
          return all_valid
      else:
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
    #
    # GUARD: require_best_files_only — some slots (like ligandfit's
    # map_coeffs_mtz) should ONLY be filled from best_files, never from
    # extension matching.  This prevents the raw Fobs MTZ from being used
    # where map coefficients are required.
    if input_def.get("require_best_files_only"):
      # Try best_files one more time (may have been skipped due to
      # uses_specific_subcategory logic above)
      best_category = self.SLOT_TO_BEST_CATEGORY.get(input_name, input_name)
      best_path = self._best_path(context.best_files.get(best_category))
      if best_path and self._file_is_available(best_path):
        if any(best_path.lower().endswith(ext) for ext in extensions):
          self._log(context, "BUILD: require_best_files_only: using best_%s for %s" %
              (best_category, input_name))
          self._record_selection(input_name, best_path, "best_files_required")
          return best_path
      self._log(context, "BUILD: require_best_files_only: no best_%s available for %s, "
          "skipping extension fallback" % (best_category, input_name))
      return None

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
      # FALLBACK: try best_files even for specific-subcategory slots.
      # When category-based lookup fails (categorized_files doesn't contain
      # the file), best_files may still have the correct file — e.g. after
      # our supplemental file discovery populates best_files["map_coeffs_mtz"]
      # from record_result/_rebuild_best_files_from_cycles.
      best_category = self.SLOT_TO_BEST_CATEGORY.get(input_name, input_name)
      best_path = self._best_path(context.best_files.get(best_category))
      if best_path and self._file_is_available(best_path):
        if any(best_path.lower().endswith(ext) for ext in extensions):
          exclude_categories_check = priorities.get("exclude_categories", [])
          excluded = self._should_exclude(
            best_path, exclude_categories_check, priority_categories,
            context.categorized_files)
          if not excluded:
            self._log(context, "BUILD: Using best_%s as fallback for %s "
                "(category lookup found no files)" % (best_category, input_name))
            self._record_selection(input_name, best_path, "best_files_fallback")
            return best_path
          else:
            self._log(context, "BUILD: best_%s is in excluded category, skipping" % best_category)

      self._log(context, "BUILD: No files found for required subcategory %s, skipping extension fallback" %
          priority_categories)
      return None

    candidates = [f for f in available_files
           if any(f.lower().endswith(ext) for ext in extensions)]

    # Apply user file_preferences.exclude (D1 fix)
    # User can say "exclude: [bad_model.pdb, old_data.mtz]" to block specific files
    user_exclude = context.excluded_files() if callable(getattr(context, 'excluded_files', None)) else []
    if user_exclude:
      user_exclude_bns = {os.path.basename(e).lower() for e in user_exclude}
      filtered = [f for f in candidates
            if os.path.basename(f).lower() not in user_exclude_bns]
      if len(filtered) < len(candidates):
        excluded_names = [os.path.basename(f) for f in candidates if f not in filtered]
        self._log(context, "BUILD: file_preference exclude removed: %s" % excluded_names)
      candidates = filtered

    # Apply exclude_categories from input_priorities (same as PRIORITY 3)
    # This prevents ligand files from being selected as search models, etc.
    if exclude_categories:
      candidates = [f for f in candidates
            if not self._should_exclude(
              f, exclude_categories, priority_categories,
              context.categorized_files)]

    # Apply exclude patterns from input_def
    exclude_patterns = input_def.get("exclude_patterns", [])
    if exclude_patterns:
      candidates = [f for f in candidates
             if not matches_exclude_pattern(
               os.path.basename(f), exclude_patterns)]

    # Content-based guard: reject small-molecule PDB files from model slots.
    # This catches files like atp.pdb or gdp.pdb that have no ligand-like
    # name pattern but are HETATM-only and should never be used as the
    # protein model.
    # Skip on server where client files aren't on disk (relies on
    # categorization having classified them correctly via ligand_hints).
    if (input_name in ("model", "protein", "pdb_file") and '.pdb' in extensions
        and context.files_local):
      filtered = []
      for f in candidates:
        if f.lower().endswith('.pdb'):
          try:
            try:
              from libtbx.langchain.agent.workflow_state import _pdb_is_small_molecule
            except ImportError:
              from agent.workflow_state import _pdb_is_small_molecule
            if _pdb_is_small_molecule(f):
              self._log(context, "BUILD: Excluded small-molecule PDB from model slot: %s" %
                  os.path.basename(f))
              continue
          except Exception:
            pass
        filtered.append(f)
      candidates = filtered

    # Inverse guard: reject protein PDB files from ligand slot.
    # The ligand slot should only contain small-molecule files.
    # This prevents the refined model from being used as the ligand.
    # Uses _pdb_is_protein_model (ratio-based protein check) rather than
    # "not _pdb_is_small_molecule" to avoid rejecting unreadable files.
    # NOTE: ratio-based check allows small files with some ATOM records
    # (e.g., nucleotide ligands like ATP, extracted binding-site fragments).
    # Skip on server where client files aren't on disk.
    if input_name == "ligand" and '.pdb' in extensions and context.files_local:
      filtered = []
      for f in candidates:
        if f.lower().endswith('.pdb'):
          try:
            try:
              from libtbx.langchain.agent.workflow_state import _pdb_is_protein_model
            except ImportError:
              from agent.workflow_state import _pdb_is_protein_model
            if _pdb_is_protein_model(f):
              self._log(context, "BUILD: Excluded protein PDB from ligand slot: %s" %
                  os.path.basename(f))
              continue
            else:
              self._log(context, "BUILD: Ligand guard passed for %s" %
                  os.path.basename(f))
          except Exception as e:
            self._log(context, "BUILD: Ligand guard skipped for %s (%s)" %
                (os.path.basename(f), e))
        filtered.append(f)
      candidates = filtered

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
            if matches_exclude_pattern(
              os.path.basename(f), prefer_patterns)]
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
      # Tertiary: path string (deterministic tiebreaker that works
      # on both client and server — server cannot stat client files)
      return (cycle, all_nums, path)

    sorted_files = sorted(file_list, key=sort_key, reverse=True)
    return sorted_files[0]

  # =========================================================================
  # STEP 2: BUILD STRATEGY
  # =========================================================================

  def _build_strategy(
      self, program, context):
    """Build strategy flags from context, LLM hints,
    AND user directives.

    Priority (highest wins):
    1. User directives (program_settings)
    2. LLM suggestions - from planning step
    3. Auto-generated defaults (output_prefix, etc.)
    """
    strategy = {}

    # Start with LLM suggestions if provided
    if context.llm_strategy:
      strategy.update(context.llm_strategy)

    # CRITICAL: Resolution must come from verified
    # sources (mtriage, xtriage, log analysis), NOT
    # from LLM guessing.  Strip LLM-provided resolution
    # if there is no verified source to confirm it.
    if ("resolution" in strategy
        and context.llm_strategy
        and "resolution" in context.llm_strategy):
      verified_resolution = context.resolution
      if verified_resolution is None:
        hallucinated = strategy.pop("resolution")
        self._log(context,
          "BUILD: Stripped LLM-hallucinated "
          "resolution=%.1f (no verified source "
          "-- must run mtriage first)"
          % float(hallucinated))

    # Override with user directive program_settings
    if context.directives:
      prog_settings = context.directives.get(
        "program_settings", {})

      # Track normalized keys already set (prevents
      # duplicates when both "cycles" and
      # "macro_cycles" map to the same key)
      directive_keys_set = set()

      # Program-specific settings first, then default
      for settings_key in [program, "default"]:
        settings = prog_settings.get(
          settings_key, {})
        if isinstance(settings, dict):
          for key, value in settings.items():
            # Map common aliases
            normalized_key = key
            if key == "cycles":
              if program == "phenix.real_space_refine":
                normalized_key = "macro_cycles"

            if normalized_key in directive_keys_set:
              continue
            directive_keys_set.add(normalized_key)

            strategy[normalized_key] = value
            self._log(context,
              "BUILD: Directive override: "
              "%s=%s (from %s)"
              % (normalized_key, value,
                 settings_key))

    # Validate space_group: LLM sometimes extracts
    # workflow descriptions (e.g. "determination")
    # or trailing prose (e.g. "P4 for refinement wit",
    # "C2 from phaser") as space group symbols.
    #
    # Strategy:
    # 1. Try cctbx.sgtbx on the full value (authoritative)
    # 2. If that fails, try progressive prefix shortening
    # 3. If cctbx unavailable, fall back to heuristic check
    if "space_group" in strategy:
      sg_val = str(strategy["space_group"]).strip()
      _sg_resolved = False
      try:
        from cctbx import sgtbx
        # Try the full value first
        try:
          sgi = sgtbx.space_group_info(sg_val)
          _ = sgi.group()
          # Valid — use canonical form
          canonical = str(sgi)
          if canonical != sg_val:
            self._log(context,
              "BUILD: Canonicalized space_group "
              "%r -> %r (cctbx.sgtbx)"
              % (sg_val, canonical))
          strategy["space_group"] = canonical
          _sg_resolved = True
        except (ValueError, RuntimeError):
          # Full value invalid — try prefix shortening
          # "P4 for refinement wit" -> "P4 for" -> "P4"
          # "C2 from phaser" -> "C2 from" -> "C2"
          salvaged = None
          tokens = sg_val.split()
          for i in range(len(tokens) - 1, 0, -1):
            candidate = " ".join(tokens[:i])
            try:
              sgi = sgtbx.space_group_info(candidate)
              _ = sgi.group()
              salvaged = str(sgi)
              break
            except (ValueError, RuntimeError):
              continue
          if salvaged:
            strategy["space_group"] = salvaged
            self._log(context,
              "BUILD: Salvaged space_group=%r "
              "from invalid %r (cctbx.sgtbx)"
              % (salvaged, sg_val))
            _sg_resolved = True
          else:
            del strategy["space_group"]
            self._log(context,
              "BUILD: Stripped invalid "
              "space_group=%r (cctbx rejects "
              "all prefixes)" % sg_val)
            _sg_resolved = True
      except ImportError:
        pass  # cctbx not available

      # Fallback: heuristic check when cctbx unavailable
      if not _sg_resolved:
        try:
          try:
            from libtbx.langchain.agent \
              .command_postprocessor \
              import _is_valid_space_group
          except ImportError:
            from agent.command_postprocessor \
              import _is_valid_space_group
          if not _is_valid_space_group(sg_val):
            del strategy["space_group"]
            self._log(context,
              "BUILD: Stripped invalid "
              "space_group=%r from strategy "
              "(heuristic, no cctbx)" % sg_val)
        except ImportError:
          pass

    # Auto-fill output_prefix for refinement programs
    if program in (
        "phenix.refine",
        "phenix.real_space_refine"):
      if "output_prefix" not in strategy:
        strategy["output_prefix"] = (
          self._generate_output_prefix(
            program, context))

    # Sanitize autosol atom_type: LLM often puts
    # multiple atoms in one field, e.g. "Se, S"
    # instead of atom_type="Se" +
    # additional_atom_types="S".  The space/comma
    # causes command-line splitting.
    if (program == "phenix.autosol"
        and "atom_type" in strategy):
      raw = str(strategy["atom_type"]).strip()
      import re as _re
      parts = [
        p.strip().rstrip(",")
        for p in _re.split(
          r'[,;\s/]+|(?:\band\b)|(?:\bor\b)',
          raw)
        if p.strip()]
      valid_atoms = []
      for p in parts:
        cleaned = p.strip().rstrip(",").strip()
        if (cleaned
            and 1 <= len(cleaned) <= 2
            and cleaned[0].isupper()):
          valid_atoms.append(cleaned)
      if valid_atoms:
        strategy["atom_type"] = valid_atoms[0]
        if len(valid_atoms) > 1:
          extras = ",".join(valid_atoms[1:])
          if "additional_atom_types" not in strategy:
            strategy["additional_atom_types"] = (
              extras)
            self._log(context,
              "BUILD: Split atom_type '%s' -> "
              "atom_type=%s, "
              "additional_atom_types=%s"
              % (raw, valid_atoms[0], extras))
          else:
            self._log(context,
              "BUILD: Trimmed atom_type "
              "'%s' -> atom_type=%s "
              "(additional_atom_types "
              "already set)"
              % (raw, valid_atoms[0]))
        elif raw != valid_atoms[0]:
          self._log(context,
            "BUILD: Cleaned atom_type "
            "'%s' -> '%s'"
            % (raw, valid_atoms[0]))

    # Sanitize autosol wavelength: for MAD the LLM
    # sometimes provides multiple wavelengths (peak,
    # inflection, high-remote).  phenix.autosol only
    # accepts a single autosol.lambda value -- extra
    # values cause duplicate or invalid keywords.
    # Keep only the first (peak) value.
    if (program == "phenix.autosol"
        and "wavelength" in strategy):
      wl = strategy["wavelength"]
      if isinstance(wl, list):
        strategy["wavelength"] = (
          wl[0] if wl else None)
        self._log(context,
          "BUILD: autosol wavelength list "
          "%s -> using first (peak) value %s"
          % (wl, strategy["wavelength"]))
        if strategy["wavelength"] is None:
          del strategy["wavelength"]

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
        using data_label_parameters from recoverable_errors.yaml (e.g., xtriage uses
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
          label_config = self._get_data_label_parameters().get(program)
          if label_config:
            param_name = label_config["parameter"]
          else:
            # Program not in data_label_parameters — it handles
            # labels internally (e.g., autosol). Skip label recovery.
            self._log(context,
                 f"BUILD: Skipping label recovery for {program} "
                 f"(not in data_label_parameters, handles labels internally)")
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

      # Skip this invariant if it is restricted to a specific experiment type
      # and the current session is a different type.
      # Example: requires_resolution for map_correlations only fires for cryoem.
      only_for = check.get("only_for_experiment_type")
      if only_for:
        session_exp_type = getattr(context, "experiment_type", None)
        if session_exp_type and session_exp_type != only_for:
          self._log(context, "BUILD: Invariant '%s' skipped (only_for_experiment_type=%s, "
               "session=%s)" % (name, only_for, session_exp_type))
          continue

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

    # X-ray specific: R-free flags — only generate when no R-free
    # MTZ has been locked in the session.  The rfree_mtz field is set
    # by the client after the first successful refinement.  Using it
    # instead of refine_count avoids two bugs:
    #   1. Input MTZ has pre-existing R-free flags (3tpp-ensemble-refine)
    #      → rfree_mtz is locked from the input, generate=True is skipped
    #   2. First refine succeeds, second refine gets generate=True
    #      (AF_POMGNT2) → rfree_mtz is locked after cycle 1, skipped
    # When rfree_mtz is NOT set, this is either the first refinement
    # or all previous refinements failed — generate=True is correct.
    if program == "phenix.refine" and context.experiment_type == "xray":
      if not context.rfree_mtz:
        if "generate_rfree_flags" not in strategy:
          strategy["generate_rfree_flags"] = True
          self._log(context,
              "BUILD: First refinement - will generate "
              "R-free flags (no rfree_mtz locked)")
      else:
        # R-free MTZ already locked — strip generate if LLM set it
        if strategy.get("generate_rfree_flags"):
          del strategy["generate_rfree_flags"]
          self._log(context,
              "BUILD: Stripped generate_rfree_flags "
              "(rfree_mtz locked: %s)"
              % os.path.basename(str(context.rfree_mtz)))

    # predict_and_build: resolution is required for building (stop_after_predict=False)
    # If no resolution available, must use stop_after_predict=True (prediction only)
    if program == "phenix.predict_and_build":
      # Check if we're trying to run the full workflow (building, not just prediction)
      wants_full_workflow = strategy.get("stop_after_predict") is False or "stop_after_predict" not in strategy

      if wants_full_workflow and not context.resolution:
        # No resolution available - must use prediction-only mode
        strategy["stop_after_predict"] = True
        self._log(context, "BUILD: No resolution available - setting stop_after_predict=True (prediction only). Run xtriage/mtriage first to get resolution for full workflow.")

    # polder: validate/fix the selection string
    # The LLM often writes invalid bare words like selection=ligand or selection=heteroatom.
    # Valid selections are PHENIX atom selection syntax, e.g. "hetero", "chain B",
    # "chain B and resseq 100", "resname ATP", etc.
    if program == "phenix.polder" and "selection" in strategy:
      strategy["selection"] = self._sanitize_polder_selection(
        strategy["selection"], files, context
      )

    return files, strategy

  # -------------------------------------------------------------------------
  # Polder selection helper
  # -------------------------------------------------------------------------

  # Bare words that the LLM writes but are NOT valid PHENIX selection syntax.
  # Any selection matching these will be replaced with an inferred/default value.
  _INVALID_POLDER_SELECTIONS = frozenset([
    "ligand", "lig", "het", "hetatm", "heteroatom", "hetero_atom",
    "ligands", "hets", "small_molecule", "smallmolecule",
  ])

  def _sanitize_polder_selection(self, selection: str,
                 files: Dict[str, Any],
                 context: "CommandContext") -> str:
    """
        Validate and fix a polder selection string.

        Rules (in priority order):
        1. If the user explicitly supplied a selection via directives/project_advice,
           pass it through unchanged (they know what they want).
        2. If selection is a known-invalid bare word (e.g. "ligand"), replace it:
           a. Try to infer a real selection from HETATM records in the model file.
           b. Fall back to "hetero" (valid PHENIX keyword for all non-water HETATM).
        3. Otherwise, pass the LLM selection through unchanged (it may be valid).

        Args:
            selection: The selection string from strategy
            files: Selected input files for this command
            context: CommandContext with directives and available file info

        Returns:
            A valid (or best-effort) PHENIX selection string
        """
    selection_stripped = selection.strip().strip("'\"").lower()

    # 1. User directive → pass through unchanged
    if context.directives:
      prog_settings = context.directives.get("program_settings", {})
      user_sel = (prog_settings.get("phenix.polder", {}) or {}).get("selection")
      if user_sel:
        self._log(context, "BUILD: Polder selection from user directive: '%s'" % user_sel)
        return user_sel

    # 2. Invalid bare word → infer from model or default to "hetero"
    if selection_stripped in self._INVALID_POLDER_SELECTIONS:
      inferred = self._infer_polder_selection_from_model(files, context)
      self._log(context,
        "BUILD: Polder selection '%s' is not valid PHENIX syntax — "
        "replaced with '%s'" % (selection, inferred))
      return inferred

    # 3. Pass through (may be valid, e.g. "chain B", "resname ATP", "hetero")
    return selection

  def _infer_polder_selection_from_model(self, files: Dict[str, Any],
                     context: "CommandContext") -> str:
    """
        Scan the model PDB file for HETATM records to build a precise selection.

        Looks for non-water HETATM residues and returns a selection like:
          "chain B and resseq 100"         (single unique residue)
          "chain B and resseq 100:105"     (residue range in one chain)
          "hetero"                          (fallback when model unreadable or complex)

        On the server (where client files are not on disk), uses pre-extracted
        HETATM data from context.model_hetatm_residues (populated by the client
        in session_info).

        Args:
            files: Selected files dict (expects "model" key)
            context: CommandContext for logging

        Returns:
            A PHENIX atom selection string
        """
    model_path = files.get("model")
    if isinstance(model_path, list):
      model_path = model_path[0] if model_path else None

    # Collect unique (chain, resname, resseq) tuples from HETATM records,
    # excluding water (HOH, WAT, DOD, H2O).
    WATER_RESNAMES = frozenset(["HOH", "WAT", "DOD", "H2O", "SOL"])
    hetatm_residues = {}  # (chain, resseq) -> resname (ordered by first seen)

    # Primary: read from file (works on client where files are on disk)
    if model_path and context.files_local and self._file_is_available(model_path):
      try:
        with open(str(model_path), 'r', errors='replace') as fh:
          for line in fh:
            if not line.startswith("HETATM"):
              continue
            if len(line) < 26:
              continue
            resname = line[17:20].strip()
            chain   = line[21:22].strip()
            resseq  = line[22:26].strip()
            if resname in WATER_RESNAMES:
              continue
            key = (chain, resseq)
            if key not in hetatm_residues:
              hetatm_residues[key] = resname
      except Exception as e:
        self._log(context, "BUILD: Could not read model for polder selection (%s)" % e)

    # Fallback: use pre-extracted HETATM data from client (works on server)
    if not hetatm_residues and context.model_hetatm_residues:
      for entry in context.model_hetatm_residues:
        if len(entry) >= 3:
          chain, resseq, resname = entry[0], entry[1], entry[2]
          key = (chain, resseq)
          if key not in hetatm_residues:
            hetatm_residues[key] = resname
      if hetatm_residues:
        self._log(context, "BUILD: Using pre-extracted HETATM data for polder selection")

    if not hetatm_residues:
      return "hetero"

    # Group by chain
    chains = {}
    for (chain, resseq), resname in hetatm_residues.items():
      chains.setdefault(chain, []).append((int(resseq) if resseq.isdigit() else resseq, resname))

    if len(chains) == 1:
      chain = list(chains.keys())[0]
      residues = sorted(chains[chain], key=lambda x: x[0] if isinstance(x[0], int) else 0)
      resname = residues[0][1]

      if len(residues) == 1:
        # Single ligand residue
        resseq_val = residues[0][0]
        sel = "chain %s and resseq %s" % (chain, resseq_val)
        self._log(context,
          "BUILD: Polder selection inferred from HETATM: '%s' (resname %s)" %
          (sel, resname))
        return sel
      else:
        # Multiple residues in same chain — use resname if consistent, else chain
        resnames = set(r[1] for r in residues)
        if len(resnames) == 1:
          sel = "resname %s" % list(resnames)[0]
        else:
          sel = "chain %s and hetero" % chain
        self._log(context,
          "BUILD: Polder selection inferred (multiple HETATM): '%s'" % sel)
        return sel
    else:
      # Multiple chains with ligands — use hetero as safe fallback
      self._log(context,
        "BUILD: Polder selection: multiple chains with HETATM, using 'hetero'")
      return "hetero"

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
