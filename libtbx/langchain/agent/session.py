"""
Agent session management - persistent tracking of agent runs.

This module handles:
- Tracking cycle-by-cycle progress
- Persisting session state to JSON
- Checking for duplicate commands
- Generating session summaries

Usage:
    from libtbx.langchain.agent import AgentSession

    session = AgentSession(session_dir='./agent_session')
    session.start_cycle(cycle_number=1)
    session.record_decision(program='phenix.xtriage', decision='...', ...)
    session.record_result(result='Success', output_files=[...])
    session.save()
"""
from __future__ import absolute_import, division, print_function

import os
import re
import json
from datetime import datetime

# Handle both PHENIX and standalone imports
try:
    from langchain_core.prompts import PromptTemplate
except Exception:
    PromptTemplate = None  # Will be mocked if needed

try:
    from libtbx.langchain.agent.best_files_tracker import BestFilesTracker
except Exception:
    from agent.best_files_tracker import BestFilesTracker


# ---------------------------------------------------------------------------
# Module-level helpers for superseded-issue provenance matching
# ---------------------------------------------------------------------------

def _strip_stage_prefix(basename):
    """Strip leading stage-prefix from a PHENIX output filename stem.

    PHENIX programs produce output files named like::

        refine_001_model.pdb
        phased_model.pdb
        ligandfit_001_model.pdb

    This removes the leading ``<word>_<digits>_`` (or ``<word>_``) prefix
    so that ``refine_001_7qz0`` and ``7qz0`` compare as equal stems.

    Args:
        basename: Filename with extension (e.g. ``"refine_001_7qz0.pdb"``).

    Returns:
        Lowercase stem with stage prefix stripped
        (e.g. ``"7qz0"`` from ``"refine_001_7qz0.pdb"``).
    """
    import re as _re
    import os as _os
    stem = _os.path.splitext(basename.lower())[0]
    # Remove one or two leading word_digits_ prefixes
    stem = _re.sub(r'^[a-z]+_\d+_', '', stem)   # e.g. refine_001_model -> model
    stem = _re.sub(r'^[a-z]+_',      '', stem)   # e.g. phased_model     -> model
    return stem


def _extract_file_basenames(command_str):
    """Extract lowercase basenames of file-like tokens from a command string.

    Splits on whitespace, strips ``key=`` prefixes and surrounding quotes,
    then returns basenames for any token whose extension is a recognised
    crystallographic file type.

    Note: uses naive ``.split()`` rather than ``shlex.split()`` deliberately.
    For a quoted path with spaces in the directory
    (e.g. ``seq_file='/data/my dir/7qz0.fa'``), ``.split()`` yields the
    trailing fragment ``7qz0.fa'`` whose quote is stripped to give the
    correct basename ``7qz0.fa``.  ``shlex.split()`` would give the wrong
    basename ``my dir/7qz0.fa``.

    Args:
        command_str: Raw command string (may be None or empty).

    Returns:
        list: Lowercase basenames of recognised file tokens (may be empty).
    """
    import os as _os
    _FILE_EXTS = {'.pdb', '.cif', '.mtz', '.sca', '.hkl', '.mrc',
                  '.ccp4', '.map', '.fa', '.fasta', '.seq', '.dat',
                  '.ncs_spec', '.eff'}
    result = []
    for token in (command_str or '').split():
        if '=' in token:
            token = token.split('=', 1)[1]
        token = token.strip('"\'')
        ext = _os.path.splitext(token)[1].lower()
        if ext in _FILE_EXTS:
            result.append(_os.path.basename(token).lower())
    return result


class AgentSession:
    """
    Manages persistent state for an agent session across multiple cycles.

    The session is stored as a JSON file and tracks:
    - Project info (advice, original files)
    - Each cycle's decision, command, and result
    - Output files from each cycle
    - A final LLM-generated summary
    """

    def __init__(self, session_dir=None, session_file=None):
        """
        Initialize or load an agent session.

        Args:
            session_dir: Directory to store session file
            session_file: Explicit path to session file (overrides session_dir)
        """
        if session_file:
            self.session_file = session_file
        elif session_dir:
            if not os.path.exists(session_dir):
                os.makedirs(session_dir)
            self.session_file = os.path.join(session_dir, "agent_session.json")
        else:
            self.session_file = "agent_session.json"

        # Load existing session or create new one
        if os.path.exists(self.session_file):
            self.load()
        else:
            self._init_new_session()

    def _init_new_session(self):
        """Initialize a new session."""
        self.data = {
            "session_id": datetime.now().strftime("%Y-%m-%d_%H-%M-%S"),
            "project_advice": "",
            "original_files": [],
            "cycles": [],
            "summary": "",
            "resolution": None,  # Single source of truth for resolution
            "resolution_source": None,  # Which program/cycle set it (e.g., "mtriage_1")
            "ignored_commands": {},  # Commands that failed due to input errors
            # Red flag detection settings
            "experiment_type": None,  # "xray" or "cryoem", locked after first detection
            "abort_on_red_flags": True,  # Abort on critical sanity check failures
            "abort_on_warnings": False,  # Also abort on warning-level issues
            # R-free MTZ tracking (X-ray only)
            # Once R-free flags are generated, we must use the same MTZ for all refinements
            "rfree_mtz": None,  # Path to the locked R-free MTZ
            "rfree_mtz_locked_at_cycle": None,  # Cycle when it was locked
            "rfree_resolution": None,  # Resolution at which R-free flags were generated (if limited)
            # User directives extracted from project_advice
            "directives": {},  # Structured directives from user advice
            "directives_extracted": False,  # Whether extraction has been attempted
            # Error recovery tracking
            # Tracks retry attempts per error type to prevent infinite loops
            "recovery_attempts": {},  # {error_type: {count: N, files_tried: {...}}}
            # File-keyed recovery strategies (persist across cycles)
            "recovery_strategies": {},  # {file_path: {flags, program_scope, reason, selected_label, selected_label_pair}}
            # Force the planner to retry a specific program (cleared after use)
            "force_retry_program": None,
        }
        # Initialize best files tracker (not stored in self.data, serialized separately)
        self.best_files = BestFilesTracker()

    def load(self):
        """Load session from file."""
        try:
            with open(self.session_file, 'r', encoding='utf-8') as f:
                self.data = json.load(f)
            # ALWAYS rebuild best_files from cycle history on load
            # This ensures consistency after cycle removal or external modification
            # The persisted best_files is just a cache that may be stale
            self._rebuild_best_files_from_cycles()
        except Exception as e:
            print(f"Warning: Could not load session file: {e}")
            self._init_new_session()

    def _get_session_dir(self):
        """Get the directory containing the session file.

        This is the agent working directory (e.g., ai_agent_directory/) that
        contains all sub_NN_program/ output directories.  It's always known
        from the session_file path and doesn't depend on any cycle data.
        """
        if hasattr(self, 'session_file') and self.session_file:
            d = os.path.dirname(os.path.abspath(self.session_file))
            if os.path.isdir(d):
                return d
        return None

    def _discover_cycle_outputs(self, cycle):
        """Discover output files for a cycle from disk.

        Tries stored output_files paths first, then falls back to scanning
        the expected output directory based on cycle number and program name.
        This ensures files are found even when output_files is empty or
        contains stale/relative paths that don't resolve.

        This is the general strategy for surviving restart scenarios where
        the agent directory exists on disk but stored paths are broken.

        Args:
            cycle: Cycle dict from session data

        Returns:
            list: Absolute paths to existing output files
        """
        import glob
        found = []
        seen = set()

        output_files = cycle.get("output_files", [])
        session_dir = self._get_session_dir()

        # Strategy 1: Try stored paths as-is
        for f in output_files:
            if f:
                abs_path = os.path.abspath(f) if not os.path.isabs(f) else f
                if os.path.exists(abs_path):
                    bn = os.path.basename(abs_path)
                    if bn not in seen:
                        found.append(abs_path)
                        seen.add(bn)

        # Strategy 2: Resolve relative paths against session directory
        # Handles: cwd changed between runs, relative paths in session JSON
        if session_dir:
            unresolved = [f for f in output_files
                          if f and os.path.basename(f) not in seen]
            for f in unresolved:
                # Try relative to session dir (ai_agent_directory/)
                for base in [session_dir, os.path.dirname(session_dir)]:
                    resolved = os.path.normpath(os.path.join(base, f))
                    if os.path.exists(resolved):
                        bn = os.path.basename(resolved)
                        if bn not in seen:
                            found.append(os.path.abspath(resolved))
                            seen.add(bn)
                        break

        # Strategy 3: Scan expected output directory
        # Uses cycle number + program name to compute the directory name,
        # then scans for all crystallography file types.  This works even
        # when output_files is completely empty.
        if session_dir:
            cycle_num = cycle.get("cycle_number", 0)
            program = cycle.get("program", "")
            result = cycle.get("result", "")

            # Strategy 2.5: Use stored output_dir if available
            # (recorded by record_result — survives cycle-number mismatches)
            stored_dir = cycle.get("output_dir")
            if stored_dir and os.path.isdir(stored_dir):
                for ext in ("*.mtz", "*.pdb", "*.cif",
                            "*.map", "*.ccp4", "*.ncs_spec"):
                    for f in glob.glob(os.path.join(stored_dir, ext)):
                        bn = os.path.basename(f)
                        if bn not in seen:
                            found.append(f)
                            seen.add(bn)

            # Only scan for successful cycles with a known program
            if (cycle_num > 0 and program
                    and "FAILED" not in result.upper()):
                shortname = program.replace("phenix.", "").replace(".", "_")

                # Try 1: Exact cycle-number-based pattern
                dir_pattern = os.path.join(
                    session_dir, "sub_%02d_%s*" % (cycle_num, shortname))
                matched_dirs = [d for d in glob.glob(dir_pattern)
                                if os.path.isdir(d)]

                # Try 2: If exact match fails, scan ALL sub_*_{program}*
                # directories.  This handles the common case where session
                # cycle numbers diverge from GUI directory numbering after
                # an agent restart (e.g., session says cycle 2 but GUI
                # created sub_04_pdbtools).
                if not matched_dirs:
                    broad_pattern = os.path.join(
                        session_dir, "sub_*_%s*" % shortname)
                    matched_dirs = sorted(
                        [d for d in glob.glob(broad_pattern)
                         if os.path.isdir(d)])

                for output_dir in matched_dirs:
                    for ext in ("*.mtz", "*.pdb", "*.cif",
                                "*.map", "*.ccp4", "*.ncs_spec"):
                        for f in glob.glob(
                                os.path.join(output_dir, ext)):
                            bn = os.path.basename(f)
                            if bn not in seen:
                                found.append(f)
                                seen.add(bn)

        return found

    def _rebuild_best_files_from_cycles(self):
        """
        Rebuild best_files tracker from cycle history.

        This is the single source of truth - we replay all successful cycles
        through the tracker to determine the current best files.
        """
        # Create fresh tracker
        self.best_files = BestFilesTracker()

        # First, evaluate original input files (cycle 0)
        for f in self.data.get("original_files", []):
            if f and os.path.exists(f):
                self.best_files.evaluate_file(
                    path=f,
                    cycle=0,
                    metrics=None,
                    stage=None
                )

        # Replay each cycle's output files with their metrics
        for cycle in self.data.get("cycles", []):
            cycle_num = cycle.get("cycle_number", 1)
            program = cycle.get("program", "")
            result = cycle.get("result", "")
            metrics = cycle.get("metrics", {})

            # Skip failed cycles
            if "FAILED" in result.upper():
                continue

            # Discover output files — tries stored paths first, then falls
            # back to scanning the expected output directory on disk.
            output_files = self._discover_cycle_outputs(cycle)

            # Determine stage from program
            stage = self._infer_stage_from_program(program)

            for f in output_files:
                if f:
                    # Build file-specific metrics
                    file_metrics = dict(metrics) if metrics else {}

                    # Determine file-specific stage based on file AND program
                    # Only apply program-specific stage to files that match expected output patterns
                    file_stage = None  # Default: let tracker infer from filename
                    basename = os.path.basename(f).lower()

                    if f.lower().endswith('.mtz'):
                        # Check for R-free patterns
                        has_rfree_patterns = (
                            'refine' in basename or
                            'refinement_data' in basename
                        )
                        # Check for phased MTZ (not suitable for refinement)
                        is_phased_mtz = (
                            'phased' in basename or
                            'phases' in basename or
                            'solve' in basename
                        )

                        # Classify as data_mtz or map_coeffs_mtz
                        # Refine output MTZ files contain map coefficients
                        # (2mFo-DFc, Fo-Fc) AND R-free flags.  They should be
                        # classified as map_coeffs, NOT data_mtz.
                        # Patterns from file_categories.yaml refine_map_coeffs:
                        #   refine_001.mtz, 7qz0_refine_001.mtz,
                        #   refine_001_001.mtz, 7qz0_refine_001_001.mtz
                        is_map_coeffs = (
                            'map_coeffs' in basename or
                            'denmod' in basename or
                            bool(re.match(
                                r'(?:.*_)?refine_\d{3}(?:_\d{3})?\.mtz$',
                                basename))
                        )

                        if is_map_coeffs:
                            # Map coefficients - for visualization/ligand fitting
                            if 'denmod' in basename:
                                file_stage = "denmod_map_coeffs"
                            else:
                                file_stage = "refine_map_coeffs"
                        elif has_rfree_patterns and not is_phased_mtz:
                            file_stage = "original_data_mtz"
                            file_metrics["has_rfree_flags"] = True
                        elif is_phased_mtz:
                            file_stage = "phased_data_mtz"
                        else:
                            file_stage = "data_mtz"  # Generic data MTZ

                    elif f.lower().endswith(('.ccp4', '.mrc', '.map')):
                        # Map files - apply stage based on filename patterns
                        if 'initial' in basename:
                            file_stage = "intermediate_map"
                        elif 'denmod' in basename or 'density_mod' in basename:
                            file_stage = "optimized_full_map"
                        elif 'sharp' in basename:
                            file_stage = "sharpened"
                        elif 'half' in basename or any(p in basename for p in ['_1.', '_2.', '_a.', '_b.']):
                            file_stage = "half_map"
                        # else: file_stage stays None, tracker will infer

                    elif f.lower().endswith('.pdb') or f.lower().endswith('.cif'):
                        # Only apply program stage if file basename matches expected output
                        if stage == "refined" and 'refine' in basename:
                            file_stage = "refined"
                        elif stage == "phaser_output" and ('phaser' in basename or basename.startswith('mr_')):
                            file_stage = "phaser_output"
                        elif stage == "predicted" and ('predict' in basename or 'alphafold' in basename):
                            file_stage = "predicted"
                        elif stage == "processed_predicted" and ('processed' in basename or 'trimmed' in basename):
                            file_stage = "processed_predicted"
                        elif stage == "docked" and ('dock' in basename or 'placed' in basename):
                            file_stage = "docked"
                        elif stage == "autobuild_output" and ('autobuild' in basename or 'overall_best' in basename):
                            file_stage = "autobuild_output"
                        elif stage == "ligand_fit_output" and ('ligand_fit' in basename or 'lig_fit' in basename):
                            file_stage = "ligand_fit_output"
                        elif stage == "with_ligand" and 'with_ligand' in basename:
                            file_stage = "with_ligand"
                        elif stage == "with_ligand" and '_modified' in basename:
                            file_stage = "with_ligand"
                        elif stage == "rsr_output" and ('real_space' in basename or 'rsr_' in basename):
                            file_stage = "rsr_output"
                        # else: file_stage stays None, tracker will infer from filename

                    self.best_files.evaluate_file(
                        path=f,
                        cycle=cycle_num,
                        metrics=file_metrics,
                        stage=file_stage
                    )

            # Also evaluate supplemental files discovered by _find_missing_outputs.
            # These are companion outputs (e.g. map coefficients MTZ from refine)
            # that the client may not have tracked in output_files.
            # Without this, best_files["map_coeffs_mtz"] won't be populated
            # and programs with require_best_files_only (like ligandfit) can't build.
            already_seen = {os.path.basename(f) for f in output_files if f}
            supplemental = self._find_missing_outputs(cycle, already_seen)
            for sf in supplemental:
                if sf and os.path.exists(sf):
                    sf_basename = os.path.basename(sf).lower()
                    sf_stage = None
                    sf_metrics = dict(metrics) if metrics else {}

                    # Classify supplemental MTZ files
                    if sf.lower().endswith('.mtz'):
                        is_map_coeffs = (
                            'map_coeffs' in sf_basename or
                            'denmod' in sf_basename or
                            bool(re.match(
                                r'(?:.*_)?refine_\d{3}(?:_\d{3})?\.mtz$',
                                sf_basename))
                        )
                        if is_map_coeffs:
                            sf_stage = "refine_map_coeffs"

                    self.best_files.evaluate_file(
                        path=sf,
                        cycle=cycle_num,
                        metrics=sf_metrics,
                        stage=sf_stage
                    )

    def save(self):
        """Save session to file."""
        try:
            # Include best files tracker in saved data
            self.data["best_files"] = self.best_files.to_dict()
            self.data["best_files_history"] = [h.to_dict() for h in self.best_files.get_history()]
            with open(self.session_file, 'w', encoding='utf-8') as f:
                json.dump(self.data, f, indent=2)
                f.flush()
                import os
                os.fsync(f.fileno())  # Ensure data is written to disk
        except Exception as e:
            print(f"Warning: Could not save session file: {e}")

    def set_project_info(self, project_advice=None, original_files=None):
        """Set project-level information.

        On resume (session already has original_files), new files are MERGED
        into the existing list rather than replacing it.  This preserves files
        from the first run (e.g. a ligand CIF) that the user may not re-supply
        when resuming with only the refinement outputs.
        """
        if project_advice:
            self.data["project_advice"] = project_advice
        if original_files:
            if isinstance(original_files, str):
                new_files = original_files.split()
            else:
                new_files = list(original_files)

            existing = self.data.get("original_files", [])
            if existing:
                # MERGE: add new files not already present (by basename)
                existing_bns = {os.path.basename(f).lower() for f in existing}
                for f in new_files:
                    if os.path.basename(f).lower() not in existing_bns:
                        existing.append(f)
                        existing_bns.add(os.path.basename(f).lower())
                self.data["original_files"] = existing
            else:
                self.data["original_files"] = new_files

        # Store run_name if not already set (compute it once on the client)
        if not self.data.get("run_name"):
            self.data["run_name"] = self._get_run_name()

        self.save()

    # =========================================================================
    # USER DIRECTIVES
    # =========================================================================

    def extract_directives(self, provider="google", model=None, log_func=None):
        """
        Extract structured directives from project_advice using LLM.

        This should be called once at the start of a session, after
        set_project_info() has been called with the user's advice.

        Args:
            provider: LLM provider ("google", "openai", "anthropic")
            model: Specific model to use (optional)
            log_func: Optional logging function

        Returns:
            dict: Extracted directives
        """
        def log(msg):
            if log_func:
                log_func(msg)
            else:
                print(msg)

        # Don't re-extract if already done
        if self.data.get("directives_extracted"):
            log("DIRECTIVES: Already extracted, skipping")
            return self.data.get("directives", {})

        project_advice = self.data.get("project_advice", "")

        if not project_advice or project_advice.strip().lower() in ("", "none", "none provided"):
            log("DIRECTIVES: No project advice to extract from")
            self.data["directives_extracted"] = True
            self.data["directives"] = {}
            self.save()
            return {}

        try:
            try:
                from libtbx.langchain.agent.directive_extractor import (
                    extract_directives,
                    extract_directives_simple,
                    format_directives_for_display,
                    merge_directives
                )
            except ImportError:
                from agent.directive_extractor import (
                    extract_directives,
                    extract_directives_simple,
                    format_directives_for_display,
                    merge_directives
                )

            # Try LLM extraction first
            log("DIRECTIVES: Extracting from project advice...")
            llm_directives = extract_directives(
                user_advice=project_advice,
                provider=provider,
                model=model,
                log_func=log
            )

            # Also run simple extraction - it's more reliable for specific patterns
            simple_directives = extract_directives_simple(project_advice)

            # Merge results: simple extraction takes precedence for stop_conditions
            # because it's more reliable for tutorial/procedure detection
            if llm_directives and simple_directives:
                # If simple found an after_program but LLM got it wrong, use simple
                simple_stop = simple_directives.get("stop_conditions", {})
                llm_stop = llm_directives.get("stop_conditions", {})

                if simple_stop.get("after_program") and not llm_stop.get("after_program"):
                    # Simple found program target, LLM missed it - use simple
                    log("DIRECTIVES: Simple extraction found after_program=%s, merging" %
                        simple_stop["after_program"])
                    llm_directives = merge_directives(llm_directives, simple_directives)
                elif simple_stop.get("after_program") and llm_stop.get("after_cycle"):
                    # LLM got after_cycle but simple found after_program - simple is probably right
                    log("DIRECTIVES: LLM got after_cycle but simple found after_program=%s, using simple" %
                        simple_stop["after_program"])
                    # Replace LLM stop_conditions with simple's
                    llm_directives["stop_conditions"] = simple_stop

                directives = llm_directives
            elif llm_directives:
                directives = llm_directives
            elif simple_directives:
                log("DIRECTIVES: LLM extraction returned empty, using simple patterns")
                directives = simple_directives
            else:
                directives = {}

            # Store and log results
            # Sanitize: strip invalid space_group values from program_settings.
            # The LLM sometimes extracts workflow descriptions as space group
            # symbols (e.g. "determination" from "space group determination").
            directives = self._sanitize_directives(directives, log)
            directives = self._fix_autosol_atom_type_order(directives, log)
            self.data["directives"] = directives
            self.data["directives_extracted"] = True
            self.save()

            if directives:
                log(format_directives_for_display(directives))
            else:
                log("DIRECTIVES: No actionable directives found")

            # Validate directives against available capabilities
            validation_result = self._validate_directives(project_advice, directives, log)

            if not validation_result.get("valid", True):
                # Store validation issues for later reference
                self.data["directive_validation"] = validation_result
                self.save()
                log("\n" + "=" * 60)
                log("DIRECTIVE VALIDATION FAILED")
                log("=" * 60)
                log(validation_result.get("message", "Unknown validation error"))
                log("=" * 60 + "\n")

            return directives

        except ImportError as e:
            log("DIRECTIVES: Could not import directive_extractor: %s" % str(e))
            # Try simple extraction as fallback
            try:
                from libtbx.langchain.agent.directive_extractor import extract_directives_simple
                directives = extract_directives_simple(project_advice)
                directives = self._fix_autosol_atom_type_order(directives, log)
                self.data["directives"] = directives
                self.data["directives_extracted"] = True
                self.save()
                return directives
            except ImportError:
                pass

            self.data["directives_extracted"] = True
            self.data["directives"] = {}
            self.save()
            return {}

        except Exception as e:
            log("DIRECTIVES: Extraction failed: %s" % str(e))
            self.data["directives_extracted"] = True
            self.data["directives"] = {}
            self.save()
            return {}

    def get_directives(self):
        """
        Get the extracted directives.

        Returns:
            dict: Directives dict, or empty dict if not extracted
        """
        return self.data.get("directives", {})

    @staticmethod
    def _sanitize_directives(directives, log=None):
        """Remove invalid crystal symmetry values from extracted directives.

        The LLM sometimes extracts workflow descriptions as parameter values.
        For example, user advice mentioning "space group determination" can
        produce program_settings.default.space_group = "determination".
        This method strips such values before they reach the command builder.

        Args:
            directives: Extracted directives dict (modified in place)
            log: Optional logging callable

        Returns:
            The (possibly modified) directives dict
        """
        if not directives or not isinstance(directives, dict):
            return directives

        prog_settings = directives.get("program_settings")
        if not prog_settings or not isinstance(prog_settings, dict):
            return directives

        try:
            try:
                from libtbx.langchain.agent.command_postprocessor import _is_valid_space_group
            except ImportError:
                from agent.command_postprocessor import _is_valid_space_group
        except ImportError:
            return directives  # Can't validate without the function

        for scope_name in list(prog_settings.keys()):
            scope = prog_settings.get(scope_name)
            if not isinstance(scope, dict):
                continue
            sg = scope.get("space_group")
            if sg and not _is_valid_space_group(sg):
                del scope["space_group"]
                if log:
                    log("DIRECTIVES: Removed invalid space_group=%r from "
                        "program_settings.%s (not a valid space group symbol)"
                        % (str(sg), scope_name))
                # Clean up empty scope
                if not scope:
                    del prog_settings[scope_name]

        return directives

    @staticmethod
    def _fix_autosol_atom_type_order(directives, log=None):
        """Ensure atom_type is the heavier element in autosol directives.

        The LLM directive extractor non-deterministically assigns which element
        goes to atom_type vs additional_atom_types.  In SAD/MAD, the heavier
        element (higher Z) should always be the primary scatterer (atom_type).

        This mirrors _ensure_primary_scatterer_is_heavier() in
        command_postprocessor.py but operates at the directive level — before
        command assembly — so the command builder starts with correct values.

        Args:
            directives: Extracted directives dict (modified in place)
            log: Optional logging callable

        Returns:
            The (possibly modified) directives dict
        """
        if not directives or not isinstance(directives, dict):
            return directives

        prog_settings = directives.get("program_settings")
        if not prog_settings or not isinstance(prog_settings, dict):
            return directives

        # Import the Z table from command_postprocessor
        try:
            try:
                from libtbx.langchain.agent.command_postprocessor import _ANOMALOUS_Z
            except ImportError:
                from agent.command_postprocessor import _ANOMALOUS_Z
        except ImportError:
            return directives  # Can't validate without the table

        # Check autosol-specific settings
        autosol_settings = prog_settings.get("phenix.autosol")
        if not autosol_settings or not isinstance(autosol_settings, dict):
            return directives

        at = autosol_settings.get("atom_type")
        ha = autosol_settings.get("additional_atom_types")
        if not at or not ha:
            return directives

        at_z = _ANOMALOUS_Z.get(str(at), 0)
        ha_z = _ANOMALOUS_Z.get(str(ha), 0)

        if at_z > 0 and ha_z > 0 and at_z < ha_z:
            # Swap: heavier element should be atom_type
            autosol_settings["atom_type"] = ha
            autosol_settings["additional_atom_types"] = at
            if log:
                log("DIRECTIVES: Swapped atom_type=%s with "
                    "additional_atom_types=%s (Z=%d < Z=%d, heavier "
                    "element is primary scatterer)" % (at, ha, at_z, ha_z))

        return directives

    def _validate_directives(self, user_advice, directives, log_func=None):
        """
        Validate directives against available capabilities.

        Checks if the user is requesting programs or parameters that
        the agent cannot deliver.

        Args:
            user_advice: Raw user advice text
            directives: Extracted directives dict
            log_func: Optional logging function

        Returns:
            dict: Validation result with 'valid', 'message', 'issues', 'warnings'
        """
        def log(msg):
            if log_func:
                log_func(msg)

        try:
            try:
                from libtbx.langchain.agent.directive_validator import validate_directives
            except ImportError:
                from agent.directive_validator import validate_directives

            result = validate_directives(user_advice, directives)

            # Convert dataclass to dict for JSON serialization
            return {
                "valid": result.valid,
                "message": result.message,
                "issues": result.issues,
                "warnings": result.warnings,
                "unavailable_programs": result.unavailable_programs,
                "unsupported_parameters": result.unsupported_parameters,
            }

        except ImportError as e:
            log("DIRECTIVES: Could not import directive_validator: %s" % str(e))
            return {"valid": True, "message": "Validation skipped (module not available)"}
        except Exception as e:
            log("DIRECTIVES: Validation failed: %s" % str(e))
            return {"valid": True, "message": "Validation skipped (error: %s)" % str(e)}

    def get_directive_validation(self):
        """
        Get the directive validation result.

        Returns:
            dict: Validation result, or None if not validated
        """
        return self.data.get("directive_validation")

    def has_directive_validation_issues(self):
        """
        Check if there are blocking validation issues.

        Returns:
            bool: True if there are issues that should prevent workflow execution
        """
        validation = self.get_directive_validation()
        if validation is None:
            return False
        return not validation.get("valid", True)

    def get_program_directive_settings(self, program_name):
        """
        Get directive settings for a specific program.

        Merges default settings with program-specific settings.

        Args:
            program_name: Name of program (e.g., "phenix.refine")

        Returns:
            dict: Settings for the program
        """
        directives = self.get_directives()
        if not directives:
            return {}

        try:
            from libtbx.langchain.agent.directive_extractor import get_program_settings
            return get_program_settings(directives, program_name)
        except ImportError:
            # Manual implementation if import fails
            prog_settings = directives.get("program_settings", {})
            default = dict(prog_settings.get("default", {}))
            specific = prog_settings.get(program_name, {})
            return {**default, **specific}

    def check_directive_stop_conditions(self, cycle_number, last_program, metrics=None):
        """
        Check if any directive stop conditions are met.

        Args:
            cycle_number: Current cycle number
            last_program: Last program that was run
            metrics: Optional dict with r_free, map_cc, etc.

        Returns:
            tuple: (should_stop: bool, reason: str or None)
        """
        directives = self.get_directives()
        if not directives:
            return False, None

        try:
            from libtbx.langchain.agent.directive_extractor import check_stop_conditions
            should_stop, reason = check_stop_conditions(
                directives, cycle_number, last_program, metrics)
        except ImportError:
            # Manual implementation if import fails
            stop_cond = directives.get("stop_conditions", {})
            if not stop_cond:
                return False, None

            should_stop = False
            reason = None

            # Check after_cycle
            if "after_cycle" in stop_cond:
                if cycle_number >= stop_cond["after_cycle"]:
                    should_stop = True
                    reason = "Reached cycle %d (directive)" % cycle_number

            # after_program — intentionally NOT a hard stop (v112.78, Bug 7)
            # See perceive_checks.py for full rationale.  after_program is
            # now a minimum-run guarantee used by PLAN to suppress auto-stop,
            # not a trigger to force-stop.  The LLM decides when to stop.

        if not should_stop:
            return False, None

        return should_stop, reason

    def _get_last_cycle_command(self):
        """Get the command string from the most recent cycle."""
        cycles = self.data.get("cycles", [])
        if cycles:
            return cycles[-1].get("command", "")
        return ""

    def should_skip_validation(self):
        """
        Check if directives say to skip validation before stopping.

        Returns:
            bool: True if validation should be skipped
        """
        directives = self.get_directives()
        stop_cond = directives.get("stop_conditions", {})
        return stop_cond.get("skip_validation", False)

    def should_allow_stop(self):
        """
        Check if STOP should be added to valid_programs based on directives.

        This is used by workflow_engine to determine if the user's directives
        allow stopping without completing the normal workflow.

        Returns:
            bool: True if STOP should be allowed in valid_programs
        """
        directives = self.get_directives()
        if not directives:
            return False  # No directives, use normal workflow

        stop_cond = directives.get("stop_conditions", {})
        return stop_cond.get("skip_validation", False)

    def get_directive_required_programs(self):
        """
        Get programs that must be in valid_programs due to directives.

        This is used by workflow_engine to ensure tutorial targets are runnable.
        For example, if the user says "stop after resolve_cryo_em", that program
        must be available even if workflow conditions wouldn't normally include it.

        Returns:
            list: Program names that directives require to be available
        """
        directives = self.get_directives()
        if not directives:
            return []

        required = []
        stop_cond = directives.get("stop_conditions", {})

        # If after_program is set, that program must be runnable
        after_program = stop_cond.get("after_program")
        if after_program:
            required.append(after_program)

        return required

    def count_program_runs(self, program_name):
        """
        Count how many times a program has been run.

        Args:
            program_name: Name of program to count

        Returns:
            int: Number of times the program has run
        """
        count = 0
        for cycle in self.data.get("cycles", []):
            if cycle.get("program") == program_name:
                count += 1
        return count

    def check_max_program_cycles(self, program_name):
        """
        Check if max cycles for a program have been reached per directives.

        Args:
            program_name: Name of program to check

        Returns:
            tuple: (limit_reached: bool, current_count: int, max_allowed: int or None)
        """
        directives = self.get_directives()
        stop_cond = directives.get("stop_conditions", {})

        # Check for program-specific max
        if program_name == "phenix.refine" and "max_refine_cycles" in stop_cond:
            max_allowed = stop_cond["max_refine_cycles"]
            current = self.count_program_runs(program_name)
            return current >= max_allowed, current, max_allowed

        return False, self.count_program_runs(program_name), None

    def get_resolution(self):
        """
        Get the session resolution (single source of truth).

        Returns:
            float or None: Resolution in Angstroms
        """
        return self.data.get("resolution")

    def set_resolution(self, resolution, source, force=False):
        """
        Set the session resolution.

        Args:
            resolution: Resolution value in Angstroms
            source: String describing source (e.g., "mtriage_1", "resolve_cryo_em_3")
            force: If True, always update. If False, only update if better (lower) or not set.

        Returns:
            tuple: (was_updated, message)
        """
        current = self.data.get("resolution")

        if current is None:
            # First time setting resolution
            self.data["resolution"] = round(resolution, 2)
            self.data["resolution_source"] = source
            self.save()
            return True, f"Resolution set to {resolution:.2f}Å from {source}"

        if force:
            old = current
            self.data["resolution"] = round(resolution, 2)
            self.data["resolution_source"] = source
            self.save()
            return True, f"Resolution updated from {old:.2f}Å to {resolution:.2f}Å (forced by {source})"

        # Only update if better (lower resolution value = higher resolution data)
        if resolution < current:
            old = current
            self.data["resolution"] = round(resolution, 2)
            self.data["resolution_source"] = source
            self.save()
            return True, f"Resolution improved from {old:.2f}Å to {resolution:.2f}Å by {source}"
        elif resolution > current + 0.5:
            # Significantly worse - warn but don't update
            return False, f"WARNING: {source} reports {resolution:.2f}Å but session resolution is {current:.2f}Å (not updating)"
        else:
            # Similar value, no update needed
            return False, None

    # =========================================================================
    # EXPERIMENT TYPE (locked after first detection)
    # =========================================================================

    def get_experiment_type(self):
        """
        Get the locked experiment type.

        Returns:
            str or None: "xray" or "cryoem", or None if not yet set
        """
        return self.data.get("experiment_type")

    def set_experiment_type(self, exp_type):
        """
        Set the experiment type (only if not already set).

        The experiment type is locked after first detection to prevent
        spurious files from changing the workflow mid-session.

        Args:
            exp_type: "xray" or "cryoem"

        Returns:
            tuple: (was_set, message)
        """
        current = self.data.get("experiment_type")

        if current is None:
            self.data["experiment_type"] = exp_type
            self.save()
            return True, f"Experiment type set to {exp_type}"

        if current == exp_type:
            return False, None  # Same type, no change needed

        # Type mismatch - this is a red flag situation
        return False, f"Experiment type already locked as {current} (attempted to set {exp_type})"

    # =========================================================================
    # R-FREE MTZ TRACKING (X-ray only)
    # =========================================================================
    # In X-ray crystallography, R-free flags are used for cross-validation.
    # Once generated, the SAME flags must be used for all subsequent refinements.
    # Using different R-free flags invalidates the R-free statistic.
    #
    # This tracking ensures:
    # 1. First refinement that generates R-free flags locks the MTZ
    # 2. All subsequent refinements MUST use this locked MTZ
    # 3. The locked MTZ is visible in session.json for debugging

    def get_rfree_mtz(self):
        """
        Get the locked R-free MTZ path.

        Returns:
            str or None: Path to the MTZ with R-free flags, or None if not yet set
        """
        return self.data.get("rfree_mtz")

    def get_rfree_resolution(self):
        """
        Get the resolution at which R-free flags were generated (if limited).

        When refinement is run at a limited resolution (e.g., 3Å), the R-free
        flags are only valid for that resolution range. Programs like polder
        must either match this resolution or generate new R-free flags.

        Returns:
            float or None: Resolution limit in Angstroms, or None if not limited
        """
        return self.data.get("rfree_resolution")

    def set_rfree_mtz(self, path, cycle):
        """
        Lock the R-free MTZ. Can only be set once per session.

        This should be called when phenix.refine first generates or uses
        R-free flags. Once locked, all subsequent refinements must use
        this MTZ to maintain valid R-free statistics.

        Args:
            path: Full path to the MTZ file with R-free flags
            cycle: Cycle number when the lock occurred

        Returns:
            tuple: (was_locked, message)
        """
        current = self.data.get("rfree_mtz")

        if current is not None:
            # Already locked - this is expected, not an error
            if current == path:
                return False, None  # Same path, no change
            return False, f"R-free MTZ already locked to {current} (cycle {self.data.get('rfree_mtz_locked_at_cycle')})"

        # Lock it
        self.data["rfree_mtz"] = path
        self.data["rfree_mtz_locked_at_cycle"] = cycle
        self.save()
        return True, f"R-free MTZ locked to {os.path.basename(path)} at cycle {cycle}"

    def is_rfree_mtz_locked(self):
        """Check if an R-free MTZ has been locked."""
        return self.data.get("rfree_mtz") is not None

    # =========================================================================
    # ERROR RECOVERY TRACKING
    # =========================================================================
    # The agent can automatically recover from certain well-defined errors.
    # These methods manage recovery state:
    # - recovery_strategies: File-keyed flags to add to commands
    # - force_retry_program: Forces the planner to retry a specific program
    # - recovery_attempts: Tracks retries to prevent infinite loops

    def set_recovery_strategy(self, file_path, flags, program, reason,
                              selected_label="", selected_label_pair=""):
        """
        Set a recovery strategy for a specific file.

        When the CommandBuilder sees this file being used, it will add
        the specified flags to the command.

        Args:
            file_path: Path to the file that triggered the error
            flags: Dict of parameter flags to add (e.g., {"obs_labels": "I(+)"})
            program: Program that triggered the error
            reason: Human-readable explanation
            selected_label: The main label (e.g., "FTOXD3")
            selected_label_pair: Full label pair (e.g., "FTOXD3,SIGFTOXD3")

        Example:
            session.set_recovery_strategy(
                "/path/to/data.mtz",
                {"scaling.input.xray_data.obs_labels": "I_CuKa(+)"},
                "phenix.autosol",
                "Selected anomalous data for MRSAD workflow"
            )
        """
        strategies = self.data.setdefault("recovery_strategies", {})
        strategies[file_path] = {
            "flags": flags,
            "program_scope": [],  # Empty = applies to all programs
            "reason": reason,
            "selected_label": selected_label,
            "selected_label_pair": selected_label_pair
        }
        self.save()

    def get_recovery_strategy(self, file_path):
        """
        Get recovery strategy for a specific file.

        Args:
            file_path: Path to check for recovery strategy

        Returns:
            dict or None: {flags: {...}, program_scope: [...], reason: "..."}
        """
        return self.data.get("recovery_strategies", {}).get(file_path)

    def get_all_recovery_strategies(self):
        """
        Get all recovery strategies.

        Returns:
            dict: {file_path: {flags, program_scope, reason}, ...}
        """
        return self.data.get("recovery_strategies", {})

    def clear_recovery_strategy(self, file_path):
        """
        Clear recovery strategy for a specific file.

        Typically called after successful recovery or when starting fresh.
        """
        strategies = self.data.get("recovery_strategies", {})
        if file_path in strategies:
            del strategies[file_path]
            self.save()

    def set_force_retry_program(self, program):
        """
        Set the program to force-retry on next cycle.

        This bypasses normal planning - the planner will immediately
        select this program without consulting the LLM.

        Args:
            program: Program name (e.g., "phenix.autosol")
        """
        self.data["force_retry_program"] = program
        self.save()

    def get_force_retry_program(self):
        """
        Get the force-retry program (if any).

        Returns:
            str or None: Program name to force, or None
        """
        return self.data.get("force_retry_program")

    def clear_force_retry_program(self):
        """
        Clear and return the force-retry program.

        This should be called by the planner after handling the forced retry.

        Returns:
            str or None: The program that was set, or None
        """
        return self.data.pop("force_retry_program", None)

    def get_recovery_attempts(self, error_type=None):
        """
        Get recovery attempt tracking info.

        Args:
            error_type: Specific error type to query, or None for all

        Returns:
            dict: Attempt info for the error type(s)
        """
        attempts = self.data.get("recovery_attempts", {})
        if error_type:
            return attempts.get(error_type, {"count": 0, "files_tried": {}})
        return attempts

    def clear_recovery_state(self):
        """
        Clear all recovery state.

        Useful when starting a completely new workflow or after
        user intervention has resolved the underlying issue.
        """
        self.data["recovery_attempts"] = {}
        self.data["recovery_strategies"] = {}
        self.data["force_retry_program"] = None
        self.save()

    # =========================================================================
    # SANITY CHECK SETTINGS
    # =========================================================================

    def get_abort_settings(self):
        """
        Get the abort settings for sanity checking.

        Returns:
            dict: {"abort_on_red_flags": bool, "abort_on_warnings": bool}
        """
        return {
            "abort_on_red_flags": self.data.get("abort_on_red_flags", True),
            "abort_on_warnings": self.data.get("abort_on_warnings", False),
        }

    def set_abort_settings(self, abort_on_red_flags=None, abort_on_warnings=None):
        """
        Set the abort settings for sanity checking.

        Args:
            abort_on_red_flags: If True, abort on critical issues
            abort_on_warnings: If True, also abort on warning-level issues
        """
        if abort_on_red_flags is not None:
            self.data["abort_on_red_flags"] = abort_on_red_flags
        if abort_on_warnings is not None:
            self.data["abort_on_warnings"] = abort_on_warnings
        self.save()

    # =========================================================================
    # IGNORED COMMANDS (Input errors that should not be repeated)
    # =========================================================================

    def record_ignored_command(self, command, reason):
        """
        Record a command that failed due to input error (wrong parameters, etc).

        These commands should not be repeated - they indicate agent mistakes
        that need to be avoided on retry.

        Args:
            command: The full command string that failed
            reason: Brief description of why it failed
        """
        if "ignored_commands" not in self.data:
            self.data["ignored_commands"] = {}
        self.data["ignored_commands"][command] = {
            "reason": reason,
            "timestamp": datetime.now().isoformat()
        }
        self.save()

    def get_ignored_commands(self):
        """
        Get dict of commands that failed due to input errors.

        Returns:
            dict: {command: {"reason": str, "timestamp": str}, ...}
        """
        return self.data.get("ignored_commands", {})

    def is_command_ignored(self, command):
        """
        Check if a command was previously ignored.

        Args:
            command: Command string to check

        Returns:
            tuple: (is_ignored, reason) - reason is None if not ignored
        """
        ignored = self.data.get("ignored_commands", {})
        if command in ignored:
            return True, ignored[command].get("reason", "Unknown")
        return False, None

    def clear_ignored_commands(self):
        """Clear all ignored commands (use when starting fresh)."""
        self.data["ignored_commands"] = {}
        self.save()

    def record_bad_inject_param(self, program_name, param_key):
        """Record a parameter name that caused an 'Unknown command line parameter'
        failure for a specific program.

        Once recorded, inject_user_params (in command_postprocessor.py) will
        skip this key for the same program on every future cycle, breaking
        the inject→crash→re-inject loop.

        Args:
            program_name: e.g. "phenix.refine"
            param_key:    The bare or dotted parameter name that was rejected,
                          e.g. "ignore_symmetry_conflicts" or "bad.scope.key"
        """
        if "bad_inject_params" not in self.data:
            self.data["bad_inject_params"] = {}
        prog_bad = self.data["bad_inject_params"].setdefault(program_name, [])
        if param_key not in prog_bad:
            prog_bad.append(param_key)
        self.save()

    def get_bad_inject_params(self, program_name):
        """Return the set of parameter keys blacklisted for *program_name*.

        Args:
            program_name: e.g. "phenix.refine"

        Returns:
            set of str — parameter keys to skip during injection
        """
        bad = self.data.get("bad_inject_params", {})
        return set(bad.get(program_name, []))

    # =========================================================================
    # BEST FILES TRACKING
    # =========================================================================

    def get_best_model(self):
        """
        Get the path to the best model file.

        Returns:
            str or None: Path to best model
        """
        return self.best_files.get_best_path("model")

    def get_best_map(self):
        """
        Get the path to the best map file.

        Returns:
            str or None: Path to best map
        """
        return self.best_files.get_best_path("map")

    def get_best_data_mtz(self):
        """
        Get the path to the best data MTZ file.

        Returns:
            str or None: Path to best data MTZ (with Fobs/R-free flags)
        """
        return self.best_files.get_best_path("data_mtz")

    def get_best_map_coeffs_mtz(self):
        """
        Get the path to the best map coefficients MTZ file.

        Returns:
            str or None: Path to best map coefficients MTZ (for ligand fitting)
        """
        return self.best_files.get_best_path("map_coeffs_mtz")

    def get_best_mtz(self):
        """
        Get the path to the best MTZ file (alias for get_best_data_mtz).

        Returns:
            str or None: Path to best data MTZ (with R-free flags if available)
        """
        return self.best_files.get_best_path("data_mtz")

    def get_best_files_dict(self):
        """
        Get all best files as a dict of category -> path.

        Returns:
            dict: {category: path, ...}
        """
        return self.best_files.get_best_dict()

    def get_best_files_summary(self):
        """
        Get a human-readable summary of best files.

        Returns:
            str: Multi-line summary
        """
        return self.best_files.get_summary()

    def evaluate_file_for_best(self, path, cycle, metrics=None, stage=None):
        """
        Manually evaluate a file for the best files tracker.

        Useful for evaluating user-provided input files at session start.

        Args:
            path: File path
            cycle: Cycle number (use 0 for input files)
            metrics: Optional metrics dict
            stage: Optional stage override

        Returns:
            bool: True if file became new best
        """
        return self.best_files.evaluate_file(path, cycle, metrics, stage)

    def update_best_model_metrics(self, metrics, cycle=None):
        """
        Update metrics for the best model and recalculate its score.

        This should be called when validation provides new metrics for
        the current best model (e.g., after molprobity validation).

        Args:
            metrics: Dict with model metrics (r_free, map_cc, clashscore, etc.)
            cycle: Optional cycle number for logging

        Returns:
            tuple: (was_updated, old_score, new_score)
        """
        result = self.best_files.update_best_model_metrics(metrics, cycle)
        if result[0]:  # was_updated
            self.save()
        return result

    def write_history_file(self, output_path=None):
        """
        Write session history as a JSON file for the agent to read.

        Writes the full cycle data so the graph can properly analyze
        what programs have been run.

        Args:
            output_path: Path to write to. If None, creates temp file.

        Returns:
            str: Path to the written file
        """
        import tempfile

        if output_path is None:
            fd, output_path = tempfile.mkstemp(suffix='.json', prefix='session_history_')
            os.close(fd)

        # Write the full session data including all cycles
        # The graph's history loader will extract the cycles
        history_data = {
            "session_id": self.data.get("session_id", ""),
            "experiment_type": self.data.get("experiment_type"),
            "cycles": self.data.get("cycles", []),
        }

        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(history_data, f, indent=2)

        return output_path

    def start_cycle(self, cycle_number):
        """
        Start a new cycle. Creates the cycle entry if it doesn't exist.

        Args:
            cycle_number: The cycle number (1-indexed)
        """
        # Ensure we have enough cycle entries
        while len(self.data["cycles"]) < cycle_number:
            self.data["cycles"].append({
                "cycle_number": len(self.data["cycles"]) + 1,
                "program": "",
                "decision": "",
                "reasoning": "",
                "explanation": "",
                "command": "",
                "result": "",
                "output_files": [],  # NEW: Track output files
                "timestamp": ""
            })
        self.save()

    def record_decision(self, cycle_number, program=None, decision=None,
                        reasoning=None, explanation=None, command=None,
                        expert_assessment=None):
        """
        Record the decision made for a cycle.

        Args:
            cycle_number: Which cycle to update
            program: The selected program
            decision: Short decision text
            reasoning: Detailed reasoning
            explanation: Full explanation including plan
            command: The generated command
            expert_assessment: Dict from thinking agent
        """
        self.start_cycle(cycle_number)
        cycle = self.data["cycles"][cycle_number - 1]

        if program is not None:
            cycle["program"] = program
        if decision is not None:
            cycle["decision"] = decision
        if reasoning is not None:
            cycle["reasoning"] = reasoning
        if explanation is not None:
            cycle["explanation"] = explanation
        if command is not None:
            cycle["command"] = command
        if expert_assessment and isinstance(
          expert_assessment, dict
        ):
            cycle["expert_assessment"] = (
              expert_assessment)

        cycle["timestamp"] = datetime.now().isoformat()
        self.save()

    def record_result(self, cycle_number, result, output_files=None):
        """
        Record the result of running a cycle's command.

        Args:
            cycle_number: Which cycle to update
            result: Result text (success/failure description)
            output_files: List of output file paths produced by the command
        """
        self.start_cycle(cycle_number)
        cycle = self.data["cycles"][cycle_number - 1]
        cycle["result"] = result

        # Store output files
        if output_files:
            # Ensure output_files is a list
            if isinstance(output_files, str):
                output_files = [output_files]
            cycle["output_files"] = list(output_files)

            # Store the output directory for reliable future lookups.
            # Session cycle_number may not match GUI directory numbering
            # after restarts, so storing the actual directory is essential.
            for f in output_files:
                if f and os.path.isabs(f) and os.path.exists(f):
                    cycle["output_dir"] = os.path.dirname(f)
                    break

        # Extract and store metrics from result for use by workflow_state
        program = cycle.get("program", "")
        metrics = self._extract_metrics_from_result(result, program)
        if metrics:
            cycle["metrics"] = metrics

            # Update session resolution from key programs
            if metrics.get("resolution"):
                resolution = metrics["resolution"]
                program_lower = program.lower()
                source = f"{program}_{cycle_number}"

                # Programs that establish resolution (first time)
                if "xtriage" in program_lower or "mtriage" in program_lower:
                    updated, msg = self.set_resolution(resolution, source)
                    if msg:
                        cycle["resolution_note"] = msg

                # Programs that might improve resolution (cryo-EM map optimization)
                elif "resolve_cryo_em" in program_lower or "map_sharpening" in program_lower:
                    # These can improve resolution - update if better
                    updated, msg = self.set_resolution(resolution, source)
                    if msg:
                        cycle["resolution_note"] = msg

                # Other programs - just note if significantly different
                else:
                    current = self.get_resolution()
                    if current and abs(resolution - current) > 0.5:
                        cycle["resolution_note"] = f"Program reports {resolution:.2f}Å (session: {current:.2f}Å)"

        # Update best files tracker with output files from successful cycles
        if output_files and "FAILED" not in result.upper():
            model_stage = self._infer_stage_from_program(program)

            # Check if this is an X-ray refinement program (produces MTZ with R-free flags)
            is_xray_refinement = (
                program and
                "refine" in program.lower() and
                "real_space" not in program.lower()  # real_space_refine is cryo-EM
            )

            for f in output_files:
                # Build file-specific metrics
                file_metrics = dict(metrics) if metrics else {}

                # Determine stage based on file type
                # MTZ files from refinement get "refined_mtz" stage
                # PDB/CIF files get the model stage (e.g., "refined")
                if f.lower().endswith('.mtz'):
                    basename = os.path.basename(f).lower()

                    # Check if this MTZ has R-free flags (from refinement or autosol)
                    # Patterns that indicate R-free flags are present:
                    # - *refine* (from phenix.refine)
                    # - *refinement_data* (from autosol)
                    has_rfree_patterns = (
                        'refine' in basename or
                        'refinement_data' in basename
                    )

                    # Patterns that indicate this is NOT suitable for refinement
                    # (phased MTZ files have multiple arrays)
                    is_phased_mtz = (
                        'phased' in basename or
                        'phases' in basename or
                        'solve' in basename
                    )

                    # Check if this is map coefficients (for visualization/ligand fitting)
                    # Matches: refine_001.mtz, refine_001_001.mtz,
                    #          7qz0_refine_001.mtz, 7qz0_refine_001_001.mtz
                    is_map_coeffs = (
                        'map_coeffs' in basename or
                        'denmod' in basename or
                        re.match(
                            r'(?:.*_)?refine_\d{3}(?:_\d{3})?\.mtz$',
                            basename)
                    )

                    if is_map_coeffs:
                        # Map coefficients - for visualization/ligand fitting
                        if 'denmod' in basename:
                            file_stage = "denmod_map_coeffs"
                        else:
                            file_stage = "refine_map_coeffs"
                    elif has_rfree_patterns and not is_phased_mtz:
                        file_stage = "original_data_mtz"
                        file_metrics["has_rfree_flags"] = True

                        # Check if this refinement was resolution-limited
                        # If so, the R-free flags only cover that resolution range
                        # and we should NOT lock this MTZ for future use
                        command = cycle.get("command", "")
                        is_resolution_limited = (
                            "high_resolution=" in command or
                            "d_min=" in command or
                            "xray_data.high_resolution=" in command
                        )

                        if is_resolution_limited:
                            # Extract the resolution value from command
                            res_match = re.search(
                                r'(?:high_resolution|d_min|xray_data\.high_resolution)\s*=\s*([\d.]+)',
                                command
                            )
                            if res_match:
                                rfree_resolution = float(res_match.group(1))
                                file_metrics["rfree_resolution"] = rfree_resolution
                                # Store in session for use by command builder
                                if self.data.get("rfree_resolution") is None:
                                    self.data["rfree_resolution"] = rfree_resolution
                                cycle["rfree_mtz_note"] = (
                                    f"R-free flags limited to {rfree_resolution}Å - "
                                    "programs using this MTZ must match this resolution"
                                )
                            else:
                                cycle["rfree_mtz_note"] = (
                                    "R-free flags limited to refinement resolution - not locking"
                                )
                            file_metrics["rfree_resolution_limited"] = True
                        else:
                            # Lock the R-free MTZ (first one wins)
                            was_locked, lock_msg = self.set_rfree_mtz(f, cycle_number)
                            if was_locked:
                                cycle["rfree_mtz_locked"] = f
                                if lock_msg:
                                    cycle["rfree_mtz_note"] = lock_msg
                    elif is_phased_mtz:
                        # Phased data MTZ - has experimental phases
                        file_stage = "phased_data_mtz"
                    else:
                        # Generic data MTZ
                        file_stage = "data_mtz"
                elif f.lower().endswith(('.ccp4', '.mrc', '.map')):
                    # Map files - apply stage only to matching patterns
                    basename = os.path.basename(f).lower()

                    # Check for intermediate maps (should be deprioritized)
                    if 'initial' in basename:
                        file_stage = "intermediate_map"
                    # Check for optimized/density-modified maps
                    elif 'denmod' in basename or 'density_mod' in basename:
                        file_stage = "optimized_full_map"
                    elif 'sharp' in basename:
                        file_stage = "sharpened"
                    # Check for half-maps
                    elif 'half' in basename or any(p in basename for p in ['_1.', '_2.', '_a.', '_b.']):
                        file_stage = "half_map"
                    else:
                        # Generic map - let tracker infer
                        file_stage = None
                elif f.lower().endswith('.pdb') or f.lower().endswith('.cif'):
                    # PDB/CIF files - only apply program stage if filename matches expected output
                    # This prevents unrelated files in output_files from getting wrong stage
                    basename = os.path.basename(f).lower()

                    if model_stage == "refined" and 'refine' in basename:
                        file_stage = "refined"
                    elif model_stage == "phaser_output" and ('phaser' in basename or basename.startswith('mr_')):
                        file_stage = "phaser_output"
                    elif model_stage == "predicted" and ('predict' in basename or 'alphafold' in basename):
                        file_stage = "predicted"
                    elif model_stage == "processed_predicted" and ('processed' in basename or 'trimmed' in basename):
                        file_stage = "processed_predicted"
                    elif model_stage == "docked" and ('dock' in basename or 'placed' in basename):
                        file_stage = "docked"
                    elif model_stage == "autobuild_output" and ('autobuild' in basename or 'overall_best' in basename):
                        file_stage = "autobuild_output"
                    elif model_stage == "ligand_fit_output" and ('ligand_fit' in basename or 'lig_fit' in basename):
                        file_stage = "ligand_fit_output"
                    elif model_stage == "with_ligand" and 'with_ligand' in basename:
                        file_stage = "with_ligand"
                    elif model_stage == "with_ligand" and '_modified' in basename:
                        # pdbtools default output: {input}_modified.pdb
                        # This IS a with_ligand model even though the
                        # filename doesn't contain "with_ligand"
                        file_stage = "with_ligand"
                    elif model_stage == "rsr_output" and ('real_space' in basename or 'rsr_' in basename):
                        file_stage = "rsr_output"
                    else:
                        # Filename doesn't match program output pattern - let tracker infer
                        # This prevents PredictAndBuild_0_predicted_model_processed.pdb from
                        # getting stage="refined" just because it was in refine's output_files
                        file_stage = None
                else:
                    # Unknown file type - let tracker handle
                    file_stage = None

                # Evaluate each output file - tracker will classify and score
                updated = self.best_files.evaluate_file(
                    path=f,
                    cycle=cycle_number,
                    metrics=file_metrics,
                    stage=file_stage
                )
                if updated:
                    category = self.best_files._classify_category(f)
                    if category:
                        cycle.setdefault("best_file_updates", []).append(category)

            # Discover and evaluate supplemental output files that the client
            # may not have tracked.  Same logic as _rebuild_best_files_from_cycles
            # but on the live path.  E.g., refine produces refine_001.mtz (map
            # coefficients) alongside refine_001_data.mtz, but the client may
            # only report the _data.mtz.  Without this, best_files["map_coeffs_mtz"]
            # stays empty and ligandfit can't build.
            already_seen = {os.path.basename(f) for f in output_files if f}
            supplemental = self._find_missing_outputs(cycle, already_seen)
            for sf in supplemental:
                if sf and os.path.exists(sf):
                    sf_basename = os.path.basename(sf).lower()
                    sf_stage = None
                    sf_metrics = dict(metrics) if metrics else {}

                    if sf.lower().endswith('.mtz'):
                        is_map_coeffs = (
                            'map_coeffs' in sf_basename or
                            'denmod' in sf_basename or
                            bool(re.match(
                                r'(?:.*_)?refine_\d{3}(?:_\d{3})?\.mtz$',
                                sf_basename))
                        )
                        if is_map_coeffs:
                            sf_stage = "refine_map_coeffs"

                    updated = self.best_files.evaluate_file(
                        path=sf,
                        cycle=cycle_number,
                        metrics=sf_metrics,
                        stage=sf_stage
                    )
                    if updated:
                        # Also add to cycle's output_files so get_available_files
                        # picks them up without needing _find_missing_outputs again
                        cycle.setdefault("output_files", []).append(sf)

        # If this was a validation program, update best model metrics
        # (validation metrics apply to the current best model, not a new file)
        if metrics and program:
            program_lower = program.lower()
            validation_programs = ["molprobity", "validation", "clashscore"]
            is_validation = any(vp in program_lower for vp in validation_programs)

            if is_validation:
                # Extract model-relevant metrics
                model_metrics = {}
                for key in ["clashscore", "rama_favored", "rama_outliers",
                           "rotamer_outliers", "cbeta_outliers"]:
                    if key in metrics:
                        model_metrics[key] = metrics[key]

                if model_metrics:
                    was_updated, old_score, new_score = self.best_files.update_best_model_metrics(
                        model_metrics, cycle_number
                    )
                    if was_updated:
                        cycle["best_model_score_update"] = {
                            "old": old_score,
                            "new": new_score,
                            "reason": "validation metrics"
                        }

        self.save()

    def _infer_stage_from_program(self, program):
        """
        Infer output file stage from program name.

        Args:
            program: Program name (e.g., "phenix.refine", "phenix.real_space_refine")

        Returns:
            str: Stage name for best files tracker
        """
        if not program:
            return None

        program_lower = program.lower()

        # Model stages - order matters for patterns that overlap
        if "refine" in program_lower and "real_space" not in program_lower:
            return "refined"
        if "real_space_refine" in program_lower:
            return "rsr_output"
        if "autobuild" in program_lower:
            return "autobuild_output"
        if "map_to_model" in program_lower:
            return "autobuild_output"
        if "dock_in_map" in program_lower:
            return "docked"
        if "process_predicted_model" in program_lower:
            return "processed_predicted"
        if "predict_and_build" in program_lower:
            return "predicted"
        if "phaser" in program_lower:
            return "phaser_output"
        if "ligandfit" in program_lower:
            return "ligand_fit_output"
        if "pdbtools" in program_lower:
            return "with_ligand"

        # Map stages
        if "resolve_cryo_em" in program_lower:
            return "optimized_full_map"
        if "map_sharpening" in program_lower or "local_aniso" in program_lower:
            return "sharpened"
        if "density_modification" in program_lower or "dm" in program_lower:
            return "density_modified"

        # MTZ stages
        if "xtriage" in program_lower:
            return "original"  # xtriage doesn't modify data

        return None

    def get_cycle(self, cycle_number):
        """Get data for a specific cycle."""
        if cycle_number <= len(self.data["cycles"]):
            return self.data["cycles"][cycle_number - 1]
        return None

    def get_all_commands(self):
        """
        Get all commands that have been run successfully (for duplicate detection).
        Only includes cycles where the command actually executed successfully.

        Returns:
            list of tuples: [(cycle_number, normalized_command), ...]
        """
        commands = []
        for cycle in self.data["cycles"]:
            # Skip CRASH cycles
            program = cycle.get("program", "")
            if program == "CRASH" or program == "":
                continue

            # Only include successful runs
            result = cycle.get("result", "")
            if not result.startswith("SUCCESS"):
                continue

            cmd = cycle.get("command", "")
            if cmd:
                # Normalize command (remove extra whitespace)
                norm_cmd = " ".join(cmd.strip().split())
                commands.append((cycle.get("cycle_number", 0), norm_cmd))

        return commands

    def cleanup_failed_cycles(self):
        """
        Remove failed/incomplete cycles from session.
        Called at session start to clean up previous crashed runs.

        Returns:
            int: Number of cycles removed
        """
        original_count = len(self.data["cycles"])

        # Keep only cycles that have both a command and a successful result
        self.data["cycles"] = [
            c for c in self.data["cycles"]
            if c.get("program") not in ("CRASH", "")
            and c.get("result", "").startswith("SUCCESS")
        ]

        # Renumber remaining cycles
        for i, cycle in enumerate(self.data["cycles"]):
            cycle["cycle_number"] = i + 1

        removed = original_count - len(self.data["cycles"])
        if removed > 0:
            print(f"Cleaned up {removed} failed/incomplete cycles from previous session")
            self.save()

        return removed

    def get_all_failed_commands(self):
        """
        Get all commands that were attempted and FAILED (for exact-duplicate loop detection).

        Only exact matches are meaningful here — we only want to catch the case where
        the same command is retried verbatim after failing.  The 80%-overlap heuristic
        is NOT applied to failed commands because partial similarity in failed attempts
        is too noisy (e.g., two different recovery strategies that share the same data
        file would incorrectly be flagged as duplicates).

        Returns:
            list of tuples: [(cycle_number, normalized_command), ...]
        """
        commands = []
        for cycle in self.data["cycles"]:
            program = cycle.get("program", "")
            if program in ("CRASH", ""):
                continue
            result = cycle.get("result", "")
            # Only include cycles that definitely failed (not incomplete / in-progress)
            if not result:
                continue
            if result.startswith("SUCCESS"):
                continue
            cmd = cycle.get("command", "")
            if cmd:
                norm_cmd = " ".join(cmd.strip().split())
                commands.append((cycle.get("cycle_number", 0), norm_cmd))
        return commands

    def get_cycle_result(self, cycle_number):
        """Return the result text of a specific cycle, or None if not found.

        Used by _build_duplicate_feedback to include the original error text
        in the feedback message when a command is rejected as a failed-duplicate.

        Args:
            cycle_number: 1-based cycle number.

        Returns:
            str or None: The result/error text stored on the cycle, or None.
        """
        for cycle in self.data.get("cycles", []):
            if cycle.get("cycle_number") == cycle_number:
                return cycle.get("result") or None
        return None


    def is_duplicate_command(self, command):
        """
        Check if a command (or very similar) has already been run.

        Two-tier check:
          1. Exact match against ALL previous commands — both successful AND failed.
             An exact repeated failure is always a loop; stop it immediately.
          2. 80%-token-overlap heuristic against SUCCESSFUL commands only.
             This catches semantically equivalent successful runs (same program,
             same core params, minor path differences) without generating false
             positives from failed recovery attempts.

        Args:
            command: The command to check

        Returns:
            tuple: (is_duplicate: bool, previous_cycle: int or None,
                    was_failure: bool)
                   was_failure is True when the matched prior cycle failed —
                   callers can use this to generate context-appropriate feedback.
        """
        if not command or command == "No command generated.":
            return False, None

        norm_new = " ".join(command.strip().split())

        # ── Tier 1: exact match against ALL commands (success AND failure) ───
        # This catches the common loop pattern where the LLM retries an identical
        # command after a failure (e.g., missing required file, unrecognised param).
        all_failed = self.get_all_failed_commands()
        failed_nums = {cn for cn, _ in all_failed}
        for cycle_num, norm_cmd in self.get_all_commands() + all_failed:
            if not norm_cmd:
                continue
            if norm_cmd == norm_new:
                return True, cycle_num, (cycle_num in failed_nums)

        # ── Tier 2: 80%-overlap heuristic against SUCCESSFUL commands only ───
        # Extract program name from new command
        new_program = norm_new.split()[0] if norm_new else ""

        for cycle_num, norm_cmd in self.get_all_commands():
            # Skip if this cycle had no real command
            if not norm_cmd:
                continue

            # Extract program name from old command
            old_program = norm_cmd.split()[0] if norm_cmd else ""

            # Different programs are never duplicates
            if os.path.basename(new_program) != os.path.basename(old_program):
                continue

            # (Exact match already handled in Tier 1 above)

            # For predict_and_build, stop_after_predict is a critical parameter
            # Commands with different stop_after_predict values are NOT duplicates
            if "predict_and_build" in new_program:
                new_has_stop = "stop_after_predict" in norm_new
                old_has_stop = "stop_after_predict" in norm_cmd
                if new_has_stop != old_has_stop:
                    continue  # Different mode, not a duplicate

            # Check if same program with same core parameters
            # (ignore path differences)
            new_parts = set(os.path.basename(p) for p in norm_new.split())
            old_parts = set(os.path.basename(p) for p in norm_cmd.split())

            # Extract file tokens (basenames with crystallographic extensions)
            _FILE_EXTS = {'.pdb', '.cif', '.mtz', '.sca', '.hkl', '.mrc',
                          '.ccp4', '.map', '.fa', '.seq', '.dat'}
            def _file_tokens(parts):
                return {p for p in parts
                        if os.path.splitext(p)[1].lower() in _FILE_EXTS}

            new_files = _file_tokens(new_parts)
            old_files = _file_tokens(old_parts)

            # If the input file basenames differ, this is NOT a duplicate.
            # Running refine with a different model is a different computation,
            # even if all other params match.
            if new_files != old_files:
                continue

            # If >80% overlap in tokens, consider it a duplicate
            if len(new_parts) > 0 and len(old_parts) > 0:
                overlap = len(new_parts & old_parts) / max(len(new_parts), len(old_parts))
                if overlap > 0.8:
                    return True, cycle_num, False  # matched successful cycle

        return False, None, False

    def get_consecutive_program_count(self, program_name):
        """Count how many times a program ran consecutively at the end of history.

        Looks at the tail of completed (SUCCESS) cycles and counts how many
        consecutive entries match *program_name* (compared by base name,
        e.g. "phenix.refine" matches "phenix.refine").

        Returns:
            int: Number of consecutive runs (0 if last cycle was different).
        """
        count = 0
        target = os.path.basename(program_name) if program_name else ""
        # Walk cycles in reverse
        for cycle in reversed(self.data.get("cycles", [])):
            prog = cycle.get("program", "")
            result = cycle.get("result", "")
            if not prog or prog in ("CRASH", ""):
                continue
            if not result.startswith("SUCCESS"):
                continue
            if os.path.basename(prog) == target:
                count += 1
            else:
                break  # different program — stop counting
        return count

    def get_history_for_agent(self):
        """
        Get cycle history in format expected by the agent.

        IMPORTANT: This now includes output_files so the agent knows
        what files were produced by each cycle, and metrics (as 'analysis')
        so workflow_state can access resolution, r_free, etc.

        Returns:
            list: List of dicts compatible with agent history format
        """
        history = []
        for cycle in self.data["cycles"]:
            if cycle.get("command") and cycle.get("result"):
                history_entry = {
                    "cycle_number": cycle.get("cycle_number", 0),
                    "job_id": str(cycle.get("cycle_number", "")),
                    "program": cycle.get("program", ""),
                    "command": cycle.get("command", ""),
                    "result": cycle.get("result", ""),
                    "summary": f"Decision: {cycle.get('decision', '')}\n"
                               f"Command: {cycle.get('command', '')}\n"
                               f"Result: {cycle.get('result', '')}",
                    # FIXED: Return metrics dict (not reasoning string) so workflow_state
                    # can access resolution, r_free, map_cc, etc.
                    "analysis": cycle.get("metrics", {}),
                    # NEW: Include output files so agent knows what was produced
                    "output_files": cycle.get("output_files", []),
                    "next_move": {
                        "command": cycle.get("command", ""),
                        "program": cycle.get("program", "")
                    }
                }
                history.append(history_entry)
        return history

    def get_all_output_files(self):
        """
        Get all output files from all successful cycles.

        Returns:
            list: List of file paths
        """
        all_files = []
        for cycle in self.data["cycles"]:
            if cycle.get("result", "").startswith("SUCCESS"):
                output_files = cycle.get("output_files", [])
                for f in output_files:
                    if f and f not in all_files:
                        all_files.append(f)
        return all_files

    def get_num_cycles(self):
        """Get the number of cycles recorded."""
        return len(self.data["cycles"])

    def format_cycle_summary(self, cycle_number):
        """
        Format a single cycle's information for display.

        Args:
            cycle_number: Which cycle to format

        Returns:
            str: Formatted cycle summary
        """
        cycle = self.get_cycle(cycle_number)
        if not cycle:
            return f"Cycle {cycle_number}: No data"

        lines = [
            f"",
            f"{'='*60}",
            f"CYCLE {cycle_number}",
            f"{'='*60}",
            f"Program: {cycle.get('program', 'N/A')}",
            f"Decision: {cycle.get('decision', 'N/A')}",
        ]

        # Truncate long reasoning
        reasoning = cycle.get('reasoning', 'N/A')
        if len(reasoning) > 600:
            reasoning = reasoning[:600] + "..."
        lines.append(f"Reasoning: {reasoning}")

        lines.extend([
            f"Command: {cycle.get('command', 'N/A')}",
            f"Result: {cycle.get('result', 'N/A')}",
        ])

        # Show output files if any
        output_files = cycle.get('output_files', [])
        if output_files:
            lines.append(f"Output files: {len(output_files)} files")
            for f in output_files[:5]:  # Show first 5
                lines.append(f"  - {os.path.basename(f)}")
            if len(output_files) > 5:
                lines.append(f"  ... and {len(output_files) - 5} more")

        return "\n".join(lines)

    def format_all_cycles(self):
        """
        Format all cycles for display.

        Returns:
            str: Formatted summary of all cycles
        """
        if not self.data["cycles"]:
            return "No cycles recorded."

        lines = [
            "",
            f"{'#'*60}",
            f"AGENT SESSION SUMMARY",
            f"{'#'*60}",
            f"Session ID: {self.data.get('session_id', 'N/A')}",
            f"Project Advice: {self.data.get('project_advice', 'None')}",
            f"Total Cycles: {sum(1 for c in self.data['cycles'] if c.get('program') not in ['STOP', None, 'unknown'])}",
        ]

        for i in range(1, len(self.data["cycles"]) + 1):
            lines.append(self.format_cycle_summary(i))

        if self.data.get("summary"):
            lines.extend([
                "",
                f"{'='*60}",
                "AI SUMMARY OF SESSION",
                f"{'='*60}",
                self.data["summary"]
            ])

        lines.append(f"{'#'*60}")
        return "\n".join(lines)

    def generate_log_for_summary(self):
        """
        Generate a log-like text suitable for LLM summarization.
        Consolidates metrics to show only final/best values.

        Returns:
            str: Log text for summarization
        """
        # Determine experiment type from files and history
        original_files = self.data.get('original_files', [])
        experiment_type = self._detect_experiment_type_for_summary(original_files)

        lines = [
            f"WORKING DIRECTORY:{self._get_working_directory_path()}",
            "COMMAND THAT WAS RUN: phenix.run_agent",
            "PHENIX AGENT SESSION LOG",
            f"Experiment Type: {experiment_type}",
            f"Project Advice: {self.data.get('project_advice', 'None')}",
            f"Original Files: {', '.join(original_files)}",
            f"Total Cycles: {sum(1 for c in self.data['cycles'] if c.get('program') not in ['STOP', None, 'unknown'])}",
        ]

        # Add resolution if set
        resolution = self.data.get("resolution")
        if resolution:
            lines.append(f"Resolution: {resolution:.2f} Å")

        lines.append("")

        # Collect key metrics by program (keep only final values)
        metrics_by_program = {}
        final_output_files = []
        programs_run = []

        for cycle in self.data["cycles"]:
            program = cycle.get('program', 'unknown')
            if program and program not in ['STOP', 'unknown', None]:
                if program not in programs_run:
                    programs_run.append(program)

            # Extract metrics from result text
            result = cycle.get('result', '')
            if isinstance(result, str):
                # Parse key metrics from result
                metrics = self._extract_metrics_from_result(result, program)
                if metrics:
                    metrics_by_program[program] = metrics  # Overwrite with latest

            # Collect output files (final cycle's files are most relevant)
            output_files = cycle.get('output_files', [])
            if output_files:
                final_output_files = output_files  # Keep updating to get latest

        # Programs run
        lines.extend([
            "="*40,
            "PROGRAMS RUN",
            "="*40,
            ", ".join(programs_run),
            "",
        ])

        # Consolidated metrics table (final values only)
        lines.extend([
            "="*40,
            "FINAL METRICS (most recent values)",
            "="*40,
        ])

        for program, metrics in metrics_by_program.items():
            if metrics:
                lines.append(f"\n{program}:")
                for key, value in metrics.items():
                    lines.append(f"  {key}: {value}")

        # Brief cycle summary (just program + outcome)
        lines.extend([
            "",
            "="*40,
            "CYCLE HISTORY (brief)",
            "="*40,
        ])

        for cycle in self.data["cycles"]:
            cycle_num = cycle.get('cycle_number', '?')
            program = cycle.get('program', 'N/A')
            result = cycle.get('result', 'N/A')
            # Truncate result to first line / 80 chars
            if isinstance(result, str):
                result_brief = result.split('\n')[0][:80]
            else:
                result_brief = str(result)[:80]
            lines.append(f"Cycle {cycle_num}: {program} - {result_brief}")

        # Final output files
        if final_output_files:
            lines.extend([
                "",
                "="*40,
                "KEY OUTPUT FILES",
                "="*40,
            ])
            for f in final_output_files[:10]:
                lines.append(f"  {os.path.basename(f)}")

        return "\n".join(lines)

    def _detect_experiment_type_for_summary(self, files):
        """
        Detect experiment type from files and history for summary generation.

        Args:
            files: List of file paths

        Returns:
            str: "Cryo-EM" or "X-ray Crystallography"
        """
        # Check if mtriage or real_space_refine was run (definitive cryo-EM indicators)
        for cycle in self.data.get("cycles", []):
            program = (cycle.get("program") or "").lower()
            if "mtriage" in program or "real_space_refine" in program:
                return "Cryo-EM"

        # Check file extensions
        has_mtz = False
        has_map = False

        for f in files:
            f_lower = f.lower()
            basename = os.path.basename(f_lower)

            # X-ray data files
            if f_lower.endswith(('.mtz', '.sca', '.hkl')):
                has_mtz = True

            # Map files (cryo-EM)
            if f_lower.endswith(('.ccp4', '.mrc', '.map')):
                has_map = True

        if has_map and not has_mtz:
            return "Cryo-EM"
        elif has_mtz:
            return "X-ray Crystallography"
        else:
            return "Unknown"

    def _extract_metrics_from_result(self, result_text, program):
        """Extract key metrics from result text using YAML patterns."""
        import re
        metrics = {}

        if not result_text:
            return metrics

        # =====================================================================
        # Use YAML-based extraction (centralized patterns)
        # =====================================================================
        try:
            from libtbx.langchain.knowledge.metric_patterns import extract_metrics_for_program
        except Exception:
            from knowledge.metric_patterns import extract_metrics_for_program
        yaml_metrics = extract_metrics_for_program(result_text, program)
        metrics.update(yaml_metrics)

        # =====================================================================
        # Fall back to hardcoded patterns for metrics not found via YAML
        # This provides backward compatibility and catches cases where
        # the YAML patterns might not be complete
        # =====================================================================
        # IMPORTANT: Use lowercase keys to match workflow_state expectations

        # Metrics that should use LAST match (can appear multiple times in refinement logs)
        # R-free and R-work appear after each macro cycle - we want the final values
        last_match_patterns = {
            'r_free': r'R.?free[:\s=]+([0-9.]+)',
            'r_work': r'R.?work[:\s=]+([0-9.]+)',
        }

        for name, pattern in last_match_patterns.items():
            if name in metrics:
                continue  # Already extracted from YAML
            matches = re.findall(pattern, result_text, re.IGNORECASE)
            if matches:
                try:
                    # Use the LAST match (final value from refinement)
                    metrics[name] = float(matches[-1])
                except ValueError:
                    pass

        # Metrics that use first match (typically appear once or first occurrence is correct)
        first_match_patterns = {
            'resolution': r'(?<!nomalous )(?<!nomalous  )[Rr]esolution\s*[=:]\s*([0-9.]+)',
            'clashscore': r'[Cc]lashscore[:\s=]+([0-9.]+)',
            'bonds_rmsd': r'bonds.?rmsd[:\s=]+([0-9.]+)',
            'angles_rmsd': r'angles.?rmsd[:\s=]+([0-9.]+)',
            'tfz': r'TFZ[:\s=]+([0-9.]+)',
            'llg': r'LLG[:\s=]+([0-9.]+)',
            # CC_mask is the primary format from real_space_refine
            # Also matches report format "Cc Mask" (from format_metrics_report)
            'map_cc': r'(?:CC[_ ]?mask|[Mm]ap[-_\s]?CC|[Mm]odel[-_\s]?vs[-_\s]?map\s+CC)\s*[:=]?\s*([0-9.]+)',
            'completeness': r'[Cc]ompleteness[:\s=]+([0-9.]+)',
            # NCS metrics from map_symmetry
            'ncs_cc': r'[Nn]cs\s*[Cc][Cc][:\s=]+([0-9.]+)',
            'ncs_copies': r'[Nn]cs\s*[Cc]opies[:\s=]+(\d+)',
        }

        for name, pattern in first_match_patterns.items():
            if name in metrics:
                continue  # Already extracted from YAML
            match = re.search(pattern, result_text, re.IGNORECASE)
            if match:
                try:
                    value = float(match.group(1)) if name != 'ncs_copies' else int(match.group(1))
                    # Skip obviously wrong values
                    if name == 'resolution' and value > 10:
                        continue  # Likely wrong (e.g., 47.31 is not resolution)
                    metrics[name] = value
                except ValueError:
                    pass

        # Symmetry type (string, not float) - only if not already extracted
        if 'symmetry_type' not in metrics:
            sym_match = re.search(r'[Ss]ymmetry\s*[Tt]ype[:\s=]+([A-Z0-9]+(?:\s*\([a-z]\))?)', result_text)
            if sym_match:
                metrics['symmetry_type'] = sym_match.group(1).strip()
            elif re.search(r'No.?[Ss]ymmetry|[Ss]ymmetry.?[Tt]ype[:\s=]+None', result_text, re.IGNORECASE):
                metrics['symmetry_type'] = "None"

        # Anomalous signal metrics (important for SAD/MAD workflow detection)
        # Pattern for anomalous measurability from xtriage output
        if 'anomalous_measurability' not in metrics:
            anom_match = re.search(r'[Aa]nomalous\s*[Mm]easurability[:\s=]+([0-9.]+)', result_text)
            if anom_match:
                try:
                    metrics['anomalous_measurability'] = float(anom_match.group(1))
                except ValueError:
                    pass

        # Anomalous resolution
        if 'anomalous_resolution' not in metrics:
            # Match both "Anomalous Resolution: 3.50" and "anomalous signal extends to 3.5 A"
            anom_res_match = re.search(
                r'[Aa]nomalous\s+(?:signal.*?(?:extends?\s+to|usable\s+to)\s+(?:about\s+)?)?[Rr]esolution[:\s=]+([0-9.]+)',
                result_text, re.IGNORECASE)
            if anom_res_match:
                try:
                    metrics['anomalous_resolution'] = float(anom_res_match.group(1))
                except ValueError:
                    pass

        # has_anomalous flag - check multiple patterns
        if 'has_anomalous' not in metrics:
            # Pattern 1: "Has Anomalous: True" from metrics report
            has_anom_match = re.search(r'[Hh]as\s+[Aa]nomalous[:\s=]+([Tt]rue|[Yy]es|1)', result_text)
            if has_anom_match:
                metrics['has_anomalous'] = True
            # Pattern 2: Strong measurability indicates usable anomalous signal
            elif metrics.get('anomalous_measurability', 0) > 0.05:
                metrics['has_anomalous'] = True
            # Pattern 3: Explicit statements in the output
            elif re.search(r'[Aa]nomalous\s+(?:data|signal).*?(?:present|significant|usable)', result_text, re.IGNORECASE):
                metrics['has_anomalous'] = True
            # Pattern 4: Has anomalous resolution implies usable anomalous
            elif metrics.get('anomalous_resolution'):
                metrics['has_anomalous'] = True

        return metrics

    def generate_summary(self, llm, use_rules_only=False):
        """
        Use LLM to generate a summary of the session (synchronous version).

        Args:
            llm: Language model to use
            use_rules_only: If True, skip LLM and return basic summary

        Returns:
            str: Generated summary
        """
        log_text = self.generate_log_for_summary()

        # Skip LLM if use_rules_only is set
        if use_rules_only:
            return "Summary generation skipped (use_rules_only=True). See session log for details."


        # Define the template using a standard string
        template = """You are an expert crystallographer reviewing an automated Phenix structure determination session.

Create a concise final report from the following session log.

IMPORTANT FORMATTING RULES:
- Do NOT repeat metrics multiple times - show only the FINAL values
- Keep the report brief and focused on outcomes
- Use a clean format without excessive duplication

Include these sections:

1. **Summary**: What was accomplished in 2-3 sentences

2. **Key Metrics**: List the most important FINAL metrics with their values. Format as a simple list:
   - For X-ray: R-work, R-free, Resolution, Clashscore
   - For Cryo-EM: Map-model CC, Resolution, Clashscore
   - Include any other relevant quality metrics (Ramachandran outliers, bonds/angles RMSD)
   Example format:
   - R-work: 0.182
   - R-free: 0.215
   - Resolution: 2.1 Å
   - Clashscore: 3.2

3. **Key Output Files**: List the most important output files with brief descriptions:
   - Final refined model (PDB filename)
   - Final data file with phases (MTZ filename)
   - Any other key outputs (maps, ligand files, etc.)
   Example format:
   - model_refine_001.pdb - Final refined atomic model
   - model_refine_001.mtz - Structure factors with R-free flags and map coefficients

4. **Status**: Whether the workflow completed successfully and any recommendations

SESSION LOG:
<session_log>
{log_text}
</session_log>

FINAL REPORT:"""

        # Create the safe prompt object
        prompt_template = PromptTemplate.from_template(template)

        # Generate the string safely
        # .format handles escaping if log_text contains { or }
        prompt = prompt_template.format(log_text=log_text)

        try:
            # Try to use rate limit handler
            try:
                from libtbx.langchain.agent.rate_limit_handler import get_google_handler
                handler = get_google_handler()
            except ImportError:
                try:
                    from agent.rate_limit_handler import get_google_handler
                    handler = get_google_handler()
                except ImportError:
                    handler = None

            if handler:
                def make_call():
                    return llm.invoke(prompt)

                response = handler.call_with_retry(make_call, lambda msg: print(msg))
            else:
                response = llm.invoke(prompt)

            summary = response.content.strip()
            self.data["summary"] = summary
            self.save()
            return summary
        except Exception as e:
            error_msg = f"Could not generate summary: {e}"
            print(error_msg)
            return error_msg

    def to_dict(self):
        """Return the session data as a dictionary."""
        return self.data.copy()

    def get_available_files(self):
        """
        Compute the list of available files from session data.

        This is the SINGLE SOURCE OF TRUTH for what files the agent can use.
        It combines:
        - original_files: The user's initial input files
        - output_files from each cycle: Files produced by completed cycles

        When cycles are removed (e.g., via session_utils --remove-last),
        this automatically returns the correct reduced file list.

        Returns:
            list: Absolute paths to all available files
        """
        import os
        files = []
        seen = set()  # Track by basename to avoid duplicates

        # 1. Start with original input files
        for f in self.data.get("original_files", []):
            if f:
                abs_path = os.path.abspath(f) if not os.path.isabs(f) else f
                basename = os.path.basename(abs_path)
                if basename not in seen:
                    files.append(abs_path)
                    seen.add(basename)

        # 2. Add output files from each cycle (in order)
        #    Uses _discover_cycle_outputs which tries stored paths first, then
        #    scans the expected output directory — handles restart scenarios
        #    where output_files is empty or paths are stale.
        for cycle in self.data.get("cycles", []):
            discovered = self._discover_cycle_outputs(cycle)
            for f in discovered:
                basename = os.path.basename(f)
                if basename not in seen:
                    files.append(f)
                    seen.add(basename)

            # Supplement: look for expected outputs that client may have missed
            # (derives companion files from known output file names/patterns)
            supplemental = self._find_missing_outputs(cycle, seen)
            for f in supplemental:
                abs_path = os.path.abspath(f) if not os.path.isabs(f) else f
                basename = os.path.basename(abs_path)
                if basename not in seen and os.path.exists(abs_path):
                    files.append(abs_path)
                    seen.add(basename)

        # 3. Catch-all: Scan agent sub-directories for output files not yet
        #    found.  Step 2's _discover_cycle_outputs handles known cycles;
        #    this catches files from cycles that may not be in session data
        #    (e.g., partially written session, manual directory copies).
        import glob
        pre_scan_count = len(files)  # Track boundary for Step 3.5
        agent_dirs = set()

        # PRIMARY: Use the session directory directly.
        session_dir = self._get_session_dir()
        if session_dir:
            agent_dirs.add(session_dir)

        # FALLBACK: Also infer agent dirs from any tracked files that happen
        # to be inside sub_* directories (original heuristic, kept for cases
        # where session_file isn't in the agent directory).
        for f in files:
            abs_f = os.path.abspath(f)
            parts = abs_f.replace("\\", "/").split("/")
            for i, part in enumerate(parts):
                if part.startswith("sub_") and "_" in part[4:]:
                    candidate = os.sep.join(parts[:i])
                    if os.path.isdir(candidate):
                        agent_dirs.add(candidate)
                    break
        for agent_dir in agent_dirs:
            # Refine outputs: map coefficients MTZ and refined models
            for entry in sorted(glob.glob(os.path.join(agent_dir, "sub_*_refine*"))):
                if not os.path.isdir(entry):
                    continue
                for mtz in glob.glob(os.path.join(entry, "refine_*.mtz")):
                    bn = os.path.basename(mtz)
                    if bn not in seen and not bn.endswith("_data.mtz"):
                        files.append(mtz)
                        seen.add(bn)
                for pdb in glob.glob(os.path.join(entry, "refine_*.pdb")):
                    bn = os.path.basename(pdb)
                    if bn not in seen:
                        files.append(pdb)
                        seen.add(bn)
            # Pdbtools outputs: *_with_ligand.pdb
            for entry in glob.glob(os.path.join(agent_dir, "sub_*_pdbtools")):
                if not os.path.isdir(entry):
                    continue
                for pdb in glob.glob(os.path.join(entry, "*_with_ligand.pdb")):
                    bn = os.path.basename(pdb)
                    if bn not in seen:
                        files.append(pdb)
                        seen.add(bn)
            # Autobuild outputs: overall_best files
            for entry in glob.glob(os.path.join(agent_dir, "sub_*_autobuild*")):
                if not os.path.isdir(entry):
                    continue
                for pattern in ["overall_best*.pdb", "overall_best*.mtz"]:
                    for match in glob.glob(os.path.join(entry, pattern)):
                        bn = os.path.basename(match)
                        if bn not in seen:
                            files.append(match)
                            seen.add(bn)

        # 3.5. Evaluate Step-3 discoveries through best_files.
        #      _rebuild_best_files_from_cycles() only processes files from known
        #      cycles.  Files discovered by the directory scan above may not be
        #      associated with any cycle (e.g., with_ligand PDBs from pdbtools
        #      when the cycle number doesn't match the directory number, or files
        #      from a previous agent run).  Without this step, best_files stays
        #      stale and the command builder selects the wrong model/MTZ.
        new_discoveries = files[pre_scan_count:]
        if new_discoveries and hasattr(self, 'best_files'):
            highest_cycle = max(
                (c.get("cycle_number", 0) for c in self.data.get("cycles", [])),
                default=0
            )
            for f in new_discoveries:
                # Let the tracker infer category and stage from filename.
                # Use highest_cycle so these files don't lose recency tiebreakers.
                self.best_files.evaluate_file(
                    path=f,
                    cycle=highest_cycle,
                    metrics=None,
                    stage=None
                )

        # 4. Filter out corrupt / zero-byte / invalid files.
        #    _is_valid_file() checks CCP4 magic bytes, PDB structural validity,
        #    etc.  Applying it here (client-side, where files are on disk) ensures
        #    the server receives only validated files, so its pass-through behavior
        #    (returning True for non-existent paths) never matters.
        try:
            try:
                from libtbx.langchain.agent.workflow_state import _is_valid_file
            except ImportError:
                from agent.workflow_state import _is_valid_file
            files = [f for f in files if _is_valid_file(f)]
        except Exception:
            pass  # Non-critical — if import fails, skip validation

        return files

    def _find_missing_outputs(self, cycle, already_seen):
        """
        Find expected output files that the client may not have tracked.

        Some clients don't discover all output files (e.g., map coefficients
        MTZ from refinement). This method uses program knowledge to look for
        expected outputs in the output directory.

        Args:
            cycle: Cycle dict with program, output_files, etc.
            already_seen: Set of basenames already tracked

        Returns:
            list: Additional file paths found on disk
        """
        import glob

        program = cycle.get("program", "")
        output_files = cycle.get("output_files", [])
        result = cycle.get("result", "")

        # Only supplement successful cycles
        if "FAILED" in result.upper():
            return []

        found = []

        # --- phenix.refine: look for map coefficients and refined model ---
        if "refine" in program.lower() and "real_space" not in program.lower():
            # Strategy 1: Derive output prefix from _data.mtz file
            # e.g., refine_001_data.mtz -> prefix = refine_001
            data_mtz_found = False
            for f in output_files:
                basename = os.path.basename(f)
                if basename.endswith("_data.mtz"):
                    data_mtz_found = True
                    prefix = basename[:-9]  # Strip "_data.mtz"
                    output_dir = os.path.dirname(os.path.abspath(f))

                    # Look for map coefficients MTZ: prefix.mtz or prefix_NNN.mtz
                    map_coeffs_patterns = [
                        os.path.join(output_dir, prefix + ".mtz"),
                        os.path.join(output_dir, prefix + "_001.mtz"),
                    ]
                    for pattern in map_coeffs_patterns:
                        if os.path.exists(pattern):
                            bn = os.path.basename(pattern)
                            if bn not in already_seen:
                                found.append(pattern)
                                break  # Only need one

                    # Look for refined model PDB: prefix.pdb or prefix_NNN.pdb
                    model_patterns = [
                        os.path.join(output_dir, prefix + ".pdb"),
                        os.path.join(output_dir, prefix + "_001.pdb"),
                    ]
                    for pattern in model_patterns:
                        if os.path.exists(pattern):
                            bn = os.path.basename(pattern)
                            if bn not in already_seen:
                                found.append(pattern)
                                break

                    # Also scan for any MTZ that's not _data.mtz
                    # (handles unexpected naming like prefix_002.mtz)
                    found_basenames = {os.path.basename(x) for x in found}
                    for mtz in glob.glob(os.path.join(output_dir, prefix + "*.mtz")):
                        bn = os.path.basename(mtz)
                        if (bn not in already_seen and
                            bn not in found_basenames and
                            not bn.endswith("_data.mtz")):
                            found.append(mtz)

                    break  # Only process one _data.mtz

            # Strategy 2: If no _data.mtz in output_files, try to derive the
            # output directory from ANY output file, or from the command itself.
            # This handles the case where the client tracked a PDB or .cif but
            # not the _data.mtz, or where output_files is empty entirely.
            if not data_mtz_found:
                output_dir = None
                # Try to get output dir from any output file
                for f in output_files:
                    if f and os.path.exists(f):
                        output_dir = os.path.dirname(os.path.abspath(f))
                        break
                # Try to derive from command's output_prefix
                if not output_dir:
                    import re
                    command = cycle.get("command", "")
                    prefix_match = re.search(r'output\.prefix\s*=\s*(\S+)', command)
                    if not prefix_match:
                        prefix_match = re.search(r'output_prefix\s*=\s*(\S+)', command)
                    if prefix_match:
                        prefix_path = prefix_match.group(1)
                        candidate_dir = os.path.dirname(os.path.abspath(prefix_path))
                        if os.path.isdir(candidate_dir):
                            output_dir = candidate_dir
                if output_dir and os.path.isdir(output_dir):
                    # Scan for refine output files
                    for mtz in glob.glob(os.path.join(output_dir, "refine_*.mtz")):
                        bn = os.path.basename(mtz)
                        if bn not in already_seen and not bn.endswith("_data.mtz"):
                            found.append(mtz)
                    for pdb in glob.glob(os.path.join(output_dir, "refine_*.pdb")):
                        bn = os.path.basename(pdb)
                        if bn not in already_seen:
                            found.append(pdb)

        # --- phenix.autobuild: look for overall_best model and map ---
        elif "autobuild" in program.lower() and "denmod" not in program.lower():
            for f in output_files:
                output_dir = os.path.dirname(os.path.abspath(f))
                # Look for overall_best files
                for pattern in ["overall_best*.pdb", "overall_best*.mtz"]:
                    for match in glob.glob(os.path.join(output_dir, pattern)):
                        bn = os.path.basename(match)
                        if bn not in already_seen:
                            found.append(match)
                break  # Only process first output file for directory

        return found

    def reset(self):
        """Reset the session (start fresh)."""
        self._init_new_session()
        self.save()

    # =========================================================================
    # AGENT SESSION SUMMARY
    # =========================================================================

    def generate_agent_session_summary(self, include_llm_assessment=True):
        """
        Generate a structured summary of the agent session.

        This produces a clean Markdown summary with:
        - Input files and user advice
        - Workflow path taken
        - Steps performed with key metrics
        - Final quality assessment
        - Output files

        Args:
            include_llm_assessment: If True, include placeholder for LLM assessment

        Returns:
            dict with:
                - markdown: The formatted Markdown summary
                - data: Structured data for LLM assessment
        """
        # Extract all structured data
        summary_data = self._extract_summary_data()

        # Generate Markdown
        markdown = self._format_summary_markdown(summary_data, include_llm_assessment)

        return {
            "markdown": markdown,
            "data": summary_data,
        }

    # ── Superseded-issue authority table ──────────────────────────────────────
    # Maps diagnosable error_type (from diagnosable_errors.yaml) to the set of
    # programs whose *success* supersedes the concern raised by that error.
    # Refinement carries the highest evidential weight because it is the most
    # sensitive test of model-data compatibility.
    _SUPERSESSION_AUTHORITY = {
        "crystal_symmetry_mismatch": {
            "phenix.refine", "phenix.real_space_refine",
        },
        "model_outside_map": {
            "phenix.dock_in_map", "phenix.real_space_refine",
        },
        # Default: any program that positions, builds, or refines a model
        # is authoritative for superseding generic probe/recoverable failures.
        # Expanded from {refine, real_space_refine} to include:
        #   phenix.phaser      — MR placement (supersedes needs_mr probe)
        #   phenix.autobuild   — model building (supersedes needs_build probe)
        #   phenix.dock_in_map — cryo-EM docking (supersedes needs_dock probe)
        # The _files_overlap check still guards each supersession — expanding
        # the authority set broadens candidacy, not provenance matching.
        "_default": {
            "phenix.refine", "phenix.real_space_refine",
            "phenix.phaser",
            "phenix.autobuild",
            "phenix.dock_in_map",
        },
    }

    @staticmethod
    def _stage_prefix_strip(basename):
        """Strip leading stage-prefix from a PHENIX output filename stem.
        Delegates to the module-level _strip_stage_prefix() function.
        """
        return _strip_stage_prefix(basename)

    @staticmethod
    def _build_lineage_graph(cycles):
        """Build a file provenance map from session cycle history.

        Maps each output file basename (lowercase) to the set of input
        file basenames that produced it.  Used by ``_files_overlap``
        Strategy 4 to trace transitive file ancestry across cycles.

        Example::

            Cycle 2: phenix.phaser beta.pdb data.mtz → PHASER.1.pdb
            Result:  {"phaser.1.pdb": {"beta.pdb", "data.mtz"}}

        Cycles whose ``output_files`` list is empty (e.g. older sessions
        recorded before output tracking was introduced) are silently skipped
        — the graph will simply have fewer entries and Strategies 1–3 still
        apply unchanged.

        Args:
            cycles: List of cycle dicts from ``session.data["cycles"]``.

        Returns:
            dict: ``{output_basename_lower: set_of_input_basename_lower}``
        """
        import os as _os
        produced_by = {}
        for cycle in cycles:
            inputs = set(_extract_file_basenames(cycle.get("command", "")))
            for f in cycle.get("output_files", []):
                bn = _os.path.basename(f).lower()
                produced_by[bn] = inputs
        return produced_by

    @staticmethod
    def _ancestors(bn, produced_by, visited=None):
        """Transitively resolve all ancestor file basenames for a given file.

        Walks the lineage graph built by ``_build_lineage_graph`` to find
        every input file that contributed — directly or indirectly — to
        producing ``bn``.

        Example::

            graph = {"phaser.1.pdb": {"beta.pdb"}, "refined.pdb": {"phaser.1.pdb"}}
            _ancestors("refined.pdb", graph)
            → {"phaser.1.pdb", "beta.pdb"}

        Cycle detection via ``visited`` prevents infinite recursion on
        degenerate or cyclic graphs.

        Args:
            bn:          Lowercase basename to resolve.
            produced_by: Lineage graph from ``_build_lineage_graph``.
            visited:     Internal set for cycle detection (callers pass None).

        Returns:
            set: All ancestor basenames (may be empty for leaf/unknown files).
        """
        if visited is None:
            visited = set()
        if bn in visited:
            return set()
        visited.add(bn)
        result = set()
        for parent in produced_by.get(bn, set()):
            result.add(parent)
            result |= AgentSession._ancestors(parent, produced_by, visited)
        return result

    @staticmethod
    def _files_overlap(probe_files, command_str, original_files, scope=None,
                       lineage_graph=None):
        """Return True if the probe and a later command share a semantic input file.

        Handles four matching strategies in order:
          1. Exact basename match
          2. Stage-prefix-stripped stem match
          3. Both trace back to the same session original file
          4. Transitive file lineage — a command file's ancestor (via the
             lineage graph) appears in the probe's input files.  Closes the
             case where an intermediate program (e.g. phaser) produces a
             derived file (PHASER.1.pdb) that a later authoritative program
             (e.g. refine) consumes, while the probe ran on the original
             input (beta.pdb) that phaser also used.

        Args:
            probe_files:    List of lowercase basenames from the probe command.
            command_str:    Raw command string of the later (authoritative) cycle.
            original_files: List of original session input file paths.
            scope:          Failure scope ("probe" | "recoverable" | "terminal").
                            When scope="probe" and either file list is empty,
                            returns True (probe programs run before derived files
                            exist and have no output history to compare against).
                            For all other scopes, empty file lists return False.
            lineage_graph:  Optional dict from ``_build_lineage_graph``.
                            When None, Strategy 4 is skipped.

        Returns:
            bool
        """
        import os as _os

        # Extract basenames from the later command
        cmd_files = _extract_file_basenames(command_str)

        if not probe_files or not cmd_files:
            # Empty file lists: only treat as overlapping for scope="probe".
            # Probe programs (e.g. model_vs_data) run early and may not yet
            # have derived output files, so we give them the benefit of the doubt.
            # For recoverable failures without parseable file provenance we
            # conservatively return False — don't auto-supersede without evidence.
            return scope == "probe"

        # Strategy 1: exact basename
        for pf in probe_files:
            if pf in cmd_files:
                return True

        # Strategy 2: stage-prefix-stripped stem
        probe_stems = {_strip_stage_prefix(f) for f in probe_files}
        cmd_stems   = {_strip_stage_prefix(f) for f in cmd_files}
        if probe_stems & cmd_stems:
            return True

        # Strategy 3: shared original file
        orig_basenames = {_os.path.basename(f).lower()
                          for f in (original_files or [])}
        if set(probe_files) & orig_basenames and set(cmd_files) & orig_basenames:
            return True

        # Strategy 4: transitive file lineage
        # Resolves cases where the later command operates on a file produced
        # from one of the probe's inputs in an intermediate cycle.
        # Example: probe ran on beta.pdb → phaser produced PHASER.1.pdb →
        #          refine ran on PHASER.1.pdb.  Ancestors of PHASER.1.pdb
        #          include beta.pdb, so overlap is confirmed.
        if lineage_graph:
            probe_set = set(probe_files)
            for cf in cmd_files:
                if AgentSession._ancestors(cf, lineage_graph) & probe_set:
                    return True

        return False

    def _annotate_superseded_issues(self, cycles):
        """Backward-scan cycles to mark probe/recoverable failures as superseded.

        For each failed cycle carrying a failure_context with scope in
        ("probe", "recoverable"), look forward for a successful authoritative
        program that ran using overlapping inputs.  If found, write
        failure_context["superseded_by"] = "<program> (cycle N)".

        File overlap is checked via four strategies (see ``_files_overlap``);
        Strategy 4 uses the lineage graph built here from ``output_files``
        to handle cases where an intermediate cycle produces a derived file
        that the authoritative cycle then consumes.

        Args:
            cycles: List of cycle dicts from session.data["cycles"].
                    Modified in-place.
        """
        lineage_graph = self._build_lineage_graph(cycles)
        original_files = self.data.get("original_files", [])

        for i, cycle in enumerate(cycles):
            fc = cycle.get("failure_context")
            if not fc:
                continue
            if fc.get("scope") not in ("probe", "recoverable"):
                continue
            if fc.get("superseded_by"):
                continue  # already annotated

            error_type  = fc.get("error_type") or "_default"
            probe_files = fc.get("input_files", [])
            authority   = (
                self._SUPERSESSION_AUTHORITY.get(error_type)
                or self._SUPERSESSION_AUTHORITY["_default"]
            )

            for j in range(i + 1, len(cycles)):
                later = cycles[j]
                if "SUCCESS" not in str(later.get("result", "")).upper():
                    continue
                later_prog = later.get("program", "")
                if later_prog not in authority:
                    continue
                later_cmd = later.get("command", "")
                if self._files_overlap(probe_files, later_cmd, original_files,
                                        scope=fc.get("scope"),
                                        lineage_graph=lineage_graph):
                    fc["superseded_by"] = "%s (cycle %s)" % (
                        later_prog, later.get("cycle_number", j + 1))
                    break

    def _extract_summary_data(self):
        """
        Extract structured summary data from the session.

        Returns:
            dict with all summary components
        """
        cycles = self.data.get("cycles", [])
        original_files = self.data.get("original_files", [])

        # Annotate superseded failures before building summary data.
        # This writes superseded_by onto failure_context dicts in-place.
        self._annotate_superseded_issues(cycles)

        # Detect experiment type
        experiment_type = self.data.get("experiment_type")
        if not experiment_type:
            experiment_type = self._detect_experiment_type_for_summary(original_files)

        # Determine workflow path
        workflow_path = self._determine_workflow_path(cycles, experiment_type)

        # Extract steps with metrics
        steps = self._extract_steps(cycles)

        # Get final metrics (from last successful refinement or analysis)
        final_metrics = self._extract_final_metrics(cycles, experiment_type)

        # Get final output files
        final_files = self._get_final_output_files(cycles)

        # Get input data quality metrics (from xtriage/mtriage)
        input_quality = self._extract_input_quality(cycles, experiment_type)

        # Get run name - use stored value first (set on client), fall back to computing it
        run_name = self.data.get("run_name") or self._get_run_name()

        # Check if this is a tutorial/focused task
        directives = self.get_directives()
        is_tutorial = False
        if directives:
            stop_cond = directives.get("stop_conditions", {})
            # It's a tutorial if there's an after_program stop condition with skip_validation
            if stop_cond.get("after_program") and stop_cond.get("skip_validation"):
                is_tutorial = True

        # Get preflight check results (v113)
        preflight_check = self.data.get("preflight_check", {})

        # Get stop info (captured during _extract_steps above)
        stop_info = getattr(self, '_last_stop_info', None)

        # Collect expert assessments across all steps
        expert_assessments = []
        for step in steps:
            ea = step.get("expert_assessment", {})
            if ea and isinstance(ea, dict) and ea.get("analysis"):
                expert_assessments.append({
                    "cycle": step["cycle"],
                    "program": step["program"],
                    "action": ea.get("action", ""),
                    "confidence": ea.get("confidence", ""),
                    "analysis": ea.get("analysis", ""),
                    "guidance": ea.get("guidance", ""),
                })

        return {
            "session_id": self.data.get("session_id", "unknown"),
            "experiment_type": experiment_type,
            "original_files": [os.path.basename(f) for f in original_files],
            "user_advice": self.data.get("project_advice", "None provided"),
            "resolution": self.data.get("resolution"),
            "workflow_path": workflow_path,
            "steps": steps,
            "final_metrics": final_metrics,
            "final_files": final_files,
            "input_quality": input_quality,
            # Exclude STOP cycles from counts - STOP is a decision, not a real program run
            "total_cycles": sum(1 for c in cycles if c.get("program") not in ["STOP", None, "unknown"]),
            "successful_cycles": sum(1 for c in cycles
                                    if c.get("program") not in ["STOP", None, "unknown"]
                                    and "SUCCESS" in str(c.get("result", ""))),
            "run_name": run_name,
            "is_tutorial": is_tutorial,
            # Job context
            "working_dir": self._get_working_directory_path(),
            "failure_diagnosis_path": self.data.get("failure_diagnosis_path"),
            # v113: thinking agent data
            "preflight_check": preflight_check,
            "stop_info": stop_info,
            "expert_assessments": expert_assessments,
        }

    def _get_run_name(self):
        """Get a descriptive name for this run based on directory.

        Returns the directory name where the session is running, which is
        typically descriptive (e.g., 'p9-xtriage', 'lysozyme-sad').
        """
        # Try input_directory first (from README location)
        input_dir = self.data.get("input_directory")
        if input_dir:
            # Convert to absolute path to handle "." properly
            abs_input_dir = os.path.abspath(input_dir)
            if os.path.isdir(abs_input_dir):
                return os.path.basename(abs_input_dir)

        # Try to get from session file path
        if hasattr(self, 'session_file') and self.session_file:
            session_dir = os.path.dirname(os.path.abspath(self.session_file))
            # Go up one level if we're in an 'agent_session' subdirectory
            dir_name = os.path.basename(session_dir)
            if dir_name in ['agent_session', 'session', '.phenix_ai', 'ai_agent_directory']:
                parent_dir = os.path.dirname(session_dir)
                if parent_dir:
                    return os.path.basename(parent_dir)
            return dir_name

        # Fall back to client_working_directory then current working directory
        client_wd = self.data.get("client_working_directory")
        if client_wd:
            return os.path.basename(client_wd)
        return os.path.basename(os.path.abspath(os.getcwd()))

    def _determine_workflow_path(self, cycles, experiment_type):
        """Determine the workflow path taken."""
        programs = [c.get("program", "") for c in cycles if c.get("program")]

        if experiment_type == "cryoem" or "Cryo-EM" in str(experiment_type):
            if "phenix.predict_and_build" in programs:
                return "Cryo-EM model building with AlphaFold prediction"
            elif "phenix.dock_in_map" in programs:
                return "Cryo-EM docking of existing model"
            else:
                return "Cryo-EM refinement"
        else:
            # X-ray
            if "phenix.autosol" in programs:
                return "X-ray experimental phasing (SAD/MAD)"
            elif "phenix.predict_and_build" in programs and "phenix.phaser" in programs:
                return "X-ray molecular replacement with AlphaFold prediction"
            elif "phenix.phaser" in programs:
                return "X-ray molecular replacement"
            elif "phenix.autobuild" in programs:
                return "X-ray automated model building"
            else:
                return "X-ray structure determination"

    def _extract_steps(self, cycles):
        """Extract step-by-step summary including reasoning and expert."""
        steps = []
        stop_info = None
        for cycle in cycles:
            program = cycle.get("program", "unknown")

            # Capture STOP cycle separately
            if program == "STOP":
                stop_info = {
                    "cycle": cycle.get("cycle_number", "?"),
                    "reasoning": cycle.get("reasoning", ""),
                    "explanation": cycle.get("explanation", ""),
                    "expert_assessment": cycle.get(
                        "expert_assessment", {}),
                }
                continue

            if program in [None, "unknown"]:
                continue

            result = cycle.get("result", "")
            success = "SUCCESS" in str(result)

            # Extract key metric for this step
            key_metric = self._get_key_metric_for_step(cycle, program)

            # Get reasoning (full, not just first sentence)
            reasoning = cycle.get("reasoning", "")
            decision = cycle.get("decision", "")

            # Get brief decision (first sentence) for table
            if decision:
                decision_brief = decision.split(".")[0] + "."
                if len(decision_brief) > 100:
                    decision_brief = decision_brief[:97] + "..."
            else:
                decision_brief = ""

            # Expert assessment
            expert = cycle.get("expert_assessment", {})

            fc = cycle.get("failure_context") or {}
            steps.append({
                "cycle": cycle.get("cycle_number", "?"),
                "program": program,
                "success": success,
                "key_metric": key_metric,
                "decision_brief": decision_brief,
                "reasoning": reasoning,
                "expert_assessment": expert if isinstance(
                    expert, dict) else {},
                "superseded_by": fc.get("superseded_by"),
            })

        # Attach stop_info to the data (will be picked up by
        # _extract_summary_data)
        self._last_stop_info = stop_info

        return steps

    def _get_key_metric_for_step(self, cycle, program):
        """
        Get the most important metric for a given step.

        Uses YAML configuration from metrics.yaml to determine what metric
        to display for each program type.
        """
        result = str(cycle.get("result", ""))

        # Check for failed runs FIRST - don't show misleading metrics
        if "FAILED" in result.upper() or "IGNORED" in result.upper():
            # Extract the error message if possible
            if ":" in result:
                error_part = result.split(":", 1)[1].strip()
                # Truncate long error messages
                if len(error_part) > 40:
                    error_part = error_part[:37] + "..."
                return f"FAILED: {error_part}"
            return "FAILED"

        try:
            from libtbx.langchain.knowledge.summary_display import format_step_metric
        except Exception:
            from knowledge.summary_display import format_step_metric
        try:
            from libtbx.langchain.knowledge.metric_patterns import extract_metrics_for_program
        except Exception:
            from knowledge.metric_patterns import extract_metrics_for_program

        # Use pre-parsed cycle metrics as primary source.
        # These are extracted by record_result -> _extract_metrics_from_result
        # using flexible hardcoded patterns that match both raw log output and
        # the reformatted metrics report (e.g., "R Free: 0.258" from
        # format_metrics_report). YAML patterns only match raw log format
        # (e.g., "R-free = 0.258") and fail on the reformatted text.
        metrics = dict(cycle.get("metrics", {}))

        # Supplement with YAML pattern extraction from result text
        # (catches anything the hardcoded patterns missed)
        yaml_metrics = extract_metrics_for_program(result, program)
        for k, v in yaml_metrics.items():
            if k not in metrics:
                metrics[k] = v

        # Format using YAML config
        formatted = format_step_metric(program, metrics)
        return formatted if formatted else ""

    def _extract_final_metrics(self, cycles, experiment_type):
        """Extract final quality metrics from the last relevant cycles."""
        metrics = {}

        # Work backwards through cycles to get most recent values
        for cycle in reversed(cycles):
            result = str(cycle.get("result", ""))
            program = (cycle.get("program") or "").lower()
            # Also check pre-parsed metrics stored in cycle
            cycle_metrics = cycle.get("metrics", {})

            # R-factors (X-ray)
            if "r_free" not in metrics:
                match = re.search(r'R.?[Ff]ree[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["r_free"] = float(match.group(1))
                elif cycle_metrics.get("r_free"):
                    metrics["r_free"] = float(cycle_metrics["r_free"])

            if "r_work" not in metrics:
                match = re.search(r'R.?[Ww]ork[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["r_work"] = float(match.group(1))
                elif cycle_metrics.get("r_work"):
                    metrics["r_work"] = float(cycle_metrics["r_work"])

            # Map CC (Cryo-EM)
            if "map_cc" not in metrics:
                match = re.search(r'[Mm]ap.?[Cc][Cc][:\s=]+([0-9.]+)', result)
                if match:
                    metrics["map_cc"] = float(match.group(1))
                elif cycle_metrics.get("map_cc"):
                    metrics["map_cc"] = float(cycle_metrics["map_cc"])

            # Map-model correlations (from phenix.map_correlations)
            # Matches both raw log "CC_mask  : 0.849" and report "Cc Mask: 0.849"
            if "cc_mask" not in metrics:
                match = re.search(r'CC[_ ]?mask\s*[=:]\s*([0-9.]+)', result, re.IGNORECASE)
                if match:
                    metrics["cc_mask"] = float(match.group(1))
                elif cycle_metrics.get("cc_mask"):
                    metrics["cc_mask"] = float(cycle_metrics["cc_mask"])

            if "cc_volume" not in metrics:
                match = re.search(r'CC[_ ]?volume\s*[=:]\s*([0-9.]+)', result, re.IGNORECASE)
                if match:
                    metrics["cc_volume"] = float(match.group(1))
                elif cycle_metrics.get("cc_volume"):
                    metrics["cc_volume"] = float(cycle_metrics["cc_volume"])

            if "cc_peaks" not in metrics:
                match = re.search(r'CC[_ ]?peaks\s*[=:]\s*([0-9.]+)', result, re.IGNORECASE)
                if match:
                    metrics["cc_peaks"] = float(match.group(1))
                elif cycle_metrics.get("cc_peaks"):
                    metrics["cc_peaks"] = float(cycle_metrics["cc_peaks"])

            if "cc_box" not in metrics:
                match = re.search(r'CC[_ ]?box\s*[=:]\s*([0-9.]+)', result, re.IGNORECASE)
                if match:
                    metrics["cc_box"] = float(match.group(1))
                elif cycle_metrics.get("cc_box"):
                    metrics["cc_box"] = float(cycle_metrics["cc_box"])

            # Geometry
            if "clashscore" not in metrics:
                match = re.search(r'[Cc]lashscore[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["clashscore"] = float(match.group(1))
                elif cycle_metrics.get("clashscore"):
                    metrics["clashscore"] = float(cycle_metrics["clashscore"])

            if "bonds_rmsd" not in metrics:
                match = re.search(r'[Bb]onds.?[Rr]msd[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["bonds_rmsd"] = float(match.group(1))

            if "angles_rmsd" not in metrics:
                match = re.search(r'[Aa]ngles.?[Rr]msd[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["angles_rmsd"] = float(match.group(1))
                elif cycle_metrics.get("angles_rmsd"):
                    metrics["angles_rmsd"] = float(cycle_metrics["angles_rmsd"])

            # Ramachandran statistics
            if "ramachandran_outliers" not in metrics:
                # Try common patterns: "Ramachandran outliers: 0.12%" or "rama_outliers=0.12"
                match = re.search(r'[Rr]ama(?:chandran)?.?[Oo]utliers?[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["ramachandran_outliers"] = float(match.group(1))
                elif cycle_metrics.get("ramachandran_outliers"):
                    metrics["ramachandran_outliers"] = float(cycle_metrics["ramachandran_outliers"])

            if "ramachandran_favored" not in metrics:
                match = re.search(r'[Rr]ama(?:chandran)?.?[Ff]avore?d?[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["ramachandran_favored"] = float(match.group(1))
                elif cycle_metrics.get("ramachandran_favored"):
                    metrics["ramachandran_favored"] = float(cycle_metrics["ramachandran_favored"])

            # Rotamer outliers
            if "rotamer_outliers" not in metrics:
                match = re.search(r'[Rr]otamer.?[Oo]utliers?[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["rotamer_outliers"] = float(match.group(1))
                elif cycle_metrics.get("rotamer_outliers"):
                    metrics["rotamer_outliers"] = float(cycle_metrics["rotamer_outliers"])

            # MolProbity score
            if "molprobity_score" not in metrics:
                match = re.search(r'[Mm]ol[Pp]robity.?[Ss]core[:\s=]+([0-9.]+)', result)
                if match:
                    metrics["molprobity_score"] = float(match.group(1))
                elif cycle_metrics.get("molprobity_score"):
                    metrics["molprobity_score"] = float(cycle_metrics["molprobity_score"])

            # Symmetry (from map_symmetry) - check both result text and pre-parsed metrics
            if "symmetry_type" not in metrics:
                match = re.search(r'[Ss]ymmetry\s*[Tt]ype[:\s=]+([A-Z0-9]+(?:\s*\([a-z]\))?)', result)
                if match:
                    metrics["symmetry_type"] = match.group(1).strip()
                elif re.search(r'No.?[Ss]ymmetry|[Ss]ymmetry.?[Tt]ype[:\s=]+None', result, re.IGNORECASE):
                    metrics["symmetry_type"] = "None"
                elif cycle_metrics.get("symmetry_type"):
                    metrics["symmetry_type"] = cycle_metrics["symmetry_type"]

            if "ncs_copies" not in metrics:
                match = re.search(r'[Nn]cs\s*[Cc]opies[:\s=]+(\d+)', result)
                if match:
                    metrics["ncs_copies"] = int(match.group(1))
                elif cycle_metrics.get("ncs_copies"):
                    metrics["ncs_copies"] = int(cycle_metrics["ncs_copies"])

            if "ncs_cc" not in metrics:
                match = re.search(r'[Nn]cs\s*[Cc][Cc][:\s=]+([0-9.]+)', result)
                if match:
                    metrics["ncs_cc"] = float(match.group(1))
                elif cycle_metrics.get("ncs_cc"):
                    metrics["ncs_cc"] = float(cycle_metrics["ncs_cc"])

        # Add quality assessments
        if "r_free" in metrics:
            metrics["r_free_assessment"] = self._assess_r_free(
                metrics["r_free"], self.data.get("resolution"))

        if "map_cc" in metrics:
            metrics["map_cc_assessment"] = self._assess_map_cc(metrics["map_cc"])

        if "cc_mask" in metrics:
            metrics["cc_mask_assessment"] = self._assess_map_cc(metrics["cc_mask"])

        if "clashscore" in metrics:
            metrics["clashscore_assessment"] = self._assess_clashscore(metrics["clashscore"])

        return metrics

    def _assess_r_free(self, r_free, resolution):
        """Assess R-free quality based on resolution."""
        if resolution is None:
            resolution = 2.5  # Default assumption

        # Resolution-dependent thresholds
        if resolution < 1.5:
            good, acceptable = 0.20, 0.25
        elif resolution < 2.5:
            good, acceptable = 0.25, 0.30
        elif resolution < 3.5:
            good, acceptable = 0.30, 0.35
        else:
            good, acceptable = 0.35, 0.40

        if r_free <= good:
            return "Good"
        elif r_free <= acceptable:
            return "Acceptable"
        else:
            return "Needs improvement"

    def _assess_map_cc(self, map_cc):
        """Assess map-model correlation."""
        if map_cc >= 0.80:
            return "Good"
        elif map_cc >= 0.70:
            return "Acceptable"
        else:
            return "Needs improvement"

    def _assess_clashscore(self, clashscore):
        """Assess geometry clashscore."""
        if clashscore <= 5:
            return "Excellent"
        elif clashscore <= 10:
            return "Good"
        elif clashscore <= 20:
            return "Acceptable"
        else:
            return "Needs improvement"

    def _get_final_output_files(self, cycles):
        """
        Get the most relevant final output files with descriptions.

        Shows ONLY the best files for each category - the actual deliverables
        the user should use. Not all intermediate outputs.
        """
        final_files = []
        seen_types = set()  # Track which types we've added

        # Use best_files from session - these are the definitive outputs
        # Handle both old format (best_files.category) and new format (best_files.best.category)
        best_files_raw = self.data.get("best_files", {})
        if "best" in best_files_raw:
            # New format: best_files.best.category.path
            best_files = {}
            for category, info in best_files_raw.get("best", {}).items():
                if isinstance(info, dict) and "path" in info:
                    best_files[category] = info["path"]
        else:
            # Old format: best_files.category = path
            best_files = best_files_raw

        if best_files:
            # Define what we want to show and in what order
            # Each entry: (category, display_type, description_template)
            categories_to_show = [
                ("model", "model", "Best refined model"),
                ("data_mtz", "data", "Reflection data with R-free flags"),
                ("map_coeffs_mtz", "map_coeffs", "Map coefficients for visualization"),
                ("map", "map", "Best electron density map"),
                ("full_map", "map", "Full reconstructed map"),
                # Note: ligandfit output is found via cycle-scanning below,
                # because "ligand_fit_output" is a stage under "model" in
                # best_files, not a separate top-level category key.
            ]

            for category, display_type, description in categories_to_show:
                filepath = best_files.get(category)
                if filepath:
                    basename = os.path.basename(filepath)
                    # Skip intermediate files
                    if self._is_intermediate_file(basename):
                        continue
                    # Note: We don't check file existence here because:
                    # 1. The session may have been created on a different machine
                    # 2. The paths in best_files are from the client's filesystem
                    final_files.append({
                        "name": basename,
                        "path": filepath,
                        "type": display_type,
                        "description": description,
                    })
                    seen_types.add(display_type)

            # Ligandfit outputs need special handling:
            # ligand_fit_output is a subcategory/stage under "model" in best_files,
            # NOT a separate top-level key. Since the ligandfit output PDB scores
            # lower than a refined model, it won't be the best "model" entry.
            # Scan cycles to find ligandfit output files.
            if "ligand" not in seen_types:
                for cycle in reversed(cycles):
                    program = (cycle.get("program") or "").lower()
                    result = str(cycle.get("result", ""))
                    if "ligandfit" in program and "SUCCESS" in result:
                        # First pass: look for ligand_fit pattern files
                        ligand_file = None
                        ligand_path = None
                        fallback_pdb = None
                        fallback_path = None
                        for f in cycle.get("output_files", []):
                            basename = os.path.basename(f)
                            if self._is_intermediate_file(basename):
                                continue
                            if basename.endswith('.pdb'):
                                if 'ligand_fit' in basename or 'lig_fit' in basename:
                                    ligand_file = basename
                                    ligand_path = f
                                    break
                                elif fallback_pdb is None:
                                    fallback_pdb = basename
                                    fallback_path = f
                        chosen = ligand_file or fallback_pdb
                        chosen_path = ligand_path or fallback_path
                        if chosen:
                            final_files.append({
                                "name": chosen,
                                "path": chosen_path,
                                "type": "ligand",
                                "description": "Fitted ligand coordinates",
                            })
                            seen_types.add("ligand")
                        if "ligand" in seen_types:
                            break

            # If we have best files, return them
            if final_files:
                return final_files

        # FALLBACK: If no best_files in session, scan all cycle outputs
        # Collect best file for each type from remaining cycles
        best_by_type = {}  # type -> (basename, file_info, cycle_num)

        for i, cycle in enumerate(cycles):
            result = str(cycle.get("result", ""))
            if "SUCCESS" not in result:
                continue
            output_files = cycle.get("output_files", [])
            program = cycle.get("program", "")

            for f in output_files:
                basename = os.path.basename(f)
                if self._is_intermediate_file(basename):
                    continue

                file_info = self._describe_output_file(basename, program)
                if file_info:
                    file_info["path"] = f  # Store full path
                    file_type = file_info.get("type", "file")
                    # Keep the latest (highest cycle number) for each type
                    if file_type not in best_by_type or i > best_by_type[file_type][2]:
                        best_by_type[file_type] = (basename, file_info, i)

        # Build final list from best_by_type, prioritizing certain types
        priority_order = ["model", "data", "map_coeffs", "map", "ligand"]
        for ptype in priority_order:
            if ptype in best_by_type:
                _, file_info, _ = best_by_type[ptype]
                if not any(existing["name"] == file_info["name"] for existing in final_files):
                    final_files.append(file_info)

        # Add any remaining types not in priority order
        for ftype, (basename, file_info, _) in best_by_type.items():
            if ftype not in priority_order:
                if not any(existing["name"] == file_info["name"] for existing in final_files):
                    final_files.append(file_info)

        return final_files[:5]  # Limit to 5 files

    def _is_intermediate_file(self, basename):
        """Check if a file is an intermediate/temporary file that shouldn't be shown."""
        basename_lower = basename.lower()

        # Skip patterns - intermediate or debug files
        # NOTE: All patterns must be lowercase since they're matched against basename_lower
        skip_patterns = [
            "_debug", "_check", ".geo", ".eff", ".def", ".log",
            "superposed_predicted", "untrimmed", "intermediate",
            "_cycle_", "carryon",
        ]
        if any(pat in basename_lower for pat in skip_patterns):
            return True

        # Skip files that end with specific intermediate suffixes
        if basename_lower.endswith(("_check.pdb", "_debug.pdb")):
            return True

        return False

    def _get_working_directory_path(self):
        """
        Get the absolute path of the working directory where the agent was run.

        Priority:
        1. client_working_directory (stored at client startup — reliable on server)
        2. input_directory (from user param or README discovery)
        3. Infer from session file path
        4. Fall back to os.getcwd() (only correct on the client machine)
        """
        # Try client_working_directory first (set by client, survives server round-trip)
        client_wd = self.data.get("client_working_directory")
        if client_wd:
            return client_wd

        # Try input_directory (set when agent starts)
        input_dir = self.data.get("input_directory")
        if input_dir:
            abs_dir = os.path.abspath(input_dir)
            if os.path.isdir(abs_dir):
                return abs_dir

        # Infer from session file path
        if hasattr(self, 'session_file') and self.session_file:
            session_dir = os.path.dirname(os.path.abspath(self.session_file))
            dir_name = os.path.basename(session_dir)
            if dir_name in ['agent_session', 'session', '.phenix_ai',
                            'ai_agent_directory']:
                parent_dir = os.path.dirname(session_dir)
                if parent_dir:
                    return parent_dir
            return session_dir

        return os.path.abspath(os.getcwd())

    def _get_display_path(self, file_entry):
        """
        Get a display path for a file entry that is relative to the working
        directory, so the user can type it directly at the command line.

        Args:
            file_entry: dict with 'name' (basename) and optional 'path' (full path)

        Returns:
            str: Relative path from working directory, or basename as fallback
        """
        full_path = file_entry.get("path")
        if not full_path:
            return file_entry.get("name", "unknown")

        try:
            working_dir = self._get_working_directory_path()
            rel_path = os.path.relpath(full_path, working_dir)
            # Sanity check: if relpath starts with many "../" it's not useful
            if rel_path.count("..") > 3:
                return file_entry.get("name", "unknown")
            return rel_path
        except (ValueError, TypeError):
            # os.path.relpath can fail on Windows with different drives
            return file_entry.get("name", "unknown")

    def _describe_output_file(self, basename, program):
        """
        Generate a description for an output file based on its name and source program.

        Returns dict with name, type, and description, or None if file should be skipped.
        """
        basename_lower = basename.lower()

        # Skip intermediate/temporary files
        skip_patterns = ["_debug", "_check", ".geo", ".eff", ".def", ".log"]
        if any(pat in basename_lower for pat in skip_patterns):
            return None

        # PDB files
        if basename.endswith(".pdb") or basename.endswith(".cif"):
            if "refine" in basename_lower and "real_space" not in basename_lower:
                return {
                    "name": basename,
                    "type": "model",
                    "description": "Refined atomic model (X-ray)",
                    "priority": True
                }
            elif "real_space_refined" in basename_lower or "rsr" in basename_lower:
                return {
                    "name": basename,
                    "type": "model",
                    "description": "Refined atomic model (cryo-EM)",
                    "priority": True
                }
            elif "overall_best" in basename_lower or "autobuild" in basename_lower:
                return {
                    "name": basename,
                    "type": "model",
                    "description": "Best model from automated building",
                    "priority": True
                }
            elif "phaser" in basename_lower:
                return {
                    "name": basename,
                    "type": "model",
                    "description": "Molecular replacement solution"
                }
            elif "autosol" in basename_lower:
                return {
                    "name": basename,
                    "type": "model",
                    "description": "Experimental phasing solution"
                }
            elif "predict" in basename_lower:
                return {
                    "name": basename,
                    "type": "model",
                    "description": "AlphaFold predicted model"
                }
            elif "dock" in basename_lower or "placed" in basename_lower:
                return {
                    "name": basename,
                    "type": "model",
                    "description": "Docked model in map"
                }
            elif "with_ligand" in basename_lower or (
                    "_modified" in basename_lower and
                    basename_lower.endswith(".pdb")):
                return {
                    "name": basename,
                    "type": "model",
                    "description": "Model with fitted ligand"
                }
            elif "ligand_fit" in basename_lower:
                return {
                    "name": basename,
                    "type": "ligand",
                    "description": "Fitted ligand coordinates"
                }
            else:
                return {
                    "name": basename,
                    "type": "model",
                    "description": "Atomic model"
                }

        # MTZ files
        elif basename.endswith(".mtz"):
            if "refine" in basename_lower:
                return {
                    "name": basename,
                    "type": "data",
                    "description": "Refined structure factors with R-free flags and map coefficients",
                    "priority": True
                }
            elif "phaser" in basename_lower:
                return {
                    "name": basename,
                    "type": "data",
                    "description": "Phased data from molecular replacement"
                }
            elif "autosol" in basename_lower:
                return {
                    "name": basename,
                    "type": "data",
                    "description": "Phased data from experimental phasing"
                }
            elif "polder" in basename_lower:
                return {
                    "name": basename,
                    "type": "data",
                    "description": "Polder omit map coefficients"
                }
            else:
                return {
                    "name": basename,
                    "type": "data",
                    "description": "Reflection data"
                }

        # Map files
        elif basename.endswith((".mrc", ".ccp4", ".map")):
            if "sharpened" in basename_lower:
                return {
                    "name": basename,
                    "type": "map",
                    "description": "Sharpened/optimized cryo-EM map"
                }
            elif "polder" in basename_lower:
                return {
                    "name": basename,
                    "type": "map",
                    "description": "Polder omit map"
                }
            else:
                return {
                    "name": basename,
                    "type": "map",
                    "description": "Electron density map"
                }

        # NCS spec files
        elif basename.endswith(".ncs_spec"):
            return {
                "name": basename,
                "type": "symmetry",
                "description": "NCS operators for map symmetry"
            }

        # CIF restraint files
        elif basename.endswith(".cif") and "ligand" in basename_lower:
            return {
                "name": basename,
                "type": "restraints",
                "description": "Ligand restraint dictionary"
            }

        return None

    def _extract_input_quality(self, cycles, experiment_type):
        """Extract input data quality metrics from xtriage/mtriage.

        Uses the pre-parsed metrics stored in cycle["metrics"] from the log parsers,
        falling back to parsing the result string if metrics are not available.
        """
        quality = {}

        for cycle in cycles:
            program = (cycle.get("program") or "").lower()
            metrics = cycle.get("metrics", {})
            result = str(cycle.get("result", ""))

            if "xtriage" in program:
                # First try pre-parsed metrics (from log_parsers._extract_xtriage_metrics)
                if metrics.get("resolution"):
                    quality["resolution"] = metrics["resolution"]

                # Completeness - handle both string "88.7%" and float formats
                if metrics.get("completeness"):
                    comp = metrics["completeness"]
                    if isinstance(comp, str):
                        # Parse "88.7%" format
                        comp_match = re.match(r'([\d.]+)', comp)
                        if comp_match:
                            quality["completeness"] = float(comp_match.group(1))
                    else:
                        quality["completeness"] = float(comp)

                # Twin fraction and twinning detection
                if metrics.get("twin_fraction") is not None:
                    tf = metrics["twin_fraction"]
                    if tf > 0.1:
                        quality["twinning"] = f"Possible (fraction: {tf:.3f})"
                    else:
                        quality["twinning"] = f"Low (fraction: {tf:.3f})"
                    quality["twin_fraction"] = tf
                elif metrics.get("twinning"):
                    quality["twinning"] = metrics["twinning"]
                elif metrics.get("no_twinning_suspected"):
                    quality["twinning"] = "None detected"

                # Twin law
                if metrics.get("twin_law"):
                    quality["twin_law"] = metrics["twin_law"]

                # Wilson B-factor
                if metrics.get("wilson_b"):
                    quality["wilson_b"] = metrics["wilson_b"]

                # I/sigI
                if metrics.get("i_over_sigma"):
                    quality["i_over_sigma"] = metrics["i_over_sigma"]

                # Anomalous signal info
                if metrics.get("anomalous_measurability"):
                    quality["anomalous_measurability"] = metrics["anomalous_measurability"]
                if metrics.get("has_anomalous"):
                    quality["has_anomalous"] = True
                if metrics.get("anomalous_resolution"):
                    quality["anomalous_resolution"] = metrics["anomalous_resolution"]

                # Fallback to parsing result string if metrics not found
                if not quality.get("resolution"):
                    match = re.search(r'(?<!nomalous )(?<!nomalous  )Resolution[:\s=]+([0-9.]+)', result, re.IGNORECASE)
                    if match:
                        res = float(match.group(1))
                        if res < 15:  # d_min is typically < 10 Å
                            quality["resolution"] = res

                if not quality.get("completeness"):
                    match = re.search(r'Completeness[:\s=]+([0-9.]+)', result, re.IGNORECASE)
                    if match:
                        quality["completeness"] = float(match.group(1))

                if not quality.get("twinning"):
                    if "No Twinning" in result or "twin" not in result.lower():
                        quality["twinning"] = "None detected"

                break  # Only need first xtriage

            elif "mtriage" in program:
                # First try pre-parsed metrics
                if metrics.get("d_fsc"):
                    quality["resolution"] = metrics["d_fsc"]
                elif metrics.get("resolution"):
                    quality["resolution"] = metrics["resolution"]

                if metrics.get("d99"):
                    quality["d99"] = metrics["d99"]

                if metrics.get("map_cc"):
                    quality["map_cc"] = metrics["map_cc"]

                # Fallback to parsing result string
                if not quality.get("resolution"):
                    match = re.search(r'd_fsc[:\s=]+([0-9.]+)', result, re.IGNORECASE)
                    if match:
                        quality["resolution"] = float(match.group(1))

                if not quality.get("d99"):
                    match = re.search(r'd99[:\s=]+([0-9.]+)', result, re.IGNORECASE)
                    if match:
                        quality["d99"] = float(match.group(1))

                break  # Only need first mtriage

        return quality

    def _format_summary_markdown(self, data, include_llm_assessment):
        """Format the summary data as Markdown.

        Produces a structured report including:
        - Input files, experiment type, user goal
        - Pre-flight check results (if thinking agent was active)
        - Input data quality metrics
        - Workflow path and detailed step-by-step narrative
        - Expert crystallographer assessments (if thinking agent)
        - Stop decision rationale (if workflow stopped early)
        - Final quality metrics
        - Output files
        - LLM assessment placeholder (optional)
        """
        lines = []

        # Header
        run_name = data.get('run_name', 'unknown')
        is_tutorial = data.get('is_tutorial', False)
        if is_tutorial:
            title = "Phenix AI Tutorial: %s" % run_name
        else:
            title = "Phenix AI Run: %s" % run_name
        lines.append("## %s" % title)
        lines.append("")

        # Session info
        session_id = data['session_id']
        if '_' in session_id:
            date_part, time_part = session_id.split('_', 1)
            time_display = time_part.replace('-', ':')
            session_display = "%s %s" % (date_part, time_display)
        else:
            session_display = session_id
        lines.append(
            "**Session:** %s | **Cycles:** %d (%d successful)"
            % (session_display, data['total_cycles'],
               data['successful_cycles']))
        working_dir = data.get("working_dir")
        if working_dir:
            lines.append("")
            lines.append("**Working directory:** `%s`" % working_dir)
        if include_llm_assessment:
            lines.append("")
            lines.append("_[STATUS_PLACEHOLDER]_")
        lines.append("")

        # ── 1. Input ────────────────────────────────────────────
        lines.append("## Input")
        lines.append("")
        lines.append(
            "- **Files:** %s" % ', '.join(data['original_files']))
        lines.append("- **User Advice:** %s" % data['user_advice'])
        lines.append(
            "- **Experiment Type:** %s" % data['experiment_type'])
        if data.get("resolution"):
            lines.append(
                "- **Resolution:** %.2f \u00C5" % data['resolution'])

        directives = self.get_directives()
        if directives:
            stop_cond = directives.get("stop_conditions", {})
            if stop_cond:
                stop_parts = []
                if stop_cond.get("after_program"):
                    stop_parts.append(
                        "stop after %s" % stop_cond['after_program'])
                if stop_cond.get("after_cycle"):
                    stop_parts.append(
                        "stop after cycle %s" % stop_cond['after_cycle'])
                if stop_parts:
                    lines.append(
                        "- **Stop Condition (Focused Task):** %s"
                        % ', '.join(stop_parts))
        lines.append("")

        # ── 2. Pre-flight Check ─────────────────────────────────
        preflight = data.get("preflight_check", {})
        if preflight and isinstance(preflight, dict):
            pf_ready = preflight.get("ready", True)
            pf_analysis = preflight.get("analysis", "")
            pf_missing = preflight.get("missing_files", [])
            pf_warnings = preflight.get("warnings", [])
            pf_recommendation = preflight.get("recommendation", "")

            if not pf_ready or pf_missing or pf_warnings:
                if not pf_ready:
                    lines.append(
                        "## \u26D4 Pre-flight Check: Missing Inputs")
                else:
                    lines.append(
                        "## \u26A0 Pre-flight Check: Warnings")
                lines.append("")
                if pf_analysis:
                    lines.append(pf_analysis.strip())
                    lines.append("")
                if pf_missing:
                    lines.append("**Missing files:**")
                    for m in pf_missing[:5]:
                        lines.append("- %s" % str(m))
                    lines.append("")
                if pf_warnings:
                    for w in pf_warnings[:3]:
                        lines.append("> **Note:** %s" % str(w))
                    lines.append("")
                if pf_recommendation:
                    lines.append(
                        "**Recommendation:** %s"
                        % pf_recommendation.strip())
                    lines.append("")
                if not pf_ready:
                    lines.append(
                        "*The workflow was stopped because "
                        "required input files are missing.*")
                    lines.append("")

        # ── 3. Input Data Quality ───────────────────────────────
        if data.get("input_quality"):
            lines.append("## Input Data Quality")
            lines.append("")
            iq = data["input_quality"]
            if "resolution" in iq:
                lines.append(
                    "- **Resolution:** %.2f \u00C5" % iq['resolution'])
            if "completeness" in iq:
                lines.append(
                    "- **Completeness:** %.1f%%" % iq['completeness'])
            if "i_over_sigma" in iq:
                lines.append(
                    "- **I/\u03C3(I):** %.1f" % iq['i_over_sigma'])
            if "wilson_b" in iq:
                lines.append(
                    "- **Wilson B-factor:** %.1f \u00C5\u00B2"
                    % iq['wilson_b'])
            if "twinning" in iq:
                lines.append(
                    "- **Twinning:** %s" % iq['twinning'])
            if "twin_law" in iq:
                lines.append(
                    "- **Twin Law:** %s" % iq['twin_law'])
            if "twin_fraction" in iq:
                lines.append(
                    "- **Twin Fraction:** %.4f" % iq['twin_fraction'])
            if iq.get("has_anomalous"):
                anom_info = "Yes"
                if "anomalous_measurability" in iq:
                    anom_info += (
                        " (measurability: %.3f)"
                        % iq['anomalous_measurability'])
                if "anomalous_resolution" in iq:
                    anom_info += (
                        ", extends to %.2f \u00C5"
                        % iq['anomalous_resolution'])
                lines.append(
                    "- **Anomalous Signal:** %s" % anom_info)
            if "d99" in iq:
                lines.append("- **d99:** %.2f \u00C5" % iq['d99'])
            if "map_cc" in iq:
                lines.append(
                    "- **Map CC:** %.3f" % iq['map_cc'])
            lines.append("")

        # ── 4. Workflow ─────────────────────────────────────────
        lines.append("## Workflow Path")
        lines.append("")
        lines.append(data["workflow_path"])
        lines.append("")

        # ── 5. Steps Performed (detailed) ───────────────────────
        lines.append("## Steps Performed")
        lines.append("")

        _has_superseded_step = any(
            step.get("superseded_by") for step in data["steps"])

        # Summary table first
        lines.append("| Cycle | Program | Result | Key Metric |")
        lines.append("|-------|---------|--------|------------|")
        for step in data["steps"]:
            if step["success"]:
                result_symbol = "\u2713"
            elif step.get("superseded_by"):
                result_symbol = "\u2717*"
            else:
                result_symbol = "\u2717"
            program_short = step["program"].replace("phenix.", "")
            key_metric = step['key_metric']
            key_metric = key_metric.replace('\n', ' ').replace(
                '\r', ' ')
            key_metric = key_metric.replace('|', '-')
            key_metric = ' '.join(key_metric.split())
            if len(key_metric) > 60:
                key_metric = key_metric[:57] + "..."
            lines.append(
                "| %s | %s | %s | %s |"
                % (step['cycle'], program_short,
                   result_symbol, key_metric))
        if _has_superseded_step:
            lines.append("")
            lines.append(
                "_\\* Failure resolved by a later successful "
                "step \u2014 does not affect the final result._")
        lines.append("")

        # Detailed narrative per step (reasoning + expert)
        _has_detail = any(
            step.get("reasoning") or step.get("expert_assessment")
            for step in data["steps"])
        if _has_detail:
            lines.append("### Step Details")
            lines.append("")
            for step in data["steps"]:
                program_short = step["program"].replace("phenix.", "")
                if step["success"]:
                    status_str = "\u2713"
                else:
                    status_str = "\u2717"
                lines.append(
                    "**Cycle %s: %s** %s"
                    % (step['cycle'], program_short, status_str))
                lines.append("")

                # Expert assessment (shown before reasoning since
                # it evaluates the PREVIOUS step's output)
                ea = step.get("expert_assessment", {})
                if ea and isinstance(ea, dict):
                    analysis = ea.get("analysis", "").strip()
                    guidance = ea.get("guidance", "").strip()
                    action = ea.get("action", "")
                    confidence = ea.get("confidence", "")
                    if analysis:
                        tag = ""
                        if action and confidence:
                            tag = " (%s, %s)" % (action, confidence)
                        lines.append(
                            "> **Expert Assessment%s:** %s" % (
                                tag, analysis[:600]))
                        if guidance:
                            lines.append(">")
                            lines.append(
                                "> **Guidance:** %s"
                                % guidance[:300])
                        lines.append("")

                # Agent reasoning
                reasoning = step.get("reasoning", "").strip()
                if reasoning:
                    display_reasoning = reasoning
                    if len(display_reasoning) > 500:
                        display_reasoning = (
                            display_reasoning[:497] + "...")
                    lines.append(
                        "**Reasoning:** %s" % display_reasoning)
                    lines.append("")

        # ── 6. Stop Decision ────────────────────────────────────
        stop_info = data.get("stop_info")
        if stop_info and isinstance(stop_info, dict):
            lines.append("## Stop Decision")
            lines.append("")
            stop_reasoning = (
                stop_info.get("reasoning")
                or stop_info.get("explanation", "")).strip()
            if stop_reasoning:
                if len(stop_reasoning) > 600:
                    stop_reasoning = stop_reasoning[:597] + "..."
                lines.append(
                    "**Why the workflow stopped:** %s"
                    % stop_reasoning)
                lines.append("")

            # Expert stop review
            stop_expert = stop_info.get("expert_assessment", {})
            if stop_expert and isinstance(stop_expert, dict):
                stop_review = stop_expert.get("stop_review", {})
                stop_analysis = stop_expert.get(
                    "analysis", "").strip()

                if stop_analysis:
                    lines.append(
                        "> **Expert analysis of stop:** %s"
                        % stop_analysis[:500])
                    lines.append("")

                if stop_review and isinstance(stop_review, dict):
                    agrees = stop_review.get(
                        "agree_with_stop", True)
                    review_analysis = stop_review.get(
                        "analysis", "").strip()
                    alternatives = stop_review.get(
                        "alternatives", [])
                    review_guidance = stop_review.get(
                        "guidance", "").strip()

                    if agrees:
                        verdict = "\u2705 Expert **agrees** " \
                                  "with stopping."
                    else:
                        verdict = "\u26A0 Expert **disagrees** " \
                                  "with stopping."
                    lines.append(verdict)
                    lines.append("")

                    if review_analysis:
                        lines.append(review_analysis[:500])
                        lines.append("")
                    if alternatives:
                        lines.append("**Suggested alternatives:**")
                        for alt in alternatives[:4]:
                            lines.append(
                                "- %s" % str(alt)[:120])
                        lines.append("")
                    if review_guidance:
                        lines.append(
                            "**Guidance:** %s"
                            % review_guidance[:300])
                        lines.append("")

        # ── 7. Final Quality ────────────────────────────────────
        if data.get("final_metrics"):
            fm = data["final_metrics"]

            try:
                from libtbx.langchain.knowledge.summary_display \
                    import format_quality_table_rows
            except Exception:
                from knowledge.summary_display \
                    import format_quality_table_rows

            rows = format_quality_table_rows(
                fm, data.get('experiment_type'))
            if rows:
                lines.append("## Final Quality")
                lines.append("")
                lines.append("| Metric | Value | Assessment |")
                lines.append("|--------|-------|------------|")
                for row in rows:
                    lines.append(
                        "| %s | %s | %s |"
                        % (row['label'], row['value'],
                           row['detail']))
                lines.append("")

        # ── 8. Output Files ─────────────────────────────────────
        if data.get("final_files"):
            lines.append("## Key Output Files")
            lines.append("")
            lines.append("| File | Type | Description |")
            lines.append("|------|------|-------------|")
            for f in data["final_files"]:
                file_type = f.get('type', 'file')
                description = f.get('description', '')
                display_path = self._get_display_path(f)
                lines.append(
                    "| %s | %s | %s |"
                    % (display_path, file_type, description))
            lines.append("")
            lines.append("")

        # ── 9. Failure Diagnosis ───────────────────────────────
        failure_diagnosis_path = data.get("failure_diagnosis_path")
        if failure_diagnosis_path:
            lines.append("## Failure Diagnosis")
            lines.append("")
            lines.append(
                "A terminal error was encountered during this "
                "run. An AI-generated diagnosis report was saved:")
            lines.append("`%s`" % failure_diagnosis_path)
            lines.append("")

        # ── 10. LLM Assessment ──────────────────────────────────
        if include_llm_assessment:
            lines.append("## Assessment")
            lines.append("")
            lines.append("_[LLM assessment will be inserted here]_")
            lines.append("")

        return "\n".join(lines)

    def get_summary_for_llm_assessment(self):
        """
        Get a concise summary suitable for LLM assessment.

        Returns a text block with key information for the LLM to evaluate:
        - Input data quality
        - Goal/strategy
        - Pre-flight check results
        - Steps taken with reasoning and expert assessments
        - Stop decision with expert review
        - Current status
        """
        data = self._extract_summary_data()

        lines = []

        # Header with run context
        run_name = data.get('run_name', 'unknown')
        is_tutorial = data.get('is_tutorial', False)
        run_type = ("TUTORIAL/FOCUSED TASK"
                    if is_tutorial else "STRUCTURE DETERMINATION")
        lines.append(
            "=== AGENT SESSION: %s (%s) ==="
            % (run_name, run_type))
        lines.append("")

        # Input
        lines.append(
            "EXPERIMENT TYPE: %s" % data['experiment_type'])
        lines.append(
            "INPUT FILES: %s" % ', '.join(data['original_files']))
        lines.append("USER GOAL: %s" % data['user_advice'])
        lines.append("")

        # Stop condition / tutorial detection
        directives = self.get_directives()
        if directives:
            stop_cond = directives.get("stop_conditions", {})
            if stop_cond:
                lines.append("STOP CONDITION (FOCUSED TASK):")
                if stop_cond.get("after_program"):
                    lines.append(
                        "  Stop after: %s"
                        % stop_cond['after_program'])
                if stop_cond.get("after_cycle"):
                    lines.append(
                        "  Stop after cycle: %s"
                        % stop_cond['after_cycle'])
                if stop_cond.get("skip_validation"):
                    lines.append(
                        "  Skip validation: Yes (user-directed stop)")
                lines.append(
                    "  NOTE: This was a FOCUSED TASK/TUTORIAL, "
                    "not full structure determination.")
                lines.append("")

        # Pre-flight check (v113)
        preflight = data.get("preflight_check", {})
        if preflight and isinstance(preflight, dict):
            pf_ready = preflight.get("ready", True)
            pf_analysis = preflight.get("analysis", "")
            pf_missing = preflight.get("missing_files", [])
            if not pf_ready or pf_missing:
                lines.append("PRE-FLIGHT CHECK: NOT READY")
                if pf_analysis:
                    lines.append("  %s" % pf_analysis.strip()[:300])
                if pf_missing:
                    lines.append(
                        "  Missing: %s" % ", ".join(
                            str(m) for m in pf_missing[:5]))
                lines.append("")

        # Input quality
        if data.get("input_quality"):
            lines.append("INPUT DATA QUALITY:")
            for k, v in data["input_quality"].items():
                lines.append("  %s: %s" % (k, v))
            lines.append("")

        # Workflow
        lines.append(
            "WORKFLOW PATH: %s" % data['workflow_path'])
        lines.append("")

        # Steps — annotate superseded failures
        raw_cycles = self.data.get("cycles", [])
        _fc_by_cycle = {}
        for _rc in raw_cycles:
            _cn = _rc.get("cycle_number")
            _fc = _rc.get("failure_context")
            if _cn is not None and _fc:
                _fc_by_cycle[_cn] = _fc

        _any_superseded = any(
            fc.get("superseded_by")
            for fc in _fc_by_cycle.values()
        )
        if _any_superseded:
            lines.append(
                "NOTE: Steps marked [SUPERSEDED] failed but "
                "were resolved by later successful steps.")
            lines.append("")

        lines.append("STEPS TAKEN:")
        for step in data["steps"]:
            status = "OK" if step["success"] else "FAILED"
            fc = _fc_by_cycle.get(step["cycle"])
            if fc and fc.get("superseded_by"):
                tag = " [SUPERSEDED by %s]" % fc["superseded_by"]
            else:
                tag = ""
            lines.append(
                "  %s. %s - %s%s %s"
                % (step['cycle'], step['program'],
                   status, tag, step['key_metric']))

            # Include expert assessment for LLM context (v113)
            ea = step.get("expert_assessment", {})
            if ea and isinstance(ea, dict):
                ea_analysis = ea.get("analysis", "").strip()
                if ea_analysis:
                    if len(ea_analysis) > 200:
                        ea_analysis = ea_analysis[:197] + "..."
                    lines.append(
                        "     [Expert: %s] %s"
                        % (ea.get("action", ""), ea_analysis))
        lines.append("")

        # Stop decision (v113)
        stop_info = data.get("stop_info")
        if stop_info and isinstance(stop_info, dict):
            stop_reasoning = (
                stop_info.get("reasoning")
                or stop_info.get("explanation", "")).strip()
            if stop_reasoning:
                if len(stop_reasoning) > 400:
                    stop_reasoning = stop_reasoning[:397] + "..."
                lines.append(
                    "STOP DECISION: %s" % stop_reasoning)

            stop_expert = stop_info.get("expert_assessment", {})
            if stop_expert and isinstance(stop_expert, dict):
                stop_review = stop_expert.get("stop_review", {})
                if stop_review:
                    agrees = stop_review.get(
                        "agree_with_stop", True)
                    lines.append(
                        "  Expert %s with stopping."
                        % ("agrees" if agrees else "DISAGREES"))
                    ra = stop_review.get("analysis", "").strip()
                    if ra:
                        lines.append("  %s" % ra[:200])
                    alts = stop_review.get("alternatives", [])
                    if alts:
                        lines.append(
                            "  Alternatives: %s" % "; ".join(
                                str(a)[:60] for a in alts[:3]))
            lines.append("")

        # Final metrics
        if data.get("final_metrics"):
            lines.append("FINAL METRICS:")
            fm = data["final_metrics"]
            if "r_free" in fm:
                lines.append(
                    "  R-free: %.4f (%s)"
                    % (fm['r_free'],
                       fm.get('r_free_assessment', '')))
            if "r_work" in fm:
                lines.append(
                    "  R-work: %.4f" % fm['r_work'])
            if "map_cc" in fm:
                lines.append(
                    "  Map CC: %.3f (%s)"
                    % (fm['map_cc'],
                       fm.get('map_cc_assessment', '')))
            if "clashscore" in fm:
                lines.append(
                    "  Clashscore: %.2f (%s)"
                    % (fm['clashscore'],
                       fm.get('clashscore_assessment', '')))
            if "bonds_rmsd" in fm:
                lines.append(
                    "  Bonds RMSD: %.4f" % fm['bonds_rmsd'])
            if "angles_rmsd" in fm:
                lines.append(
                    "  Angles RMSD: %.3f" % fm['angles_rmsd'])
            if "ramachandran_outliers" in fm:
                lines.append(
                    "  Ramachandran outliers: %.2f%%"
                    % fm['ramachandran_outliers'])
            lines.append("")

        # Key output files
        if data.get("final_files"):
            lines.append("KEY OUTPUT FILES:")
            for f in data["final_files"]:
                desc = f.get('description', f.get('type', 'file'))
                display_path = self._get_display_path(f)
                lines.append("  %s - %s" % (display_path, desc))
            lines.append("")

        lines.append(
            "TOTAL CYCLES: %d (%d successful)"
            % (data['total_cycles'], data['successful_cycles']))

        return "\n".join(lines)
