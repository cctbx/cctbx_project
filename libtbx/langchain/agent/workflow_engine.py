"""
Workflow Engine for PHENIX AI Agent.

This module interprets workflows.yaml to determine:
- Current step in workflow
- Valid programs for current step
- Transitions to next step
- Quality targets

The engine provides a higher-level abstraction than the hardcoded
state detection in workflow_state.py.

Usage:
    from libtbx.langchain.agent.workflow_engine import WorkflowEngine

    engine = WorkflowEngine()
    step_info = engine.get_current_phase("xray", context)
    valid_programs = engine.get_valid_programs("xray", step_info)
"""

from __future__ import absolute_import, division, print_function

import logging
import os
import re

logger = logging.getLogger(__name__)

# Import YAML loader
from libtbx.langchain.knowledge.yaml_loader import (
    get_workflow_steps,
    get_workflow_targets,
    get_metric_threshold,
    get_program,
)

# Import program registry for program info
from libtbx.langchain.agent.program_registry import ProgramRegistry


class WorkflowEngine:
    """
    Interprets workflow definitions from YAML.

    Provides step detection and program validation based on
    workflow configuration rather than hardcoded logic.
    """

    def __init__(self):
        """Initialize the engine."""
        self._registry = ProgramRegistry()

    # =========================================================================
    # CONTEXT BUILDING
    # =========================================================================

    def build_context(self, files, history_info, analysis=None, directives=None,
                      session_info=None, files_local=True):
        """
        Build a context dict from files, history, and analysis.

        This context is used for evaluating conditions in workflows.

        Args:
            files: Categorized files dict
            history_info: Analyzed history dict
            analysis: Current log analysis dict
            directives: Optional user directives dict

        Returns:
            dict: Context for condition evaluation
        """
        # Pre-compute user_wants_ligandfit (needed by
        # _has_placed_model before the full context
        # is built).
        _d = directives or {}
        _d_sc = _d.get("stop_conditions", {})
        _d_wf = _d.get("workflow_preferences", {})
        _d_after = (_d_sc.get("after_program")
                    or "").lower()
        _d_prefer = _d_wf.get(
            "prefer_programs", [])
        _d_explicit = (
            (session_info or {}).get(
                "explicit_program") or "").lower()
        _early_wants_ligandfit = bool(
            "ligandfit" in _d_after or
            any("ligandfit" in p.lower()
                for p in _d_prefer) or
            "ligandfit" in _d_explicit
        )

        context = {
            # File availability - using semantic categories
            "has_data_mtz": bool(files.get("data_mtz")),
            "has_phased_data_mtz": bool(files.get("phased_data_mtz")),
            "has_map_coeffs_mtz": bool(files.get("map_coeffs_mtz")),
            "has_sequence": bool(files.get("sequence")),
            # SEMANTIC: 'model' = positioned models (phaser_output, refined, docked)
            # SEMANTIC: 'search_model' = templates NOT yet positioned (predicted, pdb_template)
            "has_model": bool(files.get("model")),  # Any model file (positioned or generic PDB)
            "has_search_model": bool(files.get("search_model")),  # Templates/predictions
            "has_model_for_mr": bool(files.get("model") or files.get("search_model")),  # Either works for phaser
            "has_full_map": bool(files.get("full_map")),
            "has_half_map": len(files.get("half_map", [])) >= 2,
            "has_map": bool(files.get("map") or files.get("full_map")),
            # Composite: has a map that is NOT a half map (full, sharpened, optimized, etc.)
            # Used by map_symmetry which needs a full map but can accept processed maps
            # that may not be in the strict full_map subcategory
            "has_non_half_map": bool(
                set(files.get("map", [])) - set(files.get("half_map", []))
            ),
            "has_ligand": bool(files.get("ligand") or files.get("ligand_cif") or files.get("ligand_pdb")),
            "has_ligand_file": bool(files.get("ligand") or files.get("ligand_cif") or files.get("ligand_pdb")),  # Alias for workflow conditions
            "has_predicted_model": bool(files.get("predicted")) or history_info.get("predict_done", False),
            "has_processed_model": bool(files.get("processed_predicted")) or history_info.get("process_predicted_done", False),
            "has_placed_model": self._has_placed_model(
                files, history_info, directives,
                user_wants_ligandfit=_early_wants_ligandfit),
            # Placement confirmed by history/files only — does NOT reflect the
            # model_is_placed directive.  Used by Tier 1 routing and S1 short-circuit
            # so that a mis-set directive cannot suppress a real cell-mismatch signal.
            "has_placed_model_from_history": self._has_placed_model_from_history(files, history_info),
            # Explicit after_program directive naming a program that requires a placed model
            # (e.g. after_program=phenix.model_vs_data).  Unlike model_is_placed=True in
            # workflow_preferences, this is a deliberate user request — reliable evidence
            # that the user knows the model is placed.  Used to suppress placement_uncertain.
            "has_placed_model_from_after_program": self._has_placed_model_from_after_program(files, directives),
            "has_refined_model": self._has_refined_model(files, history_info),
            "has_ligand_fit": (
              bool(files.get("ligand_fit_output"))
              or history_info.get("ligandfit_done", False)
              or bool((session_info or {}).get(
                "input_has_ligand", False))
            ),
            "has_optimized_full_map": self._has_optimized_map(files, history_info),

            # Complex program flags (need special logic beyond simple detection)
            "phaser_done": history_info.get("phaser_done", False),
            "predict_done": history_info.get("predict_done", False),
            "predict_full_done": history_info.get("predict_full_done", False),
            "autobuild_done": history_info.get("autobuild_done", False),
            "autosol_done": history_info.get("autosol_done", False),
            "autosol_success": history_info.get("autosol_success", False),
            # autosol_attempted: True whenever autosol has been tried (success or
            # failure).  Guards step 2d from looping when autosol fails post-phaser.
            "autosol_attempted": history_info.get("autosol_attempted", False),
            "refine_done": history_info.get("refine_done", False),
            "rsr_done": history_info.get("rsr_done", False),
            "validation_done": history_info.get("validation_done", False),
            "ligandfit_done": history_info.get("ligandfit_done", False),
            "pdbtools_done": history_info.get("pdbtools_done", False),
            "needs_post_ligandfit_refine": history_info.get("needs_post_ligandfit_refine", False),
            "dock_done": history_info.get("dock_done", False),
            "process_predicted_model_done": history_info.get("process_predicted_done", False),

            # Counts
            "refine_count": history_info.get("refine_count", 0),
            "rsr_count": history_info.get("rsr_count", 0),

            # Last program that ran (for after_program directive check)
            "last_program": history_info.get("last_program"),

            # Metrics
            "r_free": self._get_metric(analysis, history_info, "r_free", "last_r_free"),
            "r_work": self._get_metric(analysis, history_info, "r_work", "last_r_work"),
            "map_cc": self._get_metric(analysis, history_info, "map_cc", "last_map_cc"),
            "clashscore": self._get_metric(analysis, history_info, "clashscore", "last_clashscore"),
            "resolution": self._get_metric(analysis, history_info, "resolution", "resolution"),
            "tfz": self._get_metric(analysis, history_info, "tfz", "last_tfz"),

            # Special conditions - check both current analysis and history
            "has_anomalous": self._get_bool(analysis, history_info, "has_anomalous"),
            "strong_anomalous": self._get_bool(analysis, history_info, "strong_anomalous"),
            "anomalous_measurability": self._get_metric(analysis, history_info, "anomalous_measurability", "anomalous_measurability"),
            "has_twinning": self._get_bool(analysis, history_info, "has_twinning"),
            "twin_law": self._get_metric(analysis, history_info, "twin_law", "twin_law"),
            "twin_fraction": self._get_metric(analysis, history_info, "twin_fraction", "twin_fraction"),
            "anomalous_resolution": self._get_metric(analysis, history_info, "anomalous_resolution", "anomalous_resolution"),
            "has_ncs": self._get_bool(analysis, history_info, "has_ncs"),

            # ── Placement detection (Tier 1 / Tier 3) ──────────────────────
            # Tier 1: unit cell comparison — fast, free, no extra cycle.
            # Computed here so routing logic can use it without a method call.
            # Note: short-circuit against history_info is applied after the
            # dict is fully built (see below) — context self-reference is not
            # possible inside the literal.
            "cell_mismatch": self._check_cell_mismatch(
                files,
                model_cell=(session_info or {}).get("unplaced_model_cell"),
                files_local=files_local),

            # Tier 3: diagnostic probe results (set by history after
            # phenix.model_vs_data / phenix.map_correlations ran as a probe).
            "placement_probed": history_info.get("placement_probed", False),
            "placement_probe_result": history_info.get("placement_probe_result", None),

            # Filled in below once has_placed_model and cell_mismatch are known.
            "placement_uncertain": False,
        }

        # =====================================================================
        # Auto-include done flags from program_registration.
        #
        # Uses get_program_done_flag_map() (all strategies) rather than the
        # old get_all_done_flags() (only strategy: run_once, 3 programs).
        # The old approach silently dropped resolve_cryo_em_done,
        # map_sharpening_done, map_to_model_done, autobuild_denmod_done,
        # and polder_done — those flags were set correctly by _analyze_history()
        # but never transferred into context.
        #
        # The "don't override" guard is intentional: flags that are already
        # set explicitly above (e.g., autobuild_done, phaser_done) keep their
        # values; this loop only fills gaps for flags not yet in context.
        # =====================================================================
        from libtbx.langchain.knowledge.program_registration import get_program_done_flag_map

        for flag_name in set(get_program_done_flag_map().values()):
            if flag_name not in context:
                context[flag_name] = history_info.get(flag_name, False)

        # Derive has_twinning from twin_fraction if not explicitly set
        if not context["has_twinning"] and context.get("twin_fraction"):
            # Twinning threshold is 0.20 (from workflows.yaml)
            if context["twin_fraction"] > 0.20:
                context["has_twinning"] = True

        # Derive strong_anomalous from measurability if not explicitly set
        if not context["strong_anomalous"] and context.get("anomalous_measurability"):
            if context["anomalous_measurability"] > 0.10:
                context["strong_anomalous"] = True
                context["has_anomalous"] = True

        # Compute derived conditions
        context["model_is_good"] = self._is_model_good(context)

        # Model-vs-Data Gate (v114.1): session-level
        # placement confirmation from evidence
        # (model_vs_data CC, refine R-free). This is
        # a hard gate that prevents destructive
        # programs from being offered.
        context["model_is_placed_confirmed"] = bool(
          (session_info or {}).get(
            "model_is_placed", False))

        # Propagate MR-SAD workflow preference from directives
        if directives:
            workflow_prefs = directives.get("workflow_preferences", {})
            context["use_mr_sad"] = workflow_prefs.get("use_mr_sad", False)

            # v115.09 Fix 3: validation-only intent from directives
            context["wants_validation_only"] = workflow_prefs.get(
                "wants_validation_only", False)

            # CRITICAL: Treat skipped programs as "done" so the workflow can
            # advance past steps that require them.  Without this, skipping
            # xtriage causes the workflow to stay stuck in "analyze" with no
            # valid programs, because _detect_xray_step checks xtriage_done
            # before allowing progression to later steps.
            skip_programs = workflow_prefs.get("skip_programs", [])
            if skip_programs:
                # get_program_done_flag_map already imported above
                program_to_done_flag = get_program_done_flag_map()
                for prog in skip_programs:
                    done_flag = program_to_done_flag.get(prog)
                    if done_flag and not context.get(done_flag):
                        context[done_flag] = True
        else:
            context["use_mr_sad"] = False
            context["wants_validation_only"] = False

        # v115.09 Fix 4: MR-SAD with model but no search_model.
        # The model is likely a search model that the categorizer
        # couldn't distinguish from an already-placed model.
        # Set force_mr so phaser is offered without reclassifying
        # the file between categories.
        if (context.get("use_mr_sad") and
            context.get("has_model") and
            not context.get("has_search_model") and
            not context.get("phaser_done")):
            context["force_mr"] = True
            context["has_model_for_mr"] = True

        # ── Tier 1 short-circuit: skip cell check when placement is resolved ────
        # _check_cell_mismatch may run phenix.show_map_info as a subprocess for
        # cryo-EM maps; skip it on subsequent cycles once placement is already known.
        # Conditions that render the check unnecessary:
        #   - history/files confirm placement (has_placed_model_from_history=True)
        #   - probe already ran this session (placement_probed=True)
        # IMPORTANT: we use has_placed_model_FROM_HISTORY, not has_placed_model.
        # A directive setting model_is_placed=True does NOT suppress the check —
        # that directive could be wrong (as seen when LLM misinterprets the goal),
        # and we need the cell-mismatch check to catch that error.
        # Note: cell_mismatch is computed before these flags are set above, so we
        # override here in post-processing once the full context is available.
        if context.get("has_placed_model_from_history") or context.get("placement_probed"):
            context["cell_mismatch"] = False

        # ── Probe result overrides has_placed_model ──────────────────────────
        # When the probe ran and confirmed placement (R-free < 0.50 or CC > 0.15),
        # treat the model as placed so normal refine routing takes over.
        # This overrides the Tier 2 heuristic which could not determine placement.
        if (context["placement_probed"] and
                context.get("placement_probe_result") == "placed"):
            context["has_placed_model"] = True

        # ── Placement uncertainty (Tier 2 gate for Tier 3 probe) ────────────
        # placement_uncertain=True when all conditions hold:
        #   - No OBJECTIVE confirmation of placement from history/files
        #   - No definitive cell mismatch (that already routes to MR/dock)
        #   - Probe has not already been run this session
        #   - Has a model AND reflection data or map (probe inputs available)
        #   - Model is NOT a predicted model (those always need MR/dock — no probe)
        #   - User has NOT explicitly requested a program that presupposes placement
        #
        # CRITICAL: we use has_placed_model_FROM_HISTORY, NOT has_placed_model.
        # The directive model_is_placed=True can be set incorrectly by the LLM
        # (e.g. "solve the structure" → model_is_placed=True is a known failure
        # mode).  Using has_placed_model here would suppress the probe even when
        # the model has never actually been placed, allowing RSR to run against a
        # mismatched map and crash with a symmetry error.
        #
        # EXCEPTION: if the user explicitly requested a program that only makes
        # sense on an already-placed model (e.g. phenix.refine, phenix.ligandfit),
        # trust that request and skip the probe.  The user is asserting that
        # placement is not in question.
        _explicit_prog = (session_info or {}).get("explicit_program") or ""
        _PROBE_SKIP_PROGRAMS = {
            "phenix.refine", "phenix.ligandfit", "phenix.polder",
            "phenix.molprobity", "phenix.autobuild",
        }
        _user_asserts_placed = _explicit_prog in _PROBE_SKIP_PROGRAMS

        # Also check raw user advice for placement clues.
        # If user says "refine" or "fit ligand" without
        # mentioning MR/phaser/solve, the model is placed.
        _raw_advice = (
            (session_info or {}).get(
                "user_advice", "")
        ).lower()
        if _raw_advice and not _user_asserts_placed:
            _mr_kw = (
                "molecular replacement", "phaser",
                "solve", "mr ", "autosol", "predict",
                "sad ", "sad,", "mad ", "mad,",
                "anomalous", "phasing",
                "heavy atom", "selenium",
                # Cryo-EM docking keywords — if the user
                # mentions docking, the model has NOT been
                # placed yet.  "refine the docked model"
                # means dock-THEN-refine, not already-done.
                "dock", "docking", "dock_in_map",
                "place the model", "place model",
                "fit into map", "fit model into",
            )
            _placed_kw = (
                "refine", "refinement",
                "fit ligand", "ligandfit",
                "fit the ligand", "fit atp",
                "fit nad", "fit fad", "fit heme",
                "polder", "validate", "molprobity",
                "just refine", "only refine",
                "further refine", "continue refin",
            )
            _advice_says_mr = any(
                kw in _raw_advice for kw in _mr_kw
            )
            _advice_says_placed = any(
                kw in _raw_advice
                for kw in _placed_kw
            )
            if _advice_says_placed and not _advice_says_mr:
                _user_asserts_placed = True

        # When user advice asserts placement AND a
        # model file exists, override has_placed_model
        if (_user_asserts_placed
            and context.get("has_model")
            and not context.get("has_placed_model")):
            context["has_placed_model"] = True

        # ── User intent: does the user explicitly want ligand fitting? ───────
        # Set when directives or session_info indicate the user wants ligandfit.
        # Used in _detect_xray_step to relax the r_free < 0.35 threshold and
        # stay in the "refine" step so ligandfit remains available, rather than
        # falling through to "validate" and offering only validation programs.
        _after_prog   = (directives or {}).get(
            "stop_conditions", {}).get(
            "after_program") or ""
        _prefer_progs = (directives or {}).get(
            "workflow_preferences", {}).get(
            "prefer_programs") or []
        _explicit_lig = (session_info or {}).get(
            "explicit_program") or ""
        context["user_wants_ligandfit"] = bool(
            "ligandfit" in _after_prog.lower() or
            any("ligandfit" in p.lower()
                for p in _prefer_progs) or
            "ligandfit" in _explicit_lig.lower()
        )

        context["placement_uncertain"] = (
            not context["has_placed_model_from_history"] and
            not context["has_placed_model_from_after_program"] and
            not context["cell_mismatch"] and
            not context["placement_probed"] and
            context["has_model"] and
            (context["has_data_mtz"] or context["has_map"]) and
            not context["has_predicted_model"] and
            not _user_asserts_placed
        )

        return context

    def _get_metric(self, analysis, history_info, analysis_key, history_key):
        """Get metric from analysis or history, ensuring numeric type.

        Values may arrive as strings from log parsing (e.g. "0.204").
        Cast to float so downstream comparisons don't crash.
        """
        val = None
        if analysis and analysis.get(analysis_key):
            val = analysis[analysis_key]
        else:
            val = history_info.get(history_key)
        if val is None:
            return None
        try:
            return float(val)
        except (ValueError, TypeError):
            return val  # Non-numeric (shouldn't happen but don't crash)

    def _get_bool(self, analysis, history_info, key):
        """Get boolean from analysis or history (True if either is True)."""
        if analysis and analysis.get(key):
            return True
        return history_info.get(key, False)

    def _has_search_model(self, files):
        """
        Check if there's a search model (template for MR/docking).

        With semantic categories, this simply checks the 'search_model' category
        which contains: predicted, processed_predicted, pdb_template
        """
        return bool(files.get("search_model"))

    def _check_cell_mismatch(self, files, model_cell=None, files_local=True):
        """
        Tier 1 placement check: compare model and data unit cells.

        Returns True **only** when both cells can be read and are definitively
        incompatible (> 5% difference on any parameter).  Any read failure
        returns False (fail-safe — never wrongly prevents MR).

        For X-ray: compares PDB CRYST1 against MTZ crystal symmetry.
        For cryo-EM: compares PDB CRYST1 against both full-map and
            present-portion cells (compatible if it matches either).

        Args:
            files (dict): Categorised files dict from the agent.
            model_cell (list|tuple|None): Pre-read CRYST1 cell transmitted from
                the client (S2L).  When supplied, skips opening the PDB file —
                which would always fail on the server because the file lives on
                the client filesystem and os.path.exists() returns False there.

        Returns:
            bool: True → definitive mismatch; False → compatible or unknown.
        """
        try:
            from libtbx.langchain.agent.placement_checker import (
                check_xray_cell_mismatch,
                check_cryoem_cell_mismatch,
                cells_are_compatible,
                read_mtz_unit_cell,
                read_map_unit_cells,
            )
        except ImportError:
            try:
                from agent.placement_checker import (
                    check_xray_cell_mismatch,
                    check_cryoem_cell_mismatch,
                    cells_are_compatible,
                    read_mtz_unit_cell,
                    read_map_unit_cells,
                )
            except ImportError:
                return False   # Module not available — fail-safe

        # ── S2L: use pre-read cell from client when available ────────────────
        # The server cannot open files that only exist on the client filesystem.
        # When the client transmitted unplaced_model_cell in session_info, use it
        # directly for comparison rather than trying to open the PDB file.
        if model_cell is not None:
            try:
                _mc = tuple(float(v) for v in model_cell)
            except (TypeError, ValueError):
                _mc = None

            if _mc is not None:
                # X-ray: compare against MTZ
                mtz_files = files.get("data_mtz", [])
                if mtz_files and files_local:
                    mtz_cell = read_mtz_unit_cell(mtz_files[0])
                    if mtz_cell is not None:
                        return not cells_are_compatible(_mc, mtz_cell)
                    return False   # MTZ unreadable → fail-safe
                elif mtz_files and not files_local:
                    return False   # Server mode: can't read MTZ → fail-safe

                # Cryo-EM: compare against map file
                map_files = (files.get("full_map", []) or
                             files.get("optimized_full_map", []) or
                             files.get("map", []))
                if map_files and files_local:
                    full_cell, present_cell = read_map_unit_cells(map_files[0])
                    if full_cell is None and present_cell is None:
                        return False   # Map unreadable → fail-safe
                    if full_cell is not None and cells_are_compatible(_mc, full_cell):
                        return False
                    if present_cell is not None and cells_are_compatible(_mc, present_cell):
                        return False
                    if full_cell is not None or present_cell is not None:
                        return True  # At least one cell readable; model matched neither
                return False

        # ── Fallback: read cell from file (works for LocalAgent / test mode) ─
        # Skip entirely on the server where client files aren't on disk.
        if not files_local:
            return False

        # Get the first model file (generic PDB, not a positioned subcategory)
        model_files = files.get("model", [])
        pdb_path = model_files[0] if model_files else None
        if not pdb_path:
            return False

        # X-ray: compare against MTZ
        mtz_files = files.get("data_mtz", [])
        if mtz_files:
            return check_xray_cell_mismatch(pdb_path, mtz_files[0])

        # Cryo-EM: compare against map (full_map preferred, then optimized, then map)
        map_files = (files.get("full_map", []) or
                     files.get("optimized_full_map", []) or
                     files.get("map", []))
        if map_files:
            return check_cryoem_cell_mismatch(pdb_path, map_files[0])

        return False   # No data to compare against

    def _promote_unclassified_for_docking(self, files, context, experiment_type):
        """
        S2c: Promote unclassified PDB files to search_model when context
        unambiguously confirms that docking is needed.

        Problem being solved
        --------------------
        Filename-based categorisation runs at session start with no workflow
        context.  A plain crystal-structure PDB like ``1aew_A.pdb`` lands in
        ``unclassified_pdb → model`` because nothing in the name says "template".
        When Tier 1 (cell_mismatch) or Tier 3 (probe) routes to ``dock_model``,
        the ``phenix.dock_in_map`` condition ``has_any: [search_model, processed_model]``
        still fails because ``files["search_model"]`` is empty.

        This method detects the unambiguous "needs docking" signal from context
        and copies the unclassified PDB paths into ``files["search_model"]`` so
        that the YAML condition is satisfied and the command builder can find
        the correct input file.

        Trigger conditions (all must hold)
        -----------------------------------
        1. ``experiment_type == "cryoem"``   — X-ray paths use ``model`` key for phaser
        2. ``files["unclassified_pdb"]`` non-empty — something to promote
        3. ``not has_placed_model_from_history``  — no historical dock/refine evidence
        4. ``not has_search_model``          — search_model already populated → no-op
        5. Any one of:
           a. ``placement_uncertain``         — Tier 3 path, probe not yet run
           b. ``placement_probed AND needs_dock`` — Tier 3 path, post-probe
           c. ``cell_mismatch AND not from_history`` — Tier 1 path (production)

        Safety
        ------
        * Input dicts are never mutated.  Shallow copies are made only when
          promotion fires.
        * ``files["model"]`` and ``files["unclassified_pdb"]`` are left intact
          so no information is lost.
        * X-ray sessions are completely unaffected (condition 1).
        * Sessions where dock/refine already ran are blocked (condition 3).

        Scope
        -----
        ``files["unclassified_pdb"]`` is only populated by ``_categorize_files_yaml``
        (the normal PHENIX path).  The hardcoded fallback never sets this key, so
        promotion silently no-ops in that degraded environment — which is fine
        because ``has_model`` is also False there and the scenario doesn't arise.

        Returns
        -------
        (files, context) — promoted copies when conditions met, originals otherwise.
        """
        if experiment_type != "cryoem":
            return files, context

        unclassified = files.get("unclassified_pdb", [])
        if not unclassified:
            return files, context

        if context.get("has_placed_model_from_history"):
            return files, context

        if context.get("has_search_model"):
            return files, context

        needs_dock = (
            context.get("placement_uncertain") or
            (context.get("placement_probed") and
             context.get("placement_probe_result") == "needs_dock") or
            (context.get("cell_mismatch") and
             not context.get("has_placed_model_from_history"))
        )
        if not needs_dock:
            return files, context

        # Shallow-copy both dicts — never mutate the caller's originals
        files = dict(files)
        files["search_model"] = list(files.get("search_model", [])) + unclassified
        if "pdb" in files:
            files["pdb"] = list(files["pdb"])
            for f in unclassified:
                if f not in files["pdb"]:
                    files["pdb"].append(f)

        context = dict(context)
        context["has_search_model"] = True
        context["unclassified_promoted_to_search_model"] = True

        return files, context

    def _has_placed_model(self, files, history_info,
                          directives=None,
                          user_wants_ligandfit=False):
        """
        Check if model is placed (after MR/building) or is ready-to-refine.

        A model is considered "placed" if:
        1. History shows MR, building, or prediction was done
        2. There's a file in a SPECIFIC positioned subcategory (refined,
           phaser_output, autobuild_output, docked, with_ligand, rsr_output)
        3. User explicitly says model is placed (model_is_placed directive)
        4. User's constraints imply model is placed (e.g., "refine", "fit ligand")
        5. User explicitly requests refinement/validation via directives
        6. User wants ligand fitting (implies model is placed)

        IMPORTANT: An unclassified PDB (generic file like "1abc.pdb") in the
        'model' parent category is NOT sufficient to be "placed". It could be
        a search model for MR. Only PDBs in explicit positioned subcategories
        count, unless history or directives confirm placement.
        """
        # ---- 0. Directive: user explicitly says model is placed ----
        if directives and files.get("model"):
            workflow_prefs = directives.get("workflow_preferences", {})
            if workflow_prefs.get("model_is_placed"):
                return True

            # Infer from constraints: if user's goal implies a placed model
            # (e.g., "refine the model", "fit a ligand"), trust that inference.
            # GUARD (Bug E fix): if the constraints ALSO mention MR/phaser as a
            # required step, the refinement references are future goals not
            # current state — e.g. "do MR then refine".  In that case, do NOT
            # infer placement from refinement keywords alone.
            constraints = directives.get("constraints", [])
            _mr_keywords = [
                "molecular replacement", "phaser", " mr ", "autobuild",
                "place the model", "molecular_replacement",
            ]
            _wants_mr_first = any(
                any(kw in (c.lower() if isinstance(c, str) else "")
                    for kw in _mr_keywords)
                for c in constraints
            )
            if not _wants_mr_first:
                placement_keywords = [
                    "refine", "refinement", "ligandfit", "fit ligand",
                    "fit a ligand", "polder", "validate", "molprobity",
                ]
                for c in constraints:
                    c_lower = c.lower() if isinstance(c, str) else ""
                    if any(kw in c_lower for kw in placement_keywords):
                        return True

        # ---- 0b. User wants ligand fitting → model must be placed ----
        if user_wants_ligandfit and files.get("model"):
            return True

        # ---- 1. History-based: definitive evidence of placement ----
        if (history_info.get("phaser_done") or
            history_info.get("autobuild_done") or
            history_info.get("dock_done") or
            history_info.get("predict_full_done") or
            history_info.get("refine_done")):
            return True

        # ---- 2. File-based: check for POSITIONED subcategories ----
        # These subcategories indicate the model has been positioned in the
        # unit cell (not just a template/search model).
        positioned_subcategories = [
            "refined", "phaser_output", "autobuild_output", "docked",
            "with_ligand", "rsr_output", "ligand_fit_output",
        ]
        for subcat in positioned_subcategories:
            if files.get(subcat):
                return True

        # ---- 3. Directive-based: user explicitly says model is ready ----
        if directives:
            # Check if directives indicate user wants prediction
            # If so, don't assume a generic PDB is positioned
            wants_prediction = False
            constraints = directives.get("constraints", [])
            for c in constraints:
                c_lower = c.lower() if isinstance(c, str) else ""
                if "predict" in c_lower or "alphafold" in c_lower:
                    wants_prediction = True
                    break

            program_settings = directives.get("program_settings", {})
            if "phenix.predict_and_build" in program_settings:
                wants_prediction = True

            stop_conditions = directives.get(
                "stop_conditions") or {}
            after_program = stop_conditions.get(
                "after_program") or ""
            if "predict" in after_program.lower():
                wants_prediction = True

            if not wants_prediction:
                workflow_prefs = directives.get("workflow_preferences", {})
                skip_programs = workflow_prefs.get("skip_programs", [])

                # If skip_programs includes phaser/autosol, user says model is ready
                if any(p in skip_programs for p in ["phenix.phaser", "phenix.autosol"]):
                    if files.get("model"):
                        return True

                # If after_program implies a placed model is needed, AND we
                # have a non-search-model, non-ligand PDB, trust the user
                programs_requiring_placed = [
                    "phenix.refine", "phenix.polder", "phenix.model_vs_data",
                    "phenix.ligandfit", "phenix.molprobity",
                    "phenix.real_space_refine",
                ]
                if after_program in programs_requiring_placed:
                    for f in files.get("pdb", []):
                        basename = os.path.basename(f).lower()
                        is_ligand = (
                            basename.startswith('lig') or
                            f in files.get("ligand_pdb", []) or
                            f in files.get("ligand", [])
                        )
                        is_search = (
                            f in files.get("search_model", []) or
                            f in files.get("predicted", []) or
                            f in files.get("processed_predicted", [])
                        )
                        if not is_ligand and not is_search:
                            return True

        return False


    def _has_placed_model_from_history(self, files, history_info):
        """
        Check if model is placed based on HISTORY and FILES only — no directives.

        This is a strict subset of _has_placed_model() that intentionally ignores
        the model_is_placed directive.  It is used by:

          - Tier 1 routing (cell_mismatch check): a directive claiming the model
            is placed cannot override a definitive cell-dimension mismatch.  Only
            concrete evidence (dock_done in history, refine_done, etc.) can.

          - S1 short-circuit: cell_mismatch subprocess is only skipped when
            placement is confirmed by history, not just by a directive.

        Returns True if:
          1. History shows MR, dock, autobuild, predict_full, or refine completed.
          2. A file in a positioned subcategory (refined, docked, phaser_output,
             autobuild_output, with_ligand, rsr_output, ligand_fit_output) exists.
        """
        # History-based: definitive evidence
        if (history_info.get("phaser_done") or
            history_info.get("autobuild_done") or
            history_info.get("dock_done") or
            history_info.get("predict_full_done") or
            history_info.get("refine_done")):
            return True

        # File-based: positioned subcategories
        positioned_subcategories = [
            "refined", "phaser_output", "autobuild_output", "docked",
            "with_ligand", "rsr_output", "ligand_fit_output",
        ]
        for subcat in positioned_subcategories:
            if files.get(subcat):
                return True

        return False

    def _has_placed_model_from_after_program(self, files, directives):
        """
        Check if an explicit after_program directive reliably implies a placed model.

        This is a SEPARATE signal from the model_is_placed workflow_preference,
        which can be set by LLM hallucination (e.g. "solve the structure" →
        model_is_placed=True is a known failure mode and is intentionally ignored
        by placement_uncertain / S1 short-circuit).

        An after_program directive is different: the user explicitly named a
        concrete program to run.  Programs in programs_requiring_placed cannot
        sensibly run without a positioned model, so the request is reliable
        evidence that the user knows the model is placed.

        Returns True when:
          - directives["stop_conditions"]["after_program"] is one of the programs
            that unambiguously require a placed model, AND
          - at least one non-ligand, non-search-model PDB is present.

        Used to suppress the placement_uncertain probe in build_context(), which
        otherwise fires whenever no history evidence of placement exists.
        """
        if not directives:
            return False
        after_program = directives.get("stop_conditions", {}).get("after_program", "")
        if not after_program:
            return False

        programs_requiring_placed = [
            "phenix.refine", "phenix.polder", "phenix.model_vs_data",
            "phenix.ligandfit", "phenix.molprobity",
            "phenix.real_space_refine",
        ]
        if after_program not in programs_requiring_placed:
            return False

        for f in files.get("pdb", []):
            basename = os.path.basename(f).lower()
            is_ligand = (
                basename.startswith("lig") or
                f in files.get("ligand_pdb", []) or
                f in files.get("ligand", [])
            )
            is_search = (
                f in files.get("search_model", []) or
                f in files.get("predicted", []) or
                f in files.get("processed_predicted", [])
            )
            if not is_ligand and not is_search:
                return True

        return False

    def _has_refined_model(self, files, history_info):
        """Check if model has been refined IN THIS SESSION.

        IMPORTANT: Only trust history-based evidence, not file names.
        User-provided input files may start with 'refine_' (e.g., refine_001_model.pdb)
        without any actual refinement having been done in this session. Relying on
        files.get('refined') would incorrectly skip to validation/STOP.
        """
        return bool(
            history_info.get("refine_count", 0) > 0 or
            history_info.get("rsr_count", 0) > 0
        )

    def _has_optimized_map(self, files, history_info):
        """Check if we have an optimized/sharpened full map."""
        # Look for sharpened maps in full_map list
        for f in files.get("full_map", []):
            if "sharpen" in f.lower():
                return True
        return False

    def _is_model_good(self, context):
        """Determine if model quality is good enough."""
        r_free = context.get("r_free")
        map_cc = context.get("map_cc")
        resolution = context.get("resolution")

        if r_free is not None:
            threshold = get_metric_threshold("r_free", "acceptable", resolution)
            if threshold and r_free < threshold:
                return True

        if map_cc is not None:
            threshold = get_metric_threshold("map_cc", "acceptable")
            if threshold and map_cc > threshold:
                return True

        return False

    # =========================================================================
    # STATE NAME MAPPING
    # =========================================================================
    # Map YAML step names to original hardcoded state names for compatibility

    XRAY_STATE_MAP = {
        "analyze": "xray_initial",
        "probe_placement": "xray_analyzed",   # Probe is pre-step-2; maps to same state
        "obtain_model": "xray_analyzed",
        "molecular_replacement": "xray_has_prediction",   # H2: distinct from xray_initial
        "build_from_phases": "xray_has_phases",
        "experimental_phasing": "xray_mr_sad",            # H2: distinct from xray_initial
        "refine": "xray_refined",
        "combine_ligand": "xray_combined",
        # H3: validate intentionally maps to the same string as refine for external API
        # compatibility. Internal code uses step_info["step"] to distinguish them.
        "validate": "xray_refined",
        "complete": "complete",
    }

    CRYOEM_STATE_MAP = {
        "analyze": "cryoem_initial",
        "probe_placement": "cryoem_analyzed",  # Probe is pre-step-2; maps to same state
        "obtain_model": "cryoem_analyzed",
        "dock_model": "cryoem_has_prediction",
        # H1: check_map and optimize_map deal with half-map → full-map conversion;
        # no model exists yet at this stage. The name "cryoem_has_model" is a legacy
        # misnomer kept for external API compatibility. No behavioral code gates on it.
        "check_map": "cryoem_has_model",
        "optimize_map": "cryoem_has_model",
        "ready_to_refine": "cryoem_docked",   # Model docked, ready for first refinement
        "refine": "cryoem_refined",
        # H3: validate intentionally maps to the same string as refine for external API
        # compatibility. Internal code uses step_info["step"] to distinguish them.
        "validate": "cryoem_refined",
        "complete": "complete",
    }

    def _map_step_to_state(self, step_name, experiment_type):
        """Map YAML step name to original state name."""
        if experiment_type == "xray":
            return self.XRAY_STATE_MAP.get(step_name, step_name)
        elif experiment_type == "cryoem":
            return self.CRYOEM_STATE_MAP.get(step_name, step_name)
        return step_name

    # =========================================================================
    # PHASE DETECTION
    # =========================================================================

    def detect_step(self, experiment_type, context):
        """
        Detect current workflow step based on context.

        Args:
            experiment_type: "xray" or "cryoem"
            context: Context dict from build_context()

        Returns:
            dict: {
                phase: str,           # Phase name
                description: str,     # Human-readable description
                goal: str,            # What we're trying to achieve
                reason: str,          # Why we're in this step
            }
        """
        # === Developer diagnostics (env-var gated) ===
        # These print() calls write to stdout/log file ONLY.
        # They do NOT enter the event system (EventType/state["events"])
        # and cannot reach the LLM's context window.
        # The LLM receives context exclusively through _emit() calls
        # in graph_nodes.py.
        _diag = os.environ.get("PHENIX_AGENT_DIAG_VALID_PROGRAMS")
        if _diag:
            print("  [GATE] detect_step input (%s):" % experiment_type)
            print("  [GATE]   placement_probed=%s result=%s"
                  % (context.get("placement_probed"),
                     context.get("placement_probe_result")))
            print("  [GATE]   validation_done=%s has_refined=%s"
                  % (context.get("validation_done"),
                     context.get("has_refined_model")))
            print("  [GATE]   has_placed_model=%s r_free=%s"
                  % (context.get("has_placed_model"),
                     context.get("r_free")))
            print("  [GATE]   xtriage_done=%s phaser_done=%s"
                  % (context.get("xtriage_done"),
                     context.get("phaser_done")))
            print("  [GATE]   has_model_for_mr=%s autobuild_done=%s"
                  % (context.get("has_model_for_mr"),
                     context.get("autobuild_done")))

        steps = get_workflow_steps(experiment_type)
        if not steps:
            return {"step": "unknown", "reason": "No workflow defined"}

        if experiment_type == "xray":
            result = self._detect_xray_step(steps, context)
        elif experiment_type == "cryoem":
            result = self._detect_cryoem_step(steps, context)
        else:
            result = {"step": "unknown", "reason": "Unknown experiment type"}

        if _diag:
            print("  [GATE] detect_step → '%s': %s"
                  % (result.get("step", "?"),
                     result.get("reason", "?")))

        return result

    def _detect_xray_step(self, steps, context):
        """Detect step in X-ray workflow."""

        # Step 1: Need analysis
        if not context.get("xtriage_done"):
            return self._make_step_result(steps, "analyze",
                "Need to analyze data quality first")

        # v115.09 Fix 3: Validation-only shortcut.
        # When the user's primary goal is validation (not refinement/MR),
        # skip directly to the validate step.  Guards: must have both
        # a model AND data (either data_mtz or phased_data_mtz — completed
        # structures often have phase columns in their MTZ).
        # model_vs_data runs first in the validate step as a crystal-
        # symmetry sanity check before molprobity.
        if (context.get("wants_validation_only") and
            context.get("has_model") and
            (context.get("has_data_mtz") or
             context.get("has_phased_data_mtz"))):
            return self._make_step_result(steps, "validate",
                "Validation-only: running model_vs_data + molprobity")

        # Step 2b: Have prediction, need to place it
        if context.get("has_predicted_model") and not context.get("has_placed_model"):
            if not context.get("has_processed_model"):
                return self._make_step_result(steps, "molecular_replacement",
                    "Have prediction, need to process for MR")
            else:
                return self._make_step_result(steps, "molecular_replacement",
                    "Model processed, need phaser")

        # Step 2c: After autosol, need autobuild
        if context.get("autosol_done") and not context.get("autobuild_done") and not context.get("has_refined_model"):
            return self._make_step_result(steps, "build_from_phases",
                "Experimental phasing complete, need autobuild")

        # Phase 2d: MR-SAD - phaser placed model, anomalous data needs autosol
        # When phaser has run AND we have anomalous signal, autosol should run
        # using the phaser output as a partial model (partpdb_file) for MR-SAD.
        # Also triggered if user explicitly requested MR-SAD workflow.
        # Guard: if autosol was already attempted (even if it failed, e.g.
        # "Unable to find anomalous amplitude arrays"), do NOT loop back here.
        if (context.get("phaser_done") and
            (context.get("has_anomalous") or context.get("use_mr_sad")) and
            not context.get("autosol_done") and
            not context.get("autosol_attempted")):
            return self._make_step_result(steps, "experimental_phasing",
                "Model placed by phaser, anomalous data detected - MR-SAD with autosol")

        # ── Tier 1: unit cell mismatch → skip probe, go straight to MR ─────
        # Both cells were readable and definitively incompatible (> 5 % on any
        # parameter).  No point probing — the model cannot be placed here.
        if context.get("cell_mismatch") and not context.get("has_placed_model_from_history"):
            # Tier 1: cell mismatch overrides even model_is_placed directive.
            # Only history evidence (phaser_done, dock_done, etc.) can suppress this.
            return self._make_step_result(steps, "molecular_replacement",
                "Unit cell mismatch: model cannot be placed in this crystal — running MR")

        # ── Tier 3: probe result known → route based on R-free ───────────────
        if context.get("placement_probed"):
            if (context.get("placement_probe_result") == "needs_mr" and
                    not context.get("has_placed_model_from_history")):
                # Only route to MR when probe says model isn't placed AND no
                # subsequent event (phaser_done, refine_done, …) has confirmed
                # placement.  Without this guard a successful phaser run loops
                # back here because placement_probed/needs_mr are still in history
                # → molecular_replacement phase → phaser excluded by not_done:phaser
                # → valid_programs=[] → STUCK/STOP.
                return self._make_step_result(steps, "molecular_replacement",
                    "Placement probe (model_vs_data) confirmed model is not placed "
                    "(R-free ≥ 0.50) — running MR")
            # probe_result == "placed" OR None (probe ran but was inconclusive,
            # e.g. crystal_symmetry_mismatch from model_vs_data).
            # In either case we treat the model as tentatively placed and proceed
            # to refine — phenix.refine is more permissive than the probe program
            # and will resolve any minor discrepancy or give a clearer error.
            if not context.get("has_placed_model"):
                context["has_placed_model"] = True

        # v115.09 Fix 4: force_mr overrides placement checks.
        # When directives indicate MR-SAD and phaser hasn't run,
        # route directly to MR instead of probing placement.
        if context.get("force_mr") and not context.get("phaser_done"):
            return self._make_step_result(steps, "molecular_replacement",
                "MR-SAD: model needs placement by phaser first")

        # ── Tier 3: probe not yet run, placement ambiguous → run it ──────────
        if context.get("placement_uncertain"):
            return self._make_step_result(steps, "probe_placement",
                "Model placement uncertain — running model_vs_data to check "
                "(R-free < 0.50 → placed, ≥ 0.50 → MR)")

        # Step 2: Need model
        if not context.get("has_placed_model"):
            return self._make_step_result(steps, "obtain_model",
                "Data analyzed, need to obtain model")

        # Step 2e (Bug F fix): R-free very high after MR — route to
        # obtain_model so the agent can retry phaser.  v115.05: also
        # require autobuild_done (incomplete model, not wrong MR).
        _r_free_now = context.get("r_free")
        if (_r_free_now is not None and
                _r_free_now >= 0.45 and
                context.get("refine_count", 0) >= 1 and
                context.get("phaser_done") and
                context.get("has_model_for_mr") and
                context.get("autobuild_done")):
            return self._make_step_result(steps, "obtain_model",
                "Refinement R-free=%.3f after MR — solution likely wrong; "
                "retry with different search model" % _r_free_now)

        # Step 3b: Ligand fitted, need to combine
        if context.get("has_ligand_fit") and not context.get("pdbtools_done"):
            return self._make_step_result(steps, "combine_ligand",
                "Ligand fitted, need to combine")

        # Step 3c: Ligand combined, need refinement of model+ligand complex.
        # After ligandfit adds a ligand, the complex always needs re-refinement.
        # This takes priority over "has_refined_model" because the model has
        # changed since the last refinement.
        if context.get("needs_post_ligandfit_refine"):
            return self._make_step_result(steps, "refine",
                "Ligand fitted, need refinement of model+ligand complex")

        # Step 3: Has model, may need refinement
        if not context.get("has_refined_model"):
            return self._make_step_result(steps, "refine",
                "Have model, need initial refinement")

        # CRITICAL: Stay in "refine" step if ligand fitting is possible.
        # This allows ligandfit to be offered as a valid program.
        # Conditions: ligandfit not already done AND refinement done.
        #
        # We do NOT require has_ligand_file — phenix.ligandfit works with a
        # ligand PDB file OR a residue code (ligand_type=ATP in guidelines).
        #
        # r_free threshold: normally < 0.35 (map quality check).  When the
        # user explicitly requested ligand fitting (user_wants_ligandfit=True),
        # relax to < 0.50 — the user knows their model quality and is asserting
        # they want to proceed regardless.  This prevents the workflow from
        # falling to the "validate" step and offering only validation programs.
        #
        # Special case: when the user provides an ALREADY-REFINED model (e.g.,
        # 7qz0.pdb with R-free=0.20), refine_count is 0 because no refinement
        # happened in this session.  If the model is already at target quality,
        # allow ligandfit without requiring a refinement cycle first.
        _wants_ligandfit = context.get("user_wants_ligandfit", False)
        if (_wants_ligandfit and
            not context.get("ligandfit_done") and
            not context.get("has_ligand_fit")):
            r_free = context.get("r_free")
            _rfree_threshold = 0.50  # Relaxed because user explicitly requested it
            refine_count = context.get("refine_count", 0)
            model_at_target = self._is_at_target(context, "xray")
            if refine_count > 0 or model_at_target:
                if r_free is None or r_free < _rfree_threshold:
                    return self._make_step_result(steps, "refine",
                        "User requested ligand fitting — staying in refine step")

        # Also stay in refine for automatic ligandfit when a ligand file is present
        if (context.get("has_ligand_file") and
            not _wants_ligandfit and
            not context.get("ligandfit_done") and
            not context.get("has_ligand_fit") and
            context.get("refine_count", 0) > 0):
            r_free = context.get("r_free")
            if r_free is None or r_free < 0.35:
                return self._make_step_result(steps, "refine",
                    "Model refined, ligand fitting available")

        # Check if validation is needed
        validation_needed = self._needs_validation(context, "xray")
        if validation_needed and not context.get("validation_done"):
            return self._make_step_result(steps, "validate",
                "Model refined, need validation before stopping")

        # Phase 3 continued: Refinement in progress
        if not self._is_at_target(context, "xray"):
            return self._make_step_result(steps, "refine",
                "Continuing refinement to reach target")

        # Phase 4: At target, validate or complete
        if not context.get("validation_done"):
            return self._make_step_result(steps, "validate",
                "Target reached, need validation")

        # Phase 5: Complete
        return self._make_step_result(steps, "complete",
            "Workflow complete")

    def _detect_cryoem_step(self, steps, context):
        """Detect step in cryo-EM workflow."""

        # Step 1: Need analysis
        # mtriage is preferred, but if other cryo-EM programs have already run
        # (e.g., resolve_cryo_em, dock_in_map, map_to_model), we're clearly past
        # analysis.  This handles tutorials that skip mtriage and also prevents
        # the workflow from getting stuck in "analyze" if mtriage's done flag
        # wasn't set (e.g., history detection missed it after a restart).
        # v115.09: added map_sharpening_done, map_symmetry_done (bgal_denmod fix),
        # and has_optimized_full_map as file-presence fallback for session resumption.
        past_analysis = (
            context.get("mtriage_done") or
            context.get("resolve_cryo_em_done", False) or
            context.get("map_sharpening_done", False) or
            context.get("map_symmetry_done", False) or
            context.get("dock_done", False) or
            context.get("rsr_done", False) or
            context.get("map_to_model_done", False) or
            context.get("autobuild_done", False) or
            context.get("refine_done", False) or
            context.get("predict_done", False) or
            # File-presence fallback: if an optimized map exists,
            # some analysis/processing program must have run
            context.get("has_optimized_full_map", False)
        )
        if not past_analysis:
            return self._make_step_result(steps, "analyze",
                "Need to analyze map quality first")

        # Phase 1.5: Density modification with resolve_cryo_em.
        # Fires when half-maps are available and resolve_cryo_em has not run.
        # We deliberately do NOT gate on has_full_map: a full map may be
        # present in the original inputs (e.g. full_map_box.ccp4) but that
        # does not mean density modification has been performed.  resolve_cryo_em
        # and map_sharpening are complementary; completing one should not
        # suppress the other.  Once resolve_cryo_em_done=True the gate closes.
        # map_sharpening_done is intentionally NOT checked here for the same
        # reason — B-factor sharpening and density modification are different
        # operations (P5 fix: log-confirmed root cause 2026-03-09).
        if (context.get("has_half_map") and
                not context.get("resolve_cryo_em_done")):
            return self._make_step_result(steps, "optimize_map",
                "Have half-maps; running resolve_cryo_em density modification "                "before model building")

        # Step 2b: Stepwise - have prediction, need to dock it
        if context.get("has_predicted_model") and not context.get("has_placed_model"):
            if context.get("has_processed_model"):
                return self._make_step_result(steps, "dock_model",
                    "Have processed prediction, need to dock in map")
            elif context.get("predict_done") and not context.get("predict_full_done"):
                return self._make_step_result(steps, "dock_model",
                    "Have prediction, need to dock in map")

        # ── Tier 1: unit cell mismatch → skip probe, go straight to docking ──
        if context.get("cell_mismatch") and not context.get("has_placed_model_from_history"):
            return self._make_step_result(steps, "dock_model",
                "Unit cell mismatch: model is not placed in this map — running docking")

        # ── Tier 3: probe result known → route based on map CC ───────────────
        if context.get("placement_probed"):
            if context.get("placement_probe_result") == "needs_dock":
                return self._make_step_result(steps, "dock_model",
                    "Placement probe (map_correlations) confirmed model is not placed "
                    "(CC ≤ 0.15) — running docking")

        # ── Tier 3: probe not yet run, placement ambiguous → run it ──────────
        if context.get("placement_uncertain"):
            return self._make_step_result(steps, "probe_placement",
                "Model placement uncertain — running map_correlations to check "
                "(CC > 0.15 → placed, ≤ 0.15 → docking)")

        # Step 2: Need model
        if not context.get("has_placed_model"):
            return self._make_step_result(steps, "obtain_model",
                "Map analyzed, need to obtain model")

        # Step 2c: Check if map needs optimization
        if context.get("has_half_map") and not context.get("has_full_map"):
            return self._make_step_result(steps, "optimize_map",
                "Have model but only half-maps, need full map for refinement")

        # Phase 3a: Have docked model but not yet refined
        if context.get("has_placed_model") and not context.get("has_refined_model"):
            return self._make_step_result(steps, "ready_to_refine",
                "Model docked in map, ready for refinement")

        # Step 3b: Refinement in progress
        if not context.get("has_refined_model"):
            return self._make_step_result(steps, "refine",
                "Have model and map, need refinement")

        # Check validation
        validation_needed = self._needs_validation(context, "cryoem")
        if validation_needed and not context.get("validation_done"):
            return self._make_step_result(steps, "validate",
                "Model refined, need validation")

        # Continue refinement if not at target
        if not self._is_at_target(context, "cryoem"):
            return self._make_step_result(steps, "refine",
                "Continuing refinement")

        # Validate or complete
        if not context.get("validation_done"):
            return self._make_step_result(steps, "validate",
                "Target reached, need validation")

        return self._make_step_result(steps, "complete",
            "Workflow complete")

    def _make_step_result(self, steps, step_name, reason):
        """Build step result dict."""
        step_def = steps.get(step_name, {})
        return {
            "step": step_name,
            "description": step_def.get("description", ""),
            "goal": step_def.get("goal", ""),
            "reason": reason,
        }

    def _needs_validation(self, context, experiment_type):
        """Check if validation is needed before stopping."""
        if experiment_type == "xray":
            r_free = context.get("r_free")
            resolution = context.get("resolution")
            if r_free is not None:
                target = get_metric_threshold("r_free", "acceptable", resolution)
                if target and r_free < target + 0.02:
                    return True
            if context.get("refine_count", 0) >= 3:
                return True
        else:
            map_cc = context.get("map_cc")
            if map_cc is not None and map_cc > 0.70:
                return True
            if context.get("rsr_count", 0) >= 3:
                return True
        return False

    def _is_at_target(self, context, experiment_type):
        """Check if we've reached quality targets.

        Returns False when ligandfit is wanted but refinement hasn't run yet,
        because refinement is a *prerequisite* for ligandfit (it produces the
        map coefficients that ligandfit needs).  Without this, the workflow
        stops prematurely for pre-refined models whose R-free is already good.
        """
        # Ligandfit prerequisite: refinement must run first to produce map
        # coefficients, regardless of current R-free quality.
        if (context.get("user_wants_ligandfit") and
                not context.get("ligandfit_done") and
                not context.get("has_ligand_fit")):
            refine_key = "rsr_count" if experiment_type == "cryoem" else "refine_count"
            if context.get(refine_key, 0) == 0:
                return False

        if experiment_type == "xray":
            r_free = context.get("r_free")
            resolution = context.get("resolution")
            refine_count = context.get("refine_count", 0)

            if r_free is not None:
                target = get_metric_threshold("r_free", "acceptable", resolution)
                if target and r_free <= target:
                    return True

                # Bail out on hopeless R-free: if R-free > 0.5 after at least
                # one refinement, something is fundamentally wrong (wrong model,
                # wrong space group, severe data issues). Further refinement
                # won't help — stop and let the user investigate.
                #
                # EXCEPTION: if autobuild hasn't been tried yet, a high R-free
                # after MR may simply mean the model is incomplete (e.g.,
                # AlphaFold model covering only part of the ASU).  Autobuild
                # can rebuild and complete the model, often dramatically
                # improving R-free (e.g., AF_POMGNT2: 0.55 → 0.28).
                # Don't declare hopeless until rebuilding has been attempted.
                if r_free > 0.50 and refine_count >= 1:
                    if context.get("autobuild_done"):
                        return True
                    # else: let autobuild have a chance

            # Also check clashscore - if it's good AND R-free is reasonable,
            # we're likely done.  A low clashscore with hopeless R-free
            # (>0.45) means the model has good geometry but is incorrectly
            # placed or incomplete (common with AlphaFold predictions after
            # MR).  Don't declare at-target in that case.
            # Also require at least one refinement cycle — a pre-refined
            # input model may have good geometry but the agent hasn't
            # produced map coefficients yet (needed by polder, ligandfit).
            clashscore = context.get("clashscore")
            if clashscore is not None and clashscore < 10:
                r_free = context.get("r_free")
                if (refine_count >= 1 and
                        (r_free is None or r_free < 0.45)):
                    return True

            # Hard limit: if we've done 3+ refine cycles, consider it at target.
            # This prevents unbounded refinement when R-free plateaus above
            # the target. The LLM may still request more if it sees improvement.
            if refine_count >= 3:
                return True
        else:
            # Cryo-EM: check map_cc OR clashscore OR max cycles
            map_cc = context.get("map_cc")
            if map_cc is not None:
                target = get_metric_threshold("map_cc", "acceptable")
                if target and map_cc >= target:
                    return True

            # Also check clashscore - if it's good, we're done
            clashscore = context.get("clashscore")
            if clashscore is not None and clashscore < 12:
                return True

            # Hard limit: if we've done 3+ RSR cycles, consider it at target
            rsr_count = context.get("rsr_count", 0)
            if rsr_count >= 3:
                return True

        return False

    # =========================================================================
    # PROGRAM SELECTION
    # =========================================================================

    def get_valid_programs(self, experiment_type, step_info, context, directives=None):
        """
        Get valid programs for current step.

        Args:
            experiment_type: "xray" or "cryoem"
            step_info: Output from detect_step()
            context: Context dict
            directives: Optional user directives dict

        Returns:
            list: Valid program names
        """
        steps = get_workflow_steps(experiment_type)
        step_name = step_info.get("step", "")
        step_def = steps.get(step_name, {})

        # === DIAGNOSTIC: trace every filter stage ===
        _diag = os.environ.get("PHENIX_AGENT_DIAG_VALID_PROGRAMS")

        # Handle completion step
        if step_def.get("stop"):
            return ["STOP"]

        # Get programs from step definition
        step_programs = step_def.get("programs", [])

        valid = []
        for prog_entry in step_programs:
            if isinstance(prog_entry, str):
                # Simple program name
                valid.append(prog_entry)
            elif isinstance(prog_entry, dict):
                # Program with conditions
                prog_name = prog_entry.get("program")
                if prog_name and self._check_conditions(prog_entry, context):
                    valid.append(prog_name)
                elif prog_name and _diag:
                    print("  [DIAG] %s excluded by conditions: %s"
                          % (prog_name,
                             prog_entry.get("conditions", [])))

        if _diag:
            print("  [DIAG] step=%s, after YAML conditions: %s"
                  % (step_name, valid))

        # Filter out programs that require full_map when only half_maps available
        has_full_map = context.get("has_full_map", False)
        has_half_map = context.get("has_half_map", False)
        only_half_maps = has_half_map and not has_full_map

        if only_half_maps:
            # Remove programs that require full maps
            filtered = []
            for prog in valid:
                prog_def = get_program(prog)
                if prog_def and prog_def.get("requires_full_map"):
                    # Skip this program - needs full map
                    continue
                filtered.append(prog)
            valid = filtered

        # Filter out run_once programs that have already been run
        filtered = []
        for prog in valid:
            prog_def = get_program(prog)
            if prog_def:
                tracking = prog_def.get("done_tracking", {})
                strategy = tracking.get("strategy", "set_flag")
                if strategy == "run_once" or tracking.get("run_once"):
                    # Check if this program has already been run
                    prog_done_key = tracking.get("flag", "")
                    if prog_done_key and context.get(prog_done_key):
                        # Skip - already run
                        continue
                # Belt-and-suspenders: Also filter any non-count program whose
                # UNIQUE done flag is set.  The "count" strategy (refine, rsr,
                # phaser) intentionally repeats; all others should stop once
                # their flag is True.  This catches programs whose strategy
                # wasn't explicitly set to run_once but should have been.
                # Shared flags like "validation_done" (used by molprobity,
                # model_vs_data, etc.) are excluded — we only filter when the
                # flag name is program-specific (contains the program short
                # name).
                if strategy != "count":
                    flag = tracking.get("flag", "")
                    short = prog.replace("phenix.", "").replace(".", "_")
                    if flag and short in flag and context.get(flag):
                        continue
            filtered.append(prog)
        valid = filtered

        if _diag:
            print("  [DIAG] after run_once/done filter: %s" % valid)

        # MR-SAD guard: When we have BOTH a search model AND anomalous data,
        # phaser must run first to place the model. AutoSol will be available
        # later in the experimental_phasing phase (after phaser completes).
        # This prevents autosol from running standalone when MR-SAD is appropriate.
        # v115.09: force_mr allows this guard to fire even when has_search_model
        # is False (model categorized as 'model' but user wants MR-SAD).
        if (context and
            (context.get("has_search_model") or
             context.get("force_mr")) and
            (context.get("has_anomalous") or context.get("use_mr_sad")) and
            not context.get("phaser_done") and
            "phenix.autosol" in valid):
            valid.remove("phenix.autosol")

        # Negligible-anomalous guard (v115.05): when anomalous
        # measurability is below 0.05 and has_anomalous is
        # NOT explicitly True, remove autosol entirely.
        # The data may contain I(+)/I(-) columns but the
        # signal is too weak for SAD phasing.  When
        # has_anomalous is explicitly True (confirmed
        # anomalous data), allow the user to try autosol
        # even with weak measurability.  When has_anomalous
        # is absent/None/False, the measurability threshold
        # applies.
        # (AF_exoV_PredictAndBuild: has_anomalous absent,
        # measurability 0.032 → autosol removed.)
        if (context and
            "phenix.autosol" in valid and
            not context.get("autosol_done")):
            _am = context.get("anomalous_measurability")
            _ha = context.get("has_anomalous")
            if _am is not None and _ha is not True:
                try:
                    _am_f = float(_am)
                except (ValueError, TypeError):
                    _am_f = None
                if _am_f is not None and _am_f < 0.05:
                    valid.remove("phenix.autosol")

        # Restore: if the user explicitly wants ligandfit, has refined at least once,
        # and ligandfit hasn't run yet, add it to valid_programs even if the r_free < 0.35
        # YAML condition was not met.  phenix.ligandfit works with a ligand PDB file OR
        # a residue code (ligand_type=ATP); we do NOT require a separately-categorized
        # ligand file — that would wrongly block runs where the user provides atp.pdb,
        # gdp.pdb, or specifies the ligand by name in their guidelines.
        if (context.get("user_wants_ligandfit") and
                context.get("refine_count", 0) > 0 and
                not context.get("ligandfit_done") and
                not context.get("has_ligand_fit") and
                "phenix.ligandfit" not in valid):
            valid.append("phenix.ligandfit")

        # Fix 3 (v116): Guard against spurious ligandfit in the absence of a
        # ligand file and without explicit user request.  After a successful
        # autobuild the YAML refine-step conditions may allow ligandfit in the
        # valid list (e.g. has_refined_model=True satisfies some conditions),
        # but there is no actual ligand to fit.  When the duplicate-detection
        # fallback exhausts autobuild_denmod it then picks ligandfit as "next
        # available program", fitting amino acids as fake ligands and wrecking
        # the R-free.  Remove ligandfit here unless the user explicitly asked
        # for it OR a ligand file is actually present in the session files.
        if ("phenix.ligandfit" in valid and
                not context.get("user_wants_ligandfit") and
                not context.get("has_ligand_file")):
            valid.remove("phenix.ligandfit")

        # === APPLY USER DIRECTIVES ===
        # This is the SINGLE PLACE where directives affect valid_programs
        if directives:
            valid = self._apply_directives(valid, directives, step_name, context,
                                           experiment_type=experiment_type)

        if _diag:
            print("  [DIAG] after _apply_directives: %s" % valid)

        # Add STOP if validation done and at target
        if step_name == "validate" and context.get("validation_done"):
            if "STOP" not in valid:
                valid.append("STOP")
                if _diag:
                    print("  [GATE] STOP added: validation_done=True "
                          "step=validate r_free=%s"
                          % context.get("r_free"))

        # Special: also allow refinement during validate step (user can choose more refinement)
        # BUT respect max_refine_cycles directive AND _is_at_target (e.g., hopeless R-free)
        if step_name == "validate":
            max_refine = directives.get("stop_conditions", {}).get("max_refine_cycles") if directives else None
            # I1 fix: use rsr_count for cryoem, refine_count for xray
            if experiment_type == "cryoem":
                refine_count = context.get("rsr_count", 0)
            else:
                refine_count = context.get("refine_count", 0)
            refine_allowed = (max_refine is None) or (refine_count < max_refine)
            at_target = self._is_at_target(context, experiment_type)

            if refine_allowed and not at_target:
                if experiment_type == "xray":
                    if "phenix.refine" not in valid:
                        valid.append("phenix.refine")
                elif experiment_type == "cryoem":
                    if "phenix.real_space_refine" not in valid:
                        valid.append("phenix.real_space_refine")
            elif at_target:
                # Remove refinement programs — further refinement won't help
                for prog in ["phenix.refine", "phenix.real_space_refine"]:
                    if prog in valid:
                        valid.remove(prog)
                # Add STOP — model is at target or hopeless
                if "STOP" not in valid:
                    valid.append("STOP")
                    if _diag:
                        print("  [GATE] STOP added: at_target=True "
                              "step=validate refine_count=%s r_free=%s"
                              % (refine_count, context.get("r_free")))

            if _diag:
                print("  [DIAG] validate step filter: at_target=%s, "
                      "refine_allowed=%s, after filter: %s"
                      % (at_target, refine_allowed, valid))

        # Refinement step: remove refinement when at target or max cycles reached
        # Exception: post-ligandfit refinement is always allowed (model changed)
        # Safety net: _is_at_target() already returns False when ligandfit needs
        # refine as a prerequisite, so at_target should be False in that case.
        # The refine_is_ligandfit_prereq check below is defense-in-depth.
        if step_name == "refine":
            needs_post_ligandfit = context.get("needs_post_ligandfit_refine", False)
            max_refine = directives.get("stop_conditions", {}).get("max_refine_cycles") if directives else None
            # I1 fix: use rsr_count for cryoem, refine_count for xray
            if experiment_type == "cryoem":
                refine_count = context.get("rsr_count", 0)
            else:
                refine_count = context.get("refine_count", 0)
            refine_allowed = (max_refine is None) or (refine_count < max_refine)
            at_target = self._is_at_target(context, experiment_type)

            # Refinement is a prerequisite for ligandfit: it produces the
            # map_coeffs_mtz that ligandfit needs.  Keep refine available when:
            #   - user wants ligandfit AND it hasn't been done yet
            #   - no refinement has run in this session (refine_count == 0)
            refine_is_ligandfit_prereq = (
                context.get("user_wants_ligandfit", False) and
                not context.get("ligandfit_done") and
                not context.get("has_ligand_fit") and
                refine_count == 0
            )

            if (not refine_allowed or at_target) and not needs_post_ligandfit and not refine_is_ligandfit_prereq:
                # Remove refinement programs from the step's own list
                for prog in ["phenix.refine", "phenix.real_space_refine"]:
                    if prog in valid:
                        valid.remove(prog)

            # Add STOP if at target AND not in post-ligandfit refinement
            # AND not keeping refine as a prerequisite for ligandfit
            if at_target and not needs_post_ligandfit and not refine_is_ligandfit_prereq:
                if "STOP" not in valid:
                    valid.append("STOP")
                    if _diag:
                        print("  [GATE] STOP added: at_target=True "
                              "step=refine refine_count=%s r_free=%s"
                              % (refine_count, context.get("r_free")))

            if _diag:
                print("  [DIAG] refine step filter: at_target=%s, "
                      "refine_allowed=%s, refine_count=%s, "
                      "after filter: %s"
                      % (at_target, refine_allowed,
                         refine_count, valid))

        # NOTE: skip_validation STOP handling is now in _apply_directives()
        # No duplicate check needed here

        # ── Model-vs-Data Gate (v114.1) ──────────
        # When model_is_placed has been confirmed by
        # evidence (model_vs_data CC, refine R-free),
        # remove destructive programs that would wipe
        # out the existing model.
        if context.get("model_is_placed_confirmed"):
          _destructive = {
            "phenix.phaser",
            "phenix.autosol",
            "phenix.predict_and_build",
          }
          _before = len(valid)
          valid = [
            p for p in valid
            if p not in _destructive
          ]
          if len(valid) < _before:
            removed = _destructive & set(
              valid + list(_destructive))
            # (logging handled by caller)

        # If no valid programs available, return STOP (stuck state)
        if not valid:
            if _diag:
                print("  [GATE_FAIL] valid_programs empty → STOP")
                print("  [GATE_FAIL]   step=%s experiment=%s"
                      % (step_name, experiment_type))
                _active = sorted(k for k in context
                    if context[k] not in (None, False, 0, [], ""))
                print("  [GATE_FAIL]   active context: %s"
                      % _active[:20])
            return ["STOP"]

        if _diag:
            print("  [DIAG] FINAL valid_programs: %s" % valid)
        return valid

    def _apply_directives(self, valid_programs, directives, step_name, context=None,
                          experiment_type="xray"):
        """
        Apply user directives to modify valid programs list.

        This is the SINGLE PLACE where directives affect what programs are available.
        All directive-related valid_programs logic should be here.

        Handles:
        - Adding STOP when skip_validation=true (allows stopping without validation)
        - Adding after_program target (ensures the tutorial target is runnable)
        - PRIORITIZING after_program (move to front of list so LLM chooses it)
        - Workflow preferences: skip_programs, prefer_programs
        - REMOVING refinement programs when max_refine_cycles reached
        - Adding programs mentioned in program_settings (user wants to configure them)

        Args:
            valid_programs: List of valid program names from workflow step
            directives: User directives dict from directive_extractor
            step_name: Current workflow step name
            context: Workflow context dict (contains refine_count, etc.)

        Returns:
            list: Modified valid programs list
        """
        if not directives:
            return valid_programs

        result = list(valid_programs)
        modifications = []  # Track what we changed for logging

        # Get stop_conditions
        stop_cond = directives.get("stop_conditions", {})

        # 0. Check max_refine_cycles - REMOVE refinement programs if limit reached
        # Exception: post-ligandfit refinement is always allowed because the model
        # changed (ligand added) and must be re-refined for scientific validity.
        max_refine_cycles = stop_cond.get("max_refine_cycles")
        needs_post_ligandfit = context.get("needs_post_ligandfit_refine", False) if context else False
        if max_refine_cycles is not None and context and not needs_post_ligandfit:
            # I1 fix: use rsr_count for cryoem, refine_count for xray
            if experiment_type == "cryoem":
                refine_count = context.get("rsr_count", 0)
            else:
                refine_count = context.get("refine_count", 0)
            if refine_count >= max_refine_cycles:
                # Remove refinement programs from valid list
                refine_programs = ["phenix.refine", "phenix.real_space_refine"]
                before_count = len(result)
                result = [p for p in result if p not in refine_programs]
                if len(result) < before_count:
                    modifications.append(
                        "Removed refinement (max_refine_cycles=%d reached, count=%d)" % (
                            max_refine_cycles, refine_count))
                # I1 fix: controlled landing — inject validate programs rather
                # than bare STOP.  max_refine_cycles means "stop refining and
                # proceed to validation", NOT "abort immediately".
                # Only inject if we're currently in refine step AND validation
                # has not yet been done.
                if (step_name in ("refine", "ready_to_refine") and
                        not (context.get("validation_done"))):
                    validate_progs_xray  = ["phenix.molprobity", "phenix.model_vs_data",
                                            "phenix.map_correlations"]
                    validate_progs_cryoem = ["phenix.molprobity", "phenix.validation_cryoem",
                                             "phenix.map_correlations"]
                    vprogs = (validate_progs_cryoem if experiment_type == "cryoem"
                              else validate_progs_xray)
                    added = []
                    for vp in vprogs:
                        if vp not in result:
                            result.append(vp)
                            added.append(vp)
                    if added:
                        modifications.append(
                            "Added validate programs %s (max_refine_cycles landed)" % added)
                # Always allow STOP (user can choose to stop without validating)
                if "STOP" not in result:
                    result.append("STOP")
                    modifications.append("Added STOP (max_refine_cycles reached)")

        # 1. If skip_validation is set, add STOP to valid programs
        #    This allows the user to stop even without completing validation
        if stop_cond.get("skip_validation"):
            if "STOP" not in result:
                result.append("STOP")
                modifications.append("Added STOP (skip_validation directive)")

        # 1b. Check for programs mentioned in program_settings
        #     If user provided settings for a program, they likely want to run it
        program_settings = directives.get("program_settings", {})
        for prog_name in program_settings.keys():
            if prog_name == "default":
                continue  # Skip default settings
            if prog_name not in result:
                # Don't re-add run_once programs that have already completed
                if self._is_program_already_done(prog_name, context):
                    modifications.append("Skipped %s from program_settings (already completed)" % prog_name)
                    continue
                # Check if we should add this program based on prerequisites
                should_add = self._check_program_prerequisites(prog_name, context, step_name)
                if should_add:
                    result.insert(0, prog_name)
                    modifications.append("Added %s (has program_settings)" % prog_name)

        # 2. If after_program is set, ensure that program is available AND prioritized
        #    This handles tutorial cases where the user wants to run a specific
        #    program. We move it to the front so the LLM is more likely to choose it.
        #    BUT: Don't add programs if prerequisites aren't met!
        #    ALSO: If after_program has already been run, add STOP to valid programs
        after_program = stop_cond.get("after_program")
        if after_program:
            # Check if after_program has already been completed
            # If so, the user's workflow is done - add STOP
            after_program_done = context.get("last_program") == after_program if context else False

            # Also check program-specific done flags
            if not after_program_done and context:
                prog_key = after_program.replace("phenix.", "").replace(".", "_")
                # Check for programs that have been run (e.g., refine_count > 0 for phenix.refine)
                if after_program == "phenix.refine" and context.get("refine_count", 0) > 0:
                    after_program_done = True
                elif after_program == "phenix.real_space_refine" and context.get("rsr_count", 0) > 0:
                    after_program_done = True
                elif after_program == "phenix.ligandfit" and context.get("ligandfit_done"):
                    after_program_done = True
                elif after_program == "phenix.polder" and context.get("polder_done"):
                    after_program_done = True
                # General fallback: check run_once done flags from YAML
                elif self._is_program_already_done(after_program, context):
                    after_program_done = True

            if after_program_done:
                # Quality override: don't declare done when metrics
                # clearly indicate the result is bad.  The program
                # "ran" but didn't accomplish its purpose.
                # AF_POMGNT2: refine ran once (refine_count=1) but
                # R-free=0.52 — autobuild hasn't been tried and
                # could dramatically improve the model.
                if context:
                    _rf = context.get("r_free")
                    _cc = context.get("map_cc")
                    _exp = experiment_type
                    _ab_done = context.get(
                        "autobuild_done", False)
                    if (_exp == "xray" and _rf is not None
                            and _rf > 0.50
                            and not _ab_done):
                        after_program_done = False
                        modifications.append(
                            "Overrode after_program_done "
                            "(R-free=%.3f > 0.50, "
                            "autobuild not tried)"
                            % _rf)
                    elif (_exp == "cryoem"
                          and _cc is not None
                          and _cc < 0.70):
                        after_program_done = False
                        modifications.append(
                            "Overrode after_program_done "
                            "(CC=%.3f < 0.70)" % _cc)

            if after_program_done:
                # User's workflow is complete — replace the entire valid_programs
                # list with just STOP.  Simply appending STOP is not enough: the
                # step detector may have already added programs like
                # predict_and_build to 'result', and the LLM (or fallback) will
                # pick those instead of stopping.
                non_stop = [p for p in result if p != "STOP"]
                if non_stop:
                    modifications.append(
                        "Cleared programs %s (after_program %s completed, workflow done)" % (
                        non_stop, after_program))
                result[:] = ["STOP"]
                modifications.append("Set valid_programs=[STOP] (after_program %s completed)" % after_program)
            else:
                # after_program not yet run - add it to valid programs
                should_add = self._check_program_prerequisites(after_program, context, step_name)

                if should_add:
                    if after_program in result:
                        # Already in list - move to front to prioritize
                        result.remove(after_program)
                        result.insert(0, after_program)
                        modifications.append("Prioritized %s (after_program directive)" % after_program)
                    else:
                        # Not in list - add at front
                        result.insert(0, after_program)
                        modifications.append("Added %s (after_program directive)" % after_program)
                else:
                    modifications.append("Skipped adding %s (prerequisites not met)" % after_program)

        # 3. Apply workflow preferences
        workflow_prefs = directives.get("workflow_preferences", {})

        # Handle use_mr_sad - phaser first, then autosol (handled by experimental_phasing phase)
        # In obtain_model phase: prioritize phaser (need to place model first)
        # In experimental_phasing phase: autosol is already the only option
        if workflow_prefs.get("use_mr_sad"):
            if "phenix.phaser" in result:
                result.remove("phenix.phaser")
                result.insert(0, "phenix.phaser")
                modifications.append("Prioritized phenix.phaser (use_mr_sad: place model first)")
            # Don't run autosol until phaser has placed the model
            # This applies to ALL phases before experimental_phasing
            phaser_done = context.get("phaser_done", False) if context else False
            if not phaser_done and "phenix.autosol" in result:
                result.remove("phenix.autosol")
                modifications.append("Removed phenix.autosol (use_mr_sad: phaser must run first)")

        # Handle use_experimental_phasing - prioritize autosol over predict_and_build
        if workflow_prefs.get("use_experimental_phasing") and not workflow_prefs.get("use_mr_sad"):
            # Move autosol to front if it's in the list
            if "phenix.autosol" in result:
                result.remove("phenix.autosol")
                result.insert(0, "phenix.autosol")
                modifications.append("Prioritized phenix.autosol (use_experimental_phasing)")
            # Also deprioritize predict_and_build (move to end)
            if "phenix.predict_and_build" in result:
                result.remove("phenix.predict_and_build")
                result.append("phenix.predict_and_build")
                modifications.append("Deprioritized phenix.predict_and_build (use_experimental_phasing)")

        # Handle use_molecular_replacement - prioritize phaser/predict_and_build over autosol
        if workflow_prefs.get("use_molecular_replacement"):
            # Deprioritize autosol (move to end)
            if "phenix.autosol" in result:
                result.remove("phenix.autosol")
                result.append("phenix.autosol")
                modifications.append("Deprioritized phenix.autosol (use_molecular_replacement)")
            # Prioritize predict_and_build or phaser
            if "phenix.predict_and_build" in result:
                result.remove("phenix.predict_and_build")
                result.insert(0, "phenix.predict_and_build")
                modifications.append("Prioritized phenix.predict_and_build (use_molecular_replacement)")
            elif "phenix.phaser" in result:
                result.remove("phenix.phaser")
                result.insert(0, "phenix.phaser")
                modifications.append("Prioritized phenix.phaser (use_molecular_replacement)")

        # Bug D fix: deprioritize autosol when anomalous signal is negligible
        # and predict_and_build is also available.  When anomalous_measurability
        # < 0.05 the signal is too weak for reliable experimental phasing, so the
        # LLM should strongly prefer predict_and_build.  We do not remove autosol
        # entirely (it remains as a fallback) but move it to the end of the list.
        # Exemptions: explicit use_mr_sad / use_experimental_phasing overrides.
        if (not workflow_prefs.get("use_mr_sad") and
                not workflow_prefs.get("use_experimental_phasing") and
                "phenix.autosol" in result and
                "phenix.predict_and_build" in result and
                not context.get("autosol_done") and
                not context.get("phaser_done")):
            anomalous_meas = context.get("anomalous_measurability") or 0
            if anomalous_meas < 0.05:
                result.remove("phenix.autosol")
                result.append("phenix.autosol")
                modifications.append(
                    "Deprioritized phenix.autosol "
                    "(anomalous_measurability=%.3f < 0.05, "
                    "prefer predict_and_build)" % anomalous_meas)

        # Handle start_with_program - add to valid programs if user explicitly requested it
        # This is for multi-step requests like "run polder then refine"
        stop_cond = directives.get("stop_conditions", {})
        start_with = stop_cond.get("start_with_program")
        if start_with and start_with not in result:
            if self._is_program_already_done(start_with, context):
                modifications.append("Skipped start_with %s (already completed)" % start_with)
            else:
                result.insert(0, start_with)
                modifications.append("Added start_with_program: %s" % start_with)

        # Remove skipped programs
        skip_programs = workflow_prefs.get("skip_programs", [])
        if skip_programs:
            before_count = len(result)
            result = [p for p in result if p not in skip_programs]
            if len(result) < before_count:
                modifications.append("Removed skipped programs: %s" % skip_programs)

        # Reorder to prefer certain programs (move to front)
        prefer_programs = workflow_prefs.get("prefer_programs", [])
        if prefer_programs:
            # Special case: when user explicitly wants ligandfit but it was
            # excluded by YAML conditions, inject it — BUT only if refinement
            # has already run (refine_count > 0).  When refine_count==0,
            # refinement is a prerequisite for ligandfit (produces the map
            # coefficients MTZ).  _is_at_target() already keeps refine
            # available in this case, so the YAML workflow naturally flows:
            #   refine (produces map_coeffs) → ligandfit (YAML conditions pass)
            if ("phenix.ligandfit" in prefer_programs and
                "phenix.ligandfit" not in result and
                context and
                context.get("user_wants_ligandfit") and
                not context.get("ligandfit_done") and
                step_name in ("refine", "ready_to_refine")):
                refine_key = "rsr_count" if experiment_type == "cryoem" else "refine_count"
                has_refined = (context.get(refine_key, 0) > 0)
                r_free = context.get("r_free")
                if has_refined and (r_free is None or r_free < 0.50):
                    result.insert(0, "phenix.ligandfit")
                    modifications.append(
                        "Injected phenix.ligandfit (user requested, model quality sufficient)")
                elif not has_refined:
                    modifications.append(
                        "Deferred phenix.ligandfit injection (refine_count=0, "
                        "refinement must run first to produce map coefficients)")

            preferred = [p for p in prefer_programs
                         if p in result]
            dropped = [p for p in prefer_programs
                       if p not in result and p]
            others = [p for p in result
                      if p not in prefer_programs]
            if preferred:
                result = preferred + others
                modifications.append(
                    "Prioritized: %s" % preferred)
            if dropped:
                modifications.append(
                    "Preferred not available: %s"
                    % dropped)
                logger.debug(
                    "Preferred programs %s not in "
                    "valid list %s — skipped",
                    dropped, result,
                )

        # Log modifications for debug and store for
        # caller retrieval.
        for mod in modifications:
            logger.debug(
                "Directive modification: %s", mod)
        self._last_directive_modifications = modifications

        return result

    def _is_program_already_done(self, program, context):
        """Check if a program has already been completed.

        Returns True when:
        1. The program has a run_once done_tracking strategy and its done
           flag is set (original check).
        2. The program is a non-count program with a program-specific done
           flag that is already True in context.  This mirrors the
           belt-and-suspenders filter in get_valid_programs() and prevents
           _apply_directives from re-adding completed programs via the
           program_settings path (e.g., autosol re-running after success
           because the directive extractor found autosol settings in user
           advice).

        Args:
            program: Program name (e.g., 'phenix.map_symmetry')
            context: Workflow context dict with done flags

        Returns:
            bool: True if program is done and should not be re-run
        """
        if not context:
            return False
        prog_def = get_program(program)
        if not prog_def:
            return False
        tracking = prog_def.get("done_tracking", {})
        strategy = tracking.get("strategy", "set_flag")

        # Check 1: run_once programs (original)
        if strategy == "run_once" or tracking.get("run_once"):
            done_key = tracking.get("flag", "")
            if done_key and context.get(done_key):
                return True

        # Check 2: non-count programs with program-specific done flag
        # (mirrors belt-and-suspenders filter in get_valid_programs)
        # Count-strategy programs (refine, rsr, phaser) intentionally repeat.
        # Shared flags (e.g., "validation_done") are excluded — only filter
        # when the flag name contains the program short name.
        if strategy != "count":
            flag = tracking.get("flag", "")
            short = program.replace("phenix.", "").replace(".", "_")
            if flag and short in flag and context.get(flag):
                return True

        return False

    def _check_program_prerequisites(self, program, context, step_name):
        """
        Check if a program's prerequisites are met.

        Args:
            program: Program name to check
            context: Workflow context dict
            step_name: Current workflow step name

        Returns:
            bool: True if program can be added
        """
        if not context:
            return True  # No context to check against

        # Don't add refinement programs if we're in an early step without a model
        if program in ("phenix.refine", "phenix.real_space_refine"):
            has_model_to_refine = (
                context.get("refine_count", 0) > 0 or
                context.get("has_refined_model") or
                context.get("has_placed_model") or
                step_name in ("refine", "validate")
            )
            if not has_model_to_refine:
                return False

        # Don't add ligandfit if refinement hasn't been done yet
        # (ligandfit requires map coefficients from refined MTZ)
        if program == "phenix.ligandfit":
            has_refined = (
                context.get("refine_count", 0) > 0 or
                context.get("rsr_count", 0) > 0 or
                context.get("has_refined_model")
            )
            if not has_refined:
                return False

        # predict_and_build: for X-ray, only add if xtriage has been run (we need resolution)
        # for cryo-EM, only add if mtriage has been run (we need resolution)
        # Exception: always allow if we're in a step that includes predict_and_build
        if program == "phenix.predict_and_build":
            # Don't re-run if the full workflow already completed
            if context.get("predict_full_done"):
                return False

            # In STEPWISE mode, after prediction is done, don't allow predict_and_build again
            # The workflow should go: predict(stop_after) -> process_predicted -> phaser -> refine
            # NOT: predict(stop_after) -> predict(full)
            automation_path = context.get("automation_path", "automated")
            if automation_path == "stepwise" and context.get("predict_done"):
                return False  # Force stepwise path through process_predicted -> phaser

            # Check if we're already in a step that includes predict_and_build
            # Note: molecular_replacement is NOT included - that step is for stepwise path only
            if step_name in ("obtain_model", "dock_model"):
                return True  # Let the step conditions handle it

            # For X-ray, require xtriage to be done (to get resolution for building)
            # For cryo-EM, require mtriage to be done
            if context.get("xtriage_done") or context.get("mtriage_done"):
                return True

            # Allow prediction-only if explicitly requested, but warn
            # (This case is handled by command builder forcing stop_after_predict=True)
            return False  # Don't add to early steps - let workflow proceed normally

        # Default: allow
        # phenix.autosol prerequisites:
        # - Always needs xtriage first (to detect anomalous signal, resolution)
        # - For MR-SAD (explicit or implicit), needs phaser to place model first
        #   Implicit MR-SAD: user has search model + anomalous data
        if program == "phenix.autosol":
            if not context.get("xtriage_done"):
                return False
            # MR-SAD: phaser must run before autosol when we have a search model
            # and anomalous data (whether or not use_mr_sad directive is set)
            needs_phaser_first = (
                context.get("use_mr_sad") or
                (context.get("has_search_model") and context.get("has_anomalous"))
            )
            if needs_phaser_first and not context.get("phaser_done"):
                return False
            # PredictAndBuild guard: if a sequence file and a model are both present
            # and predict_and_build has not yet completed, the appropriate workflow is
            # predict_and_build (MR with an AF-predicted model), NOT autosol (SAD/MAD
            # experimental phasing).  Block autosol unless the user has explicitly
            # requested MR-SAD (use_mr_sad) — in that case phaser must run first
            # anyway (caught above), so this branch is only reached in the pure
            # experimental-phasing case where both signals appear together.
            if (context.get("has_sequence") and
                    context.get("has_model_for_mr") and
                    not context.get("predict_full_done") and
                    not context.get("use_mr_sad")):
                return False
        return True

    # All condition keywords recognized by _check_conditions.
    # Any keyword in workflows.yaml NOT in this set is a bug — see A6 guard below.
    _KNOWN_CONDITION_KEYS = frozenset([
        "has",          # {"has": "sequence"}         → context["has_sequence"] is truthy
        "has_any",      # {"has_any": ["a", "b"]}     → any(context["has_a"], context["has_b"])
        "not_has",      # {"not_has": "full_map"}      → not context["has_full_map"]
        "not_done",     # {"not_done": "autobuild"}    → not context["autobuild_done"]
        "not",          # {"not": "model"}              → not context["has_model"]
        "r_free",       # {"r_free": "> 0.35"}         → metric comparison
        "map_cc",       # {"map_cc": "> 0.6"}          → metric comparison
        "refine_count", # {"refine_count": "> 0"}      → metric comparison
        "rsr_count",    # {"rsr_count": "> 0"}         → metric comparison
    ])

    def _check_conditions(self, prog_entry, context):
        """Check if program conditions are met.

        Each condition in prog_entry["conditions"] is a single-key dict.  All
        keys are processed independently (implicit AND — all must pass).

        Returns True only if every condition is satisfied.
        """
        conditions = prog_entry.get("conditions", [])

        for cond in conditions:
            if not isinstance(cond, dict):
                continue

            # ── A6 guard: catch unknown keywords early ─────────────────────
            unknown = set(cond.keys()) - self._KNOWN_CONDITION_KEYS
            if unknown:
                msg = (
                    "Unknown condition keyword(s) %s in program entry %r — "
                    "add a handler to _check_conditions() or fix the YAML. "
                    "Condition is being SKIPPED (treated as satisfied)."
                    % (sorted(unknown), prog_entry.get("program", "?"))
                )
                if os.environ.get("PHENIX_AGENT_STRICT"):
                    raise ValueError(msg)
                else:
                    # Use print for now; callers may not pass a logger.
                    # Deliberately loud so it isn't buried in INFO spam.
                    print("ERROR: " + msg)
                # Skip the unknown keys but continue evaluating known ones.

            # ── has ────────────────────────────────────────────────────────
            # {"has": "sequence"} → context["has_sequence"] must be truthy
            if "has" in cond:
                key = "has_" + cond["has"]
                if not context.get(key):
                    return False

            # ── has_any ───────────────────────────────────────────────────
            # {"has_any": ["full_map", "optimized_full_map"]}
            # → at least one of context["has_full_map"] / context["has_optimized_full_map"]
            if "has_any" in cond:
                keys = ["has_" + k for k in cond["has_any"]]
                if not any(context.get(k) for k in keys):
                    return False

            # ── not_has ───────────────────────────────────────────────────
            # {"not_has": "full_map"} → context["has_full_map"] must be falsy
            # Relies on _categorize_files() already excluding zero-byte/corrupt
            # files so that has_full_map reflects genuine usable file presence.
            if "not_has" in cond:
                key = "has_" + cond["not_has"]
                if context.get(key):
                    return False

            # ── not_done ──────────────────────────────────────────────────
            # {"not_done": "autobuild"} → context["autobuild_done"] must be falsy
            if "not_done" in cond:
                key = cond["not_done"] + "_done"
                if context.get(key):
                    return False

            # ── not ───────────────────────────────────────────────────────
            # {"not": "model"} → context["has_model"] must be falsy
            # Negates a has-style check: the named resource must be absent.
            if "not" in cond:
                key = "has_" + cond["not"]
                if context.get(key):
                    return False

            # ── metric comparisons ────────────────────────────────────────
            # {"r_free": "> 0.35"}, {"map_cc": "> 0.6"},
            # {"refine_count": "> 0"}, {"rsr_count": "> 0"}
            for metric in ("r_free", "map_cc", "refine_count", "rsr_count"):
                if metric in cond:
                    if not self._check_metric_condition(context, metric, cond[metric]):
                        return False

        return True

    def explain_unavailable_program(self, experiment_type, program, context):
        """
        Explain why a specific program is not available.

        Args:
            experiment_type: "xray" or "cryoem"
            program: Program name to check
            context: Context dict

        Returns:
            str: Explanation of why the program is unavailable, or None if it's available
        """
        steps = get_workflow_steps(experiment_type)
        reasons = []
        found_program = False

        # Find the program in any step
        for step_name, step_def in steps.items():
            if not isinstance(step_def, dict):
                continue
            step_programs = step_def.get("programs", [])

            for prog_entry in step_programs:
                prog_name = None
                conditions = []

                if isinstance(prog_entry, str):
                    prog_name = prog_entry
                elif isinstance(prog_entry, dict):
                    prog_name = prog_entry.get("program")
                    conditions = prog_entry.get("conditions", [])

                if prog_name != program:
                    continue

                found_program = True

                # Found the program - check each condition
                for cond in conditions:
                    if isinstance(cond, dict):
                        # Check "has" condition
                        if "has" in cond:
                            key = "has_" + cond["has"]
                            if not context.get(key):
                                reasons.append("missing required file: %s" % cond["has"])

                        # Check "not_has" condition
                        if "not_has" in cond:
                            key = "has_" + cond["not_has"]
                            if context.get(key):
                                reasons.append(
                                    "file already present: %s (program only runs without it)"
                                    % cond["not_has"])

                        # Check "not_done" condition
                        if "not_done" in cond:
                            key = cond["not_done"] + "_done"
                            if context.get(key):
                                reasons.append("already completed: %s" % cond["not_done"])

                        # Check metric conditions
                        for metric in ["r_free", "map_cc", "refine_count", "rsr_count"]:
                            if metric in cond:
                                value = context.get(metric)
                                condition_str = cond[metric]
                                # Only report failure if we have a value and it doesn't satisfy
                                # (matching the actual _check_metric_condition which allows None)
                                if value is not None:
                                    if not self._check_metric_condition(context, metric, condition_str):
                                        reasons.append("%s=%s does not satisfy condition '%s'" % (
                                            metric, value, condition_str))
                                elif metric == "refine_count":
                                    # Special case: refine_count=0 means no successful refinements
                                    # Report this as a reason since it's a common issue
                                    reasons.append("%s=0 does not satisfy condition '%s'" % (metric, condition_str))

        # If program not found in any step, it's not relevant for this experiment type
        if not found_program:
            return None  # Not a program for this workflow

        # Check run_once
        prog_def = get_program(program)
        if prog_def:
            tracking = prog_def.get("done_tracking", {})
            if tracking.get("strategy") == "run_once" or tracking.get("run_once"):
                prog_done_key = tracking.get("flag", "")
                if prog_done_key and context.get(prog_done_key):
                    reasons.append("run_once program already executed")

        # Check MR-SAD guard (autosol blocked when search model + anomalous + no phaser)
        if program == "phenix.autosol":
            if not context.get("xtriage_done"):
                reasons.append("xtriage must run first")
            elif (context.get("has_search_model") and
                  (context.get("has_anomalous") or context.get("use_mr_sad")) and
                  not context.get("phaser_done")):
                reasons.append("MR-SAD: phaser must place model before autosol runs")
            elif (context.get("has_sequence") and
                  context.get("has_model_for_mr") and
                  not context.get("predict_full_done") and
                  not context.get("use_mr_sad")):
                reasons.append(
                    "predict_and_build workflow detected: a sequence file and model "
                    "are present and predict_and_build has not yet run. Use "
                    "phenix.predict_and_build (MR with predicted model) instead of "
                    "phenix.autosol (SAD/MAD experimental phasing). To override, set "
                    "use_mr_sad=true in workflow_preferences."
                )

        # P3 fix: autobuild requires phasing to be done in X-ray mode.
        # Without a phased MTZ (PHIB/FOM columns), autobuild will crash
        # with 'MTZ lacks phase/FOM columns'.  This explanation surfaces
        # in the PLAN node via unavailable_explanations so the LLM can
        # reason about it (and so that it appears in the session log).
        if (program == "phenix.autobuild" and
                experiment_type == "xray" and
                not context.get("phaser_done") and
                not context.get("autosol_done") and
                not context.get("has_placed_model_from_history")):
            reasons.append(
                "autobuild requires a phased MTZ (PHIB/FOM). "                "Run phenix.phaser or phenix.autosol first."
            )

        if reasons:
            return "; ".join(reasons)
        return None

    def _check_metric_condition(self, context, metric, condition):
        """Check a metric condition like "> 0.35" or "< target_r_free"."""
        value = context.get(metric)
        if value is None:
            return True  # If we don't have the metric, don't block

        # Ensure numeric type (values may arrive as strings from log analysis)
        try:
            value = float(value)
        except (ValueError, TypeError):
            return True  # Non-numeric value, don't block

        # Parse condition string
        match = re.match(r'([<>=!]+)\s*(\S+)', str(condition))
        if not match:
            return True

        op, threshold_str = match.groups()

        # Handle special threshold names
        if threshold_str == "target_r_free":
            threshold = get_metric_threshold("r_free", "acceptable", context.get("resolution"))
        elif threshold_str == "autobuild_threshold":
            threshold = 0.35
        else:
            try:
                threshold = float(threshold_str)
            except ValueError:
                return True

        if threshold is None:
            return True

        # Evaluate
        if op == ">":
            return value > threshold
        elif op == ">=":
            return value >= threshold
        elif op == "<":
            return value < threshold
        elif op == "<=":
            return value <= threshold
        elif op == "==":
            return abs(value - threshold) < 0.001

        return True

    # =========================================================================
    # TARGET ACCESS
    # =========================================================================

    def get_target(self, experiment_type, metric, resolution=None):
        """
        Get quality target for a metric.

        Args:
            experiment_type: "xray" or "cryoem"
            metric: Metric name (e.g., "r_free")
            resolution: Optional resolution for resolution-dependent targets

        Returns:
            float: Target value or None
        """
        targets = get_workflow_targets(experiment_type)
        target_def = targets.get(metric, {})

        # Check resolution-dependent targets
        if resolution and "by_resolution" in target_def:
            for entry in target_def["by_resolution"]:
                range_min, range_max = entry.get("range", [0, 999])
                if range_min <= resolution < range_max:
                    return entry.get("value")

        return target_def.get("default")

    # =========================================================================
    # HIGH-LEVEL API
    # =========================================================================

    def get_workflow_state(self, experiment_type, files, history_info, analysis=None,
                           directives=None, maximum_automation=True, session_info=None,
                           files_local=True):
        """
        Get complete workflow state (compatible with workflow_state.py output).

        Args:
            experiment_type: "xray" or "cryoem"
            files: Categorized files dict
            history_info: Analyzed history dict
            analysis: Current log analysis
            directives: Optional user directives dict
            maximum_automation: If False, use stepwise path (process_predicted -> phaser)
            session_info: Optional session state dict; used for client-supplied
                          data that the server cannot read directly (e.g.
                          unplaced_model_cell — CRYST1 from the client-side PDB)

        Returns:
            dict: Workflow state compatible with existing code
        """
        context = self.build_context(files, history_info, analysis, directives,
                                     session_info, files_local=files_local)

        # Add automation_path to context for program filtering
        context["automation_path"] = "automated" if maximum_automation else "stepwise"

        # S2c: promote unclassified PDB to search_model when docking is needed.
        # Must run after build_context (needs placement flags) and before
        # detect_step + get_valid_programs (both consume context and files).
        files, context = self._promote_unclassified_for_docking(
            files, context, experiment_type)

        step_info = self.detect_step(experiment_type, context)

        _diag = os.environ.get("PHENIX_AGENT_DIAG_VALID_PROGRAMS")
        if _diag:
            print("  [DIAG] detect_step: %s" % step_info.get("step"))
            print("  [DIAG] context: phaser_done=%s, refine_done=%s, "
                  "refine_count=%s, autobuild_done=%s, r_free=%s, "
                  "has_placed_model=%s, has_refined_model=%s, "
                  "has_sequence=%s, has_predicted_model=%s"
                  % (context.get("phaser_done"),
                     context.get("refine_done"),
                     context.get("refine_count"),
                     context.get("autobuild_done"),
                     context.get("r_free"),
                     context.get("has_placed_model"),
                     context.get("has_refined_model"),
                     context.get("has_sequence"),
                     context.get("has_predicted_model")))

        valid_programs = self.get_valid_programs(experiment_type, step_info, context, directives)

        # Get priority_when info for each valid program
        program_priorities = self._get_program_priorities(experiment_type, step_info, context)

        # Map YAML step name to original state name for compatibility
        step_name = step_info.get("step", "unknown")
        state_name = self._map_step_to_state(step_name, experiment_type)

        # Build reason string
        reason = step_info.get("reason", "")
        if step_info.get("goal"):
            reason += " - Goal: " + step_info["goal"]

        # Add metric info to reason
        if experiment_type == "xray" and context.get("r_free"):
            reason += " (R-free: %.3f)" % context["r_free"]
        elif experiment_type == "cryoem" and context.get("map_cc"):
            reason += " (CC: %.3f)" % context["map_cc"]

        # Check for stuck state (no programs available except STOP)
        if valid_programs == ["STOP"] and step_name not in ["complete", "validate"]:
            reason = "STUCK: " + reason
            if experiment_type == "xray":
                reason += " - Upload a sequence (.fa) for AlphaFold or a model (.pdb) for MR"
            else:
                reason += " - Upload a sequence (.fa) or model (.pdb)"

        # Check for forced_program from directives (after_program)
        forced_program = None
        if directives:
            stop_cond = directives.get("stop_conditions", {})
            after_program = stop_cond.get("after_program")
            if after_program and after_program in valid_programs:
                # Only force after_program when it is the last meaningful program
                # remaining. If other programs (e.g. phenix.mtriage) are still
                # available they should run first; forcing after_program here would
                # skip them entirely.
                other_valid = [p for p in valid_programs
                               if p != after_program and p != "STOP"]
                if not other_valid:
                    forced_program = after_program

        # P4 fix: remove session-blocked programs from valid_programs.
        # session_blocked_programs is populated by the PLAN node pivot handler
        # when should_pivot returns a 'session block' reason.  It persists
        # across cycles (valid_programs is reconstructed each cycle, so
        # per-cycle exclusion in the pivot handler alone does not survive).
        session_blocked = set(session_info.get(
            "session_blocked_programs", []) if session_info else [])
        if session_blocked:
            before_block = valid_programs[:]
            valid_programs = [
                p for p in valid_programs
                if p not in session_blocked
                and p.replace("phenix.", "") not in session_blocked
            ]
            if valid_programs != before_block:
                logger.debug(
                    "WORKFLOW: session_blocked filtered %s -> %s",
                    before_block, valid_programs)

        # P3 fix: remove autobuild from valid_programs when phasing prerequisites
        # are not met for X-ray experiments.  Without this, PLAN sees autobuild
        # as a valid choice, the LLM picks it, and VALIDATE has to reject it
        # after BUILD has already built the command — wasting an attempt slot
        # and confusing the LLM.  Filtering here is cleaner and mirrors the
        # session_blocked filtering above.
        # Exemptions: cryoem (different pipeline), has_placed_model_from_history
        # (refine_done etc. mean a map exists), autosol_done, phaser_done.
        if (experiment_type == "xray" and
                "phenix.autobuild" in valid_programs):
            _ab_ctx = context  # already built above
            _phasing_satisfied = (
                _ab_ctx.get("phaser_done") or
                _ab_ctx.get("autosol_done") or
                _ab_ctx.get("has_placed_model_from_history")
            )
            if not _phasing_satisfied:
                valid_programs = [
                    p for p in valid_programs
                    if p != "phenix.autobuild"
                ]
                logger.debug(
                    "WORKFLOW: autobuild excluded (P3): "
                    "phaser_done=%s autosol_done=%s "
                    "has_placed_model_from_history=%s",
                    _ab_ctx.get("phaser_done"),
                    _ab_ctx.get("autosol_done"),
                    _ab_ctx.get("has_placed_model_from_history"),
                )

        # Check for common programs that users might expect but aren't available
        # This helps explain why certain programs can't be run
        unavailable_explanations = {}
        common_programs = ["phenix.ligandfit", "phenix.autobuild", "phenix.phaser",
                          "phenix.predict_and_build", "phenix.dock_in_map",
                          "phenix.autosol"]
        for prog in common_programs:
            if prog not in valid_programs:
                explanation = self.explain_unavailable_program(experiment_type, prog, context)
                if explanation:
                    unavailable_explanations[prog] = explanation

        return {
            "state": state_name,  # Use mapped state name
            "experiment_type": experiment_type,
            "valid_programs": valid_programs,
            "program_priorities": program_priorities,  # NEW: priority_when triggers
            "forced_program": forced_program,  # NEW: from after_program directive
            "reason": reason,
            "conditions": {},
            "step_info": step_info,
            "context": context,
            "resolution": context.get("resolution"),
            "unavailable_explanations": unavailable_explanations,  # NEW: why programs aren't available
            # S2c: promoted files so detect_workflow_state doesn't overwrite with originals
            "categorized_files": files,
        }

    def _get_program_priorities(self, experiment_type, step_info, context):
        """
        Get programs that should be prioritized based on priority_when conditions.

        Args:
            experiment_type: "xray" or "cryoem"
            step_info: Phase detection result
            context: Context dict

        Returns:
            list: Programs that should be prioritized, in priority order
        """
        steps = get_workflow_steps(experiment_type)
        step_name = step_info.get("step", "")
        step_def = steps.get(step_name, {})

        prioritized = []
        step_programs = step_def.get("programs", [])

        for prog_entry in step_programs:
            if isinstance(prog_entry, dict):
                prog_name = prog_entry.get("program")
                priority_when = prog_entry.get("priority_when")

                if prog_name and priority_when:
                    # Check if priority_when condition is met
                    if self._check_priority_condition(priority_when, context):
                        prioritized.append(prog_name)

        return prioritized

    def _check_priority_condition(self, condition, context):
        """
        Check if a priority_when condition is satisfied.

        Args:
            condition: Condition string (e.g., "strong_anomalous")
            context: Context dict

        Returns:
            bool: True if condition is met
        """
        if condition == "strong_anomalous":
            return context.get("strong_anomalous", False)

        # Add more conditions as needed
        return False


# =============================================================================
# SINGLETON AND CONVENIENCE FUNCTIONS
# =============================================================================

_engine = None

def get_engine():
    """Get global WorkflowEngine instance."""
    global _engine
    if _engine is None:
        _engine = WorkflowEngine()
    return _engine


def detect_workflow_phase(experiment_type, files, history_info, analysis=None):
    """Convenience function to detect current step."""
    engine = get_engine()
    context = engine.build_context(files, history_info, analysis)
    return engine.detect_step(experiment_type, context)


# =============================================================================
# TESTING
# =============================================================================

if __name__ == "__main__":
    print("Testing WorkflowEngine...")
    print()

    engine = WorkflowEngine()

    # Test X-ray initial state
    files = {"data_mtz": ["data.mtz"], "sequence": ["seq.fa"]}
    history = {}
    context = engine.build_context(files, history)
    step = engine.detect_step("xray", context)
    programs = engine.get_valid_programs("xray", step, context)

    print("X-ray initial state:")
    print("  Step:", step["step"])
    print("  Reason:", step["reason"])
    print("  Valid programs:", programs)
    print()

    # Test X-ray after xtriage
    history = {"xtriage_done": True}
    context = engine.build_context(files, history)
    step = engine.detect_step("xray", context)
    programs = engine.get_valid_programs("xray", step, context)

    print("X-ray after xtriage:")
    print("  Step:", step["step"])
    print("  Reason:", step["reason"])
    print("  Valid programs:", programs)
    print()

    # Test cryo-EM initial state
    files = {"full_map": ["map.mrc"], "sequence": ["seq.fa"]}
    history = {}
    context = engine.build_context(files, history)
    step = engine.detect_step("cryoem", context)
    programs = engine.get_valid_programs("cryoem", step, context)

    print("Cryo-EM initial state:")
    print("  Step:", step["step"])
    print("  Reason:", step["reason"])
    print("  Valid programs:", programs)
    print()

    # Test target retrieval
    print("Quality targets:")
    print("  R-free at 2.0A:", engine.get_target("xray", "r_free", 2.0))
    print("  Map CC:", engine.get_target("cryoem", "map_cc"))
