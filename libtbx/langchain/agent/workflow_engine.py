"""
Workflow Engine for PHENIX AI Agent.

This module interprets workflows.yaml to determine:
- Current phase in workflow
- Valid programs for current phase
- Transitions to next phase
- Quality targets

The engine provides a higher-level abstraction than the hardcoded
state detection in workflow_state.py.

Usage:
    from libtbx.langchain.agent.workflow_engine import WorkflowEngine

    engine = WorkflowEngine()
    phase_info = engine.get_current_phase("xray", context)
    valid_programs = engine.get_valid_programs("xray", phase_info)
"""

from __future__ import absolute_import, division, print_function

import os
import re

# Import YAML loader
from libtbx.langchain.knowledge.yaml_loader import (
    get_workflow_phases,
    get_workflow_targets,
    get_metric_threshold,
    get_program,
)

# Import program registry for program info
from libtbx.langchain.agent.program_registry import ProgramRegistry


class WorkflowEngine:
    """
    Interprets workflow definitions from YAML.

    Provides phase detection and program validation based on
    workflow configuration rather than hardcoded logic.
    """

    def __init__(self):
        """Initialize the engine."""
        self._registry = ProgramRegistry()

    # =========================================================================
    # CONTEXT BUILDING
    # =========================================================================

    def build_context(self, files, history_info, analysis=None, directives=None):
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
        context = {
            # File availability - using semantic categories
            "has_data_mtz": bool(files.get("data_mtz")),
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
            "has_placed_model": self._has_placed_model(files, history_info, directives),
            "has_refined_model": self._has_refined_model(files, history_info),
            "has_ligand_fit": bool(files.get("ligand_fit_output")) or history_info.get("ligandfit_done", False),
            "has_optimized_full_map": self._has_optimized_map(files, history_info),

            # Complex program flags (need special logic beyond simple detection)
            "phaser_done": history_info.get("phaser_done", False),
            "predict_done": history_info.get("predict_done", False),
            "predict_full_done": history_info.get("predict_full_done", False),
            "autobuild_done": history_info.get("autobuild_done", False),
            "autosol_done": history_info.get("autosol_done", False),
            "autosol_success": history_info.get("autosol_success", False),
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
        }

        # =====================================================================
        # Auto-include done flags from program_registration
        # This automatically adds flags for all programs with run_once: true
        # =====================================================================
        from libtbx.langchain.knowledge.program_registration import get_all_done_flags

        for flag_name in get_all_done_flags():
            # Only add if not already in context (don't override complex flags)
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

        # Propagate MR-SAD workflow preference from directives
        if directives:
            workflow_prefs = directives.get("workflow_preferences", {})
            context["use_mr_sad"] = workflow_prefs.get("use_mr_sad", False)

            # CRITICAL: Treat skipped programs as "done" so the workflow can
            # advance past phases that require them.  Without this, skipping
            # xtriage causes the workflow to stay stuck in "analyze" with no
            # valid programs, because _detect_xray_phase checks xtriage_done
            # before allowing progression to later phases.
            skip_programs = workflow_prefs.get("skip_programs", [])
            if skip_programs:
                from libtbx.langchain.knowledge.program_registration import get_program_done_flag_map
                program_to_done_flag = get_program_done_flag_map()
                for prog in skip_programs:
                    done_flag = program_to_done_flag.get(prog)
                    if done_flag and not context.get(done_flag):
                        context[done_flag] = True
                for prog in skip_programs:
                    done_flag = program_to_done_flag.get(prog)
                    if done_flag and not context.get(done_flag):
                        context[done_flag] = True
        else:
            context["use_mr_sad"] = False

        return context

    def _get_metric(self, analysis, history_info, analysis_key, history_key):
        """Get metric from analysis or history."""
        if analysis and analysis.get(analysis_key):
            return analysis[analysis_key]
        return history_info.get(history_key)

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

    def _has_placed_model(self, files, history_info, directives=None):
        """
        Check if model is placed (after MR/building) or is ready-to-refine.

        A model is considered "placed" if:
        1. History shows MR, building, or prediction was done
        2. There's a file in a SPECIFIC positioned subcategory (refined,
           phaser_output, autobuild_output, docked, with_ligand, rsr_output)
        3. User explicitly says model is placed (model_is_placed directive)
        4. User's constraints imply model is placed (e.g., "refine", "fit ligand")
        5. User explicitly requests refinement/validation via directives

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
            # (e.g., "refine the model", "fit a ligand"), trust that inference
            constraints = directives.get("constraints", [])
            placement_keywords = [
                "refine", "refinement", "ligandfit", "fit ligand",
                "fit a ligand", "polder", "validate", "molprobity",
            ]
            for c in constraints:
                c_lower = c.lower() if isinstance(c, str) else ""
                if any(kw in c_lower for kw in placement_keywords):
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

            stop_conditions = directives.get("stop_conditions", {})
            after_program = stop_conditions.get("after_program", "")
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
    # Map YAML phase names to original hardcoded state names for compatibility

    XRAY_STATE_MAP = {
        "analyze": "xray_initial",
        "obtain_model": "xray_analyzed",
        "molecular_replacement": "xray_has_prediction",
        "build_from_phases": "xray_has_phases",
        "experimental_phasing": "xray_mr_sad",
        "refine": "xray_refined",
        "combine_ligand": "xray_combined",
        "validate": "xray_refined",  # Validation is part of refined state
        "complete": "complete",
    }

    CRYOEM_STATE_MAP = {
        "analyze": "cryoem_initial",
        "obtain_model": "cryoem_analyzed",
        "dock_model": "cryoem_has_prediction",
        "check_map": "cryoem_has_model",
        "optimize_map": "cryoem_has_model",
        "ready_to_refine": "cryoem_docked",  # Model docked, ready for first refinement
        "refine": "cryoem_refined",
        "validate": "cryoem_refined",  # Validation is part of refined state
        "complete": "complete",
    }

    def _map_phase_to_state(self, phase, experiment_type):
        """Map YAML phase name to original state name."""
        if experiment_type == "xray":
            return self.XRAY_STATE_MAP.get(phase, phase)
        elif experiment_type == "cryoem":
            return self.CRYOEM_STATE_MAP.get(phase, phase)
        return phase

    # =========================================================================
    # PHASE DETECTION
    # =========================================================================

    def detect_phase(self, experiment_type, context):
        """
        Detect current workflow phase based on context.

        Args:
            experiment_type: "xray" or "cryoem"
            context: Context dict from build_context()

        Returns:
            dict: {
                phase: str,           # Phase name
                description: str,     # Human-readable description
                goal: str,            # What we're trying to achieve
                reason: str,          # Why we're in this phase
            }
        """
        phases = get_workflow_phases(experiment_type)
        if not phases:
            return {"phase": "unknown", "reason": "No workflow defined"}

        if experiment_type == "xray":
            return self._detect_xray_phase(phases, context)
        elif experiment_type == "cryoem":
            return self._detect_cryoem_phase(phases, context)
        else:
            return {"phase": "unknown", "reason": "Unknown experiment type"}

    def _detect_xray_phase(self, phases, context):
        """Detect phase in X-ray workflow."""

        # Phase 1: Need analysis
        if not context["xtriage_done"]:
            return self._make_phase_result(phases, "analyze",
                "Need to analyze data quality first")

        # Phase 2b: Have prediction, need to place it
        if context["has_predicted_model"] and not context["has_placed_model"]:
            if not context["has_processed_model"]:
                return self._make_phase_result(phases, "molecular_replacement",
                    "Have prediction, need to process for MR")
            else:
                return self._make_phase_result(phases, "molecular_replacement",
                    "Model processed, need phaser")

        # Phase 2c: After autosol, need autobuild
        if context["autosol_done"] and not context["autobuild_done"] and not context["has_refined_model"]:
            return self._make_phase_result(phases, "build_from_phases",
                "Experimental phasing complete, need autobuild")

        # Phase 2d: MR-SAD - phaser placed model, anomalous data needs autosol
        # When phaser has run AND we have anomalous signal, autosol should run
        # using the phaser output as a partial model (partpdb_file) for MR-SAD.
        # Also triggered if user explicitly requested MR-SAD workflow.
        if (context["phaser_done"] and
            (context.get("has_anomalous") or context.get("use_mr_sad")) and
            not context["autosol_done"]):
            return self._make_phase_result(phases, "experimental_phasing",
                "Model placed by phaser, anomalous data detected - MR-SAD with autosol")

        # Phase 2: Need model
        if not context["has_placed_model"]:
            return self._make_phase_result(phases, "obtain_model",
                "Data analyzed, need to obtain model")

        # Phase 3b: Ligand fitted, need to combine
        if context["has_ligand_fit"] and not context["pdbtools_done"]:
            return self._make_phase_result(phases, "combine_ligand",
                "Ligand fitted, need to combine")

        # Phase 3c: Ligand combined, need refinement of model+ligand complex.
        # After ligandfit adds a ligand, the complex always needs re-refinement.
        # This takes priority over "has_refined_model" because the model has
        # changed since the last refinement.
        if context.get("needs_post_ligandfit_refine"):
            return self._make_phase_result(phases, "refine",
                "Ligand fitted, need refinement of model+ligand complex")

        # Phase 3: Has model, may need refinement
        if not context["has_refined_model"]:
            return self._make_phase_result(phases, "refine",
                "Have model, need initial refinement")

        # CRITICAL: Stay in "refine" phase if ligand fitting is possible
        # This allows ligandfit to be offered as a valid program
        # Conditions: have ligand file, r_free is good enough, not already done, AND refinement done
        # (ligandfit requires map coefficients from refined MTZ)
        if (context.get("has_ligand_file") and
            not context.get("ligandfit_done") and
            not context.get("has_ligand_fit") and
            context.get("refine_count", 0) > 0):  # Must have refined first
            # Check r_free threshold for ligandfit (< 0.35)
            r_free = context.get("r_free")
            if r_free is not None and r_free < 0.35:
                return self._make_phase_result(phases, "refine",
                    "Model refined, ligand fitting available")

        # Check if validation is needed
        validation_needed = self._needs_validation(context, "xray")
        if validation_needed and not context["validation_done"]:
            return self._make_phase_result(phases, "validate",
                "Model refined, need validation before stopping")

        # Phase 3 continued: Refinement in progress
        if not self._is_at_target(context, "xray"):
            return self._make_phase_result(phases, "refine",
                "Continuing refinement to reach target")

        # Phase 4: At target, validate or complete
        if not context["validation_done"]:
            return self._make_phase_result(phases, "validate",
                "Target reached, need validation")

        # Phase 5: Complete
        return self._make_phase_result(phases, "complete",
            "Workflow complete")

    def _detect_cryoem_phase(self, phases, context):
        """Detect phase in cryo-EM workflow."""

        # Phase 1: Need analysis
        # mtriage is preferred, but if other cryo-EM programs have already run
        # (e.g., resolve_cryo_em, dock_in_map), we're clearly past analysis.
        # This handles tutorials that skip mtriage.
        past_analysis = (
            context["mtriage_done"] or
            context.get("resolve_cryo_em_done", False) or
            context.get("dock_done", False) or
            context.get("rsr_done", False)
        )
        if not past_analysis:
            return self._make_phase_result(phases, "analyze",
                "Need to analyze map quality first")

        # Phase 1.5: Create full map from half-maps before model building.
        # When we have half-maps but no full map and resolve_cryo_em hasn't run,
        # go to optimize_map first.  Without this, phase detection falls through
        # to obtain_model and picks predict_and_build (listed first) even though
        # resolve_cryo_em is the correct next step for the tutorial workflow
        # mtriage → resolve_cryo_em → map_symmetry.
        if (context.get("has_half_map") and not context.get("has_full_map") and
                not context.get("resolve_cryo_em_done") and
                not context.get("map_sharpening_done")):
            return self._make_phase_result(phases, "optimize_map",
                "Have half-maps but no full map; creating optimized full map "
                "with resolve_cryo_em before model building")

        # Phase 2b: Stepwise - have prediction, need to dock it
        # This handles the case where predict_and_build ran with stop_after_predict=True
        if context["has_predicted_model"] and not context["has_placed_model"]:
            # For stepwise path, need to dock the predicted model
            if context["has_processed_model"]:
                # Already processed, need to dock
                return self._make_phase_result(phases, "dock_model",
                    "Have processed prediction, need to dock in map")
            elif context["predict_done"] and not context["predict_full_done"]:
                # Prediction only (stop_after_predict was used)
                # Check if we need to process first or can dock directly
                return self._make_phase_result(phases, "dock_model",
                    "Have prediction, need to dock in map")

        # Phase 2: Need model
        if not context["has_placed_model"]:
            return self._make_phase_result(phases, "obtain_model",
                "Map analyzed, need to obtain model")

        # Phase 2c: Check if map needs optimization
        if context["has_half_map"] and not context["has_full_map"]:
            return self._make_phase_result(phases, "optimize_map",
                "Have model but only half-maps, need full map for refinement")

        # Phase 3a: Have docked model but not yet refined
        if context["has_placed_model"] and not context["has_refined_model"]:
            return self._make_phase_result(phases, "ready_to_refine",
                "Model docked in map, ready for refinement")

        # Phase 3b: Refinement in progress
        if not context["has_refined_model"]:
            return self._make_phase_result(phases, "refine",
                "Have model and map, need refinement")

        # Check validation
        validation_needed = self._needs_validation(context, "cryoem")
        if validation_needed and not context["validation_done"]:
            return self._make_phase_result(phases, "validate",
                "Model refined, need validation")

        # Continue refinement if not at target
        if not self._is_at_target(context, "cryoem"):
            return self._make_phase_result(phases, "refine",
                "Continuing refinement")

        # Validate or complete
        if not context["validation_done"]:
            return self._make_phase_result(phases, "validate",
                "Target reached, need validation")

        return self._make_phase_result(phases, "complete",
            "Workflow complete")

    def _make_phase_result(self, phases, phase_name, reason):
        """Build phase result dict."""
        phase_def = phases.get(phase_name, {})
        return {
            "phase": phase_name,
            "description": phase_def.get("description", ""),
            "goal": phase_def.get("goal", ""),
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
        """Check if we've reached quality targets."""
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
                if r_free > 0.50 and refine_count >= 1:
                    return True

            # Also check clashscore - if it's good, we're likely done
            clashscore = context.get("clashscore")
            if clashscore is not None and clashscore < 10:
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

    def get_valid_programs(self, experiment_type, phase_info, context, directives=None):
        """
        Get valid programs for current phase.

        Args:
            experiment_type: "xray" or "cryoem"
            phase_info: Output from detect_phase()
            context: Context dict
            directives: Optional user directives dict

        Returns:
            list: Valid program names
        """
        phases = get_workflow_phases(experiment_type)
        phase_name = phase_info.get("phase", "")
        phase_def = phases.get(phase_name, {})

        # Handle completion phase
        if phase_def.get("stop"):
            return ["STOP"]

        # Get programs from phase definition
        phase_programs = phase_def.get("programs", [])

        valid = []
        for prog_entry in phase_programs:
            if isinstance(prog_entry, str):
                # Simple program name
                valid.append(prog_entry)
            elif isinstance(prog_entry, dict):
                # Program with conditions
                prog_name = prog_entry.get("program")
                if prog_name and self._check_conditions(prog_entry, context):
                    valid.append(prog_name)

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
                if tracking.get("strategy") == "run_once" or tracking.get("run_once"):
                    # Check if this program has already been run
                    prog_done_key = tracking.get("flag", "")
                    if prog_done_key and context.get(prog_done_key):
                        # Skip - already run
                        continue
            filtered.append(prog)
        valid = filtered

        # MR-SAD guard: When we have BOTH a search model AND anomalous data,
        # phaser must run first to place the model. AutoSol will be available
        # later in the experimental_phasing phase (after phaser completes).
        # This prevents autosol from running standalone when MR-SAD is appropriate.
        if (context and
            context.get("has_search_model") and
            (context.get("has_anomalous") or context.get("use_mr_sad")) and
            not context.get("phaser_done") and
            "phenix.autosol" in valid):
            valid.remove("phenix.autosol")

        # === APPLY USER DIRECTIVES ===
        # This is the SINGLE PLACE where directives affect valid_programs
        if directives:
            valid = self._apply_directives(valid, directives, phase_name, context)

        # Add STOP if validation done and at target
        if phase_name == "validate" and context.get("validation_done"):
            if "STOP" not in valid:
                valid.append("STOP")

        # Special: also allow refinement during validate phase (user can choose more refinement)
        # BUT respect max_refine_cycles directive AND _is_at_target (e.g., hopeless R-free)
        if phase_name == "validate":
            max_refine = directives.get("stop_conditions", {}).get("max_refine_cycles") if directives else None
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

        # Refinement phase: remove refinement when at target or max cycles reached
        # Exception: post-ligandfit refinement is always allowed (model changed)
        if phase_name == "refine":
            needs_post_ligandfit = context.get("needs_post_ligandfit_refine", False)
            max_refine = directives.get("stop_conditions", {}).get("max_refine_cycles") if directives else None
            refine_count = context.get("refine_count", 0)
            refine_allowed = (max_refine is None) or (refine_count < max_refine)
            at_target = self._is_at_target(context, experiment_type)

            if (not refine_allowed or at_target) and not needs_post_ligandfit:
                # Remove refinement programs from the phase's own list
                for prog in ["phenix.refine", "phenix.real_space_refine"]:
                    if prog in valid:
                        valid.remove(prog)

            # Add STOP if at target AND not in post-ligandfit refinement
            if at_target and not needs_post_ligandfit:
                if "STOP" not in valid:
                    valid.append("STOP")

        # NOTE: skip_validation STOP handling is now in _apply_directives()
        # No duplicate check needed here

        # If no valid programs available, return STOP (stuck state)
        if not valid:
            return ["STOP"]

        return valid

    def _apply_directives(self, valid_programs, directives, phase_name, context=None):
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
            valid_programs: List of valid program names from workflow phase
            directives: User directives dict from directive_extractor
            phase_name: Current workflow phase name
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
                # Also add STOP since we can't refine anymore
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
                should_add = self._check_program_prerequisites(prog_name, context, phase_name)
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
                # User's workflow is complete — replace the entire valid_programs
                # list with just STOP.  Simply appending STOP is not enough: the
                # phase detector may have already added programs like
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
                should_add = self._check_program_prerequisites(after_program, context, phase_name)

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
            preferred = [p for p in prefer_programs if p in result]
            others = [p for p in result if p not in prefer_programs]
            if preferred:
                result = preferred + others
                modifications.append("Prioritized: %s" % preferred)

        # Log modifications (these will show up in debug output)
        # Store modifications in a class-level list so caller can retrieve them
        self._last_directive_modifications = modifications

        return result

    def _is_program_already_done(self, program, context):
        """Check if a run_once program has already been completed.

        Returns True if the program has a run_once done_tracking strategy
        and its done flag is set in the context. This prevents directives
        from re-adding programs that the run_once filter already removed.

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
        if not (tracking.get("strategy") == "run_once" or tracking.get("run_once")):
            return False
        done_key = tracking.get("flag", "")
        if done_key and context.get(done_key):
            return True
        return False

    def _check_program_prerequisites(self, program, context, phase_name):
        """
        Check if a program's prerequisites are met.

        Args:
            program: Program name to check
            context: Workflow context dict
            phase_name: Current workflow phase name

        Returns:
            bool: True if program can be added
        """
        if not context:
            return True  # No context to check against

        # Don't add refinement programs if we're in an early phase without a model
        if program in ("phenix.refine", "phenix.real_space_refine"):
            has_model_to_refine = (
                context.get("refine_count", 0) > 0 or
                context.get("has_refined_model") or
                context.get("has_placed_model") or
                phase_name in ("refine", "validate")
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
        # Exception: always allow if we're in a phase that includes predict_and_build
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

            # Check if we're already in a phase that includes predict_and_build
            # Note: molecular_replacement is NOT included - that phase is for stepwise path only
            if phase_name in ("obtain_model", "dock_model"):
                return True  # Let the phase conditions handle it

            # For X-ray, require xtriage to be done (to get resolution for building)
            # For cryo-EM, require mtriage to be done
            if context.get("xtriage_done") or context.get("mtriage_done"):
                return True

            # Allow prediction-only if explicitly requested, but warn
            # (This case is handled by command builder forcing stop_after_predict=True)
            return False  # Don't add to early phases - let workflow proceed normally

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
        return True

    def _check_conditions(self, prog_entry, context):
        """Check if program conditions are met."""
        conditions = prog_entry.get("conditions", [])

        for cond in conditions:
            if isinstance(cond, dict):
                # Condition like {"has": "sequence"}
                if "has" in cond:
                    key = "has_" + cond["has"]
                    if not context.get(key):
                        return False

                # Condition like {"has_any": ["full_map", "optimized_full_map"]}
                # True if ANY of the listed keys is present in context
                if "has_any" in cond:
                    keys = ["has_" + k for k in cond["has_any"]]
                    if not any(context.get(k) for k in keys):
                        return False

                # Condition like {"not_done": "autobuild"}
                if "not_done" in cond:
                    key = cond["not_done"] + "_done"
                    if context.get(key):
                        return False

                # Condition like {"r_free": "> 0.35"}
                for metric in ["r_free", "map_cc", "refine_count"]:
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
        phases = get_workflow_phases(experiment_type)
        reasons = []
        found_program = False

        # Find the program in any phase
        for phase_name, phase_def in phases.items():
            if not isinstance(phase_def, dict):
                continue
            phase_programs = phase_def.get("programs", [])

            for prog_entry in phase_programs:
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

        # If program not found in any phase, it's not relevant for this experiment type
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

        if reasons:
            return "; ".join(reasons)
        return None

    def _check_metric_condition(self, context, metric, condition):
        """Check a metric condition like "> 0.35" or "< target_r_free"."""
        value = context.get(metric)
        if value is None:
            return True  # If we don't have the metric, don't block

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
                           directives=None, maximum_automation=True):
        """
        Get complete workflow state (compatible with workflow_state.py output).

        Args:
            experiment_type: "xray" or "cryoem"
            files: Categorized files dict
            history_info: Analyzed history dict
            analysis: Current log analysis
            directives: Optional user directives dict
            maximum_automation: If False, use stepwise path (process_predicted -> phaser)

        Returns:
            dict: Workflow state compatible with existing code
        """
        context = self.build_context(files, history_info, analysis, directives)

        # Add automation_path to context for program filtering
        context["automation_path"] = "automated" if maximum_automation else "stepwise"

        phase_info = self.detect_phase(experiment_type, context)
        valid_programs = self.get_valid_programs(experiment_type, phase_info, context, directives)

        # Get priority_when info for each valid program
        program_priorities = self._get_program_priorities(experiment_type, phase_info, context)

        # Map YAML phase name to original state name for compatibility
        phase_name = phase_info.get("phase", "unknown")
        state_name = self._map_phase_to_state(phase_name, experiment_type)

        # Build reason string
        reason = phase_info.get("reason", "")
        if phase_info.get("goal"):
            reason += " - Goal: " + phase_info["goal"]

        # Add metric info to reason
        if experiment_type == "xray" and context.get("r_free"):
            reason += " (R-free: %.3f)" % context["r_free"]
        elif experiment_type == "cryoem" and context.get("map_cc"):
            reason += " (CC: %.3f)" % context["map_cc"]

        # Check for stuck state (no programs available except STOP)
        if valid_programs == ["STOP"] and phase_name not in ["complete", "validate"]:
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
            "phase_info": phase_info,
            "context": context,
            "resolution": context.get("resolution"),
            "unavailable_explanations": unavailable_explanations,  # NEW: why programs aren't available
        }

    def _get_program_priorities(self, experiment_type, phase_info, context):
        """
        Get programs that should be prioritized based on priority_when conditions.

        Args:
            experiment_type: "xray" or "cryoem"
            phase_info: Phase detection result
            context: Context dict

        Returns:
            list: Programs that should be prioritized, in priority order
        """
        phases = get_workflow_phases(experiment_type)
        phase_name = phase_info.get("phase", "")
        phase_def = phases.get(phase_name, {})

        prioritized = []
        phase_programs = phase_def.get("programs", [])

        for prog_entry in phase_programs:
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
    """Convenience function to detect current phase."""
    engine = get_engine()
    context = engine.build_context(files, history_info, analysis)
    return engine.detect_phase(experiment_type, context)


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
    phase = engine.detect_phase("xray", context)
    programs = engine.get_valid_programs("xray", phase, context)

    print("X-ray initial state:")
    print("  Phase:", phase["phase"])
    print("  Reason:", phase["reason"])
    print("  Valid programs:", programs)
    print()

    # Test X-ray after xtriage
    history = {"xtriage_done": True}
    context = engine.build_context(files, history)
    phase = engine.detect_phase("xray", context)
    programs = engine.get_valid_programs("xray", phase, context)

    print("X-ray after xtriage:")
    print("  Phase:", phase["phase"])
    print("  Reason:", phase["reason"])
    print("  Valid programs:", programs)
    print()

    # Test cryo-EM initial state
    files = {"full_map": ["map.mrc"], "sequence": ["seq.fa"]}
    history = {}
    context = engine.build_context(files, history)
    phase = engine.detect_phase("cryoem", context)
    programs = engine.get_valid_programs("cryoem", phase, context)

    print("Cryo-EM initial state:")
    print("  Phase:", phase["phase"])
    print("  Reason:", phase["reason"])
    print("  Valid programs:", programs)
    print()

    # Test target retrieval
    print("Quality targets:")
    print("  R-free at 2.0A:", engine.get_target("xray", "r_free", 2.0))
    print("  Map CC:", engine.get_target("cryoem", "map_cc"))
