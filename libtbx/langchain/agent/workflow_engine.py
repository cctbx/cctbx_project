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
    from agent.workflow_engine import WorkflowEngine

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

    def build_context(self, files, history_info, analysis=None):
        """
        Build a context dict from files, history, and analysis.

        This context is used for evaluating conditions in workflows.

        Args:
            files: Categorized files dict
            history_info: Analyzed history dict
            analysis: Current log analysis dict

        Returns:
            dict: Context for condition evaluation
        """
        context = {
            # File availability
            "has_mtz": bool(files.get("mtz")),
            "has_sequence": bool(files.get("sequence")),
            "has_model": bool(files.get("pdb")),
            "has_full_map": bool(files.get("full_map")),
            "has_half_map": len(files.get("half_map", [])) >= 2,
            "has_map": bool(files.get("map") or files.get("full_map")),
            "has_ligand": bool(files.get("ligand_cif") or files.get("ligand_pdb")),
            "has_search_model": self._has_search_model(files),
            "has_predicted_model": bool(files.get("predicted")) or history_info.get("predict_done", False),
            "has_processed_model": bool(files.get("processed_predicted")) or history_info.get("process_predicted_done", False),
            "has_placed_model": self._has_placed_model(files, history_info),
            "has_refined_model": self._has_refined_model(files, history_info),
            "has_ligand_fit": bool(files.get("ligand_fit")) or history_info.get("ligandfit_done", False),
            "has_optimized_full_map": self._has_optimized_map(files, history_info),

            # Program completion
            "xtriage_done": history_info.get("xtriage_done", False),
            "mtriage_done": history_info.get("mtriage_done", False),
            "phaser_done": history_info.get("phaser_done", False),
            "predict_done": history_info.get("predict_done", False),
            "predict_full_done": history_info.get("predict_full_done", False),
            "autobuild_done": history_info.get("autobuild_done", False),
            "autosol_done": history_info.get("autosol_done", False),
            "refine_done": history_info.get("refine_done", False),
            "rsr_done": history_info.get("rsr_done", False),
            "validation_done": history_info.get("validation_done", False),
            "ligandfit_done": history_info.get("ligandfit_done", False),
            "pdbtools_done": history_info.get("pdbtools_done", False),
            "dock_done": history_info.get("dock_done", False),
            "process_predicted_model_done": history_info.get("process_predicted_done", False),

            # Counts
            "refine_count": history_info.get("refine_count", 0),
            "rsr_count": history_info.get("rsr_count", 0),

            # Metrics
            "r_free": self._get_metric(analysis, history_info, "r_free", "last_r_free"),
            "r_work": self._get_metric(analysis, history_info, "r_work", "last_r_work"),
            "map_cc": self._get_metric(analysis, history_info, "map_cc", "last_map_cc"),
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
        }

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
        """Check if there's a search model (not from refinement/prediction)."""
        for f in files.get("pdb", []):
            basename = os.path.basename(f).lower()
            if not any(x in basename for x in ['refine', 'ligand', 'phaser', 'autobuild', 'predicted', 'processed']):
                return True
        return False

    def _has_placed_model(self, files, history_info):
        """Check if model is placed (after MR/building)."""
        return bool(
            files.get("phaser_output") or
            files.get("refined") or
            files.get("autobuild_output") or
            files.get("docked") or
            history_info.get("phaser_done") or
            history_info.get("autobuild_done") or
            history_info.get("dock_done") or
            history_info.get("predict_full_done")
        )

    def _has_refined_model(self, files, history_info):
        """Check if model has been refined."""
        return bool(
            files.get("refined") or
            files.get("rsr_output") or
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

        # Phase 2: Need model
        if not context["has_placed_model"]:
            return self._make_phase_result(phases, "obtain_model",
                "Data analyzed, need to obtain model")

        # Phase 3b: Ligand fitted, need to combine
        if context["has_ligand_fit"] and not context["pdbtools_done"]:
            return self._make_phase_result(phases, "combine_ligand",
                "Ligand fitted, need to combine")

        # Phase 3: Has model, may need refinement
        if not context["has_refined_model"]:
            return self._make_phase_result(phases, "refine",
                "Have model, need initial refinement")

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
        if not context["mtriage_done"]:
            return self._make_phase_result(phases, "analyze",
                "Need to analyze map quality first")

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

        # Phase 3: Refinement
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
            if r_free is not None:
                target = get_metric_threshold("r_free", "acceptable", resolution)
                if target and r_free <= target:
                    return True
        else:
            map_cc = context.get("map_cc")
            if map_cc is not None:
                target = get_metric_threshold("map_cc", "acceptable")
                if target and map_cc >= target:
                    return True
        return False

    # =========================================================================
    # PROGRAM SELECTION
    # =========================================================================

    def get_valid_programs(self, experiment_type, phase_info, context):
        """
        Get valid programs for current phase.

        Args:
            experiment_type: "xray" or "cryoem"
            phase_info: Output from detect_phase()
            context: Context dict

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

        # Add STOP if validation done and at target
        if phase_name == "validate" and context.get("validation_done"):
            valid.append("STOP")

        # Special: also allow refinement during validate phase (user can choose more refinement)
        if phase_name == "validate":
            if experiment_type == "xray":
                if "phenix.refine" not in valid:
                    valid.append("phenix.refine")
            elif experiment_type == "cryoem":
                if "phenix.real_space_refine" not in valid:
                    valid.append("phenix.real_space_refine")

        # Special: always allow refinement to continue
        if phase_name == "refine":
            if experiment_type == "xray" and "phenix.refine" not in valid:
                valid.append("phenix.refine")
            elif experiment_type == "cryoem" and "phenix.real_space_refine" not in valid:
                valid.append("phenix.real_space_refine")

            # Add STOP if at target and validated
            if context.get("validation_done") and self._is_at_target(context, experiment_type):
                if "STOP" not in valid:
                    valid.append("STOP")

        # If no valid programs available, return STOP (stuck state)
        if not valid:
            return ["STOP"]

        return valid

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

    def get_workflow_state(self, experiment_type, files, history_info, analysis=None):
        """
        Get complete workflow state (compatible with workflow_state.py output).

        Args:
            experiment_type: "xray" or "cryoem"
            files: Categorized files dict
            history_info: Analyzed history dict
            analysis: Current log analysis

        Returns:
            dict: Workflow state compatible with existing code
        """
        context = self.build_context(files, history_info, analysis)
        phase_info = self.detect_phase(experiment_type, context)
        valid_programs = self.get_valid_programs(experiment_type, phase_info, context)

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

        return {
            "state": state_name,  # Use mapped state name
            "experiment_type": experiment_type,
            "valid_programs": valid_programs,
            "program_priorities": program_priorities,  # NEW: priority_when triggers
            "reason": reason,
            "conditions": {},
            "phase_info": phase_info,
            "context": context,
            "resolution": context.get("resolution"),
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
    files = {"mtz": ["data.mtz"], "sequence": ["seq.fa"]}
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
