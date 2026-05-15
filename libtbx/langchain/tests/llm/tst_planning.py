"""
LLM Decision Tests — Planning (D4, Phase 2)

Eight scenarios exercising the PLAN node's LLM call.  Each scenario
hand-authors a state dict that get_planning_prompt accepts, calls the
production-faithful planning chain (langchain SystemMessage +
HumanMessage), and asserts on the parsed intent.

All scenarios are reliability tests at threshold 0.8, max_runs=5.

Design philosophy:
  We assert on the LLM's RAW intent, BEFORE production's post-
  processing (validate_intent, forced_program override, etc.).  This
  isolates the LLM's decision from deterministic correction logic.
  If the LLM is unreliable on a scenario but production handles it
  correctly via post-processing, that's diagnostically useful —
  fail-count reporting surfaces this rather than hiding it.

Run from a parent directory with framework.py importable, e.g.:

    cd $PHENIX/modules/cctbx_project/libtbx/langchain
    phenix.python tests/llm/run_llm_tests.py --suite planning

See PHASE2_PLAN.md for design rationale and scenario details.
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Make the framework importable when invoked directly
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from framework import (
    Scenario, is_stop_intent, make_planning_run_fn,
    validate_planning_state)


# =====================================================================
# Assertion helpers
# =====================================================================

VALIDATION_PROGRAMS_XRAY = (
    "phenix.molprobity",
    "phenix.model_vs_data",
)

VALIDATION_PROGRAMS_CRYOEM = (
    "phenix.validation_cryoem",
    "phenix.molprobity",
    "phenix.map_correlations",
)


def _picks_validation_xray(intent):
    """Intent's program is a validation program (X-ray) AND not a STOP.

    The "AND not a STOP" guard catches contradictory outputs like
    {"program": "phenix.molprobity", "stop": true, ...} — production
    treats stop=true as the stop signal regardless of program field,
    so we should too.
    """
    if is_stop_intent(intent):
        return (False, "picked STOP (program=%r, stop=%r)"
                       % (intent.get("program"), intent.get("stop")))
    prog = (intent or {}).get("program")
    if prog in VALIDATION_PROGRAMS_XRAY:
        return (True, "validation program chosen: %s" % prog)
    return (False, "expected validation program, got %r" % prog)


def _picks_validation_cryoem(intent):
    """Intent's program is a validation program (cryo-EM) AND not a STOP."""
    if is_stop_intent(intent):
        return (False, "picked STOP (program=%r, stop=%r)"
                       % (intent.get("program"), intent.get("stop")))
    prog = (intent or {}).get("program")
    if prog in VALIDATION_PROGRAMS_CRYOEM:
        return (True, "validation program chosen: %s" % prog)
    return (False, "expected validation program, got %r" % prog)


def _picks_stop(intent):
    """Intent indicates a STOP decision (either program=='STOP' or
    stop==True).  Matches production's stop detection."""
    if is_stop_intent(intent):
        return (True, "STOP intent (program=%r, stop=%r)"
                       % (intent.get("program"), intent.get("stop")))
    return (False, "expected STOP, got program=%r stop=%r"
                   % (intent.get("program"), intent.get("stop")))


def _picks_not_program(forbidden):
    """Intent's program is NOT the forbidden program AND is not STOP.

    The "AND is not STOP" clause catches a degenerate pass mode where
    the LLM gives up entirely (picks STOP) rather than doing the
    workflow.  Picking STOP technically avoids the forbidden program
    but isn't the behavior we want to test.  If the LLM stopping is
    a legitimate response for a particular scenario, use _picks_stop
    or a custom helper instead.
    """
    def check(intent):
        prog = (intent or {}).get("program")
        if is_stop_intent(intent):
            return (False, ("picked STOP (degenerate avoid of %s)"
                            % forbidden))
        if prog == forbidden:
            return (False, "picked forbidden program: %s" % prog)
        return (True, "avoided %s (picked %r)" % (forbidden, prog))
    return check


def _picks_valid_non_stop(valid_programs):
    """Intent's program is in valid_programs AND is not STOP."""
    def check(intent):
        prog = (intent or {}).get("program")
        if is_stop_intent(intent):
            return (False, "picked STOP")
        if prog not in valid_programs:
            return (False, "picked %r, not in valid_programs %r"
                           % (prog, valid_programs))
        return (True, "valid non-STOP choice: %s" % prog)
    return check


def _picks_different_or_changes_strategy(failed_program, failed_strategy):
    """Intent picks a DIFFERENT program OR same program with DIFFERENT
    strategy.  STOP is rejected as a degenerate pass — the point of
    this assertion is that the LLM should ATTEMPT recovery, not give
    up after a single failure.
    """
    def check(intent):
        prog = (intent or {}).get("program")
        strat = (intent or {}).get("strategy") or {}
        if is_stop_intent(intent):
            return (False, "picked STOP (gave up rather than recover)")
        if prog != failed_program:
            return (True, "picked different program: %s" % prog)
        if strat != failed_strategy:
            return (True, ("same program but different strategy: %r vs %r"
                           % (strat, failed_strategy)))
        return (False, ("retried %s with identical strategy %r"
                        % (prog, strat)))
    return check


# =====================================================================
# Scenario states — hand-authored Python dicts
# =====================================================================
#
# Each state populates EVERY key in PLANNING_REQUIRED_KEYS (validated
# at scenario load time by validate_planning_state).
#
# workflow_state has the three required sub-keys: state,
# experiment_type, valid_programs.  Other sub-keys (step_info, reason,
# context, etc.) are optional but populated for realism.


# ---------------------------------------------------------------------
# Scenario 1: validate_pending_dont_stop
# ---------------------------------------------------------------------
#
# X-ray workflow, refinement reached R-free 0.22 (under 0.25 target),
# no validation has been run yet.  LLM must pick a validation program,
# NOT STOP.  Tests rule #7 of the system prompt
# ("Validate before stopping").

_STATE_VALIDATE_PENDING = {
    "history": [
        {"cycle_number": 1, "program": "phenix.xtriage",
         "command": "phenix.xtriage data.mtz",
         "result": "Twinning analysis complete; no twinning detected",
         "summary": "xtriage ok", "error": "",
         "metrics": {},
         "output_files": ["xtriage.log"],
         "decision": "continue"},
        {"cycle_number": 2, "program": "phenix.refine",
         "command": "phenix.refine data.mtz init.pdb",
         "result": "R-free 0.31, R-work 0.26",
         "summary": "refine_001", "error": "",
         "metrics": {"r_free": 0.31, "r_work": 0.26},
         "output_files": ["refine_001.pdb", "refine_001.mtz"],
         "decision": "continue"},
        {"cycle_number": 3, "program": "phenix.refine",
         "command": "phenix.refine refine_001.mtz refine_001.pdb",
         "result": "R-free 0.24, R-work 0.20",
         "summary": "refine_002", "error": "",
         "metrics": {"r_free": 0.24, "r_work": 0.20},
         "output_files": ["refine_002.pdb", "refine_002.mtz"],
         "decision": "continue"},
        {"cycle_number": 4, "program": "phenix.refine",
         "command": "phenix.refine refine_002.mtz refine_002.pdb",
         "result": "R-free 0.22, R-work 0.18",
         "summary": "refine_003", "error": "",
         "metrics": {"r_free": 0.22, "r_work": 0.18},
         "output_files": ["refine_003.pdb", "refine_003.mtz"],
         "decision": "continue"},
    ],
    "analysis": {
        "program": "phenix.refine",
        "metrics": {"r_free": 0.22, "r_work": 0.18, "rms_bonds": 0.011},
        "output_files": ["refine_003.pdb", "refine_003.mtz"],
        "resolution": 2.1,
        "errors": [],
        "warnings": [],
    },
    "available_files": [
        "data.mtz", "init.pdb", "sequence.fa",
        "refine_001.pdb", "refine_001.mtz",
        "refine_002.pdb", "refine_002.mtz",
        "refine_003.pdb", "refine_003.mtz",
        "xtriage.log",
    ],
    "previous_attempts": [],
    "user_advice": "Refine the model.",
    "metrics_trend": {
        "should_stop": False,
        "reason": "",
        "trend_summary": "R-free improving steadily: 0.31 → 0.24 → 0.22",
        "recommendation": "continue",
        "consecutive_refines": 3,
        "consecutive_rsr": 0,
    },
    "workflow_state": {
        "state": "xray_refined",
        "experiment_type": "xray",
        "step_info": {"step": "validate",
                      "description": "Refinement converged; "
                                     "validation pending"},
        "valid_programs": ["phenix.molprobity",
                            "phenix.model_vs_data",
                            "STOP"],
        "reason": ("R-free 0.22 below target 0.25, but no validation "
                   "has been run yet.  Per rule #7, validate before "
                   "stopping."),
        "context": {},
        "categorized_files": {
            "model": ["refine_003.pdb"],
            "data": ["refine_003.mtz"],
        },
        "automation_path": "stepwise",
    },
    "directives": {},
    "best_files": {
        "model": "refine_003.pdb",
        "data_mtz": "refine_003.mtz",
        "map": None,
        "map_coeffs_mtz": None,
    },
}


# ---------------------------------------------------------------------
# Scenario 2: plateau_with_validation_done
# ---------------------------------------------------------------------
#
# X-ray workflow, metrics plateaued (should_stop=True), molprobity
# already in history.  LLM must pick STOP.

_STATE_PLATEAU_VAL_DONE = {
    "history": [
        {"cycle_number": 1, "program": "phenix.xtriage",
         "result": "ok", "summary": "xtriage", "error": "",
         "metrics": {}, "output_files": ["xtriage.log"],
         "decision": "continue"},
        {"cycle_number": 2, "program": "phenix.refine",
         "result": "R-free 0.28", "summary": "refine_001", "error": "",
         "metrics": {"r_free": 0.28, "r_work": 0.23},
         "output_files": ["refine_001.pdb", "refine_001.mtz"],
         "decision": "continue"},
        {"cycle_number": 3, "program": "phenix.refine",
         "result": "R-free 0.235", "summary": "refine_002", "error": "",
         "metrics": {"r_free": 0.235, "r_work": 0.19},
         "output_files": ["refine_002.pdb", "refine_002.mtz"],
         "decision": "continue"},
        {"cycle_number": 4, "program": "phenix.molprobity",
         "result": "Clashscore 4.2, all-atom 95th percentile",
         "summary": "molprobity_ok", "error": "",
         "metrics": {"clashscore": 4.2, "all_atom_percentile": 95},
         "output_files": ["molprobity.log"],
         "decision": "continue"},
        {"cycle_number": 5, "program": "phenix.refine",
         "result": "R-free 0.232", "summary": "refine_003", "error": "",
         "metrics": {"r_free": 0.232, "r_work": 0.188},
         "output_files": ["refine_003.pdb", "refine_003.mtz"],
         "decision": "continue"},
        {"cycle_number": 6, "program": "phenix.refine",
         "result": "R-free 0.231", "summary": "refine_004", "error": "",
         "metrics": {"r_free": 0.231, "r_work": 0.187},
         "output_files": ["refine_004.pdb", "refine_004.mtz"],
         "decision": "continue"},
    ],
    "analysis": {
        "program": "phenix.refine",
        "metrics": {"r_free": 0.231, "r_work": 0.187, "rms_bonds": 0.011},
        "output_files": ["refine_004.pdb", "refine_004.mtz"],
        "resolution": 2.1,
        "errors": [],
        "warnings": [],
    },
    "available_files": [
        "data.mtz", "init.pdb",
        "refine_001.pdb", "refine_001.mtz",
        "refine_002.pdb", "refine_002.mtz",
        "refine_003.pdb", "refine_003.mtz",
        "refine_004.pdb", "refine_004.mtz",
        "molprobity.log", "xtriage.log",
    ],
    "previous_attempts": [],
    "user_advice": "Refine and validate the structure.",
    "metrics_trend": {
        "should_stop": True,
        "reason": ("Metrics plateaued for 3+ cycles "
                   "(<0.3% improvement)"),
        "trend_summary": ("R-free trend: 0.235 → 0.232 → 0.231 "
                          "(plateau detected)"),
        "recommendation": "consider_stopping",
        "consecutive_refines": 2,
        "consecutive_rsr": 0,
    },
    "workflow_state": {
        "state": "xray_refined",
        "experiment_type": "xray",
        "step_info": {"step": "complete",
                      "description": "Refinement plateaued; "
                                     "validation complete"},
        "valid_programs": ["phenix.refine", "STOP"],
        "reason": ("R-free 0.231 below target 0.25; molprobity already "
                   "run; metrics plateaued."),
        "context": {},
        "categorized_files": {
            "model": ["refine_004.pdb"],
            "data": ["refine_004.mtz"],
        },
        "automation_path": "stepwise",
    },
    "directives": {},
    "best_files": {
        "model": "refine_004.pdb",
        "data_mtz": "refine_004.mtz",
        "map": None,
        "map_coeffs_mtz": None,
    },
}


# ---------------------------------------------------------------------
# Scenario 3: dont_stop_on_plateau_when_validation_pending
# ---------------------------------------------------------------------
#
# Same plateau signal as scenario 2, but molprobity has NOT been run.
# Validate-before-stopping (rule #7) must win over the plateau signal.
# LLM must pick a validation program, NOT STOP.

_STATE_PLATEAU_VAL_PENDING = {
    "history": [
        {"cycle_number": 1, "program": "phenix.xtriage",
         "result": "ok", "summary": "xtriage", "error": "",
         "metrics": {}, "output_files": ["xtriage.log"],
         "decision": "continue"},
        {"cycle_number": 2, "program": "phenix.refine",
         "result": "R-free 0.28", "summary": "refine_001", "error": "",
         "metrics": {"r_free": 0.28, "r_work": 0.23},
         "output_files": ["refine_001.pdb", "refine_001.mtz"],
         "decision": "continue"},
        {"cycle_number": 3, "program": "phenix.refine",
         "result": "R-free 0.235", "summary": "refine_002", "error": "",
         "metrics": {"r_free": 0.235, "r_work": 0.19},
         "output_files": ["refine_002.pdb", "refine_002.mtz"],
         "decision": "continue"},
        {"cycle_number": 4, "program": "phenix.refine",
         "result": "R-free 0.232", "summary": "refine_003", "error": "",
         "metrics": {"r_free": 0.232, "r_work": 0.188},
         "output_files": ["refine_003.pdb", "refine_003.mtz"],
         "decision": "continue"},
        {"cycle_number": 5, "program": "phenix.refine",
         "result": "R-free 0.231", "summary": "refine_004", "error": "",
         "metrics": {"r_free": 0.231, "r_work": 0.187},
         "output_files": ["refine_004.pdb", "refine_004.mtz"],
         "decision": "continue"},
    ],
    "analysis": {
        "program": "phenix.refine",
        "metrics": {"r_free": 0.231, "r_work": 0.187, "rms_bonds": 0.011},
        "output_files": ["refine_004.pdb", "refine_004.mtz"],
        "resolution": 2.1,
        "errors": [],
        "warnings": [],
    },
    "available_files": [
        "data.mtz", "init.pdb",
        "refine_001.pdb", "refine_001.mtz",
        "refine_002.pdb", "refine_002.mtz",
        "refine_003.pdb", "refine_003.mtz",
        "refine_004.pdb", "refine_004.mtz",
        "xtriage.log",
    ],
    "previous_attempts": [],
    "user_advice": "Refine the structure to convergence.",
    "metrics_trend": {
        "should_stop": True,
        "reason": ("Metrics plateaued for 3+ cycles "
                   "(<0.3% improvement)"),
        "trend_summary": ("R-free trend: 0.235 → 0.232 → 0.231 "
                          "(plateau detected)"),
        "recommendation": "consider_stopping",
        "consecutive_refines": 3,
        "consecutive_rsr": 0,
    },
    "workflow_state": {
        "state": "xray_refined",
        "experiment_type": "xray",
        "step_info": {"step": "validate",
                      "description": ("Refinement plateaued; "
                                      "no validation has been run yet")},
        "valid_programs": ["phenix.molprobity",
                            "phenix.model_vs_data",
                            "STOP"],
        "reason": ("R-free below target and plateaued, BUT no "
                   "validation has been run.  Per rule #7, validate "
                   "before stopping."),
        "context": {},
        "categorized_files": {
            "model": ["refine_004.pdb"],
            "data": ["refine_004.mtz"],
        },
        "automation_path": "stepwise",
    },
    "directives": {},
    "best_files": {
        "model": "refine_004.pdb",
        "data_mtz": "refine_004.mtz",
        "map": None,
        "map_coeffs_mtz": None,
    },
}


# ---------------------------------------------------------------------
# Scenario 4: cryoem_stop_target_with_validation
# ---------------------------------------------------------------------
#
# Cryo-EM workflow, map_cc 0.72 (above 0.70 target),
# validation_cryoem already run.  LLM must pick STOP.

_STATE_CRYOEM_STOP_TARGET = {
    "history": [
        {"cycle_number": 1, "program": "phenix.predict_and_build",
         "result": "Model placed; map_cc 0.42",
         "summary": "predict_and_build", "error": "",
         "metrics": {"map_cc": 0.42},
         "output_files": ["predicted_model.pdb"],
         "decision": "continue"},
        {"cycle_number": 2, "program": "phenix.resolve_cryo_em",
         "result": "Density modification; CC 0.55",
         "summary": "denmod", "error": "",
         "metrics": {"map_cc": 0.55},
         "output_files": ["denmod_map.ccp4"],
         "decision": "continue"},
        {"cycle_number": 3, "program": "phenix.real_space_refine",
         "result": "Refined; CC 0.68",
         "summary": "rsr_001", "error": "",
         "metrics": {"map_cc": 0.68},
         "output_files": ["rsr_001.pdb"],
         "decision": "continue"},
        {"cycle_number": 4, "program": "phenix.validation_cryoem",
         "result": ("Geometry score 0.93, model-map FSC at 2.86A; "
                    "all validations pass"),
         "summary": "validation_cryoem", "error": "",
         "metrics": {"geometry_score": 0.93},
         "output_files": ["validation_cryoem.log"],
         "decision": "continue"},
        {"cycle_number": 5, "program": "phenix.real_space_refine",
         "result": "Refined; CC 0.72",
         "summary": "rsr_002", "error": "",
         "metrics": {"map_cc": 0.72},
         "output_files": ["rsr_002.pdb"],
         "decision": "continue"},
    ],
    "analysis": {
        "program": "phenix.real_space_refine",
        "metrics": {"map_cc": 0.72, "geometry_score": 0.93},
        "output_files": ["rsr_002.pdb"],
        "resolution": 2.86,
        "errors": [],
        "warnings": [],
    },
    "available_files": [
        "sequence.fa", "half_map_1.ccp4", "half_map_2.ccp4",
        "full_map.ccp4",
        "predicted_model.pdb", "denmod_map.ccp4",
        "rsr_001.pdb", "rsr_002.pdb",
        "validation_cryoem.log",
    ],
    "previous_attempts": [],
    "user_advice": ("Build a model into the cryo-EM map and refine "
                    "until validated."),
    "metrics_trend": {
        "should_stop": True,
        "reason": ("Map-model CC 0.72 exceeds cryo-EM target 0.70; "
                   "validation passed."),
        "trend_summary": "CC progression: 0.42 → 0.55 → 0.68 → 0.72",
        "recommendation": "consider_stopping",
        "consecutive_refines": 0,
        "consecutive_rsr": 2,
    },
    "workflow_state": {
        "state": "cryoem_refined",
        "experiment_type": "cryoem",
        "step_info": {"step": "complete",
                      "description": ("Map-model CC > 0.70 and "
                                      "validation passed")},
        "valid_programs": ["phenix.real_space_refine", "STOP"],
        "reason": ("Map-model CC 0.72 above 0.70 target AND "
                   "validation_cryoem passed."),
        "context": {},
        "categorized_files": {
            "model": ["rsr_002.pdb"],
            "map": ["denmod_map.ccp4"],
            "half_maps": ["half_map_1.ccp4", "half_map_2.ccp4"],
        },
        "automation_path": "stepwise",
    },
    "directives": {},
    "best_files": {
        "model": "rsr_002.pdb",
        "data_mtz": None,
        "map": "denmod_map.ccp4",
        "map_coeffs_mtz": None,
    },
}


# ---------------------------------------------------------------------
# Scenario 5: directives_force_stop_after_phaser
# ---------------------------------------------------------------------
#
# X-ray, phaser just ran, directives.stop_conditions.after_program =
# phenix.phaser.  LLM must pick STOP (honor the directive without
# help from post-processing).

_STATE_DIRECTIVES_STOP_AFTER_PHASER = {
    "history": [
        {"cycle_number": 1, "program": "phenix.phaser",
         "command": "phenix.phaser data.mtz model.pdb",
         "result": ("Molecular replacement: LLG 285, TFZ 12.3; "
                    "solution found"),
         "summary": "phaser_ok", "error": "",
         "metrics": {"phaser_llg": 285, "phaser_tfz": 12.3},
         "output_files": ["PHASER.1.pdb", "PHASER.1.mtz"],
         "decision": "continue"},
    ],
    "analysis": {
        "program": "phenix.phaser",
        "metrics": {"phaser_llg": 285, "phaser_tfz": 12.3},
        "output_files": ["PHASER.1.pdb", "PHASER.1.mtz"],
        "resolution": 2.5,
        "errors": [],
        "warnings": [],
    },
    "available_files": [
        "data.mtz", "model.pdb", "sequence.fa",
        "PHASER.1.pdb", "PHASER.1.mtz",
    ],
    "previous_attempts": [],
    "user_advice": ("Run phaser for molecular replacement.  Stop "
                    "after phaser completes; I'll look at the "
                    "output before continuing."),
    "metrics_trend": {
        "should_stop": False,
        "reason": "",
        "trend_summary": "",
        "recommendation": "continue",
        "consecutive_refines": 0,
        "consecutive_rsr": 0,
    },
    "workflow_state": {
        "state": "xray_phased",
        "experiment_type": "xray",
        "step_info": {"step": "user_stop",
                      "description": ("User requested stop after "
                                      "phaser")},
        "valid_programs": ["phenix.refine", "phenix.xtriage", "STOP"],
        "reason": ("User directive: stop_conditions.after_program = "
                   "phenix.phaser.  Phaser just ran; pick STOP."),
        "context": {},
        "categorized_files": {
            "model": ["PHASER.1.pdb"],
            "data": ["PHASER.1.mtz"],
        },
        "automation_path": "stepwise",
    },
    "directives": {
        "stop_conditions": {"after_program": "phenix.phaser"},
    },
    "best_files": {
        "model": "PHASER.1.pdb",
        "data_mtz": "PHASER.1.mtz",
        "map": None,
        "map_coeffs_mtz": None,
    },
}


# ---------------------------------------------------------------------
# Scenario 6: directives_skip_program_honored
# ---------------------------------------------------------------------
#
# X-ray workflow start (cycle 0).  directives.workflow_preferences.
# skip_programs = [phenix.mtriage].  valid_programs includes
# phenix.mtriage.  LLM must NOT pick mtriage.

_STATE_DIRECTIVES_SKIP_MTRIAGE = {
    "history": [],
    "analysis": {
        "program": None,
        "metrics": {},
        "output_files": [],
        "resolution": None,
        "errors": [],
        "warnings": [],
    },
    "available_files": ["data.mtz", "model.pdb", "sequence.fa"],
    "previous_attempts": [],
    "user_advice": ("Run refinement on this model.  Don't run "
                    "mtriage - map analysis is not needed."),
    "metrics_trend": {
        "should_stop": False,
        "reason": "",
        "trend_summary": "Workflow not yet started",
        "recommendation": "continue",
        "consecutive_refines": 0,
        "consecutive_rsr": 0,
    },
    "workflow_state": {
        "state": "xray_start",
        "experiment_type": "xray",
        "step_info": {"step": "begin",
                      "description": "Start of X-ray workflow"},
        "valid_programs": ["phenix.xtriage",
                            "phenix.refine",
                            "phenix.mtriage",
                            "STOP"],
        "reason": "Initial workflow start; data + model present.",
        "context": {},
        "categorized_files": {
            "model": ["model.pdb"],
            "data": ["data.mtz"],
            "sequence": ["sequence.fa"],
        },
        "automation_path": "stepwise",
    },
    "directives": {
        "workflow_preferences": {
            "skip_programs": ["phenix.mtriage"],
        },
    },
    "best_files": {
        "model": "model.pdb",
        "data_mtz": "data.mtz",
        "map": None,
        "map_coeffs_mtz": None,
    },
}


# ---------------------------------------------------------------------
# Scenario 7: cryoem_next_step_after_placement
# ---------------------------------------------------------------------
#
# Cryo-EM, just-finished predict_and_build, map_cc 0.42 (below 0.70
# target).  LLM must pick SOMETHING in valid_programs that is NOT
# STOP.  Capability-style check at reliability threshold.

_STATE_CRYOEM_AFTER_PLACEMENT = {
    "history": [
        {"cycle_number": 1, "program": "phenix.predict_and_build",
         "command": ("phenix.predict_and_build "
                     "seq_file=sequence.fa "
                     "map_file=half_map_1.ccp4 "
                     "alphafold.use_templates=False"),
         "result": "Model placed; map_cc 0.42",
         "summary": "predict_and_build", "error": "",
         "metrics": {"map_cc": 0.42, "cc_box": 0.45},
         "output_files": ["predicted_model.pdb"],
         "decision": "continue"},
    ],
    "analysis": {
        "program": "phenix.predict_and_build",
        "metrics": {"map_cc": 0.42, "cc_box": 0.45},
        "output_files": ["predicted_model.pdb"],
        "resolution": 2.86,
        "errors": [],
        "warnings": [],
    },
    "available_files": [
        "sequence.fa", "half_map_1.ccp4", "half_map_2.ccp4",
        "full_map.ccp4", "predicted_model.pdb",
    ],
    "previous_attempts": [],
    "user_advice": ("Run PredictAndBuild for the cryo-EM tutorial, "
                    "rebuild missing loops, and refine."),
    "metrics_trend": {
        "should_stop": False,
        "reason": "",
        "trend_summary": ("Initial placement; map_cc 0.42 below "
                          "0.70 target"),
        "recommendation": "continue",
        "consecutive_refines": 0,
        "consecutive_rsr": 0,
    },
    "workflow_state": {
        "state": "cryoem_model_placed",
        "experiment_type": "cryoem",
        "step_info": {"step": "rebuild",
                      "description": ("Predicted model placed; "
                                      "needs rebuilding and "
                                      "refinement")},
        "valid_programs": ["phenix.resolve_cryo_em",
                            "phenix.real_space_refine",
                            "phenix.validation_cryoem",
                            "STOP"],
        "reason": ("After predict_and_build placement.  map_cc 0.42 "
                   "below 0.70 target — need to improve via density "
                   "modification, rebuilding, or refinement."),
        "context": {},
        "categorized_files": {
            "model": ["predicted_model.pdb"],
            "map": ["full_map.ccp4"],
            "half_maps": ["half_map_1.ccp4", "half_map_2.ccp4"],
            "sequence": ["sequence.fa"],
        },
        "automation_path": "stepwise",
    },
    "directives": {},
    "best_files": {
        "model": "predicted_model.pdb",
        "data_mtz": None,
        "map": "full_map.ccp4",
        "map_coeffs_mtz": None,
    },
}


# ---------------------------------------------------------------------
# Scenario 8: dont_retry_same_failure
# ---------------------------------------------------------------------
#
# X-ray, last cycle of refine FAILED with clash score 80.
# previous_attempts has the failure.  LLM must pick a DIFFERENT
# program OR same program with DIFFERENT strategy.

_DONT_RETRY_FAILED_STRATEGY = {"strategy": "individual_sites"}

_STATE_DONT_RETRY = {
    "history": [
        {"cycle_number": 1, "program": "phenix.xtriage",
         "result": "ok", "summary": "xtriage", "error": "",
         "metrics": {}, "output_files": ["xtriage.log"],
         "decision": "continue"},
        {"cycle_number": 2, "program": "phenix.refine",
         "result": "R-free 0.31", "summary": "refine_001", "error": "",
         "metrics": {"r_free": 0.31, "r_work": 0.26, "clashscore": 18.5},
         "output_files": ["refine_001.pdb", "refine_001.mtz"],
         "decision": "continue"},
        {"cycle_number": 3, "program": "phenix.refine",
         "command": ("phenix.refine refine_001.mtz refine_001.pdb "
                     "strategy=individual_sites"),
         "result": ("FAILED: clash score 80 — model contains many "
                    "atom-atom clashes after individual_sites "
                    "refinement"),
         "summary": "refine_failed",
         "error": "Clash score 80 (target < 10)",
         "metrics": {"clashscore": 80, "r_free": 0.40},
         "output_files": ["refine_002_failed.log"],
         "decision": "recover"},
    ],
    "analysis": {
        "program": "phenix.refine",
        "metrics": {"clashscore": 80, "r_free": 0.40},
        "output_files": ["refine_002_failed.log"],
        "resolution": 2.5,
        "errors": ["Clash score 80 exceeds threshold 10"],
        "warnings": [],
    },
    "available_files": [
        "data.mtz", "init.pdb", "xtriage.log",
        "refine_001.pdb", "refine_001.mtz",
        "refine_002_failed.log",
    ],
    "previous_attempts": [
        # Schema matches production's _increment_attempt:
        # {intent, command, error}.  We also keep program/strategy for
        # the assertion helper to compare against.
        {"intent": {"program": "phenix.refine",
                    "strategy": _DONT_RETRY_FAILED_STRATEGY},
         "command": ("phenix.refine refine_001.mtz refine_001.pdb "
                     "strategy=individual_sites"),
         "error": ("Clash score 80 — model contains many clashes "
                   "after individual_sites refinement"),
         "program": "phenix.refine",
         "strategy": _DONT_RETRY_FAILED_STRATEGY,
         "cycle_number": 3},
    ],
    "user_advice": "Refine the model.",
    "metrics_trend": {
        "should_stop": False,
        "reason": "",
        "trend_summary": ("Refine just failed with high clash score; "
                          "recovery needed"),
        "recommendation": "change_strategy",
        "consecutive_refines": 0,
        "consecutive_rsr": 0,
    },
    "workflow_state": {
        "state": "xray_recovery",
        "experiment_type": "xray",
        "step_info": {"step": "recover",
                      "description": ("Refine failed with high "
                                      "clashes; try different "
                                      "program or strategy")},
        "valid_programs": ["phenix.refine",
                            "phenix.real_space_refine",
                            "phenix.geometry_minimization",
                            "STOP"],
        "reason": ("phenix.refine with strategy=individual_sites "
                   "produced unacceptable clashes.  Try a different "
                   "program OR same program with different strategy."),
        "context": {"last_failure": "high_clash_score"},
        "categorized_files": {
            "model": ["refine_001.pdb"],
            "data": ["refine_001.mtz"],
        },
        "automation_path": "stepwise",
    },
    "directives": {},
    "best_files": {
        "model": "refine_001.pdb",
        "data_mtz": "refine_001.mtz",
        "map": None,
        "map_coeffs_mtz": None,
    },
}


# =====================================================================
# Scenarios
# =====================================================================

def build_scenarios():
    """Build the list of D4 planning scenarios.

    Each state is validated against PLANNING_REQUIRED_KEYS at load
    time; a malformed state raises ValueError here, not silently at
    test runtime.
    """
    scenarios = [
        Scenario(
            name="validate_pending_dont_stop",
            description=(
                "X-ray refinement reached R-free 0.22, no validation "
                "run yet — LLM should pick a validation program, "
                "NOT STOP.  Tests rule #7 'validate before stopping'."),
            decision_point="planning",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_STATE_VALIDATE_PENDING,
            expected_fn=_picks_validation_xray,
        ),
        Scenario(
            name="plateau_with_validation_done",
            description=(
                "X-ray metrics plateaued (should_stop=True), "
                "molprobity already run — LLM should pick STOP."),
            decision_point="planning",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_STATE_PLATEAU_VAL_DONE,
            expected_fn=_picks_stop,
        ),
        Scenario(
            name="dont_stop_on_plateau_when_validation_pending",
            description=(
                "X-ray plateau detected BUT no validation run yet — "
                "validate-first must win over plateau signal.  LLM "
                "should pick a validation program, NOT STOP."),
            decision_point="planning",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_STATE_PLATEAU_VAL_PENDING,
            expected_fn=_picks_validation_xray,
        ),
        Scenario(
            name="cryoem_stop_target_with_validation",
            description=(
                "Cryo-EM map_cc 0.72 above 0.70 target, "
                "validation_cryoem already run, metrics plateaued — "
                "LLM should pick STOP."),
            decision_point="planning",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_STATE_CRYOEM_STOP_TARGET,
            expected_fn=_picks_stop,
        ),
        Scenario(
            name="directives_force_stop_after_phaser",
            description=(
                "Phaser just ran, "
                "directives.stop_conditions.after_program = "
                "phenix.phaser — LLM should pick STOP without help "
                "from production's forced_program override."),
            decision_point="planning",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_STATE_DIRECTIVES_STOP_AFTER_PHASER,
            expected_fn=_picks_stop,
        ),
        Scenario(
            name="directives_skip_program_honored",
            description=(
                "directives.workflow_preferences.skip_programs = "
                "[phenix.mtriage], mtriage in valid_programs — LLM "
                "should NOT pick mtriage."),
            decision_point="planning",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_STATE_DIRECTIVES_SKIP_MTRIAGE,
            expected_fn=_picks_not_program("phenix.mtriage"),
        ),
        Scenario(
            name="cryoem_next_step_after_placement",
            description=(
                "Cryo-EM after predict_and_build placement, map_cc "
                "0.42 below target — LLM should pick any valid "
                "non-STOP next step."),
            decision_point="planning",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_STATE_CRYOEM_AFTER_PLACEMENT,
            expected_fn=_picks_valid_non_stop(
                _STATE_CRYOEM_AFTER_PLACEMENT["workflow_state"][
                    "valid_programs"]),
        ),
        Scenario(
            name="dont_retry_same_failure",
            description=(
                "phenix.refine with strategy=individual_sites just "
                "failed with clash score 80 — LLM should pick a "
                "different program OR same program with different "
                "strategy."),
            decision_point="planning",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_STATE_DONT_RETRY,
            expected_fn=_picks_different_or_changes_strategy(
                "phenix.refine", _DONT_RETRY_FAILED_STRATEGY),
        ),
    ]

    # Validate every scenario's state at load time
    for sc in scenarios:
        try:
            validate_planning_state(sc.input)
        except ValueError as e:
            raise ValueError(
                "Scenario %r has malformed state: %s" % (sc.name, e))

    return scenarios


# =====================================================================
# Standalone runner — useful for debugging this file directly
# =====================================================================

def main(argv=None):
    """Run all D4 scenarios against all available providers."""
    if argv is None:
        argv = sys.argv[1:]

    from framework import (
        AVAILABLE_PROVIDERS, build_summary,
        create_log_dir, run_scenario_against_providers,
        write_run_summary)
    import time

    if not AVAILABLE_PROVIDERS:
        sys.stderr.write(
            "No LLM provider API keys configured.  Set at least one of\n"
            "  GOOGLE_API_KEY, OPENAI_API_KEY, ANTHROPIC_API_KEY\n")
        return 2

    providers = list(AVAILABLE_PROVIDERS)
    scenarios = build_scenarios()
    run_one = make_planning_run_fn()
    log_dir = create_log_dir(phase="phase2-d4")

    print()
    print("LLM Decision Tests — Planning (Phase 2, D4)")
    print("=" * 65)
    print("Available providers: %s" % ", ".join(providers))
    skipped = [p for p in ("google", "openai", "anthropic")
               if p not in providers]
    if skipped:
        print("Skipped providers:   %s (no API key)" % ", ".join(skipped))
    print("Log directory:       %s" % log_dir)
    print()

    all_verdicts = []
    t0 = time.time()
    for scenario in scenarios:
        verdicts = run_scenario_against_providers(
            scenario, providers, run_one, log_dir=log_dir)
        all_verdicts.extend(verdicts)
    wall_time = time.time() - t0

    summary = build_summary(
        all_verdicts, providers_tested=providers,
        providers_skipped=skipped, wall_time_s=wall_time,
        phase="phase2-d4")
    write_run_summary(log_dir, summary)

    print()
    print("Summary: %d PASS, %d FAIL, %d ERROR, %d SKIP"
          % (summary["totals"]["passed"], summary["totals"]["failed"],
             summary["totals"]["errored"], summary["totals"]["skipped"]))
    print("Total LLM calls: %d" % summary["totals"]["total_llm_calls"])
    print("Wall time:       %.1fs" % wall_time)
    print("Log:             %s" % log_dir)

    return 0 if (summary["totals"]["failed"] == 0
                 and summary["totals"]["errored"] == 0) else 1


if __name__ == "__main__":
    sys.exit(main())
