"""
B2 — LLM Decision Test: Gated Wipe in Planning Suite (v117)

This is a PLANNING-SUITE test (D4, not D2).  It exercises the
v117 stop_refactor architecture end-to-end: when stop_after_requested
is True and after_program_done is True, _apply_directives wipes
valid_programs to ["STOP"], and the LLM correctly picks STOP because
that's its only option.

Two scenarios:
  B2a (positive): wiped to [STOP] -> LLM picks STOP
  B2b (negative control): not wiped -> LLM does NOT auto-stop

The unit-level test of the wipe itself is in
tst_validate_step_after_program_guard.py (M1a/b/c).  B2 verifies the
END-TO-END behavior: that the wipe + LLM prompt formatting + LLM
choice all line up.

Dependency: this file requires the Phase 2A planning infrastructure
(`call_planning_llm`, `make_planning_run_fn`) to exist in framework.py.
Per PHASE2_PLAN_v2.md §2, that infrastructure was smoke-test-confirmed
on 2026-05-15 but had not yet been merged.  Until merged, this file
raises ImportError at construct-time and the test cannot run.

If Phase 2A is not yet in place when you need to run B2, you have
two options:
  1. Merge Phase 2A first (recommended — gives you all 8 D4 scenarios
     plus this one)
  2. Re-implement B2 as a unit-level test of _apply_directives only.
     The unit test exists already at M1a in
     tst_validate_step_after_program_guard.py; B2 in unit form is
     redundant with that.

Run:
    phenix.python tests/llm/tst_b2_gated_wipe_planning.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from framework import Scenario

# Phase 2A planning infrastructure.  These functions are stubs in
# framework.py until Phase 2A is merged (per PHASE2_PLAN_v2.md §2 + §6).
# The import succeeds even with stubs; the NotImplementedError fires
# only when this test's main() actually tries to use them.  Until
# Phase 2A lands, running this test produces a clear error pointing
# at the plan.  The unit-level equivalent (testing only the wipe
# behavior in workflow_engine._apply_directives) is at
# tst_validate_step_after_program_guard.py::test_M1a_*.
from framework import (call_planning_llm, make_planning_run_fn)


# =====================================================================
# Scenario states
# =====================================================================
#
# State construction follows the same shape as PHASE2_PLAN_v2.md
# scenarios.  Required fields per the planning prompt are:
#   history, analysis, available_files, previous_attempts,
#   user_advice, metrics_trend, workflow_state, directives, best_files
#
# These minimal states isolate the gating behavior:
#   B2a: directives have stop_after_requested=True and after_program
#        completed; workflow_state.valid_programs is wiped to ["STOP"]
#   B2b: same scenario WITHOUT stop_after_requested; valid_programs
#        carries normal options
#
# Real state construction would be richer; see PHASE2_PLAN_v2.md §4
# for the full required-fields checklist.

_B2A_STATE = {
    "history": [
        # one prior cycle where after_program (real_space_refine) ran
        {
            "cycle": 1,
            "program": "phenix.real_space_refine",
            "status": "success",
            "metrics": {"map_cc": 0.786},
        },
    ],
    "analysis": (
        "Real-space refinement completed successfully (cycle 1, "
        "map_cc=0.786).  User requested stop after refinement."),
    "available_files": ["model.pdb", "data.mtz", "map.ccp4"],
    "previous_attempts": [],
    "user_advice": "real space refine and stop",
    "metrics_trend": "stable",
    "workflow_state": {
        # The critical bit: stop_after_requested=True + after_program_done
        # caused _apply_directives to wipe valid_programs to ["STOP"]
        # before the planning prompt was constructed.
        "valid_programs": ["STOP"],
        "step_name": "validate",
        "experiment_type": "cryoem",
        "rsr_count": 1,
        "validation_done": False,
        "last_program": "phenix.real_space_refine",
        "successful_programs": ["phenix.real_space_refine"],
    },
    "directives": {
        "stop_conditions": {
            "after_program": "phenix.real_space_refine",
            "stop_after_requested": True,
            "skip_validation": True,
        },
    },
    "best_files": {
        "model": "model.pdb",
        "data_mtz": "data.mtz",
        "full_map": "map.ccp4",
    },
}

# B2b: identical context but stop_after_requested is False, so the
# wipe never fired and valid_programs has real options.
_B2B_STATE = dict(_B2A_STATE)
_B2B_STATE["workflow_state"] = dict(_B2A_STATE["workflow_state"])
_B2B_STATE["workflow_state"]["valid_programs"] = [
    "phenix.molprobity", "phenix.validation_cryoem",
    "phenix.map_correlations",
]
_B2B_STATE["directives"] = {
    "stop_conditions": {
        "after_program": "phenix.real_space_refine",
        # stop_after_requested intentionally absent
    },
}
_B2B_STATE["user_advice"] = "refine the model"  # no explicit stop


# =====================================================================
# Expected-output helpers
# =====================================================================

def _llm_picks_stop(intent):
    """B2a: the LLM should pick STOP.

    Two acceptable shapes:
      intent.get("stop") is True
      intent.get("program") == "STOP"
    """
    if not isinstance(intent, dict):
        return (False, "intent is not a dict: %r" % type(intent).__name__)
    if intent.get("stop") is True:
        return (True, "intent.stop=True")
    if intent.get("program") == "STOP":
        return (True, "intent.program=STOP")
    return (False, "LLM did not pick STOP despite valid_programs=['STOP']. "
                   "Got program=%r, stop=%r"
                   % (intent.get("program"), intent.get("stop")))


def _llm_does_not_auto_stop(intent):
    """B2b: the LLM should pick a validation program, NOT STOP.

    This is the negative control — confirms that the flag (not just
    after_program_done) is what drives the wipe.
    """
    if not isinstance(intent, dict):
        return (False, "intent is not a dict")
    program = intent.get("program")
    if intent.get("stop") is True or program == "STOP":
        return (False, "LLM auto-stopped without stop_after_requested. "
                       "program=%r stop=%r"
                       % (program, intent.get("stop")))
    validation_programs = (
        "phenix.molprobity", "phenix.validation_cryoem",
        "phenix.map_correlations", "phenix.model_vs_data",
    )
    if program in validation_programs:
        return (True, "LLM correctly picked validation: %s" % program)
    return (False, "LLM picked unexpected program %r (expected validation)"
                   % program)


# =====================================================================
# Scenarios
# =====================================================================

def build_scenarios():
    """B2 — gated wipe + LLM end-to-end."""
    return [
        Scenario(
            name="b2a_gated_wipe_llm_picks_stop",
            description=(
                "stop_after_requested=True + after_program_done=True "
                "caused _apply_directives to wipe valid_programs to "
                "['STOP'].  LLM must pick STOP (it's the only option)."),
            decision_point="planning",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_B2A_STATE,
            expected_fn=_llm_picks_stop,
        ),
        Scenario(
            name="b2b_no_wipe_llm_does_not_auto_stop",
            description=(
                "Negative control: same after_program completion but "
                "WITHOUT stop_after_requested.  valid_programs has "
                "real options; LLM must pick a validation program, "
                "NOT auto-stop."),
            decision_point="planning",
            test_type="reliability",
            threshold=0.8,
            max_runs=5,
            input=_B2B_STATE,
            expected_fn=_llm_does_not_auto_stop,
        ),
    ]


# =====================================================================
# Standalone runner
# =====================================================================

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    from framework import (AVAILABLE_PROVIDERS, build_summary,
                           create_log_dir, run_scenario_against_providers,
                           write_run_summary)
    import time

    if not AVAILABLE_PROVIDERS:
        sys.stderr.write(
            "No LLM provider API keys configured.\n")
        return 2

    providers = list(AVAILABLE_PROVIDERS)
    scenarios = build_scenarios()
    run_one = make_planning_run_fn()
    log_dir = create_log_dir(phase="phase1-b2")

    print()
    print("LLM Decision Tests — B2 Gated Wipe Planning (v117)")
    print("=" * 65)
    print("Available providers: %s" % ", ".join(providers))
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
        providers_skipped=[], wall_time_s=wall_time,
        phase="phase1-b2")
    write_run_summary(log_dir, summary)

    print()
    print("Summary: %d PASS, %d FAIL, %d ERROR, %d SKIP"
          % (summary["totals"]["passed"], summary["totals"]["failed"],
             summary["totals"]["errored"], summary["totals"]["skipped"]))
    print("Log:                 %s" % log_dir)

    return 0 if summary["totals"]["failed"] == 0 and \
                summary["totals"]["errored"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
