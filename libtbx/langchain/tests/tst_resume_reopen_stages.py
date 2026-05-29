"""K_H15_ITEM_2: Resume reopen targeted stages (v119.H15).

Tom's bromodomain resume failure (run 144) had three
compounding causes:
  - Bug 1 corrupted final_refinement.status to COMPLETE
  - Resume cleared gate_stop but NOT per-stage statuses
  - Gate immediately re-fired "all stages complete" →
    LLM saw STATE=complete → chose polder

H15 Item 1 prevents the corruption going forward.  Item
2 is the resume safety net: when new advice is provided
on resume, walk the directives and reopen affected
stages.

Item 2 uses a TARGETED single-stage reopen (per Gemini
critique of the original H15 plan):
  - Find the LATEST completed stage whose `programs`
    contains a program named in directives.program_settings.
  - Reset ONLY that stage to PENDING.  cycles_used → 0,
    runtime fields cleared.
  - Don't cascade-reset stages downstream of the
    reopened one.
  - Earlier stages with the same program are NOT reopened.

This keeps blast radius O(1) regardless of plan size.

7 tests:
  §A: Tom's exact scenario
  §B: No advice change → no-op
  §C: Directive for unmatched program → no reopen
  §D: Multiple programs in directives
  §E: Skipped stages stay skipped
  §F: Empty plan / empty directives → no-op, no exception
  §G: Stage's strategy already honors directive — H15 design says reopen anyway
"""
from __future__ import absolute_import, division, print_function

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.dirname(_HERE)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from knowledge.plan_schema import (
    StageDef, StructurePlan,
    STAGE_PENDING, STAGE_ACTIVE, STAGE_COMPLETE,
    STAGE_FAILED, STAGE_SKIPPED,
)
from agent.plan_generator import reopen_stages_for_directives


# =====================================================================
# Test fixtures
# =====================================================================

def _toms_corrupted_plan():
    """Build Tom's plan in the post-Bug-1 corrupted state:
    every stage marked COMPLETE despite final_refinement
    having r_free=0.272 > 0.25 target."""
    stages = [
        StageDef(
            id="data_assessment",
            programs=["phenix.xtriage"],
            max_cycles=1,
            success_criteria={"xtriage_completed": "true"},
        ),
        StageDef(
            id="model_prediction",
            programs=["phenix.predict_and_build"],
            max_cycles=1,
            success_criteria={"predict_completed": "true"},
        ),
        StageDef(
            id="initial_refinement",
            programs=["phenix.refine"],
            max_cycles=3,
            success_criteria={"r_free": "<0.35"},
        ),
        StageDef(
            id="model_rebuilding",
            programs=["phenix.autobuild", "phenix.refine"],
            max_cycles=2,
            success_criteria={"r_free": "<0.30"},
            skip_if="r_free < 0.28",
        ),
        StageDef(
            id="ligand_fitting",
            programs=["phenix.ligandfit"],
            max_cycles=1,
            success_criteria={"ligand_cc": ">0.7"},
        ),
        StageDef(
            id="final_refinement",
            programs=["phenix.refine"],
            max_cycles=3,
            success_criteria={"r_free": "<0.25"},
        ),
        StageDef(
            id="validation",
            programs=["phenix.molprobity"],
            max_cycles=1,
            success_criteria={},
        ),
    ]
    plan = StructurePlan(stages=stages)
    # Corrupted state: everything COMPLETE except model_rebuilding (SKIPPED)
    for i, s in enumerate(plan.stages):
        if s.id == "model_rebuilding":
            s.status = STAGE_SKIPPED
        else:
            s.status = STAGE_COMPLETE
        s.cycles_used = 1
    plan.current_stage_index = len(plan.stages)  # past end
    return plan


# =====================================================================
# §A: Tom's exact scenario
# =====================================================================

def test_toms_scenario_reopens_only_final_refinement():
    """Tom's resume: directives say phenix.refine needs
    generate_rfree_flags=True.  TWO completed stages
    contain phenix.refine (initial_refinement at idx 2,
    final_refinement at idx 5).  Targeted reopen picks
    ONLY the latest (final_refinement) — earlier work is
    preserved."""
    plan = _toms_corrupted_plan()
    directives = {
        "program_settings": {
            "phenix.refine": {
                "generate_rfree_flags": True,
            },
        },
    }

    reopened = reopen_stages_for_directives(plan, directives)

    assert reopened == ["final_refinement"], (
        "Should reopen ONLY final_refinement (latest with "
        "phenix.refine), got %r" % reopened)
    # Verify the right stage is reset
    final_ref = plan.stages[5]
    assert final_ref.status == STAGE_PENDING
    assert final_ref.cycles_used == 0
    assert final_ref.start_cycle is None
    assert final_ref.end_cycle is None
    assert final_ref.result_metrics == {}
    # Earlier refine stage stays COMPLETE
    initial_ref = plan.stages[2]
    assert initial_ref.status == STAGE_COMPLETE, (
        "initial_refinement should NOT be reopened (earlier "
        "than the targeted final_refinement), got %r"
        % initial_ref.status)
    # current_stage_index pulled back to the reopened stage
    assert plan.current_stage_index == 5
    print("  PASS: test_toms_scenario_reopens_only_final_refinement")


# =====================================================================
# §B: No advice change → no-op
# =====================================================================

def test_no_directives_no_reopen():
    """Empty directives → no stages reopened.  No
    exception, no mutation."""
    plan = _toms_corrupted_plan()
    snapshot = [(s.id, s.status) for s in plan.stages]

    reopened = reopen_stages_for_directives(plan, {})

    assert reopened == [], (
        "Empty directives → no reopen, got %r" % reopened)
    assert [(s.id, s.status) for s in plan.stages] == snapshot, (
        "Plan should be unchanged")
    print("  PASS: test_no_directives_no_reopen")


# =====================================================================
# §C: Directive for unmatched program → no reopen
# =====================================================================

def test_unmatched_program_no_reopen():
    """If directives mention a program not in any plan
    stage, no reopen.  No exception."""
    plan = _toms_corrupted_plan()
    directives = {
        "program_settings": {
            "phenix.imaginary_program": {"foo": "bar"},
        },
    }
    snapshot = [(s.id, s.status) for s in plan.stages]

    reopened = reopen_stages_for_directives(plan, directives)

    assert reopened == [], (
        "Unmatched program → no reopen, got %r" % reopened)
    assert [(s.id, s.status) for s in plan.stages] == snapshot
    print("  PASS: test_unmatched_program_no_reopen")


# =====================================================================
# §D: Multiple programs — each reopens its latest matching stage
# =====================================================================

def test_multiple_programs_each_reopens_latest():
    """Two programs in directives.  Each finds its OWN
    latest matching stage."""
    stages = [
        StageDef(id="a", programs=["phenix.refine"]),
        StageDef(id="b", programs=["phenix.autobuild"]),
        StageDef(id="c", programs=["phenix.refine"]),
        StageDef(id="d", programs=["phenix.autobuild"]),
        StageDef(id="e", programs=["phenix.molprobity"]),
    ]
    plan = StructurePlan(stages=stages)
    for s in plan.stages:
        s.status = STAGE_COMPLETE
        s.cycles_used = 1
    plan.current_stage_index = len(stages)
    directives = {
        "program_settings": {
            "phenix.refine": {"x": 1},
            "phenix.autobuild": {"y": 2},
        },
    }

    reopened = reopen_stages_for_directives(plan, directives)

    # Both reopens occurred — order may vary by set iteration
    assert set(reopened) == {"c", "d"}, (
        "Should reopen latest of each program (c for refine, "
        "d for autobuild), got %r" % reopened)
    assert plan.stages[0].status == STAGE_COMPLETE  # earlier refine
    assert plan.stages[2].status == STAGE_PENDING   # latest refine
    assert plan.stages[1].status == STAGE_COMPLETE  # earlier autobuild
    assert plan.stages[3].status == STAGE_PENDING   # latest autobuild
    assert plan.stages[4].status == STAGE_COMPLETE  # untouched
    # current_stage_index = min reopened index = 2
    assert plan.current_stage_index == 2
    print("  PASS: test_multiple_programs_each_reopens_latest")


# =====================================================================
# §E: Skipped stages stay skipped
# =====================================================================

def test_skipped_stages_not_reopened():
    """A SKIPPED stage matching a program in directives
    should NOT be reopened.  Skips are structural
    decisions (skip_if triggered) that the directive
    doesn't override."""
    stages = [
        StageDef(id="a", programs=["phenix.refine"]),
        StageDef(id="b", programs=["phenix.refine"]),
    ]
    plan = StructurePlan(stages=stages)
    plan.stages[0].status = STAGE_SKIPPED
    plan.stages[1].status = STAGE_SKIPPED
    plan.current_stage_index = 2
    directives = {
        "program_settings": {"phenix.refine": {"x": 1}},
    }

    reopened = reopen_stages_for_directives(plan, directives)

    assert reopened == [], (
        "Skipped stages should not be reopened, got %r"
        % reopened)
    assert plan.stages[0].status == STAGE_SKIPPED
    assert plan.stages[1].status == STAGE_SKIPPED
    print("  PASS: test_skipped_stages_not_reopened")


# =====================================================================
# §F: Empty plan / empty directives → no-op
# =====================================================================

def test_empty_plan_no_op():
    """Empty plan or None plan should be no-op, no
    exception."""
    plan = StructurePlan(stages=[])
    reopened = reopen_stages_for_directives(
        plan, {"program_settings": {"phenix.refine": {"x": 1}}},
    )
    assert reopened == []

    # None plan → returns empty, no exception
    reopened2 = reopen_stages_for_directives(
        None, {"program_settings": {"phenix.refine": {"x": 1}}},
    )
    assert reopened2 == []

    # Non-dict directives → no-op, no exception
    reopened3 = reopen_stages_for_directives(
        StructurePlan(stages=[]), None,
    )
    assert reopened3 == []
    print("  PASS: test_empty_plan_no_op")


# =====================================================================
# §G: 'default' key in program_settings is skipped
# =====================================================================

def test_default_key_in_program_settings_skipped():
    """The 'default' key in program_settings applies to
    every program — it's not a specific program name and
    should not trigger reopens."""
    stages = [
        StageDef(id="a", programs=["phenix.refine"]),
        StageDef(id="b", programs=["phenix.molprobity"]),
    ]
    plan = StructurePlan(stages=stages)
    for s in plan.stages:
        s.status = STAGE_COMPLETE
    plan.current_stage_index = len(stages)
    directives = {
        "program_settings": {
            "default": {"anisotropic_adp": True},
        },
    }

    reopened = reopen_stages_for_directives(plan, directives)

    assert reopened == [], (
        "'default' key should be skipped (applies globally, "
        "not to a specific program), got %r" % reopened)
    print("  PASS: test_default_key_in_program_settings_skipped")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    test_toms_scenario_reopens_only_final_refinement()
    test_no_directives_no_reopen()
    test_unmatched_program_no_reopen()
    test_multiple_programs_each_reopens_latest()
    test_skipped_stages_not_reopened()
    test_empty_plan_no_op()
    test_default_key_in_program_settings_skipped()


if __name__ == "__main__":
    print("K_H15_ITEM_2: Resume reopen targeted stages (v119.H15)")
    print("=" * 70)
    run_all_tests()
    print("=" * 70)
    print("K_H15_ITEM_2 complete.")
