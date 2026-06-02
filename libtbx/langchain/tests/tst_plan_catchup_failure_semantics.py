"""K_H15_ITEM_1: Plan catch-up failure semantics (v119.H15).

The H14.3 ship had a latent bug in
`StructurePlan.record_stage_cycle`'s "plan catch-up" path:
when the LLM ran a program belonging to a later stage,
the catch-up loop walked forward by calling
`self.advance()` repeatedly.  advance() unconditionally
stamped STAGE_COMPLETE on whatever was active, ignoring
the abandoned stage's success_criteria.

Tom's bromodomain session (run 135 → 144) hit this:
cycle 7's `phenix.refine` failed with no R-free flags;
cycle 8 picked `phenix.molprobity`; catch-up advanced
`final_refinement` to STAGE_COMPLETE even though its
criterion `r_free < 0.25` was unmet (r_free was stuck
at 0.272 from cycle 3).  Resume then saw "all stages
complete" and picked polder instead of re-refining.

H15 Item 1 fixes this by:
  1. Adding `force=False, structure_model=None, via=None`
     parameters to advance().  When force=False, the
     current stage's success_criteria are evaluated
     against structure_model.  Stages that meet criteria
     → STAGE_COMPLETE.  Stages that don't → STAGE_FAILED
     with failure_reason="abandoned_by_deviation".
  2. record_stage_cycle's catch-up uses strict semantics
     (force=False).  Direct callers of advance() still
     get the legacy force=True behavior unchanged.
  3. Recording structured PLAN_DEVIATION events into
     session_data so future architectural work has the
     data shape it needs.

This file pins the new semantics with 7 tests:
  - §A: Tom's exact scenario, expecting STAGE_FAILED
  - §B: Legitimate catch-up (criteria met) → COMPLETE
  - §C: Multi-stage catch-up with mixed criteria
  - §D: No structure_model → legacy COMPLETE behavior
  - §E: advance(force=True) regression — unchanged
  - §F: Empty success_criteria → COMPLETE
  - §G: STAGE_FAILED constant actually reached
"""
from __future__ import absolute_import, division, print_function

import os
import sys


# Sandbox sys.path
_HERE = os.path.dirname(os.path.abspath(__file__))
_PARENT = os.path.dirname(_HERE)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from knowledge.plan_schema import (
    StageDef, StructurePlan,
    STAGE_ACTIVE, STAGE_COMPLETE,
    STAGE_FAILED, STAGE_SKIPPED,
)


# =====================================================================
# Test fixtures
# =====================================================================

class FakeStructureModel(object):
    """Minimal StructureModel mock — only needs get_metric.

    The real StructureModel has many fields; the gate evaluator's
    _check_success uses _get_metric_value(structure_model, name)
    which itself uses getattr/dict-access fallback.  A dict-like
    that returns scalars is enough for criteria evaluation."""

    def __init__(self, metrics=None):
        self._metrics = dict(metrics or {})

    def get_metric(self, name):
        return self._metrics.get(name)

    # _get_metric_value also tries getattr, so expose as attrs
    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        m = self.__dict__.get("_metrics", {})
        if name in m:
            return m[name]
        raise AttributeError(name)


def _toms_plan():
    """Build the 7-stage plan from Tom's bromodomain session.

    Mirrors the plan in agent_session_144.json verbatim
    where it matters: stage IDs, programs lists,
    success_criteria, max_cycles."""
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
    plan = StructurePlan(
        goal="Predict + MR + refine + ligand",
        stages=stages,
    )
    return plan


def _fast_forward_to(plan, target_id, mark_complete=True):
    """Move the plan's current_stage_index to target_id.

    Marks all stages before target as COMPLETE (or SKIPPED for
    stages whose skip_if would fire).  Sets the target stage
    ACTIVE.  Doesn't touch result_metrics — tests that need
    metrics on completed stages should set them after."""
    target_idx = None
    for i, s in enumerate(plan.stages):
        if s.id == target_id:
            target_idx = i
            break
    assert target_idx is not None, target_id
    for i in range(target_idx):
        plan.stages[i].status = (
            STAGE_COMPLETE if mark_complete else STAGE_SKIPPED
        )
    plan.current_stage_index = target_idx
    plan.stages[target_idx].status = STAGE_ACTIVE
    return plan


# =====================================================================
# §A: Tom's exact scenario
# =====================================================================

def test_toms_scenario_marks_final_refinement_failed():
    """Reproduce Tom's bromodomain bug exactly.

    Setup: 7-stage plan, final_refinement is ACTIVE with
    cycles_used=1, structure_model has r_free=0.272.
    Action: record_stage_cycle('phenix.molprobity') —
    triggers catch-up to `validation`.
    Expected (H15): final_refinement.status == STAGE_FAILED
    with failure_reason='abandoned_by_deviation'.  Pre-H15:
    final_refinement was silently marked STAGE_COMPLETE."""
    plan = _toms_plan()
    # Skip model_rebuilding (Tom's session had it SKIPPED
    # via skip_if since r_free=0.272 < 0.28)
    plan.stages[3].status = STAGE_SKIPPED
    _fast_forward_to(plan, "final_refinement")
    plan.stages[5].cycles_used = 1
    sm = FakeStructureModel({"r_free": 0.272})
    session_data = {}

    plan.record_stage_cycle(
        "phenix.molprobity",
        structure_model=sm,
        session_data=session_data,
        cycle_number=8,
    )

    final_ref = plan.stages[5]
    assert final_ref.status == STAGE_FAILED, (
        "final_refinement should be FAILED (r_free=0.272 "
        ">= 0.25), got %r" % final_ref.status)
    assert final_ref.failure_reason == "abandoned_by_deviation", (
        "failure_reason should be 'abandoned_by_deviation', got %r"
        % final_ref.failure_reason)
    # Validation should now be the active stage
    validation = plan.stages[6]
    assert validation.status == STAGE_ACTIVE, (
        "validation should be ACTIVE, got %r" % validation.status)
    assert validation.entered_via == "deviation", (
        "validation.entered_via should be 'deviation', got %r"
        % validation.entered_via)
    # Deviation event recorded
    deviations = session_data.get("plan_deviations", [])
    assert len(deviations) == 1, (
        "Expected exactly 1 deviation event, got %d"
        % len(deviations))
    dev = deviations[0]
    assert dev["from_stage"] == "final_refinement"
    assert dev["to_stage"] == "validation"
    assert dev["program"] == "phenix.molprobity"
    assert dev["cycle"] == 8
    print("  PASS: test_toms_scenario_marks_final_refinement_failed")


# =====================================================================
# §B: Legitimate catch-up — criteria met → COMPLETE
# =====================================================================

def test_legitimate_catchup_marks_complete():
    """When the abandoned stage's criteria ARE met, mark
    STAGE_COMPLETE (not STAGE_FAILED).  This is the
    legitimate-catch-up case the original code was designed
    for: agent ran a later program because the current stage
    was effectively done.

    Setup: final_refinement ACTIVE, r_free=0.22 (meets <0.25).
    Action: record_stage_cycle('phenix.molprobity').
    Expected: final_refinement.status == STAGE_COMPLETE
    (NOT failed) and failure_reason stays None."""
    plan = _toms_plan()
    plan.stages[3].status = STAGE_SKIPPED
    _fast_forward_to(plan, "final_refinement")
    plan.stages[5].cycles_used = 2
    sm = FakeStructureModel({"r_free": 0.22})  # meets <0.25

    plan.record_stage_cycle(
        "phenix.molprobity",
        structure_model=sm,
        session_data={},
        cycle_number=8,
    )

    final_ref = plan.stages[5]
    assert final_ref.status == STAGE_COMPLETE, (
        "criteria met (r_free=0.22 < 0.25), expected COMPLETE, "
        "got %r" % final_ref.status)
    assert final_ref.failure_reason is None, (
        "failure_reason should stay None on legit advance, got %r"
        % final_ref.failure_reason)
    validation = plan.stages[6]
    assert validation.status == STAGE_ACTIVE
    # Even on legitimate catch-up, entered_via is still
    # "deviation" — it was an LLM choice that drove the
    # transition, even if the outcome happened to be okay.
    assert validation.entered_via == "deviation"
    print("  PASS: test_legitimate_catchup_marks_complete")


# =====================================================================
# §C: Multi-stage catch-up with mixed criteria
# =====================================================================

def test_multistage_catchup_mixed_criteria():
    """When catch-up walks past multiple intermediate stages,
    each gets its OWN criteria check.

    Setup: 3-stage chain where stage_1 meets criteria,
    stage_2 doesn't, target is stage_3.  Initial state:
    stage_1 ACTIVE.  Catch-up should:
      - advance stage_1 (criteria met → COMPLETE)
      - advance stage_2 (criteria not met → FAILED)
      - land on stage_3 ACTIVE.

    We use distinct metric names to give each stage its own
    pass/fail story."""
    stages = [
        StageDef(
            id="s1",
            programs=["prog1"],
            success_criteria={"r_free": "<0.40"},
            max_cycles=2,
        ),
        StageDef(
            id="s2",
            programs=["prog2"],
            success_criteria={"r_free": "<0.20"},
            max_cycles=2,
        ),
        StageDef(
            id="s3",
            programs=["prog3"],
            success_criteria={},
            max_cycles=1,
        ),
    ]
    plan = StructurePlan(goal="multistage", stages=stages)
    plan.stages[0].status = STAGE_ACTIVE
    sm = FakeStructureModel({"r_free": 0.30})  # meets s1, fails s2
    session_data = {}

    plan.record_stage_cycle(
        "prog3",
        structure_model=sm,
        session_data=session_data,
        cycle_number=10,
    )

    assert plan.stages[0].status == STAGE_COMPLETE, (
        "s1 criteria met (r_free=0.30 < 0.40), expected "
        "COMPLETE, got %r" % plan.stages[0].status)
    assert plan.stages[1].status == STAGE_FAILED, (
        "s2 criteria not met (r_free=0.30 NOT < 0.20), expected "
        "FAILED, got %r" % plan.stages[1].status)
    assert plan.stages[1].failure_reason == "abandoned_by_deviation"
    assert plan.stages[2].status == STAGE_ACTIVE, (
        "s3 should be ACTIVE, got %r" % plan.stages[2].status)
    # Only the first advance gets entered_via=deviation
    # (subsequent advances cascade from the same deviation event)
    assert plan.stages[1].entered_via == "deviation"
    assert plan.stages[2].entered_via is None, (
        "s3 (second advance) should NOT have entered_via set, "
        "got %r" % plan.stages[2].entered_via)
    # Single deviation event for the whole catch-up
    deviations = session_data.get("plan_deviations", [])
    assert len(deviations) == 1, (
        "Expected 1 deviation event for the catch-up, got %d"
        % len(deviations))
    print("  PASS: test_multistage_catchup_mixed_criteria")


# =====================================================================
# §D: No structure_model → legacy COMPLETE behavior
# =====================================================================

def test_no_structure_model_falls_back_to_legacy():
    """When structure_model is None (e.g., ai_agent.py is
    pre-H15 and doesn't pass one), _criteria_met returns
    True → stage marked COMPLETE.  Preserves legacy
    behavior; graceful degradation during partial deployment.

    Setup: same as Tom's bug, but pass structure_model=None.
    Expected: final_refinement marked COMPLETE (legacy)."""
    plan = _toms_plan()
    plan.stages[3].status = STAGE_SKIPPED
    _fast_forward_to(plan, "final_refinement")

    plan.record_stage_cycle(
        "phenix.molprobity",
        structure_model=None,
        session_data=None,
        cycle_number=8,
    )

    final_ref = plan.stages[5]
    assert final_ref.status == STAGE_COMPLETE, (
        "Without structure_model, should fall back to legacy "
        "COMPLETE, got %r" % final_ref.status)
    assert final_ref.failure_reason is None
    print("  PASS: test_no_structure_model_falls_back_to_legacy")


# =====================================================================
# §E: advance(force=True) regression — unchanged behavior
# =====================================================================

def test_advance_force_true_unchanged():
    """All pre-H15 callers of advance() pass no arguments,
    which defaults force=True.  That path must mark
    STAGE_COMPLETE unconditionally (no criteria check).

    Setup: stage ACTIVE with success_criteria that AREN'T
    met by a (passed-but-ignored) model.
    Expected: STAGE_COMPLETE, not STAGE_FAILED."""
    stages = [
        StageDef(
            id="a",
            programs=["p1"],
            success_criteria={"r_free": "<0.20"},
        ),
        StageDef(
            id="b",
            programs=["p2"],
        ),
    ]
    plan = StructurePlan(stages=stages)
    plan.stages[0].status = STAGE_ACTIVE
    sm = FakeStructureModel({"r_free": 0.50})  # would fail criteria

    # Default: force=True
    advanced = plan.advance()
    assert advanced is True
    assert plan.stages[0].status == STAGE_COMPLETE, (
        "advance() with default force=True should always mark "
        "COMPLETE; got %r" % plan.stages[0].status)
    assert plan.stages[0].failure_reason is None, (
        "force=True path doesn't set failure_reason; got %r"
        % plan.stages[0].failure_reason)
    print("  PASS: test_advance_force_true_unchanged")


# =====================================================================
# §F: Empty success_criteria → COMPLETE
# =====================================================================

def test_empty_criteria_marks_complete():
    """Validation stages (and other observe-only stages)
    have empty success_criteria.  _criteria_met returns
    True for empty criteria, so advance(force=False) marks
    STAGE_COMPLETE — not FAILED.

    Setup: stage ACTIVE with success_criteria={}, then
    catch up over it via record_stage_cycle for a later
    program."""
    stages = [
        StageDef(
            id="observe",
            programs=["phenix.molprobity"],
            success_criteria={},
        ),
        StageDef(
            id="target",
            programs=["phenix.refine"],
        ),
    ]
    plan = StructurePlan(stages=stages)
    plan.stages[0].status = STAGE_ACTIVE
    sm = FakeStructureModel({"r_free": 0.99})

    plan.record_stage_cycle(
        "phenix.refine",
        structure_model=sm,
        session_data={},
        cycle_number=1,
    )

    assert plan.stages[0].status == STAGE_COMPLETE, (
        "Empty criteria should advance to COMPLETE, got %r"
        % plan.stages[0].status)
    assert plan.stages[1].status == STAGE_ACTIVE
    print("  PASS: test_empty_criteria_marks_complete")


# =====================================================================
# §G: STAGE_FAILED constant actually reached
# =====================================================================

def test_stage_failed_constant_is_reachable():
    """Regression guard against future refactors that might
    accidentally drop the FAILED path from advance().

    Asserts that the STAGE_FAILED value appears as the
    status of a stage after a deviation catch-up over an
    unmet-criteria stage.  This is independent of test_A's
    assertions about Tom's specific plan; this is a
    structural guarantee about the new semantics."""
    stages = [
        StageDef(
            id="must_fail",
            programs=["p1"],
            success_criteria={"r_free": "<0.10"},
        ),
        StageDef(
            id="target",
            programs=["p2"],
        ),
    ]
    plan = StructurePlan(stages=stages)
    plan.stages[0].status = STAGE_ACTIVE
    sm = FakeStructureModel({"r_free": 0.50})

    plan.record_stage_cycle(
        "p2",
        structure_model=sm,
        session_data={},
        cycle_number=1,
    )

    # Direct assertion against the constant value
    statuses = [s.status for s in plan.stages]
    assert STAGE_FAILED in statuses, (
        "STAGE_FAILED must appear in stage statuses after a "
        "deviation over an unmet-criteria stage.  Got: %r"
        % statuses)
    # And it must be on the right stage
    assert plan.stages[0].status == STAGE_FAILED
    print("  PASS: test_stage_failed_constant_is_reachable")


# =====================================================================
# §H: H1 micro-fix — gate stop reason surfaces FAILED stages
# =====================================================================

def test_get_failed_stage_ids_basic():
    """plan.get_failed_stage_ids() returns the list of stages
    currently in STAGE_FAILED.  Pure query, no mutation."""
    stages = [
        StageDef(id="a"),
        StageDef(id="b"),
        StageDef(id="c"),
    ]
    plan = StructurePlan(stages=stages)
    plan.stages[0].status = STAGE_COMPLETE
    plan.stages[1].status = STAGE_FAILED
    plan.stages[2].status = STAGE_FAILED

    failed = plan.get_failed_stage_ids()
    assert failed == ["b", "c"], (
        "Expected ['b', 'c'], got %r" % failed)

    # No FAILED stages → empty list
    plan.stages[1].status = STAGE_COMPLETE
    plan.stages[2].status = STAGE_COMPLETE
    failed = plan.get_failed_stage_ids()
    assert failed == [], (
        "Expected [] with no FAILED stages, got %r" % failed)
    print("  PASS: test_get_failed_stage_ids_basic")


def test_gate_stop_reason_surfaces_failed_stages():
    """When is_complete() returns True AND at least one stage
    is STAGE_FAILED, the gate evaluator's stop reason must
    surface the failure rather than say 'all stages complete.'

    This is the H1 fix — without it, Tom's bromodomain
    scenario stopped with the misleading 'all stages complete'
    message even though final_refinement was FAILED."""
    # Inline import so this test still runs even if gate_evaluator
    # isn't importable (we'd want graceful skip in that case)
    try:
        from agent.gate_evaluator import GateEvaluator
    except ImportError:
        print("  SKIP: test_gate_stop_reason_surfaces_failed_stages "
              "(gate_evaluator unavailable)")
        return

    stages = [
        StageDef(id="early", success_criteria={}),
        StageDef(id="final_refinement", success_criteria={"r_free": "<0.25"}),
        StageDef(id="validation", success_criteria={}),
    ]
    plan = StructurePlan(stages=stages)
    # Mimic the post-Bug-1-fix state: final_refinement FAILED,
    # validation COMPLETE
    plan.stages[0].status = STAGE_COMPLETE
    plan.stages[1].status = STAGE_FAILED
    plan.stages[1].failure_reason = "abandoned_by_deviation"
    plan.stages[2].status = STAGE_COMPLETE
    plan.current_stage_index = len(stages)  # past end

    # Sanity: is_complete returns True even with one FAILED
    assert plan.is_complete() is True, (
        "is_complete should return True with no PENDING/ACTIVE "
        "stages even if one is FAILED — got %r"
        % plan.is_complete())

    # Run the gate
    gate = GateEvaluator()
    result = gate.evaluate(plan, None, None, 10)

    assert result.action == "stop", (
        "Expected action=stop, got %r" % result.action)
    # The new reason should mention the failure
    assert "fail" in result.reason.lower(), (
        "Stop reason should surface the FAILED stage; got %r"
        % result.reason)
    assert "final_refinement" in result.reason, (
        "Stop reason should name the failed stage; got %r"
        % result.reason)
    # And it should NOT be the misleading old message
    assert result.reason != "all stages complete", (
        "Stop reason should NOT be the legacy 'all stages "
        "complete' when stages are FAILED; got %r"
        % result.reason)
    print("  PASS: test_gate_stop_reason_surfaces_failed_stages")


def test_gate_stop_reason_unchanged_when_no_failed_stages():
    """Regression guard: when NO stages are FAILED, the gate
    must still emit the legacy 'all stages complete' reason.
    The H1 fix only changes behavior for the failure case."""
    try:
        from agent.gate_evaluator import GateEvaluator
    except ImportError:
        print("  SKIP: test_gate_stop_reason_unchanged_when_no_failed_stages")
        return

    stages = [
        StageDef(id="a", success_criteria={}),
        StageDef(id="b", success_criteria={}),
    ]
    plan = StructurePlan(stages=stages)
    plan.stages[0].status = STAGE_COMPLETE
    plan.stages[1].status = STAGE_COMPLETE
    plan.current_stage_index = len(stages)

    gate = GateEvaluator()
    result = gate.evaluate(plan, None, None, 5)

    assert result.action == "stop"
    assert result.reason == "all stages complete", (
        "When no stages are FAILED, the legacy reason must be "
        "preserved unchanged; got %r" % result.reason)
    print("  PASS: test_gate_stop_reason_unchanged_when_no_failed_stages")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    test_toms_scenario_marks_final_refinement_failed()
    test_legitimate_catchup_marks_complete()
    test_multistage_catchup_mixed_criteria()
    test_no_structure_model_falls_back_to_legacy()
    test_advance_force_true_unchanged()
    test_empty_criteria_marks_complete()
    test_stage_failed_constant_is_reachable()
    # H1 micro-fix tests
    test_get_failed_stage_ids_basic()
    test_gate_stop_reason_surfaces_failed_stages()
    test_gate_stop_reason_unchanged_when_no_failed_stages()


if __name__ == "__main__":
    print("K_H15_ITEM_1: Plan catch-up failure semantics (v119.H15)")
    print("=" * 70)
    run_all_tests()
    print("=" * 70)
    print("K_H15_ITEM_1 complete.")
