"""
Unit tests for skip_to_program (v116.10).

Tests StructurePlan.skip_to_program() and the format_plan_header()
display integration.

skip_to_program is called at plan initialization when the user
has requested a target program that isn't in Stage 1 of the plan
(e.g. "predict and stop" — target is predict_and_build but the
plan template always starts with xtriage).  Without this, the
agent would try to run xtriage first and stop with an
"after_program_not_available" error.

Coverage:
  Section A: skip_to_program core behavior (8 tests)
  Section B: serialization / resume (2 tests)
  Section C: format_plan_header display (7 tests)
  Section D: must-pass scenarios from skip_stages_plan.md (7 tests)
  Section E: backward compatibility (2 tests)
  Section F: edge cases (4 tests)
"""

from __future__ import absolute_import, division, print_function
import sys
import os

# Setup import paths
sys.path.insert(0, os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))

try:
    from libtbx.langchain.knowledge.plan_schema import (
        StructurePlan, StageDef,
        STAGE_PENDING, STAGE_ACTIVE, STAGE_COMPLETE,
        STAGE_SKIPPED, STAGE_FAILED,
    )
except ImportError:
    from knowledge.plan_schema import (
        StructurePlan, StageDef,
        STAGE_PENDING, STAGE_ACTIVE, STAGE_COMPLETE,
        STAGE_SKIPPED, STAGE_FAILED,
    )


# =====================================================================
# Helpers — plan templates matching real PHENIX workflows
# =====================================================================

def _predict_plan():
    """X-ray predict-build workflow: xtriage → predict → refine."""
    return StructurePlan(
        goal="Predict model + MR + refinement",
        stages=[
            StageDef(id="triage",
                     description="Data triage",
                     programs=["phenix.xtriage"]),
            StageDef(id="predict",
                     description="AlphaFold prediction and MR",
                     programs=["phenix.predict_and_build"]),
            StageDef(id="refine",
                     description="Reciprocal-space refinement",
                     programs=["phenix.refine"]),
            StageDef(id="validate",
                     description="Validation",
                     programs=["phenix.molprobity"]),
        ])


def _mr_refine_plan():
    """X-ray MR workflow: xtriage → phaser → refine."""
    return StructurePlan(
        goal="Molecular replacement + refinement",
        stages=[
            StageDef(id="triage",
                     description="Data triage",
                     programs=["phenix.xtriage"]),
            StageDef(id="mr",
                     description="Molecular replacement",
                     programs=["phenix.phaser"]),
            StageDef(id="refine",
                     description="Reciprocal-space refinement",
                     programs=["phenix.refine"]),
        ])


def _cryo_plan():
    """Cryo-EM workflow: mtriage → sharpen → density mod."""
    return StructurePlan(
        goal="Cryo-EM density modification",
        stages=[
            StageDef(id="mtriage",
                     description="Map triage",
                     programs=["phenix.mtriage"]),
            StageDef(id="sharpen",
                     description="Map sharpening",
                     programs=["phenix.map_sharpening"]),
            StageDef(id="denmod",
                     description="Density modification",
                     programs=["phenix.resolve_cryo_em"]),
        ])


def _multi_program_plan():
    """Plan with multiple programs per stage."""
    return StructurePlan(
        goal="Multi-program stages",
        stages=[
            StageDef(id="triage",
                     description="Data triage",
                     programs=["phenix.xtriage",
                               "phenix.mtriage"]),
            StageDef(id="solve",
                     description="Solve",
                     programs=["phenix.phaser",
                               "phenix.autosol"]),
            StageDef(id="refine",
                     description="Refinement",
                     programs=["phenix.refine"]),
        ])


# =====================================================================
# SECTION A: skip_to_program core behavior
# =====================================================================

def test_skip_to_program_forward_basic():
    """Skipping to a later stage marks Stage 1 as SKIPPED."""
    print("Test: skip_to_program_forward_basic")
    p = _predict_plan()
    assert p.current_stage_index == 0
    n = p.skip_to_program("phenix.predict_and_build")
    assert n == 1, \
        "Expected 1 stage skipped, got %d" % n
    assert p.stages[0].status == STAGE_SKIPPED, \
        "Stage 0 should be SKIPPED, got %s" % \
        p.stages[0].status
    assert p.stages[1].status == STAGE_PENDING, \
        "Stage 1 should remain PENDING, got %s" % \
        p.stages[1].status
    assert p.current_stage_index == 1, \
        "current_stage_index should be 1, got %d" % \
        p.current_stage_index
    print("  PASS")


def test_skip_to_program_target_is_current():
    """No-op when target program is in the current stage."""
    print("Test: skip_to_program_target_is_current")
    p = _predict_plan()
    # current_stage_index is 0, xtriage IS Stage 0
    n = p.skip_to_program("phenix.xtriage")
    assert n == 0, \
        "Expected 0 stages skipped, got %d" % n
    assert p.stages[0].status == STAGE_PENDING, \
        "Stage 0 should remain PENDING, got %s" % \
        p.stages[0].status
    assert p.current_stage_index == 0
    print("  PASS")


def test_skip_to_program_unknown_program():
    """No-op when program is not in any stage."""
    print("Test: skip_to_program_unknown_program")
    p = _predict_plan()
    n = p.skip_to_program("phenix.does_not_exist")
    assert n == 0
    assert p.stages[0].status == STAGE_PENDING
    assert p.current_stage_index == 0
    print("  PASS")


def test_skip_to_program_multiple_stages():
    """Skipping multiple stages marks all earlier ones SKIPPED."""
    print("Test: skip_to_program_multiple_stages")
    p = _predict_plan()  # 4 stages
    n = p.skip_to_program("phenix.molprobity")  # last stage
    assert n == 3, \
        "Expected 3 stages skipped, got %d" % n
    assert p.stages[0].status == STAGE_SKIPPED
    assert p.stages[1].status == STAGE_SKIPPED
    assert p.stages[2].status == STAGE_SKIPPED
    assert p.stages[3].status == STAGE_PENDING
    assert p.current_stage_index == 3
    print("  PASS")


def test_skip_to_program_refuses_backward():
    """Cannot skip to a stage before the current one."""
    print("Test: skip_to_program_refuses_backward")
    p = _predict_plan()
    p.current_stage_index = 2  # at refine
    n = p.skip_to_program("phenix.xtriage")  # Stage 0
    assert n == 0
    assert p.current_stage_index == 2, \
        "current_stage_index should stay at 2"
    # Stage 0 unchanged
    assert p.stages[0].status == STAGE_PENDING
    print("  PASS")


def test_skip_to_program_preserves_complete():
    """Skip does NOT mark COMPLETE stages as SKIPPED.

    Important for the resume path: if Stage 0 already
    completed in a previous session, skip_to_program
    should leave it alone (not overwrite COMPLETE
    with SKIPPED).
    """
    print("Test: skip_to_program_preserves_complete")
    p = _predict_plan()
    # Mark Stage 0 complete, agent at Stage 1
    p.stages[0].status = STAGE_COMPLETE
    p.current_stage_index = 1
    n = p.skip_to_program("phenix.refine")
    assert n == 1, "Expected 1 (just stage 1)"
    assert p.stages[0].status == STAGE_COMPLETE, \
        "Stage 0 (COMPLETE) must not be marked SKIPPED"
    assert p.stages[1].status == STAGE_SKIPPED
    assert p.current_stage_index == 2
    print("  PASS")


def test_skip_to_program_first_match_wins():
    """When target program appears in multiple stages,
    the FIRST matching stage is the target."""
    print("Test: skip_to_program_first_match_wins")
    # Build a plan with the same program in stages 1 and 3
    p = StructurePlan(
        goal="duplicate program",
        stages=[
            StageDef(id="s1", description="Stage 1",
                     programs=["phenix.xtriage"]),
            StageDef(id="s2", description="Stage 2",
                     programs=["phenix.refine"]),
            StageDef(id="s3", description="Stage 3",
                     programs=["phenix.refine"]),
        ])
    n = p.skip_to_program("phenix.refine")
    # First match is stage index 1, so 1 stage skipped
    assert n == 1
    assert p.current_stage_index == 1
    assert p.stages[1].status == STAGE_PENDING
    assert p.stages[2].status == STAGE_PENDING
    print("  PASS")


def test_skip_to_program_multi_program_stage():
    """Target can match any program in a multi-program stage."""
    print("Test: skip_to_program_multi_program_stage")
    p = _multi_program_plan()
    # autosol is in Stage 1 (index 1), second program
    n = p.skip_to_program("phenix.autosol")
    assert n == 1
    assert p.stages[0].status == STAGE_SKIPPED
    assert p.current_stage_index == 1
    print("  PASS")


# =====================================================================
# SECTION B: serialization / resume
# =====================================================================

def test_serialization_roundtrip_preserves_skipped():
    """STAGE_SKIPPED status and current_stage_index
    survive to_dict / from_dict roundtrip."""
    print("Test: serialization_roundtrip_preserves_skipped")
    p = _predict_plan()
    p.skip_to_program("phenix.refine")
    # Capture state
    expected_statuses = [s.status for s in p.stages]
    expected_idx = p.current_stage_index

    d = p.to_dict()
    p2 = StructurePlan.from_dict(d)

    got_statuses = [s.status for s in p2.stages]
    assert got_statuses == expected_statuses, \
        "Stage statuses changed across roundtrip:\n" \
        "  before: %s\n  after:  %s" % (
            expected_statuses, got_statuses)
    assert p2.current_stage_index == expected_idx, \
        "current_stage_index changed: %d → %d" % (
            expected_idx, p2.current_stage_index)
    print("  PASS")


def test_resume_skip_is_idempotent():
    """Re-running skip_to_program on a plan that has
    already been skipped is a safe no-op.

    Models the resume path: a plan with skipped stages
    is loaded from session data, then the skip logic
    runs again — should return 0 the second time.
    """
    print("Test: resume_skip_is_idempotent")
    p = _predict_plan()
    n1 = p.skip_to_program("phenix.predict_and_build")
    assert n1 == 1
    # Simulate save/restore
    p = StructurePlan.from_dict(p.to_dict())
    # Re-run skip — should be no-op
    n2 = p.skip_to_program("phenix.predict_and_build")
    assert n2 == 0, \
        "Second skip should be no-op, got %d" % n2
    assert p.current_stage_index == 1
    print("  PASS")


# =====================================================================
# SECTION C: format_plan_header display
# =====================================================================

def test_format_hides_skipped_stages():
    """STAGE_SKIPPED stages are absent from header output."""
    print("Test: format_hides_skipped_stages")
    p = _predict_plan()
    p.skip_to_program("phenix.predict_and_build")
    h = p.format_plan_header()
    assert "Stage 1" not in h, \
        "Stage 1 (skipped) should be hidden"
    assert "Data triage" not in h, \
        "Description of skipped stage should be hidden"
    assert "Stage 2" in h, "Stage 2 should be shown"
    print("  PASS")


def test_format_starting_here_when_skipped():
    """First displayed stage gets 'Starting here' annotation
    when any prior stages were skipped."""
    print("Test: format_starting_here_when_skipped")
    p = _mr_refine_plan()
    p.skip_to_program("phenix.phaser")  # skip xtriage
    h = p.format_plan_header()
    assert "Starting here" in h, \
        "Expected 'Starting here' annotation"
    # Should appear after Stage 2 (the first shown), not Stage 3
    lines = h.split("\n")
    stage2_idx = next(i for i, ln in enumerate(lines)
                       if "Stage 2" in ln)
    stage3_idx = next(i for i, ln in enumerate(lines)
                       if "Stage 3" in ln)
    starting_idx = next(i for i, ln in enumerate(lines)
                         if "Starting here" in ln)
    assert stage2_idx < starting_idx < stage3_idx, \
        "'Starting here' should appear between " \
        "Stage 2 header and Stage 3 header"
    print("  PASS")


def test_format_will_stop_here():
    """Last displayed stage gets 'Will stop here' when
    stop_after is active."""
    print("Test: format_will_stop_here")
    p = _predict_plan()
    h = p.format_plan_header(stop_after="phenix.refine")
    assert "Will stop here" in h
    # Stages 1, 2, 3 shown; Stage 4 (validate) hidden
    assert "Stage 1" in h
    assert "Stage 3" in h
    assert "Stage 4" not in h, \
        "Stage 4 (past stop) should be hidden"
    print("  PASS")


def test_format_combined_starting_and_stopping():
    """When the same stage is both first shown and the
    stop target, the annotation reads
    'Starting here — Will stop here'."""
    print("Test: format_combined_starting_and_stopping")
    p = _predict_plan()
    p.skip_to_program("phenix.predict_and_build")
    h = p.format_plan_header(
        stop_after="phenix.predict_and_build")
    assert "Starting here — Will stop here" in h, \
        "Expected combined annotation, got:\n%s" % h
    # Standalone forms should NOT appear separately
    # (the combined form replaces them)
    # Specifically, "Starting here" line on its own
    # should not exist
    lines = h.split("\n")
    standalone_starting = [
        ln for ln in lines
        if ln.strip() == "▸ Starting here"]
    standalone_stopping = [
        ln for ln in lines
        if ln.strip() == "▸ Will stop here"]
    assert not standalone_starting, \
        "Should not have standalone 'Starting here'"
    assert not standalone_stopping, \
        "Should not have standalone 'Will stop here'"
    print("  PASS")


def test_format_no_annotations_when_no_skip_no_stop():
    """Normal plan with no directives shows all stages,
    no annotations."""
    print("Test: format_no_annotations_when_no_skip_no_stop")
    p = _predict_plan()
    h = p.format_plan_header()
    assert "Starting here" not in h
    assert "Will stop here" not in h
    # All 4 stages shown
    assert "Stage 1" in h
    assert "Stage 2" in h
    assert "Stage 3" in h
    assert "Stage 4" in h
    print("  PASS")


def test_format_preserves_original_stage_numbers():
    """When stages are skipped, the displayed stage
    numbers reflect ORIGINAL positions (Stage 2, Stage 3),
    not renumbered (Stage 1, Stage 2).  This is intentional
    — the user sees where they are in the plan."""
    print("Test: format_preserves_original_stage_numbers")
    p = _predict_plan()
    p.skip_to_program("phenix.refine")  # skip 0 and 1
    h = p.format_plan_header()
    # Skipped stages 1, 2; remaining shown as 3, 4
    assert "Stage 3" in h
    assert "Stage 4" in h
    assert "Stage 1" not in h
    assert "Stage 2" not in h
    print("  PASS")


def test_format_does_not_break_on_runtime_skip():
    """The gate evaluator's existing skip_stage() call
    (used at runtime when a stage is dynamically skipped)
    still produces a sensible header.  Under the new
    rules, a stage skipped at runtime is hidden — but
    that's the only display change, no crash."""
    print("Test: format_does_not_break_on_runtime_skip")
    p = _predict_plan()
    # Simulate runtime skip via existing skip_stage()
    # (used by the gate evaluator at line ~5007 of
    # ai_agent.py)
    p.skip_stage("triage", reason="redundant")
    p.advance()
    h = p.format_plan_header()
    # Header should render without error
    assert "STRATEGY PLAN" in h
    # The runtime-skipped stage is now hidden
    assert "Stage 1" not in h, \
        "Runtime-skipped stage should also be hidden"
    # But other stages still show
    assert "Stage 2" in h
    # 'Starting here' annotation does NOT appear here
    # because the first shown stage is Stage 2 which
    # is now the current stage after advance(); the
    # _any_skipped flag is set, so 'Starting here'
    # appears on Stage 2.  This is consistent behavior.
    print("  PASS")


# =====================================================================
# SECTION D: must-pass scenarios from skip_stages_plan.md
# =====================================================================

def _simulate_init(plan, directives):
    """Mirror what _initialize_plan_inner does for skip + display.

    Returns (skipped_count, header_text).
    """
    sc = directives.get("stop_conditions", {})
    stop_after = sc.get("after_program")
    skip_target = sc.get("start_with_program") or stop_after
    skipped = 0
    if skip_target and hasattr(plan, "skip_to_program"):
        skipped = plan.skip_to_program(skip_target)
    header = plan.format_plan_header(stop_after=stop_after)
    return skipped, header


def test_scenario_predict_and_stop():
    """'predict and stop' → skip xtriage, show predict
    with combined Starting/Stop annotation."""
    print("Test: scenario_predict_and_stop")
    skipped, h = _simulate_init(
        _predict_plan(),
        {"stop_conditions": {
            "after_program": "phenix.predict_and_build"}})
    assert skipped == 1
    assert "Stage 1" not in h
    assert "Stage 2" in h
    assert "Stage 3" not in h
    assert "Starting here — Will stop here" in h
    print("  PASS")


def test_scenario_density_modify_and_stop_cryo():
    """'density modify and stop' (cryo) → skip mtriage
    and sharpen, show denmod only."""
    print("Test: scenario_density_modify_and_stop_cryo")
    skipped, h = _simulate_init(
        _cryo_plan(),
        {"stop_conditions": {
            "after_program": "phenix.resolve_cryo_em"}})
    assert skipped == 2
    assert "Stage 1" not in h
    assert "Stage 2" not in h
    assert "Stage 3" in h
    assert "Starting here — Will stop here" in h
    print("  PASS")


def test_scenario_xtriage_and_stop():
    """'run xtriage and stop' → no skip (xtriage IS Stage 1),
    show Stage 1 with 'Will stop here' only."""
    print("Test: scenario_xtriage_and_stop")
    skipped, h = _simulate_init(
        _predict_plan(),
        {"stop_conditions": {
            "after_program": "phenix.xtriage"}})
    assert skipped == 0
    assert "Stage 1" in h
    assert "Stage 2" not in h
    assert "Will stop here" in h
    assert "Starting here" not in h
    print("  PASS")


def test_scenario_solve_the_structure():
    """'solve the structure' → no skip, no stop, full plan."""
    print("Test: scenario_solve_the_structure")
    skipped, h = _simulate_init(
        _predict_plan(),
        {"stop_conditions": {}})
    assert skipped == 0
    assert "Stage 1" in h
    assert "Stage 2" in h
    assert "Stage 3" in h
    assert "Stage 4" in h
    assert "Starting here" not in h
    assert "Will stop here" not in h
    print("  PASS")


def test_scenario_phaser_refine_and_stop():
    """'phaser, refine, and stop' → start_with=phaser,
    after=refine.  Skip xtriage, show stages 2 & 3."""
    print("Test: scenario_phaser_refine_and_stop")
    skipped, h = _simulate_init(
        _mr_refine_plan(),
        {"stop_conditions": {
            "start_with_program": "phenix.phaser",
            "after_program": "phenix.refine"}})
    assert skipped == 1
    assert "Stage 1" not in h
    assert "Stage 2" in h
    assert "Stage 3" in h
    assert "Starting here" in h
    assert "Will stop here" in h
    # Combined form does NOT apply here — Stage 2 and
    # Stage 3 are different stages
    assert "Starting here — Will stop here" not in h
    print("  PASS")


def test_scenario_sharpen_and_density_modify_cryo():
    """'sharpen and density modify' (cryo) → start_with
    =sharpen, no stop.  Skip mtriage, show stages 2 & 3."""
    print("Test: scenario_sharpen_and_density_modify_cryo")
    skipped, h = _simulate_init(
        _cryo_plan(),
        {"stop_conditions": {
            "start_with_program": "phenix.map_sharpening"}})
    assert skipped == 1
    assert "Stage 1" not in h
    assert "Stage 2" in h
    assert "Stage 3" in h
    assert "Starting here" in h
    assert "Will stop here" not in h
    print("  PASS")


def test_gemini_concern_xtriage_refine_stop():
    """Gemini concern #1: 'xtriage, refine, stop' must
    NOT skip xtriage.  The resolver sets start_with=
    xtriage (first action) and after=refine (last action);
    skip_target prioritizes start_with → skip to xtriage
    = Stage 1 → no skip."""
    print("Test: gemini_concern_xtriage_refine_stop")
    skipped, h = _simulate_init(
        _predict_plan(),
        {"stop_conditions": {
            "start_with_program": "phenix.xtriage",
            "after_program": "phenix.refine"}})
    assert skipped == 0, \
        "Expected NO skip when start_with=xtriage"
    assert "Stage 1" in h, "xtriage must still be shown"
    assert "Starting here" not in h, \
        "No skip → no 'Starting here' annotation"
    assert "Will stop here" in h
    print("  PASS")


# =====================================================================
# SECTION E: backward compatibility
# =====================================================================

def test_backward_compat_hasattr_guard():
    """The hasattr guard pattern correctly bypasses the
    call when plan_schema.py is the old version that
    lacks skip_to_program."""
    print("Test: backward_compat_hasattr_guard")

    class OldStylePlan(object):
        """Simulates an older StructurePlan."""
        def __init__(self):
            self.goal = "old plan"
            self.stages = []

        def format_plan_header(self, stop_after=None):
            return "(old format)"

    old = OldStylePlan()
    assert not hasattr(old, "skip_to_program"), \
        "Old plan should lack skip_to_program"

    # Replicate the guarded call from ai_agent.py
    directives = {"stop_conditions": {
        "after_program": "phenix.predict_and_build"}}
    sc = directives.get("stop_conditions", {})
    skip_target = (
        sc.get("start_with_program")
        or sc.get("after_program"))
    skipped = 0
    if skip_target and hasattr(old, "skip_to_program"):
        # Should NOT reach here
        skipped = old.skip_to_program(skip_target)
        raise AssertionError("guard failed!")

    assert skipped == 0
    # Old plan can still render its (old-style) header
    header = old.format_plan_header()
    assert header == "(old format)"
    print("  PASS")


def test_backward_compat_stage_skipped_constant():
    """STAGE_SKIPPED constant existed before this change
    (line 32 of plan_schema.py).  Confirming this means
    a new plan_schema.py with STAGE_SKIPPED stages can be
    serialized and deserialized correctly by old code
    paths (the constant is just the string 'skipped')."""
    print("Test: backward_compat_stage_skipped_constant")
    assert STAGE_SKIPPED == "skipped", \
        "STAGE_SKIPPED value changed; would break " \
        "cross-version session compatibility"
    print("  PASS")


# =====================================================================
# SECTION F: edge cases
# =====================================================================

def test_edge_empty_plan():
    """skip_to_program on an empty plan returns 0
    gracefully."""
    print("Test: edge_empty_plan")
    p = StructurePlan(goal="empty", stages=[])
    n = p.skip_to_program("phenix.anything")
    assert n == 0
    print("  PASS")


def test_edge_single_stage_plan():
    """skip_to_program on a single-stage plan returns 0
    (target is already current or absent)."""
    print("Test: edge_single_stage_plan")
    p = StructurePlan(
        goal="single",
        stages=[StageDef(id="only", description="Only stage",
                         programs=["phenix.refine"])])
    # Target present → no-op (target IS current)
    n = p.skip_to_program("phenix.refine")
    assert n == 0
    assert p.current_stage_index == 0
    # Target absent → no-op
    n = p.skip_to_program("phenix.does_not_exist")
    assert n == 0
    print("  PASS")


def test_edge_advance_after_skip_works():
    """advance() still operates correctly when stages
    have been marked SKIPPED by skip_to_program.

    Verifies that the gate evaluator's normal flow
    is unaffected — after running the target stage,
    advance() should move to the NEXT pending stage
    (skipping over any SKIPPED stage in between if
    present)."""
    print("Test: edge_advance_after_skip_works")
    p = _predict_plan()
    p.skip_to_program("phenix.predict_and_build")
    # We're at stage 1 (predict)
    assert p.current_stage_index == 1
    # Mark current as complete and advance
    p.stages[1].status = STAGE_COMPLETE
    p.advance()
    # Should move to stage 2 (refine), NOT back to
    # stage 0 (which is SKIPPED)
    assert p.current_stage_index == 2, \
        "advance() should move to stage 2, got %d" % \
        p.current_stage_index
    # advance() marks the new current stage ACTIVE
    assert p.stages[2].status == STAGE_ACTIVE, \
        "advance() should mark stage 2 ACTIVE, got %s" % \
        p.stages[2].status
    # SKIPPED stage 0 unchanged
    assert p.stages[0].status == STAGE_SKIPPED
    print("  PASS")


def test_edge_skip_with_no_stop_conditions():
    """When stop_conditions is empty/missing, simulate_init
    correctly performs no skip and no annotations.

    This is the most common case — most user requests
    do not set start_with_program or after_program."""
    print("Test: edge_skip_with_no_stop_conditions")
    skipped, h = _simulate_init(
        _predict_plan(), {})
    assert skipped == 0
    assert "Starting here" not in h
    assert "Will stop here" not in h
    print("  PASS")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    """Run all tests using cctbx fail-fast convention."""
    try:
        from libtbx.langchain.tests.tst_utils import (
            run_tests_with_fail_fast)
    except ImportError:
        try:
            from tests.tst_utils import (
                run_tests_with_fail_fast)
        except ImportError:
            # Fall through to standalone runner
            _standalone_runner()
            return
    run_tests_with_fail_fast()


def _standalone_runner():
    """Simple runner when tst_utils is not available."""
    test_fns = [v for k, v in sorted(globals().items())
                if k.startswith("test_") and callable(v)]
    passed = 0
    failed = 0
    for fn in test_fns:
        try:
            fn()
            passed += 1
        except Exception as e:
            print("  FAIL: %s" % e)
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    _standalone_runner()
