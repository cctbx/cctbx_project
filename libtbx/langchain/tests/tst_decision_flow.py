"""
Tests for the clean decision flow architecture.

These tests verify that:
1. Directives properly modify valid_programs (workflow_engine)
2. Validation gate is a simple "in list?" check (graph_nodes)
3. Stop conditions are checked only post-execution (ai_agent)

Run with: python tests/tst_decision_flow.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_false, assert_in, assert_not_in,
    assert_greater,
    run_tests_with_fail_fast
)

# Try to import workflow_engine - may fail without libtbx
try:
    from agent.workflow_engine import WorkflowEngine
    HAS_WORKFLOW_ENGINE = True
except ImportError:
    HAS_WORKFLOW_ENGINE = False
    WorkflowEngine = None

# Try to import directive functions
try:
    from agent.directive_extractor import check_stop_conditions
    HAS_DIRECTIVE_EXTRACTOR = True
except ImportError:
    HAS_DIRECTIVE_EXTRACTOR = False
    check_stop_conditions = None

# Import validate_intent
try:
    from agent.directive_validator import validate_intent
    HAS_VALIDATE_INTENT = True
except ImportError:
    HAS_VALIDATE_INTENT = False
    validate_intent = None


# =============================================================================
# WORKFLOW ENGINE DIRECTIVE TESTS (require libtbx)
# =============================================================================

def test_skip_validation_adds_stop():
    """skip_validation=true should add STOP to valid_programs."""
    if not HAS_WORKFLOW_ENGINE:
        print("  SKIPPED (requires libtbx)")
        return

    engine = WorkflowEngine()
    valid_programs = ["phenix.refine", "phenix.molprobity"]
    directives = {
        "stop_conditions": {"skip_validation": True}
    }

    result = engine._apply_directives(valid_programs, directives, "refine")
    assert_in("STOP", result)


def test_skip_validation_no_duplicate_stop():
    """Should not add duplicate STOP if already present."""
    if not HAS_WORKFLOW_ENGINE:
        print("  SKIPPED (requires libtbx)")
        return

    engine = WorkflowEngine()
    valid_programs = ["phenix.refine", "STOP"]
    directives = {
        "stop_conditions": {"skip_validation": True}
    }

    result = engine._apply_directives(valid_programs, directives, "refine")
    assert_equal(result.count("STOP"), 1)


def test_after_program_adds_target():
    """after_program should add the target program to valid_programs."""
    if not HAS_WORKFLOW_ENGINE:
        print("  SKIPPED (requires libtbx)")
        return

    engine = WorkflowEngine()
    valid_programs = ["phenix.predict_and_build"]
    directives = {
        "stop_conditions": {
            "after_program": "phenix.resolve_cryo_em",
            "skip_validation": True
        }
    }

    result = engine._apply_directives(valid_programs, directives, "obtain_model")
    assert_in("phenix.resolve_cryo_em", result)
    # Should be at front (high priority)
    assert_equal(result[0], "phenix.resolve_cryo_em")


def test_after_program_no_duplicate():
    """Should not add duplicate if target already in list."""
    if not HAS_WORKFLOW_ENGINE:
        print("  SKIPPED (requires libtbx)")
        return

    engine = WorkflowEngine()
    valid_programs = ["phenix.resolve_cryo_em", "phenix.predict_and_build"]
    directives = {
        "stop_conditions": {
            "after_program": "phenix.resolve_cryo_em"
        }
    }

    result = engine._apply_directives(valid_programs, directives, "obtain_model")
    assert_equal(result.count("phenix.resolve_cryo_em"), 1)


def test_skip_programs():
    """skip_programs should remove programs from list."""
    if not HAS_WORKFLOW_ENGINE:
        print("  SKIPPED (requires libtbx)")
        return

    engine = WorkflowEngine()
    valid_programs = ["phenix.refine", "phenix.molprobity", "phenix.autobuild"]
    directives = {
        "workflow_preferences": {
            "skip_programs": ["phenix.autobuild"]
        }
    }

    result = engine._apply_directives(valid_programs, directives, "refine")
    assert_not_in("phenix.autobuild", result)
    assert_in("phenix.refine", result)
    assert_in("phenix.molprobity", result)


def test_prefer_programs():
    """prefer_programs should move programs to front."""
    if not HAS_WORKFLOW_ENGINE:
        print("  SKIPPED (requires libtbx)")
        return

    engine = WorkflowEngine()
    valid_programs = ["phenix.refine", "phenix.molprobity", "phenix.autobuild"]
    directives = {
        "workflow_preferences": {
            "prefer_programs": ["phenix.autobuild"]
        }
    }

    result = engine._apply_directives(valid_programs, directives, "refine")
    assert_equal(result[0], "phenix.autobuild")


def test_empty_directives():
    """Empty directives should return list unchanged."""
    if not HAS_WORKFLOW_ENGINE:
        print("  SKIPPED (requires libtbx)")
        return

    engine = WorkflowEngine()
    valid_programs = ["phenix.refine", "phenix.molprobity"]

    result = engine._apply_directives(valid_programs, {}, "refine")
    assert_equal(result, valid_programs)


def test_none_directives():
    """None directives should return list unchanged."""
    if not HAS_WORKFLOW_ENGINE:
        print("  SKIPPED (requires libtbx)")
        return

    engine = WorkflowEngine()
    valid_programs = ["phenix.refine", "phenix.molprobity"]

    result = engine._apply_directives(valid_programs, None, "refine")
    assert_equal(result, valid_programs)


# =============================================================================
# VALIDATE INTENT TESTS
# =============================================================================

def test_no_should_stop_in_result():
    """validate_intent should not return should_stop (stop conditions checked post-execution)."""
    if not HAS_VALIDATE_INTENT:
        print("  SKIPPED (validate_intent not available)")
        return

    intent = {"program": "phenix.refine", "strategy": {}}
    directives = {
        "stop_conditions": {"after_cycle": 1}
    }

    result = validate_intent(intent, directives, cycle_number=1)
    assert_not_in("should_stop", result)
    assert_not_in("stop_reason", result)


def test_still_applies_program_settings():
    """validate_intent should still apply program settings."""
    if not HAS_VALIDATE_INTENT:
        print("  SKIPPED (validate_intent not available)")
        return

    intent = {"program": "phenix.refine", "strategy": {}}
    directives = {
        "program_settings": {
            "default": {"resolution": 2.5}
        }
    }

    result = validate_intent(intent, directives, cycle_number=1)
    assert_equal(result["validated_intent"]["strategy"]["resolution"], 2.5)


def test_returns_modifications():
    """validate_intent should still return modifications list."""
    if not HAS_VALIDATE_INTENT:
        print("  SKIPPED (validate_intent not available)")
        return

    intent = {"program": "phenix.refine", "strategy": {}}
    directives = {
        "program_settings": {
            "default": {"resolution": 2.5}
        }
    }

    result = validate_intent(intent, directives, cycle_number=1)
    assert_in("modifications", result)
    assert_greater(len(result["modifications"]), 0)


# =============================================================================
# POST-EXECUTION STOP CHECK TESTS
# =============================================================================

def test_after_program_triggers_stop():
    """after_program should trigger stop after program completes."""
    if not HAS_DIRECTIVE_EXTRACTOR:
        print("  SKIPPED (check_stop_conditions not available)")
        return

    directives = {
        "stop_conditions": {"after_program": "phenix.resolve_cryo_em"}
    }

    should_stop, reason = check_stop_conditions(
        directives,
        cycle_number=2,
        last_program="phenix.resolve_cryo_em"
    )
    assert_true(should_stop)
    assert_in("resolve_cryo_em", reason)


def test_after_program_normalized_match():
    """Should match with or without phenix. prefix."""
    if not HAS_DIRECTIVE_EXTRACTOR:
        print("  SKIPPED (check_stop_conditions not available)")
        return

    directives = {
        "stop_conditions": {"after_program": "phenix.xtriage"}
    }

    should_stop, reason = check_stop_conditions(
        directives,
        cycle_number=1,
        last_program="xtriage"
    )
    assert_true(should_stop)


def test_after_program_no_match():
    """Should not stop if different program ran."""
    if not HAS_DIRECTIVE_EXTRACTOR:
        print("  SKIPPED (check_stop_conditions not available)")
        return

    directives = {
        "stop_conditions": {"after_program": "phenix.resolve_cryo_em"}
    }

    should_stop, reason = check_stop_conditions(
        directives,
        cycle_number=2,
        last_program="phenix.mtriage"
    )
    assert_false(should_stop)


def test_after_cycle_triggers_stop():
    """after_cycle should trigger stop at specified cycle."""
    if not HAS_DIRECTIVE_EXTRACTOR:
        print("  SKIPPED (check_stop_conditions not available)")
        return

    directives = {
        "stop_conditions": {"after_cycle": 3}
    }

    should_stop, reason = check_stop_conditions(
        directives,
        cycle_number=3,
        last_program="phenix.refine"
    )
    assert_true(should_stop)
    assert_in("cycle", reason.lower())


def test_after_cycle_before_target():
    """Should not stop before target cycle."""
    if not HAS_DIRECTIVE_EXTRACTOR:
        print("  SKIPPED (check_stop_conditions not available)")
        return

    directives = {
        "stop_conditions": {"after_cycle": 3}
    }

    should_stop, reason = check_stop_conditions(
        directives,
        cycle_number=2,
        last_program="phenix.refine"
    )
    assert_false(should_stop)


def test_empty_directives_no_stop():
    """Empty directives should not trigger stop."""
    if not HAS_DIRECTIVE_EXTRACTOR:
        print("  SKIPPED (check_stop_conditions not available)")
        return

    should_stop, reason = check_stop_conditions(
        {},
        cycle_number=5,
        last_program="phenix.refine"
    )
    assert_false(should_stop)


def test_no_stop_conditions_no_stop():
    """Directives without stop_conditions should not trigger stop."""
    if not HAS_DIRECTIVE_EXTRACTOR:
        print("  SKIPPED (check_stop_conditions not available)")
        return

    directives = {
        "program_settings": {"default": {"resolution": 2.5}}
    }

    should_stop, reason = check_stop_conditions(
        directives,
        cycle_number=5,
        last_program="phenix.refine"
    )
    assert_false(should_stop)


# =============================================================================
# INTEGRATION TESTS (require libtbx)
# =============================================================================

def test_tutorial_flow_resolve_cryo_em():
    """
    Tutorial: "Run density modification and stop"
    Should add resolve_cryo_em to valid_programs and allow stop after.
    """
    if not HAS_WORKFLOW_ENGINE or not HAS_DIRECTIVE_EXTRACTOR:
        print("  SKIPPED (requires libtbx)")
        return

    engine = WorkflowEngine()

    # Initial valid programs (from workflow phase - might not include resolve_cryo_em)
    valid_programs = ["phenix.predict_and_build"]

    # User directive for tutorial
    directives = {
        "stop_conditions": {
            "after_program": "phenix.resolve_cryo_em",
            "skip_validation": True
        }
    }

    # Step 1: Workflow engine adds required program
    result = engine._apply_directives(valid_programs, directives, "obtain_model")
    assert_in("phenix.resolve_cryo_em", result)
    assert_in("STOP", result)

    # Step 2: LLM can now choose resolve_cryo_em (it's in valid_programs)
    # (graph_nodes validation would pass)

    # Step 3: After resolve_cryo_em completes, check stop condition
    should_stop, reason = check_stop_conditions(
        directives,
        cycle_number=2,
        last_program="phenix.resolve_cryo_em"
    )
    assert_true(should_stop)


def test_normal_workflow_no_directives():
    """Normal workflow without directives should work unchanged."""
    if not HAS_WORKFLOW_ENGINE or not HAS_DIRECTIVE_EXTRACTOR:
        print("  SKIPPED (requires libtbx)")
        return

    engine = WorkflowEngine()
    valid_programs = ["phenix.refine", "phenix.molprobity"]

    # No directives
    result = engine._apply_directives(valid_programs, {}, "refine")

    # Should be unchanged
    assert_equal(result, valid_programs)
    assert_not_in("STOP", result)

    # No stop condition triggered
    should_stop, _ = check_stop_conditions({}, 1, "phenix.refine")
    assert_false(should_stop)


def test_use_experimental_phasing_prioritizes_autosol():
    """Test that use_experimental_phasing prioritizes autosol over predict_and_build."""
    print("Test: use_experimental_phasing_prioritizes_autosol")

    if not HAS_WORKFLOW_ENGINE:
        print("  SKIPPED (requires workflow_engine)")
        return

    engine = WorkflowEngine()

    # Simulate obtain_model phase with both programs available
    valid_programs = ["phenix.predict_and_build", "phenix.autosol", "phenix.phaser"]

    directives = {
        "workflow_preferences": {
            "use_experimental_phasing": True
        }
    }

    result = engine._apply_directives(valid_programs, directives, "obtain_model")

    # autosol should be first (prioritized)
    assert_equal(result[0], "phenix.autosol")
    # predict_and_build should be last (deprioritized)
    assert_equal(result[-1], "phenix.predict_and_build")

    print("  PASSED")


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
