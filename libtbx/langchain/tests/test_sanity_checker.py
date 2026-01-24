#!/usr/bin/env python
"""
Unit tests for sanity_checker.py

Tests:
- Critical checks: experiment type change, no model for refine, no data, repeated failures
- Warning checks: resolution unknown, R-free spike, map CC drop
- Abort settings: abort_on_red_flags, abort_on_warnings
- Message formatting
- Edge cases
"""

from __future__ import absolute_import, division, print_function

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from libtbx.langchain.agent.sanity_checker import (
        SanityChecker,
        SanityIssue,
        SanityResult,
    )
except ImportError:
    from agent.sanity_checker import (
        SanityChecker,
        SanityIssue,
        SanityResult,
    )


# =============================================================================
# TEST FIXTURES
# =============================================================================

def make_context(experiment_type="xray", has_model=True, has_mtz=True,
                 has_map=False, state="xray_refine", history=None,
                 metrics_history=None, resolution=2.5):
    """Create a test context dict."""
    return {
        "experiment_type": experiment_type,
        "has_model": has_model,
        "has_mtz": has_mtz,
        "has_map": has_map,
        "state": state,
        "history": history or [],
        "metrics_history": metrics_history or [],
        "resolution": resolution,
    }


def make_session_info(experiment_type=None):
    """Create a test session_info dict."""
    return {"experiment_type": experiment_type}


# =============================================================================
# CRITICAL CHECK TESTS
# =============================================================================

def test_experiment_type_stable_no_change():
    """Test that no issue is raised when experiment type is stable."""
    print("Test: experiment_type_stable_no_change")

    checker = SanityChecker()
    context = make_context(experiment_type="xray")
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    assert result.ok, "Should be ok when experiment type is stable"
    assert not result.should_abort, "Should not abort"
    assert not any(i.code == "experiment_type_changed" for i in result.issues)

    print("  PASSED")


def test_experiment_type_changed():
    """Test that experiment type change is detected."""
    print("Test: experiment_type_changed")

    checker = SanityChecker()
    context = make_context(experiment_type="cryoem", has_mtz=False, has_map=True)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info, abort_on_red_flags=True)

    assert not result.ok, "Should not be ok when experiment type changed"
    assert result.should_abort, "Should abort on experiment type change"
    assert any(i.code == "experiment_type_changed" for i in result.issues)

    # Check the issue details
    issue = next(i for i in result.issues if i.code == "experiment_type_changed")
    assert issue.severity == "critical"
    assert "xray" in issue.message
    assert "cryoem" in issue.message

    print("  PASSED")


def test_experiment_type_not_set_yet():
    """Test that no issue when experiment type not yet locked."""
    print("Test: experiment_type_not_set_yet")

    checker = SanityChecker()
    context = make_context(experiment_type="xray")
    session_info = make_session_info(experiment_type=None)  # Not set yet

    result = checker.check(context, session_info)

    assert result.ok, "Should be ok when experiment type not yet locked"
    assert not any(i.code == "experiment_type_changed" for i in result.issues)

    print("  PASSED")


def test_no_model_for_refine():
    """Test that missing model for refinement is detected."""
    print("Test: no_model_for_refine")

    checker = SanityChecker()
    context = make_context(has_model=False, state="xray_refine")
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    assert not result.ok, "Should not be ok when no model for refine"
    assert result.should_abort
    assert any(i.code == "no_model_for_refine" for i in result.issues)

    print("  PASSED")


def test_model_exists_for_refine():
    """Test that no issue when model exists for refinement."""
    print("Test: model_exists_for_refine")

    checker = SanityChecker()
    context = make_context(has_model=True, state="xray_refine")
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    assert not any(i.code == "no_model_for_refine" for i in result.issues)

    print("  PASSED")


def test_no_model_non_refine_state():
    """Test that missing model is ok in non-refine states."""
    print("Test: no_model_non_refine_state")

    checker = SanityChecker()
    context = make_context(has_model=False, state="xray_initial")
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    assert not any(i.code == "no_model_for_refine" for i in result.issues)

    print("  PASSED")


def test_no_mtz_for_xray():
    """Test that missing MTZ for X-ray workflow is detected."""
    print("Test: no_mtz_for_xray")

    checker = SanityChecker()
    context = make_context(experiment_type="xray", has_mtz=False, state="xray_refine")
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    assert not result.ok
    assert any(i.code == "no_data_for_workflow" for i in result.issues)

    print("  PASSED")


def test_no_map_for_cryoem():
    """Test that missing map for cryo-EM workflow is detected."""
    print("Test: no_map_for_cryoem")

    checker = SanityChecker()
    context = make_context(
        experiment_type="cryoem",
        has_mtz=False,
        has_map=False,
        state="cryoem_refine"
    )
    session_info = make_session_info(experiment_type="cryoem")

    result = checker.check(context, session_info)

    assert not result.ok
    assert any(i.code == "no_data_for_workflow" for i in result.issues)

    print("  PASSED")


def test_repeated_failures():
    """Test that repeated failures are detected."""
    print("Test: repeated_failures")

    checker = SanityChecker()

    # Create history with 3 identical failures
    history = [
        {"program": "phenix.refine", "result": "FAILED: Some error message", "cycle": 1},
        {"program": "phenix.refine", "result": "FAILED: Some error message", "cycle": 2},
        {"program": "phenix.refine", "result": "FAILED: Some error message", "cycle": 3},
    ]

    context = make_context(history=history)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    assert not result.ok
    assert any(i.code == "repeated_failures" for i in result.issues)

    # Check the issue details
    issue = next(i for i in result.issues if i.code == "repeated_failures")
    assert "phenix.refine" in issue.message
    assert "3" in issue.message or "three" in issue.message.lower()

    print("  PASSED")


def test_no_repeated_failures_success_with_errors_word():
    """Test that 'SUCCESS: Command completed without errors' is not a failure."""
    print("Test: no_repeated_failures_success_with_errors_word")

    checker = SanityChecker()

    # Create history with successful results that contain the word "errors"
    history = [
        {"program": "phenix.refine", "result": "SUCCESS: Command completed without errors", "cycle": 1},
        {"program": "phenix.refine", "result": "SUCCESS: Command completed without errors", "cycle": 2},
        {"program": "phenix.refine", "result": "SUCCESS: Command completed without errors", "cycle": 3},
    ]

    context = make_context(history=history)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    # Should NOT detect repeated failures
    assert not any(i.code == "repeated_failures" for i in result.issues)

    print("  PASSED")


def test_no_repeated_failures_different_errors():
    """Test that different errors don't trigger repeated failures."""
    print("Test: no_repeated_failures_different_errors")

    checker = SanityChecker()

    history = [
        {"program": "phenix.refine", "result": "FAILED: Error A", "cycle": 1},
        {"program": "phenix.refine", "result": "FAILED: Error B", "cycle": 2},
        {"program": "phenix.refine", "result": "FAILED: Error C", "cycle": 3},
    ]

    context = make_context(history=history)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    assert not any(i.code == "repeated_failures" for i in result.issues)

    print("  PASSED")


def test_no_repeated_failures_two_times():
    """Test that 2 failures don't trigger (need 3+)."""
    print("Test: no_repeated_failures_two_times")

    checker = SanityChecker()

    history = [
        {"program": "phenix.refine", "result": "FAILED: Same error", "cycle": 1},
        {"program": "phenix.refine", "result": "FAILED: Same error", "cycle": 2},
    ]

    context = make_context(history=history)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    assert not any(i.code == "repeated_failures" for i in result.issues)

    print("  PASSED")


# =============================================================================
# WARNING CHECK TESTS
# =============================================================================

def test_resolution_unknown_warning():
    """Test that missing resolution in refine state triggers warning."""
    print("Test: resolution_unknown_warning")

    checker = SanityChecker()
    context = make_context(state="xray_refine", resolution=None)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    # Should be a warning, not critical
    assert result.ok, "Should be ok (warning, not critical)"
    assert any(i.code == "resolution_unknown" for i in result.issues)

    issue = next(i for i in result.issues if i.code == "resolution_unknown")
    assert issue.severity == "warning"

    print("  PASSED")


def test_r_free_spike_warning():
    """Test that R-free spike triggers warning."""
    print("Test: r_free_spike_warning")

    checker = SanityChecker()

    metrics_history = [
        {"cycle": 1, "r_free": 0.30},
        {"cycle": 2, "r_free": 0.50},  # Jumped by 0.20!
    ]

    context = make_context(metrics_history=metrics_history)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    assert any(i.code == "r_free_spike" for i in result.issues)

    issue = next(i for i in result.issues if i.code == "r_free_spike")
    assert issue.severity == "warning"
    assert "0.30" in issue.message or "0.3" in issue.message
    assert "0.50" in issue.message or "0.5" in issue.message

    print("  PASSED")


def test_no_r_free_spike_normal_change():
    """Test that normal R-free change doesn't trigger warning."""
    print("Test: no_r_free_spike_normal_change")

    checker = SanityChecker()

    metrics_history = [
        {"cycle": 1, "r_free": 0.30},
        {"cycle": 2, "r_free": 0.32},  # Small increase (0.02)
    ]

    context = make_context(metrics_history=metrics_history)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    assert not any(i.code == "r_free_spike" for i in result.issues)

    print("  PASSED")


def test_map_cc_drop_warning():
    """Test that map CC drop triggers warning."""
    print("Test: map_cc_drop_warning")

    checker = SanityChecker()

    metrics_history = [
        {"cycle": 1, "map_cc": 0.85},
        {"cycle": 2, "map_cc": 0.45},  # Dropped by 0.40!
    ]

    context = make_context(
        experiment_type="cryoem",
        has_mtz=False,
        has_map=True,
        metrics_history=metrics_history
    )
    session_info = make_session_info(experiment_type="cryoem")

    result = checker.check(context, session_info)

    assert any(i.code == "map_cc_drop" for i in result.issues)

    issue = next(i for i in result.issues if i.code == "map_cc_drop")
    assert issue.severity == "warning"

    print("  PASSED")


# =============================================================================
# ABORT SETTINGS TESTS
# =============================================================================

def test_abort_on_red_flags_true():
    """Test that critical issues cause abort when abort_on_red_flags=True."""
    print("Test: abort_on_red_flags_true")

    checker = SanityChecker()
    context = make_context(experiment_type="cryoem", has_mtz=False, has_map=True)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info, abort_on_red_flags=True)

    assert result.should_abort, "Should abort when abort_on_red_flags=True"

    print("  PASSED")


def test_abort_on_red_flags_false():
    """Test that critical issues don't abort when abort_on_red_flags=False."""
    print("Test: abort_on_red_flags_false")

    checker = SanityChecker()
    context = make_context(experiment_type="cryoem", has_mtz=False, has_map=True)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info, abort_on_red_flags=False)

    assert not result.should_abort, "Should not abort when abort_on_red_flags=False"
    # But issues should still be detected
    assert len(result.issues) > 0

    print("  PASSED")


def test_abort_on_warnings_false():
    """Test that warnings don't abort when abort_on_warnings=False."""
    print("Test: abort_on_warnings_false")

    checker = SanityChecker()

    # Create context with only warning-level issue
    metrics_history = [
        {"cycle": 1, "r_free": 0.30},
        {"cycle": 2, "r_free": 0.50},  # R-free spike
    ]
    context = make_context(metrics_history=metrics_history)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info,
                          abort_on_red_flags=True,
                          abort_on_warnings=False)

    # Should have warning but not abort
    assert any(i.code == "r_free_spike" for i in result.issues)
    assert not result.should_abort, "Should not abort on warnings when abort_on_warnings=False"

    print("  PASSED")


def test_abort_on_warnings_true():
    """Test that warnings cause abort when abort_on_warnings=True."""
    print("Test: abort_on_warnings_true")

    checker = SanityChecker()

    metrics_history = [
        {"cycle": 1, "r_free": 0.30},
        {"cycle": 2, "r_free": 0.50},  # R-free spike
    ]
    context = make_context(metrics_history=metrics_history)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info,
                          abort_on_red_flags=True,
                          abort_on_warnings=True)

    assert result.should_abort, "Should abort on warnings when abort_on_warnings=True"

    print("  PASSED")


# =============================================================================
# MESSAGE FORMATTING TESTS
# =============================================================================

def test_abort_message_formatting():
    """Test that abort message is properly formatted."""
    print("Test: abort_message_formatting")

    checker = SanityChecker()
    context = make_context(experiment_type="cryoem", has_mtz=False, has_map=True)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info, abort_on_red_flags=True)

    assert result.abort_message is not None
    assert "AGENT STOPPED" in result.abort_message
    assert "Sanity check failed" in result.abort_message
    assert "CRITICAL ISSUES" in result.abort_message
    assert "experiment_type_changed" in result.abort_message
    assert "TO RESUME" in result.abort_message

    print("  PASSED")


def test_abort_message_includes_suggestions():
    """Test that abort message includes suggestions."""
    print("Test: abort_message_includes_suggestions")

    checker = SanityChecker()
    context = make_context(has_model=False, state="xray_refine")
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    assert "→" in result.abort_message, "Should include suggestion indicator"

    print("  PASSED")


# =============================================================================
# EDGE CASES
# =============================================================================

def test_empty_context():
    """Test with minimal/empty context."""
    print("Test: empty_context")

    checker = SanityChecker()
    context = {}
    session_info = {}

    # Should not crash
    result = checker.check(context, session_info)
    assert isinstance(result, SanityResult)

    print("  PASSED")


def test_empty_history():
    """Test with empty history."""
    print("Test: empty_history")

    checker = SanityChecker()
    context = make_context(history=[])
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    # Should not have repeated_failures issue
    assert not any(i.code == "repeated_failures" for i in result.issues)

    print("  PASSED")


def test_empty_metrics_history():
    """Test with empty metrics history."""
    print("Test: empty_metrics_history")

    checker = SanityChecker()
    context = make_context(metrics_history=[])
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    # Should not have metric anomaly issues
    assert not any(i.code == "r_free_spike" for i in result.issues)
    assert not any(i.code == "map_cc_drop" for i in result.issues)

    print("  PASSED")


def test_single_metrics_entry():
    """Test with single metrics entry (no previous to compare)."""
    print("Test: single_metrics_entry")

    checker = SanityChecker()

    metrics_history = [
        {"cycle": 1, "r_free": 0.50},  # High but no previous to compare
    ]
    context = make_context(metrics_history=metrics_history)
    session_info = make_session_info(experiment_type="xray")

    result = checker.check(context, session_info)

    # Should not have spike/drop issues (need 2 entries to compare)
    assert not any(i.code == "r_free_spike" for i in result.issues)
    assert not any(i.code == "map_cc_drop" for i in result.issues)

    print("  PASSED")


def test_sanity_issue_to_dict():
    """Test SanityIssue.to_dict() method."""
    print("Test: sanity_issue_to_dict")

    issue = SanityIssue(
        severity="critical",
        code="test_code",
        message="Test message",
        suggestion="Test suggestion",
        details={"key": "value"}
    )

    d = issue.to_dict()

    assert d["severity"] == "critical"
    assert d["code"] == "test_code"
    assert d["message"] == "Test message"
    assert d["suggestion"] == "Test suggestion"
    assert d["details"] == {"key": "value"}

    print("  PASSED")


def test_sanity_result_to_dict():
    """Test SanityResult.to_dict() method."""
    print("Test: sanity_result_to_dict")

    issue = SanityIssue(
        severity="critical",
        code="test_code",
        message="Test message",
        suggestion="Test suggestion"
    )

    result = SanityResult(
        ok=False,
        issues=[issue],
        should_abort=True,
        abort_message="Test abort message"
    )

    d = result.to_dict()

    assert d["ok"] == False
    assert d["should_abort"] == True
    assert d["abort_message"] == "Test abort message"
    assert len(d["issues"]) == 1
    assert d["issues"][0]["code"] == "test_code"

    print("  PASSED")


def test_multiple_issues():
    """Test that multiple issues are all captured."""
    print("Test: multiple_issues")

    checker = SanityChecker()

    # Create context with multiple issues
    metrics_history = [
        {"cycle": 1, "r_free": 0.30},
        {"cycle": 2, "r_free": 0.50},  # R-free spike
    ]

    context = make_context(
        experiment_type="cryoem",
        has_mtz=False,
        has_map=True,
        state="cryoem_refine",
        resolution=None,  # Resolution unknown
        metrics_history=metrics_history
    )
    session_info = make_session_info(experiment_type="xray")  # Type changed

    result = checker.check(context, session_info)

    # Should have multiple issues
    codes = [i.code for i in result.issues]
    assert "experiment_type_changed" in codes
    assert "resolution_unknown" in codes

    print("  PASSED")


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all sanity checker tests."""
    print("\n" + "=" * 60)
    print("SANITY CHECKER TESTS")
    print("=" * 60 + "\n")

    # Critical check tests
    print("\n--- Critical Check Tests ---")
    test_experiment_type_stable_no_change()
    test_experiment_type_changed()
    test_experiment_type_not_set_yet()
    test_no_model_for_refine()
    test_model_exists_for_refine()
    test_no_model_non_refine_state()
    test_no_mtz_for_xray()
    test_no_map_for_cryoem()
    test_repeated_failures()
    test_no_repeated_failures_different_errors()
    test_no_repeated_failures_two_times()
    test_no_repeated_failures_success_with_errors_word()

    # Warning check tests
    print("\n--- Warning Check Tests ---")
    test_resolution_unknown_warning()
    test_r_free_spike_warning()
    test_no_r_free_spike_normal_change()
    test_map_cc_drop_warning()

    # Abort settings tests
    print("\n--- Abort Settings Tests ---")
    test_abort_on_red_flags_true()
    test_abort_on_red_flags_false()
    test_abort_on_warnings_false()
    test_abort_on_warnings_true()

    # Message formatting tests
    print("\n--- Message Formatting Tests ---")
    test_abort_message_formatting()
    test_abort_message_includes_suggestions()

    # Edge cases
    print("\n--- Edge Case Tests ---")
    test_empty_context()
    test_empty_history()
    test_empty_metrics_history()
    test_single_metrics_entry()
    test_sanity_issue_to_dict()
    test_sanity_result_to_dict()
    test_multiple_issues()

    print("\n" + "=" * 60)
    print("ALL SANITY CHECKER TESTS PASSED")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    run_all_tests()
