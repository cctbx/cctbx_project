"""
Tests for error_analyzer.py - Recoverable Error Handling.

This module tests the automatic error recovery system that detects
structured errors in PHENIX output and determines appropriate recovery
strategies.

Run with: python tests/test_error_analyzer.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.test_utils import (
    assert_equal, assert_true, assert_false, assert_none, assert_not_none,
    assert_in,
    run_tests_with_fail_fast
)

from agent.error_analyzer import ErrorAnalyzer, ErrorRecovery


# =============================================================================
# HELPER CLASSES
# =============================================================================

class MockSession:
    """Mock session for testing."""

    def __init__(self, data=None):
        self.data = data or {}

    def save(self):
        pass


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def create_analyzer():
    """Create an ErrorAnalyzer instance."""
    return ErrorAnalyzer()


def get_error_def(analyzer, error_type="ambiguous_data_labels"):
    """Get error definition from config."""
    return analyzer._config.get("errors", {}).get(error_type, {})


# =============================================================================
# ERROR DETECTION TESTS
# =============================================================================

def test_detect_ambiguous_data_error_standard():
    """Test detection of standard ambiguous data labels error."""
    analyzer = create_analyzer()
    log = """
    Multiple equally suitable arrays of observed xray data found.

    Possible choices:
      /path/data.mtz:IMEAN_CuKa,SIGIMEAN_CuKa
      /path/data.mtz:I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)

    Please use scaling.input.xray_data.obs_labels
    to specify an unambiguous substring of the target label.
    """

    error_type = analyzer._detect_error_type(log)
    assert_equal(error_type, "ambiguous_data_labels")


def test_detect_ambiguous_data_error_multiline_keyword():
    """Test detection when keyword is on separate line."""
    analyzer = create_analyzer()
    log = """
    Multiple equally suitable arrays found.

    Please use
    obs_labels
    to specify an unambiguous substring.
    """

    error_type = analyzer._detect_error_type(log)
    assert_equal(error_type, "ambiguous_data_labels")


def test_no_error_detected_success_log():
    """Test that successful logs don't trigger false positives."""
    analyzer = create_analyzer()
    log = """
    phenix.autosol completed successfully.
    R-factor = 0.20
    Final model contains 150 residues.
    """

    error_type = analyzer._detect_error_type(log)
    assert_none(error_type)


def test_no_error_detected_unrelated_failure():
    """Test that unrelated failures don't match."""
    analyzer = create_analyzer()
    log = """
    FATAL ERROR: File not found: model.pdb
    Program terminated with error code 1.
    """

    error_type = analyzer._detect_error_type(log)
    assert_none(error_type)


# =============================================================================
# KEYWORD EXTRACTION TESTS
# =============================================================================

def test_extract_full_keyword():
    """Test extraction of full parameter path."""
    analyzer = create_analyzer()
    error_def = get_error_def(analyzer)
    log = """
    Multiple equally suitable arrays found.
    Possible choices:
      /data/test.mtz:IMEAN,SIGIMEAN
    Please use scaling.input.xray_data.obs_labels
    to specify an unambiguous substring.
    """

    info = analyzer._extract_ambiguous_labels_info(log, error_def)

    assert_not_none(info)
    assert_equal(info.get("keyword"), "scaling.input.xray_data.obs_labels")


def test_extract_simple_keyword():
    """Test extraction of simple keyword."""
    analyzer = create_analyzer()
    error_def = get_error_def(analyzer)
    log = """
    Ambiguous data labels detected.
    Possible choices:
      /data/test.mtz:IMEAN,SIGIMEAN
    Please use obs_labels to specify.
    """

    info = analyzer._extract_ambiguous_labels_info(log, error_def)

    assert_not_none(info)
    assert_in("keyword", info)


# =============================================================================
# CHOICE EXTRACTION TESTS
# =============================================================================

def test_extract_two_choices():
    """Test extraction of two label choices."""
    analyzer = create_analyzer()
    error_def = get_error_def(analyzer)
    log = """
    Possible choices:
      /data/test.mtz:IMEAN_CuKa,SIGIMEAN_CuKa
      /data/test.mtz:I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)

    Please use obs_labels to specify.
    """

    info = analyzer._extract_ambiguous_labels_info(log, error_def)

    assert_not_none(info)
    assert_equal(len(info["choices"]), 2)
    assert_in("IMEAN_CuKa,SIGIMEAN_CuKa", info["choices"])
    assert_in("I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)", info["choices"])


def test_extract_affected_file():
    """Test extraction of affected file path."""
    analyzer = create_analyzer()
    error_def = get_error_def(analyzer)
    log = """
    Possible choices:
      /data/experiment/lyso.mtz:IMEAN,SIGIMEAN
      /data/experiment/lyso.mtz:I(+),SIGI(+)
    """

    info = analyzer._extract_ambiguous_labels_info(log, error_def)

    assert_not_none(info)
    assert_equal(info["affected_file"], "/data/experiment/lyso.mtz")


def test_extract_no_choices_returns_none():
    """Test that no choices returns None."""
    analyzer = create_analyzer()
    error_def = get_error_def(analyzer)
    log = """
    Some error without choice lines.
    Please use obs_labels to specify.
    """

    info = analyzer._extract_ambiguous_labels_info(log, error_def)
    assert_none(info)


# =============================================================================
# LABEL CLASSIFICATION TESTS
# =============================================================================

def test_anomalous_label_plus_minus():
    """Test detection of (+)/(-) anomalous labels."""
    analyzer = create_analyzer()
    labels = [
        "I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)",
        "F(+),SIGF(+),F(-),SIGF(-)",
        "I(+),SIGI(+)",
    ]

    for label in labels:
        assert_true(analyzer._is_anomalous_label(label),
            f"Expected {label} to be detected as anomalous")


def test_anomalous_label_dano():
    """Test detection of DANO labels."""
    analyzer = create_analyzer()
    assert_true(analyzer._is_anomalous_label("DANO,SIGDANO"))


def test_anomalous_label_anom_suffix():
    """Test detection of *anom labels."""
    analyzer = create_analyzer()
    assert_true(analyzer._is_anomalous_label("I_anom,SIGI_anom"))
    assert_true(analyzer._is_anomalous_label("Fanom,SIGFanom"))


def test_merged_label_imean():
    """Test detection of IMEAN merged labels."""
    analyzer = create_analyzer()
    merged_labels = [
        "IMEAN_CuKa,SIGIMEAN_CuKa",
        "IMEAN,SIGIMEAN",
    ]

    for label in merged_labels:
        assert_true(analyzer._is_merged_label(label),
            f"Expected {label} to be detected as merged")


def test_merged_label_fmean():
    """Test detection of FMEAN merged labels."""
    analyzer = create_analyzer()
    assert_true(analyzer._is_merged_label("FMEAN,SIGFMEAN"))


def test_merged_label_fobs():
    """Test detection of F_obs merged labels."""
    analyzer = create_analyzer()
    assert_true(analyzer._is_merged_label("F_obs,SIGF_obs"))


def test_label_neither_returns_false():
    """Test that unclassified labels return False for both."""
    analyzer = create_analyzer()
    weird_label = "CUSTOM_COLUMN_XYZ"

    assert_false(analyzer._is_anomalous_label(weird_label))
    assert_false(analyzer._is_merged_label(weird_label))


# =============================================================================
# PROGRAM DATA PREFERENCE TESTS
# =============================================================================

def test_autosol_needs_anomalous():
    """Test that autosol prefers anomalous data."""
    analyzer = create_analyzer()
    needs_anom = analyzer._needs_anomalous_data("phenix.autosol", {})
    assert_true(needs_anom)


def test_hyss_needs_anomalous():
    """Test that hyss prefers anomalous data."""
    analyzer = create_analyzer()
    needs_anom = analyzer._needs_anomalous_data("phenix.hyss", {})
    assert_true(needs_anom)


def test_refine_prefers_merged():
    """Test that refine prefers merged data (doesn't need anomalous)."""
    analyzer = create_analyzer()
    needs_anom = analyzer._needs_anomalous_data("phenix.refine", {})
    assert_false(needs_anom)


def test_phaser_prefers_merged():
    """Test that phaser prefers merged data (doesn't need anomalous)."""
    analyzer = create_analyzer()
    needs_anom = analyzer._needs_anomalous_data("phenix.phaser", {})
    assert_false(needs_anom)


# =============================================================================
# LABEL RESOLUTION TESTS
# =============================================================================

def test_select_anomalous_for_autosol():
    """Test that autosol selects anomalous labels."""
    analyzer = create_analyzer()
    error_info = {
        "keyword": "obs_labels",
        "choices": ["IMEAN,SIGIMEAN", "I(+),SIGI(+),I(-),SIGI(-)"],
        "affected_file": "/path/data.mtz"
    }

    recovery = analyzer._resolve_ambiguous_labels(
        error_info, "phenix.autosol", {}
    )

    assert_not_none(recovery)
    assert_in("(+)", recovery.selected_choice)


def test_select_merged_for_refine():
    """Test that refine selects merged labels."""
    analyzer = create_analyzer()
    error_info = {
        "keyword": "obs_labels",
        "choices": ["IMEAN,SIGIMEAN", "I(+),SIGI(+),I(-),SIGI(-)"],
        "affected_file": "/path/data.mtz"
    }

    recovery = analyzer._resolve_ambiguous_labels(
        error_info, "phenix.refine", {}
    )

    assert_not_none(recovery)
    assert_in("IMEAN", recovery.selected_choice)


def test_fallback_when_no_matching_type():
    """Test fallback to first choice when no matching type."""
    analyzer = create_analyzer()
    error_info = {
        "keyword": "obs_labels",
        "choices": ["CUSTOM_LABEL1", "CUSTOM_LABEL2"],
        "affected_file": "/path/data.mtz"
    }

    recovery = analyzer._resolve_ambiguous_labels(
        error_info, "phenix.refine", {}
    )

    assert_not_none(recovery)
    assert_equal(recovery.selected_choice, "CUSTOM_LABEL1")


# =============================================================================
# RETRY LIMITS TESTS
# =============================================================================

def test_allow_first_retry():
    """Test that first retry is allowed."""
    analyzer = create_analyzer()
    session = MockSession({"recovery_attempts": {}})

    can_retry, reason = analyzer._check_retry_limits(
        session, "ambiguous_data_labels"
    )

    assert_true(can_retry)
    assert_none(reason)


def test_allow_second_retry():
    """Test that second retry is allowed."""
    analyzer = create_analyzer()
    session = MockSession({
        "recovery_attempts": {
            "ambiguous_data_labels": {"count": 1}
        }
    })

    can_retry, reason = analyzer._check_retry_limits(
        session, "ambiguous_data_labels"
    )

    assert_true(can_retry)


def test_deny_after_max_retries():
    """Test that we deny after max retries."""
    analyzer = create_analyzer()
    session = MockSession({
        "recovery_attempts": {
            "ambiguous_data_labels": {"count": 3}
        }
    })

    can_retry, reason = analyzer._check_retry_limits(
        session, "ambiguous_data_labels"
    )

    assert_false(can_retry)
    assert_in("Max", reason)


def test_retry_tracking_increments():
    """Test that retry tracking increments correctly."""
    analyzer = create_analyzer()
    session = MockSession({})

    analyzer._update_retry_tracking(
        session, "ambiguous_data_labels", "/path/data.mtz", "I_CuKa(+)"
    )

    attempts = session.data["recovery_attempts"]["ambiguous_data_labels"]
    assert_equal(attempts["count"], 1)
    assert_in("/path/data.mtz", attempts["files_tried"])


# =============================================================================
# FULL RECOVERY FLOW TESTS
# =============================================================================

def test_full_analyze_flow_autosol():
    """Test complete analyze() call for autosol."""
    analyzer = create_analyzer()
    log = """
    Multiple equally suitable arrays of observed xray data found.

    Possible choices:
      /data/lyso.mtz:IMEAN_CuKa,SIGIMEAN_CuKa
      /data/lyso.mtz:I_CuKa(+),SIGI_CuKa(+),I_CuKa(-),SIGI_CuKa(-)

    Please use scaling.input.xray_data.obs_labels
    to specify an unambiguous substring of the target label.
    """

    session = MockSession({})

    recovery = analyzer.analyze(
        log_text=log,
        program="phenix.autosol",
        context={"project_advice": "MRSAD phasing"},
        session=session
    )

    assert_not_none(recovery)
    assert_equal(recovery.error_type, "ambiguous_data_labels")
    assert_equal(recovery.retry_program, "phenix.autosol")
    assert_in("I_CuKa(+)", recovery.selected_choice)

    # Check flags contain the right parameter
    assert_in("scaling.input.xray_data.obs_labels", recovery.flags)
    assert_equal(recovery.flags["scaling.input.xray_data.obs_labels"], "I_CuKa(+)")


def test_full_analyze_flow_refine():
    """Test complete analyze() call for refine."""
    analyzer = create_analyzer()
    log = """
    Multiple equally suitable arrays found.

    Possible choices:
      /data/data.mtz:IMEAN,SIGIMEAN
      /data/data.mtz:I(+),SIGI(+),I(-),SIGI(-)

    Please use obs_labels to specify.
    """

    session = MockSession({})

    recovery = analyzer.analyze(
        log_text=log,
        program="phenix.refine",
        context={},
        session=session
    )

    assert_not_none(recovery)
    assert_in("IMEAN", recovery.selected_choice)


def test_returns_none_for_non_recoverable():
    """Test that non-recoverable errors return None."""
    analyzer = create_analyzer()
    log = "FATAL ERROR: Model file not found"

    session = MockSession({})

    recovery = analyzer.analyze(
        log_text=log,
        program="phenix.refine",
        context={},
        session=session
    )

    assert_none(recovery)


def test_returns_none_after_max_retries():
    """Test that analyze returns None after max retries."""
    analyzer = create_analyzer()
    log = """
    Multiple equally suitable arrays found.
    Possible choices:
      /data/data.mtz:IMEAN,SIGIMEAN
    Please use obs_labels to specify.
    """

    session = MockSession({
        "recovery_attempts": {
            "ambiguous_data_labels": {"count": 3}
        }
    })

    recovery = analyzer.analyze(
        log_text=log,
        program="phenix.refine",
        context={},
        session=session
    )

    assert_none(recovery)


# =============================================================================
# GET SUGGESTION TESTS
# =============================================================================

def test_suggestion_for_ambiguous_labels():
    """Test that suggestion is generated for ambiguous labels."""
    analyzer = create_analyzer()
    log = """
    Multiple equally suitable arrays found.
    Possible choices:
      /data/data.mtz:IMEAN,SIGIMEAN
      /data/data.mtz:I(+),SIGI(+)
    Please use obs_labels to specify.
    """

    suggestion = analyzer.get_suggestion(log, "phenix.refine")

    assert_not_none(suggestion)
    assert_in("obs_labels", suggestion)
    assert_in("IMEAN", suggestion)


def test_no_suggestion_for_unrecoverable():
    """Test that no suggestion is generated for unrecoverable errors."""
    analyzer = create_analyzer()
    log = "FATAL ERROR: Something went wrong"

    suggestion = analyzer.get_suggestion(log, "phenix.refine")
    assert_none(suggestion)


# =============================================================================
# ERROR RECOVERY DATACLASS TESTS
# =============================================================================

def test_create_error_recovery():
    """Test creating an ErrorRecovery instance."""
    recovery = ErrorRecovery(
        error_type="ambiguous_data_labels",
        affected_file="/path/to/data.mtz",
        flags={"obs_labels": "I(+)"},
        reason="Selected anomalous data",
        retry_program="phenix.autosol",
        selected_choice="I(+),SIGI(+)",
        all_choices=["IMEAN,SIGIMEAN", "I(+),SIGI(+)"]
    )

    assert_equal(recovery.error_type, "ambiguous_data_labels")
    assert_equal(recovery.flags["obs_labels"], "I(+)")
    assert_equal(len(recovery.all_choices), 2)


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
