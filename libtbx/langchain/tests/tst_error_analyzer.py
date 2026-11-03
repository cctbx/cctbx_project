"""
Tests for error_analyzer.py - Recoverable Error Handling.

This module tests the automatic error recovery system that detects
structured errors in PHENIX output and determines appropriate recovery
strategies.

Run with: python tests/tst_error_analyzer.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
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
        all_choices=["IMEAN,SIGIMEAN", "I(+),SIGI(+)"],
        selected_label="I(+)",
        selected_label_pair="I(+),SIGI(+)"
    )

    assert_equal(recovery.error_type, "ambiguous_data_labels")
    assert_equal(recovery.flags["obs_labels"], "I(+)")
    assert_equal(len(recovery.all_choices), 2)
    assert_equal(recovery.selected_label, "I(+)")
    assert_equal(recovery.selected_label_pair, "I(+),SIGI(+)")


def test_error_recovery_label_fields_default_empty():
    """Test that selected_label fields default to empty strings."""
    recovery = ErrorRecovery(
        error_type="ambiguous_data_labels",
        affected_file="/path/to/data.mtz",
        flags={"obs_labels": "I(+)"},
        reason="test",
        retry_program="phenix.refine",
        selected_choice="I(+),SIGI(+)",
    )
    assert_equal(recovery.selected_label, "")
    assert_equal(recovery.selected_label_pair, "")


def test_full_analyze_populates_label_fields():
    """Test that analyze() populates selected_label and selected_label_pair."""
    analyzer = create_analyzer()
    log = """
    Multiple equally suitable arrays found.

    Possible choices:
      /data/toxd.mtz:FTOXD3,SIGFTOXD3
      /data/toxd.mtz:FAU20,SIGFAU20

    Please use scaling.input.xray_data.obs_labels to specify.
    """

    session = MockSession({})
    recovery = analyzer.analyze(
        log_text=log,
        program="phenix.xtriage",
        context={},
        session=session
    )

    assert_not_none(recovery)
    assert_equal(recovery.selected_label, "FTOXD3")
    assert_equal(recovery.selected_label_pair, "FTOXD3,SIGFTOXD3")


# =============================================================================
# EXPERIMENTAL PHASES TESTS
# =============================================================================

def test_detect_ambiguous_experimental_phases():
    """Test detection of ambiguous experimental phases error."""
    analyzer = create_analyzer()
    log = """
    Multiple equally suitable arrays of experimental phases found.

    Possible choices:
      /Users/terwill/unix/run_examples/nsf-d2-sad/AutoSol_run_1_/overall_best_refine_data.mtz:HLAM,HLBM,HLCM,HLDM
      /Users/terwill/unix/run_examples/nsf-d2-sad/AutoSol_run_1_/overall_best_refine_data.mtz:HLanomA,HLanomB,HLanomC,HLanomD

    Please use miller_array.labels.name.
    to specify an unambiguous substring of the target label.

    Sorry: Multiple equally suitable arrays of experimental phases found.
    """

    error_type = analyzer._detect_error_type(log)
    assert_equal(error_type, "ambiguous_experimental_phases")


def test_extract_experimental_phases_choices():
    """Test extraction of HL coefficient choices."""
    analyzer = create_analyzer()
    log = """
    Multiple equally suitable arrays of experimental phases found.

    Possible choices:
      /path/to/data.mtz:HLAM,HLBM,HLCM,HLDM
      /path/to/data.mtz:HLanomA,HLanomB,HLanomC,HLanomD

    Please use miller_array.labels.name.
    to specify an unambiguous substring of the target label.
    """

    error_def = get_error_def(analyzer, "ambiguous_experimental_phases")
    error_info = analyzer._extract_ambiguous_labels_info(log, error_def, "ambiguous_experimental_phases")

    assert_not_none(error_info)
    assert_equal(len(error_info["choices"]), 2)
    assert_in("HLAM,HLBM,HLCM,HLDM", error_info["choices"])
    assert_in("HLanomA,HLanomB,HLanomC,HLanomD", error_info["choices"])


def test_extract_experimental_phases_keyword():
    """Test extraction of miller_array.labels.name keyword."""
    analyzer = create_analyzer()
    log = """
    Multiple equally suitable arrays of experimental phases found.

    Possible choices:
      /path/to/data.mtz:HLAM,HLBM,HLCM,HLDM

    Please use miller_array.labels.name.
    to specify an unambiguous substring of the target label.
    """

    error_def = get_error_def(analyzer, "ambiguous_experimental_phases")
    error_info = analyzer._extract_ambiguous_labels_info(log, error_def, "ambiguous_experimental_phases")

    assert_not_none(error_info)
    assert_equal(error_info["keyword"], "miller_array.labels.name")


def test_resolve_experimental_phases_for_refine():
    """Test that phenix.refine selects standard (non-anomalous) HL coefficients."""
    analyzer = create_analyzer()

    error_info = {
        "keyword": "miller_array.labels.name",
        "choices": ["HLAM,HLBM,HLCM,HLDM", "HLanomA,HLanomB,HLanomC,HLanomD"],
        "affected_file": "/path/to/data.mtz"
    }

    recovery = analyzer._resolve_ambiguous_phases(
        error_info,
        program="phenix.refine",
        context={}
    )

    assert_not_none(recovery)
    assert_equal(recovery.error_type, "ambiguous_experimental_phases")
    assert_equal(recovery.flags["miller_array.labels.name"], "HLAM")
    assert_in("standard", recovery.reason.lower())


def test_resolve_experimental_phases_for_autosol():
    """Test that phenix.autosol selects anomalous HL coefficients."""
    analyzer = create_analyzer()

    error_info = {
        "keyword": "miller_array.labels.name",
        "choices": ["HLAM,HLBM,HLCM,HLDM", "HLanomA,HLanomB,HLanomC,HLanomD"],
        "affected_file": "/path/to/data.mtz"
    }

    recovery = analyzer._resolve_ambiguous_phases(
        error_info,
        program="phenix.autosol",
        context={}
    )

    assert_not_none(recovery)
    assert_equal(recovery.flags["miller_array.labels.name"], "HLanomA")
    assert_in("anomalous", recovery.reason.lower())


def test_is_anomalous_hl():
    """Test anomalous HL coefficient detection."""
    analyzer = create_analyzer()

    assert_true(analyzer._is_anomalous_hl("HLanomA,HLanomB,HLanomC,HLanomD"))
    assert_true(analyzer._is_anomalous_hl("HLanom"))
    assert_false(analyzer._is_anomalous_hl("HLAM,HLBM,HLCM,HLDM"))
    assert_false(analyzer._is_anomalous_hl("HLA,HLB,HLC,HLD"))


def test_full_experimental_phases_recovery():
    """Test full recovery flow for experimental phases error."""
    analyzer = create_analyzer()
    session = MockSession()

    log = """
    Multiple equally suitable arrays of experimental phases found.

    Possible choices:
      /path/to/data.mtz:HLAM,HLBM,HLCM,HLDM
      /path/to/data.mtz:HLanomA,HLanomB,HLanomC,HLanomD

    Please use miller_array.labels.name.
    to specify an unambiguous substring of the target label.
    """

    recovery = analyzer.analyze(
        log_text=log,
        program="phenix.refine",
        context={},
        session=session
    )

    assert_not_none(recovery)
    assert_equal(recovery.error_type, "ambiguous_experimental_phases")
    assert_equal(recovery.retry_program, "phenix.refine")
    assert_equal(recovery.flags["miller_array.labels.name"], "HLAM")
    assert_equal(recovery.affected_file, "/path/to/data.mtz")


# =============================================================================
# v119.H17: strip_parameter resolution tests
# =============================================================================
#
# H17 adds the `missing_phib_input_map_file` recovery (autobuild given
# a map_file= MTZ lacking PHIB phase columns) and wires the generic
# `strip_parameter` resolution that was declared in YAML for
# `rfree_flags_mismatch` but had no code dispatch before H17.
#
# Tests below cover:
#   - Detection of the exact lysozyme-MRSAD error message
#   - ErrorRecovery.strip_flags content
#   - ErrorRecovery.flags stays empty (no spurious additions)
#   - rfree_flags_mismatch retroactive resolution (bonus side benefit)
#   - Pattern variants
#   - False-positive guard: "input_map_file" mention in unrelated
#     contexts (without PHIB) must NOT trigger the H17 recovery
# =============================================================================

H17_LYSOZYME_LOG = """
Running: phenix.autobuild data=.../overall_best_refine_data.mtz \
seq_file=.../hewl.seq model=.../1fkq_prot.pdb \
map_file=.../lyso2001_scala1.mtz rebuild_in_place=False nproc=4

Sorry:
Sorry, PHIB is required for input_map_file
( ['F_CuKa', 'None', 'None'] supplied for ['FP', 'PHIB', 'FOM'])
"""


def test_h17_detect_phib_missing_from_lysozyme_log():
    """Detection of the exact PHIB-missing error from Tom's lysozyme-MRSAD run."""
    analyzer = create_analyzer()
    error_type = analyzer._detect_error_type(H17_LYSOZYME_LOG)
    assert_equal(error_type, "missing_phib_input_map_file")


def test_h17_recovery_strips_map_file_flags():
    """Recovery for PHIB error includes the expected flag prefixes to strip."""
    analyzer = create_analyzer()
    session = MockSession({"recovery_attempts": {}})

    recovery = analyzer.analyze(
        log_text=H17_LYSOZYME_LOG,
        program="phenix.autobuild",
        context={},
        session=session,
    )

    assert_not_none(recovery)
    assert_equal(recovery.error_type, "missing_phib_input_map_file")
    assert_equal(recovery.retry_program, "phenix.autobuild")
    # All three PHIL variants for the slot must be in strip_flags
    assert_in("map_file", recovery.strip_flags)
    assert_in("input_map_file", recovery.strip_flags)
    assert_in("input_files.map_file", recovery.strip_flags)


def test_h17_recovery_does_not_add_flags():
    """strip_parameter recoveries leave .flags empty — they only remove,
    never add."""
    analyzer = create_analyzer()
    session = MockSession({"recovery_attempts": {}})

    recovery = analyzer.analyze(
        log_text=H17_LYSOZYME_LOG,
        program="phenix.autobuild",
        context={},
        session=session,
    )

    assert_not_none(recovery)
    assert_equal(recovery.flags, {},
                 "strip_parameter recovery must not populate .flags")


def test_h17_false_positive_input_map_file_without_phib():
    """A log mentioning 'input_map_file' in an unrelated context
    (WITHOUT 'PHIB') must NOT trigger missing_phib_input_map_file.

    This is the false-positive guard Gemini specifically asked for:
    a loose match on just 'input_map_file' or 'required for
    input_map_file' (without the PHIB anchor) would false-positive
    on unrelated PHENIX validation errors that mention the slot in
    other contexts.
    """
    analyzer = create_analyzer()
    # Realistic unrelated error that mentions input_map_file but NOT PHIB
    log = """
    Sorry, please check your input_map_file path — the file
    appears to be missing or unreadable.  Verify the path exists.
    """

    error_type = analyzer._detect_error_type(log)
    assert_true(error_type != "missing_phib_input_map_file",
                "Unrelated 'input_map_file' message should NOT trigger "
                "missing_phib_input_map_file; got error_type=%r" % error_type)


def test_h17_rfree_flags_mismatch_now_resolves():
    """Regression for the bonus retroactive fix: H17 wires the generic
    strip_parameter dispatch, which retroactively resolves
    rfree_flags_mismatch (previously declared in YAML but unimplemented).
    """
    analyzer = create_analyzer()
    session = MockSession({"recovery_attempts": {}})

    log = """
    Sorry, please resolve the R-free flags mismatch.
    The input MTZ already contains R-free flags.
    """

    recovery = analyzer.analyze(
        log_text=log,
        program="phenix.refine",
        context={},
        session=session,
    )

    assert_not_none(recovery)
    assert_equal(recovery.error_type, "rfree_flags_mismatch")
    # Should strip the generate flag (both variants from YAML)
    assert_in("xray_data.r_free_flags.generate=True", recovery.strip_flags)
    assert_in("xray_data.r_free_flags.generate", recovery.strip_flags)
    assert_equal(recovery.flags, {})


def test_h17_phib_pattern_variant_with_spacing():
    """The second YAML pattern 'PHIB.*required.*input_map_file' allows
    for wording variants (different word ordering, extra spaces).
    Ensures the detector isn't overly literal."""
    analyzer = create_analyzer()
    # Variant with different spacing — DOTALL means .* spans newlines
    log = """
    Sorry, PHIB column is required for the input_map_file parameter
    in autobuild.  Got Fobs only.
    """

    error_type = analyzer._detect_error_type(log)
    assert_equal(error_type, "missing_phib_input_map_file")


def test_h17_strip_flags_field_default_empty_list():
    """ErrorRecovery.strip_flags defaults to an empty list for
    backward compatibility with existing recoveries
    (ambiguous_data_labels, ambiguous_experimental_phases) that
    don't strip anything."""
    analyzer = create_analyzer()
    session = MockSession({"recovery_attempts": {}})

    log = """
    Multiple equally suitable arrays of observed xray data found.
    Possible choices:
      /path/data.mtz:IMEAN,SIGIMEAN
      /path/data.mtz:I(+),SIGI(+),I(-),SIGI(-)
    Please use scaling.input.xray_data.obs_labels
    to specify an unambiguous substring.
    """

    recovery = analyzer.analyze(
        log_text=log,
        program="phenix.xtriage",
        context={},
        session=session,
    )

    assert_not_none(recovery)
    assert_equal(recovery.error_type, "ambiguous_data_labels")
    # Existing recoveries should have empty strip_flags (no behavior change)
    assert_equal(recovery.strip_flags, [],
                 "Pre-H17 recoveries should have empty strip_flags")


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
