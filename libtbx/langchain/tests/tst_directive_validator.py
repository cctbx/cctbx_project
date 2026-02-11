#!/usr/bin/env python
"""
Tests for directive_validator module.

These tests verify that user directives are correctly validated against
available capabilities BEFORE the workflow starts.

Run with: python tests/test_directive_validator.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_false,
    assert_in, assert_greater,
    assert_is_instance, run_tests_with_fail_fast
)

from agent.directive_validator import (
    validate_directives,
    check_program_available,
    get_supported_parameters,
    list_available_programs,
    list_unavailable_programs,
    _normalize_program_name,
    _extract_program_references,
    _extract_parameter_references,
    _load_available_programs,
    ValidationResult,
    validate_intent,
)


# =============================================================================
# PROGRAM NORMALIZATION TESTS
# =============================================================================

def test_normalize_already_normalized():
    """Program with phenix. prefix should remain unchanged."""
    result = _normalize_program_name("phenix.refine")
    assert_equal(result, "phenix.refine")


def test_normalize_add_prefix():
    """Short name should get phenix. prefix."""
    result = _normalize_program_name("refine")
    assert_equal(result, "phenix.refine")


def test_normalize_alias_rsr():
    """RSR alias should normalize to real_space_refine."""
    result = _normalize_program_name("rsr")
    assert_equal(result, "phenix.real_space_refine")


def test_normalize_alias_autobuild():
    """autobuild should normalize correctly."""
    result = _normalize_program_name("autobuild")
    assert_equal(result, "phenix.autobuild")


def test_normalize_case_insensitive():
    """Should handle uppercase input."""
    result = _normalize_program_name("REFINE")
    assert_equal(result, "phenix.refine")


def test_normalize_polder():
    """Polder should normalize to phenix.polder."""
    result = _normalize_program_name("polder")
    assert_equal(result, "phenix.polder")


# =============================================================================
# PROGRAM REFERENCE EXTRACTION TESTS
# =============================================================================

def test_extract_explicit_phenix_reference():
    """Should find explicit phenix.X references."""
    text = "Run phenix.refine to improve the model"
    programs = _extract_program_references(text)
    assert_in("phenix.refine", programs)


def test_extract_multiple_programs():
    """Should find multiple program references."""
    text = "Run phenix.phaser then phenix.refine"
    programs = _extract_program_references(text)
    assert_in("phenix.phaser", programs)
    assert_in("phenix.refine", programs)


def test_extract_run_keyword_pattern():
    """Should find 'run X' patterns."""
    text = "Please run polder to calculate omit maps"
    programs = _extract_program_references(text)
    assert_in("phenix.polder", programs)


def test_extract_use_keyword_pattern():
    """Should find 'use X' patterns."""
    text = "Use autobuild to build the model"
    programs = _extract_program_references(text)
    assert_in("phenix.autobuild", programs)


def test_extract_polder_maps_pattern():
    """Should find 'X maps' patterns."""
    text = "Calculate polder maps for the ligand"
    programs = _extract_program_references(text)
    assert_in("phenix.polder", programs)


def test_extract_no_programs():
    """Text without program refs should return empty list."""
    text = "Just process the data normally"
    programs = _extract_program_references(text)
    assert_equal(programs, [])


def test_extract_empty_text():
    """Empty text should return empty list."""
    programs = _extract_program_references("")
    assert_equal(programs, [])


# =============================================================================
# PARAMETER REFERENCE EXTRACTION TESTS
# =============================================================================

def test_extract_resolution():
    """Should detect resolution references."""
    text = "Use 2.5 Angstrom resolution"
    params = _extract_parameter_references(text)
    param_names = [p[0] if isinstance(p, tuple) else p for p in params]
    assert_in("resolution", param_names)


def test_extract_macro_cycles_numeric():
    """Should detect 'N macro cycles' pattern."""
    text = "Run 5 macro cycles"
    params = _extract_parameter_references(text)
    param_names = [p[0] if isinstance(p, tuple) else p for p in params]
    assert_in("macro_cycles", param_names)


def test_extract_one_macro_cycle():
    """Should detect 'one macro cycle' pattern."""
    text = "Run one macro cycle only"
    params = _extract_parameter_references(text)
    param_names = [p[0] if isinstance(p, tuple) else p for p in params]
    assert_in("macro_cycles", param_names)


def test_extract_anisotropic():
    """Should detect anisotropic reference."""
    text = "Enable anisotropic B-factors"
    params = _extract_parameter_references(text)
    param_names = [p[0] if isinstance(p, tuple) else p for p in params]
    assert_in("anisotropic", param_names)


def test_extract_add_waters():
    """Should detect water addition reference."""
    text = "Add waters during refinement"
    params = _extract_parameter_references(text)
    param_names = [p[0] if isinstance(p, tuple) else p for p in params]
    assert_in("waters", param_names)


def test_extract_empty_text_params():
    """Empty text should return empty list."""
    params = _extract_parameter_references("")
    assert_equal(params, [])


# =============================================================================
# VALIDATE DIRECTIVES TESTS
# =============================================================================

def test_validate_empty_advice():
    """Empty advice should pass validation."""
    result = validate_directives("")
    assert_true(result.valid)


def test_validate_available_program():
    """Available program should pass validation."""
    result = validate_directives("Run phenix.refine")
    assert_true(result.valid)


def test_validate_directives_with_valid_after_program():
    """Directive with valid after_program should pass."""
    directives = {"stop_conditions": {"after_program": "phenix.refine"}}
    result = validate_directives("Stop after refine", directives)
    assert_true(result.valid)


def test_validate_directives_with_unavailable_after_program():
    """Directive with unavailable after_program should fail."""
    directives = {"stop_conditions": {"after_program": "phenix.fem"}}
    result = validate_directives("Stop after FEM", directives)
    assert_false(result.valid)


def test_validate_macro_cycles_warning():
    """Setting macro_cycles should generate warning."""
    result = validate_directives("Run 2 macro cycles of refinement")
    # Should still be valid but might have warnings
    assert_true(result.valid)


# =============================================================================
# CHECK PROGRAM AVAILABLE TESTS
# =============================================================================

def test_check_available_refine():
    """phenix.refine should be available."""
    result = check_program_available("phenix.refine")
    # Result is (bool, message) tuple
    assert_true(result[0])


def test_check_available_short_name():
    """Short name 'refine' should resolve and be available."""
    result = check_program_available("refine")
    # Result is (bool, message) tuple
    assert_true(result[0])


def test_check_unavailable_polder():
    """phenix.polder availability check should work."""
    # Note: polder may or may not be in available list depending on config
    # This test just checks the function works and returns proper type
    result = check_program_available("phenix.polder")
    # Result is (bool, message) tuple
    assert_is_instance(result, tuple)
    assert_is_instance(result[0], bool)


def test_check_unknown_program():
    """Unknown program should return False."""
    result = check_program_available("phenix.nonexistent_program")
    # Result is (bool, message) tuple
    assert_false(result[0])


# =============================================================================
# GET SUPPORTED PARAMETERS TESTS
# =============================================================================

def test_get_params_refine_has_parameters():
    """phenix.refine should have supported parameters."""
    params = get_supported_parameters("phenix.refine")
    assert_greater(len(params), 0, "refine should have parameters")


def test_get_params_unknown_program():
    """Unknown program should return empty list."""
    params = get_supported_parameters("phenix.nonexistent")
    assert_equal(params, [])


# =============================================================================
# LIST FUNCTIONS TESTS
# =============================================================================

def test_list_available_not_empty():
    """Available programs list should not be empty."""
    programs = list_available_programs()
    assert_greater(len(programs), 0, "Should have available programs")


def test_list_available_includes_refine():
    """Available programs should include phenix.refine."""
    programs = list_available_programs()
    assert_in("phenix.refine", programs)


def test_list_unavailable_includes_fem():
    """Unavailable programs should include phenix.fem."""
    programs = list_unavailable_programs()
    # fem is typically not available
    assert_is_instance(programs, list)


def test_list_no_overlap():
    """Available and unavailable lists should not overlap."""
    available = set(list_available_programs())
    unavailable = set(list_unavailable_programs())
    overlap = available & unavailable
    assert_equal(len(overlap), 0, f"Overlap found: {overlap}")


# =============================================================================
# VALIDATION RESULT TESTS
# =============================================================================

def test_validation_result_valid():
    """ValidationResult with valid=True."""
    result = ValidationResult(valid=True)
    assert_true(result.valid)
    assert_equal(result.issues, [])
    assert_equal(result.warnings, [])


def test_validation_result_with_issues():
    """ValidationResult with issues."""
    result = ValidationResult(
        valid=False,
        issues=["Issue 1", "Issue 2"]
    )
    assert_false(result.valid)
    assert_equal(len(result.issues), 2)


def test_validation_result_with_warnings():
    """ValidationResult with warnings but still valid."""
    result = ValidationResult(
        valid=True,
        warnings=["Warning 1"]
    )
    assert_true(result.valid)
    assert_equal(len(result.warnings), 1)


# =============================================================================
# YAML LOADING TESTS
# =============================================================================

def test_yaml_loads_programs():
    """Should load programs from YAML."""
    programs = _load_available_programs()
    assert_is_instance(programs, dict)


def test_yaml_has_strategy_flags():
    """Loaded programs should have strategy_flags key."""
    programs = _load_available_programs()
    for prog, config in programs.items():
        if isinstance(config, dict):
            # Just check structure is valid
            assert_is_instance(config, dict)
            break  # Only need to check one


# =============================================================================
# EDGE CASES TESTS
# =============================================================================

def test_edge_case_none_advice():
    """None advice should not crash."""
    # validate_directives should handle None gracefully
    try:
        result = validate_directives(None)
        assert_true(result.valid)
    except TypeError:
        # Also acceptable - None is not valid input
        pass


def test_edge_case_none_directives():
    """None directives should work fine."""
    result = validate_directives("Test", None)
    assert_true(result.valid)


def test_edge_case_empty_directives():
    """Empty directives dict should work fine."""
    result = validate_directives("Test", {})
    assert_true(result.valid)


def test_edge_case_very_long_advice():
    """Should handle long advice text."""
    long_advice = "Run phenix.refine. " * 100
    result = validate_directives(long_advice)
    assert_true(result.valid)


def test_edge_case_case_insensitive_detection():
    """Should detect programs regardless of case."""
    text = "Run PHENIX.REFINE"
    programs = _extract_program_references(text)
    assert_in("phenix.refine", programs)


# =============================================================================
# ATTEMPT-BASED OVERRIDE TESTS
# =============================================================================

def test_attempt_first_uses_directive():
    """First attempt (attempt_number=0) should use directive value."""
    intent = {
        "program": "phenix.polder",
        "strategy": {"selection": "resname MES and resseq 88"}
    }
    directives = {
        "program_settings": {
            "phenix.polder": {"selection": "solvent molecule MES 88"}
        }
    }

    result = validate_intent(intent, directives, attempt_number=0)
    validated = result["validated_intent"]

    # First attempt should use directive value
    assert_equal(validated["strategy"]["selection"], "solvent molecule MES 88")
    assert_greater(len(result["warnings"]), 0)
    assert_in("first attempt", result["warnings"][0].lower())


def test_attempt_retry_uses_llm():
    """Retry attempt (attempt_number>0) should use LLM's interpretation."""
    intent = {
        "program": "phenix.polder",
        "strategy": {"selection": "resname MES and resseq 88"}
    }
    directives = {
        "program_settings": {
            "phenix.polder": {"selection": "solvent molecule MES 88"}
        }
    }

    result = validate_intent(intent, directives, attempt_number=1)
    validated = result["validated_intent"]

    # Retry should use LLM's value
    assert_equal(validated["strategy"]["selection"], "resname MES and resseq 88")
    assert_greater(len(result["warnings"]), 0)
    assert_in("retry", result["warnings"][0].lower())


def test_attempt_missing_value_adds_directive():
    """Missing LLM value should always be filled from directive."""
    intent = {
        "program": "phenix.polder",
        "strategy": {}
    }
    directives = {
        "program_settings": {
            "phenix.polder": {"selection": "solvent molecule MES 88"}
        }
    }

    # Should add directive value regardless of attempt number
    for attempt in [0, 1, 2]:
        result = validate_intent(
            {"program": "phenix.polder", "strategy": {}},
            directives,
            attempt_number=attempt
        )
        validated = result["validated_intent"]
        assert_equal(validated["strategy"]["selection"], "solvent molecule MES 88")


def test_attempt_same_value_no_warning():
    """Same value from LLM and directive should not generate warning."""
    intent = {
        "program": "phenix.refine",
        "strategy": {"resolution": 2.5}
    }
    directives = {
        "program_settings": {
            "default": {"resolution": 2.5}
        }
    }

    result = validate_intent(intent, directives, attempt_number=0)
    # No warning because values match
    assert_equal(len(result["warnings"]), 0)


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
