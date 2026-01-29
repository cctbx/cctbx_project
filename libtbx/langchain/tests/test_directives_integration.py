"""
Integration tests for the complete directives system.

These tests verify that directives flow correctly through the entire system,
from extraction to validation to workflow state modification.

Run with: python tests/test_directives_integration.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.test_utils import (
    assert_equal, assert_true, assert_false, assert_in, assert_not_in,
    assert_not_none, assert_greater,
    run_tests_with_fail_fast
)

from agent.directive_extractor import (
    extract_directives_simple,
    get_program_settings,
    check_stop_conditions,
    validate_directives,
)

from agent.directive_validator import (
    validate_intent,
    augment_intent_with_directives,
    _check_stop_conditions as validator_check_stop,
)


# =============================================================================
# END-TO-END TESTS
# =============================================================================

def test_resolution_flow_simple_extraction():
    """Test resolution directive from extraction to validation."""
    # Step 1: Extract directives from user advice
    user_advice = "use resolution 3.0 in autosol but 2.5 in other programs"

    # Using simple extraction (no LLM) - construct directives manually
    directives = {
        "program_settings": {
            "phenix.autosol": {"resolution": 3.0},
            "default": {"resolution": 2.5}
        }
    }

    # Step 2: Validate directives
    validated = validate_directives(directives)
    assert_equal(validated["program_settings"]["phenix.autosol"]["resolution"], 3.0)
    assert_equal(validated["program_settings"]["default"]["resolution"], 2.5)

    # Step 3: Get program settings
    autosol_settings = get_program_settings(validated, "phenix.autosol")
    refine_settings = get_program_settings(validated, "phenix.refine")

    assert_equal(autosol_settings["resolution"], 3.0)
    assert_equal(refine_settings["resolution"], 2.5)

    # Step 4: Validate LLM intent
    # Case A: LLM correctly uses autosol resolution
    intent_a = {
        "program": "phenix.autosol",
        "strategy": {"resolution": 3.0}
    }
    result_a = validate_intent(intent_a, validated, cycle_number=1)
    assert_equal(result_a["validated_intent"]["strategy"]["resolution"], 3.0)
    assert_equal(len(result_a["warnings"]), 0)

    # Case B: LLM forgets resolution for refine
    intent_b = {
        "program": "phenix.refine",
        "strategy": {}
    }
    result_b = validate_intent(intent_b, validated, cycle_number=2)
    assert_equal(result_b["validated_intent"]["strategy"]["resolution"], 2.5)
    assert_greater(len(result_b["modifications"]), 0)

    # Case C: LLM uses wrong resolution for autosol
    intent_c = {
        "program": "phenix.autosol",
        "strategy": {"resolution": 2.5}  # Wrong!
    }
    result_c = validate_intent(intent_c, validated, cycle_number=1)
    assert_equal(result_c["validated_intent"]["strategy"]["resolution"], 3.0)  # Corrected
    assert_greater(len(result_c["warnings"]), 0)


def test_stop_after_first_refinement():
    """Test 'stop after first refinement' directive."""
    directives = extract_directives_simple("stop after the first refinement")

    assert_in("stop_conditions", directives)
    assert_equal(directives["stop_conditions"]["after_program"], "phenix.refine")
    assert_equal(directives["stop_conditions"]["max_refine_cycles"], 1)

    validated = validate_directives(directives)

    # Before refine runs
    history_before = [{"program": "phenix.phaser"}]
    should_stop, reason = validator_check_stop(
        validated, cycle_number=2, last_program="phenix.phaser",
        history=history_before, log=lambda x: None
    )
    assert_false(should_stop)

    # After refine completes
    history_after = [
        {"program": "phenix.phaser"},
        {"program": "phenix.refine"}
    ]
    should_stop, reason = validator_check_stop(
        validated, cycle_number=3, last_program="phenix.refine",
        history=history_after, log=lambda x: None
    )
    assert_true(should_stop)
    assert_in("refine", reason.lower())


def test_stop_at_cycle_n():
    """Test 'stop at cycle N' directive."""
    directives = extract_directives_simple("stop at cycle 4")

    assert_in("stop_conditions", directives)
    assert_equal(directives["stop_conditions"]["after_cycle"], 4)

    validated = validate_directives(directives)

    # Cycle 3 - should not stop
    should_stop, _ = check_stop_conditions(validated, cycle_number=3, last_program="phenix.refine")
    assert_false(should_stop)

    # Cycle 4 - should stop
    should_stop, reason = check_stop_conditions(validated, cycle_number=4, last_program="phenix.refine")
    assert_true(should_stop)
    assert_in("4", reason)


def test_anisotropic_refinement():
    """Test 'run anisotropic refinement' directive."""
    directives = extract_directives_simple("run anisotropic refinement")

    assert_in("program_settings", directives)
    assert_true(directives["program_settings"]["phenix.refine"]["anisotropic_adp"])

    validated = validate_directives(directives)

    # LLM intent without anisotropic
    intent = {
        "program": "phenix.refine",
        "strategy": {"resolution": 2.0}
    }

    result = validate_intent(intent, validated, cycle_number=1)

    # Should have anisotropic_adp added
    assert_true(result["validated_intent"]["strategy"]["anisotropic_adp"])
    assert_true(any("anisotropic" in m.lower() for m in result["modifications"]))


def test_skip_validation_directive():
    """Test skip_validation directive allows stopping."""
    directives = extract_directives_simple("stop after the first refinement, skip validation")
    directives["stop_conditions"]["skip_validation"] = True

    validated = validate_directives(directives)

    intent = {
        "program": None,
        "stop": True,
        "stop_reason": "User requested"
    }

    result = validate_intent(intent, validated, cycle_number=2, history=[])
    assert_in("validated_intent", result)


def test_multiple_directives_combined():
    """Test multiple directives working together."""
    directives = {
        "program_settings": {
            "phenix.autosol": {"resolution": 3.0, "atom_type": "Se"},
            "phenix.refine": {"anisotropic_adp": True},
            "default": {"resolution": 2.5}
        },
        "stop_conditions": {
            "r_free_target": 0.25,
            "skip_validation": True
        },
        "workflow_preferences": {
            "skip_programs": ["phenix.autobuild"]
        }
    }

    validated = validate_directives(directives)

    # Check autosol settings
    autosol_settings = get_program_settings(validated, "phenix.autosol")
    assert_equal(autosol_settings["resolution"], 3.0)
    assert_equal(autosol_settings["atom_type"], "Se")

    # Check refine settings (should merge default + specific)
    refine_settings = get_program_settings(validated, "phenix.refine")
    assert_equal(refine_settings["resolution"], 2.5)  # From default
    assert_true(refine_settings["anisotropic_adp"])  # From specific

    # Check r_free target stop
    history = [{"program": "phenix.refine", "metrics": {"r_free": 0.24}}]
    should_stop, reason = validator_check_stop(
        validated, cycle_number=3, last_program="phenix.refine",
        history=history, log=lambda x: None
    )
    assert_true(should_stop)
    assert_in("R-free", reason)


def test_empty_directives_passthrough():
    """Test that empty directives don't affect normal operation."""
    intent = {
        "program": "phenix.refine",
        "strategy": {"resolution": 2.5},
        "files": {"data_mtz": "data.mtz"}
    }

    result = validate_intent(intent, {}, cycle_number=1)

    assert_equal(result["validated_intent"]["program"], "phenix.refine")
    assert_equal(result["validated_intent"]["strategy"]["resolution"], 2.5)
    assert_equal(result["modifications"], [])
    assert_equal(result["warnings"], [])


def test_directive_validation_with_history():
    """Test directive validation applies program settings correctly."""
    directives = {
        "program_settings": {
            "phenix.refine": {"resolution": 2.0}
        }
    }

    validated = validate_directives(directives)

    history = [
        {"program": "phenix.phaser"},
        {"program": "phenix.refine"},
    ]

    intent = {"program": "phenix.refine", "strategy": {}}
    result = validate_intent(intent, validated, cycle_number=3, history=history)

    assert_equal(result["validated_intent"]["strategy"]["resolution"], 2.0)


# =============================================================================
# EDGE CASES
# =============================================================================

def test_invalid_resolution_ignored():
    """Invalid resolution values should be ignored."""
    directives = {
        "program_settings": {
            "default": {"resolution": "invalid"}
        }
    }

    validated = validate_directives(directives)
    assert_not_in("program_settings", validated)


def test_unknown_program_fixed():
    """Unknown program names should be fixed if possible."""
    directives = {
        "program_settings": {
            "refine": {"resolution": 2.5}  # Missing "phenix." prefix
        }
    }

    validated = validate_directives(directives)
    assert_in("phenix.refine", validated.get("program_settings", {}))


def test_contradictory_stop_conditions():
    """Test handling of potentially contradictory stop conditions."""
    directives = {
        "stop_conditions": {
            "after_cycle": 3,
            "after_program": "phenix.refine",
            "r_free_target": 0.20
        }
    }

    validated = validate_directives(directives)

    # after_cycle should trigger
    should_stop, reason = check_stop_conditions(
        validated, cycle_number=3, last_program="phenix.refine"
    )
    assert_true(should_stop)


def test_none_values_handled():
    """Test that None values don't cause crashes."""
    intent = {
        "program": "phenix.refine",
        "strategy": None,
        "files": None
    }

    directives = {
        "program_settings": {
            "default": {"resolution": 2.5}
        }
    }

    result = validate_intent(intent, directives, cycle_number=1)
    assert_not_none(result["validated_intent"])
    assert_equal(result["validated_intent"]["strategy"]["resolution"], 2.5)


def test_augment_preserves_existing():
    """augment_intent_with_directives should not override existing values."""
    intent = {
        "program": "phenix.refine",
        "strategy": {"resolution": 3.0}  # Already set
    }

    directives = {
        "program_settings": {
            "default": {"resolution": 2.5, "cycles": 5}
        }
    }

    result = augment_intent_with_directives(intent, directives)

    # Resolution should NOT be changed (already set)
    assert_equal(result["strategy"]["resolution"], 3.0)
    # Cycles should be added (was missing)
    assert_equal(result["strategy"]["cycles"], 5)


# =============================================================================
# WORKFLOW PREFERENCES
# =============================================================================

def test_skip_programs_in_directives():
    """Test that skip_programs is properly validated."""
    directives = {
        "workflow_preferences": {
            "skip_programs": ["phenix.autobuild", "phenix.ligandfit"]
        }
    }

    validated = validate_directives(directives)

    assert_in("workflow_preferences", validated)
    assert_equal(len(validated["workflow_preferences"]["skip_programs"]), 2)
    assert_in("phenix.autobuild", validated["workflow_preferences"]["skip_programs"])


def test_prefer_programs_in_directives():
    """Test that prefer_programs is properly validated."""
    directives = {
        "workflow_preferences": {
            "prefer_programs": ["phenix.phaser", "phenix.refine"]
        }
    }

    validated = validate_directives(directives)

    assert_in("workflow_preferences", validated)
    assert_equal(len(validated["workflow_preferences"]["prefer_programs"]), 2)


def test_experimental_phasing_preference():
    """Test experimental phasing preference."""
    directives = {
        "workflow_preferences": {
            "use_experimental_phasing": True
        }
    }

    validated = validate_directives(directives)
    assert_true(validated["workflow_preferences"]["use_experimental_phasing"])


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
