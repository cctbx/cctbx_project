#!/usr/bin/env python
"""
Tests for directive_validator module.

These tests verify that LLM intents are correctly validated and
augmented against stored directives.

Run with: python tests/test_directive_validator.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from agent.directive_validator import (
    validate_intent,
    augment_intent_with_directives,
    get_stop_reason_from_directives,
    format_validation_result,
    _apply_program_settings,
    _check_stop_conditions,
    _check_max_program_cycles,
)


# =============================================================================
# Tests for validate_intent
# =============================================================================

def test_empty_directives():
    """Empty directives should return intent unchanged."""
    print("Test: empty_directives")

    intent = {
        "program": "phenix.refine",
        "strategy": {"resolution": 2.5},
        "files": {"mtz": "data.mtz"}
    }

    result = validate_intent(intent, {}, cycle_number=1)

    assert result["validated_intent"]["program"] == "phenix.refine"
    assert result["validated_intent"]["strategy"]["resolution"] == 2.5
    assert result["modifications"] == []
    assert result["warnings"] == []

    print("  PASSED")


def test_none_directives():
    """None directives should return intent unchanged."""
    print("Test: none_directives")

    intent = {"program": "phenix.refine", "strategy": {}}
    result = validate_intent(intent, None, cycle_number=1)

    assert result["validated_intent"]["program"] == "phenix.refine"

    print("  PASSED")


def test_adds_missing_strategy_from_default():
    """Should add missing strategy params from default directive."""
    print("Test: adds_missing_strategy_from_default")

    intent = {
        "program": "phenix.refine",
        "strategy": {}
    }
    directives = {
        "program_settings": {
            "default": {"resolution": 2.5}
        }
    }

    result = validate_intent(intent, directives, cycle_number=1)

    assert result["validated_intent"]["strategy"]["resolution"] == 2.5
    assert len(result["modifications"]) > 0
    assert "resolution" in result["modifications"][0]

    print("  PASSED")


def test_adds_missing_strategy_from_specific():
    """Should add missing strategy params from program-specific directive."""
    print("Test: adds_missing_strategy_from_specific")

    intent = {
        "program": "phenix.refine",
        "strategy": {}
    }
    directives = {
        "program_settings": {
            "phenix.refine": {"anisotropic_adp": True}
        }
    }

    result = validate_intent(intent, directives, cycle_number=1)

    assert result["validated_intent"]["strategy"]["anisotropic_adp"] == True

    print("  PASSED")


def test_specific_overrides_default():
    """Program-specific settings should override defaults."""
    print("Test: specific_overrides_default")

    intent = {
        "program": "phenix.autosol",
        "strategy": {}
    }
    directives = {
        "program_settings": {
            "default": {"resolution": 2.5},
            "phenix.autosol": {"resolution": 3.0}
        }
    }

    result = validate_intent(intent, directives, cycle_number=1)

    assert result["validated_intent"]["strategy"]["resolution"] == 3.0

    print("  PASSED")


def test_warns_and_overrides_llm_contradiction():
    """Should warn and override when LLM contradicts directive."""
    print("Test: warns_and_overrides_llm_contradiction")

    intent = {
        "program": "phenix.refine",
        "strategy": {"resolution": 3.5}  # LLM says 3.5
    }
    directives = {
        "program_settings": {
            "default": {"resolution": 2.5}  # Directive says 2.5
        }
    }

    result = validate_intent(intent, directives, cycle_number=1)

    # Should override to directive value
    assert result["validated_intent"]["strategy"]["resolution"] == 2.5
    # Should have a warning
    assert len(result["warnings"]) > 0
    assert "3.5" in result["warnings"][0]

    print("  PASSED")


# =============================================================================
# Tests for _apply_program_settings
# =============================================================================

def test_empty_settings():
    """Empty settings should make no changes."""
    print("Test: empty_settings")

    intent = {"strategy": {"existing": "value"}}
    directives = {"program_settings": {}}

    mods, warns = _apply_program_settings(
        intent, directives, "phenix.refine", lambda x: None
    )

    assert mods == []
    assert warns == []
    assert intent["strategy"]["existing"] == "value"

    print("  PASSED")


def test_merge_multiple_settings():
    """Should merge multiple settings correctly."""
    print("Test: merge_multiple_settings")

    intent = {"strategy": {}}
    directives = {
        "program_settings": {
            "default": {"resolution": 2.5, "cycles": 5},
            "phenix.refine": {"anisotropic_adp": True}
        }
    }

    mods, warns = _apply_program_settings(
        intent, directives, "phenix.refine", lambda x: None
    )

    assert intent["strategy"]["resolution"] == 2.5
    assert intent["strategy"]["cycles"] == 5
    assert intent["strategy"]["anisotropic_adp"] == True
    assert len(mods) == 3

    print("  PASSED")


# =============================================================================
# Tests for _check_stop_conditions
# =============================================================================

def test_no_conditions():
    """No conditions should not trigger stop."""
    print("Test: no_conditions")

    should_stop, reason = _check_stop_conditions(
        {}, 5, "phenix.refine", [], lambda x: None
    )
    assert should_stop == False
    assert reason is None

    print("  PASSED")


def test_r_free_target_reached():
    """Should stop when R-free target is reached."""
    print("Test: r_free_target_reached")

    directives = {
        "stop_conditions": {"r_free_target": 0.25}
    }
    history = [
        {"program": "phenix.refine", "metrics": {"r_free": 0.30}},
        {"program": "phenix.refine", "metrics": {"r_free": 0.24}}  # Below target
    ]

    should_stop, reason = _check_stop_conditions(
        directives, 3, "phenix.refine", history, lambda x: None
    )

    assert should_stop == True
    assert "R-free" in reason

    print("  PASSED")


def test_r_free_target_not_reached():
    """Should not stop when R-free target not reached."""
    print("Test: r_free_target_not_reached")

    directives = {
        "stop_conditions": {"r_free_target": 0.20}
    }
    history = [
        {"program": "phenix.refine", "metrics": {"r_free": 0.25}}
    ]

    should_stop, reason = _check_stop_conditions(
        directives, 2, "phenix.refine", history, lambda x: None
    )

    assert should_stop == False

    print("  PASSED")


def test_map_cc_target():
    """Should stop when map CC target is reached."""
    print("Test: map_cc_target")

    directives = {
        "stop_conditions": {"map_cc_target": 0.70}
    }
    history = [
        {"program": "phenix.real_space_refine", "metrics": {"map_cc": 0.75}}
    ]

    should_stop, reason = _check_stop_conditions(
        directives, 2, "phenix.real_space_refine", history, lambda x: None
    )

    assert should_stop == True
    assert "CC" in reason

    print("  PASSED")


# =============================================================================
# Tests for _check_max_program_cycles
# =============================================================================

def test_no_max_set():
    """Should not stop when no max is set."""
    print("Test: no_max_set")

    directives = {"stop_conditions": {}}
    history = [
        {"program": "phenix.refine"},
        {"program": "phenix.refine"},
        {"program": "phenix.refine"},
    ]

    should_stop, reason = _check_max_program_cycles(
        directives, "phenix.refine", history, lambda x: None
    )

    assert should_stop == False

    print("  PASSED")


def test_max_refine_cycles_reached():
    """Should stop when max refine cycles reached."""
    print("Test: max_refine_cycles_reached")

    directives = {
        "stop_conditions": {"max_refine_cycles": 2}
    }
    history = [
        {"program": "phenix.phaser"},
        {"program": "phenix.refine"},
        {"program": "phenix.refine"},  # 2nd refine
    ]

    should_stop, reason = _check_max_program_cycles(
        directives, "phenix.refine", history, lambda x: None
    )

    assert should_stop == True
    assert "max" in reason.lower()

    print("  PASSED")


def test_max_refine_cycles_not_reached():
    """Should not stop when max refine cycles not reached."""
    print("Test: max_refine_cycles_not_reached")

    directives = {
        "stop_conditions": {"max_refine_cycles": 3}
    }
    history = [
        {"program": "phenix.refine"},
    ]

    should_stop, reason = _check_max_program_cycles(
        directives, "phenix.refine", history, lambda x: None
    )

    assert should_stop == False

    print("  PASSED")


# =============================================================================
# Tests for augment_intent_with_directives
# =============================================================================

def test_augments_missing_only():
    """Should only add missing settings, not override existing."""
    print("Test: augments_missing_only")

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

    # Should keep existing resolution
    assert result["strategy"]["resolution"] == 3.0
    # Should add missing cycles
    assert result["strategy"]["cycles"] == 5

    print("  PASSED")


def test_handles_stop_program():
    """Should handle STOP program gracefully."""
    print("Test: handles_stop_program")

    intent = {"program": "STOP", "strategy": {}}
    directives = {
        "program_settings": {"default": {"resolution": 2.5}}
    }

    result = augment_intent_with_directives(intent, directives)

    # Should not add settings to STOP
    assert "resolution" not in result.get("strategy", {})

    print("  PASSED")


# =============================================================================
# Tests for get_stop_reason_from_directives
# =============================================================================

def test_returns_reason_when_should_stop():
    """Should return reason when stop condition met."""
    print("Test: returns_reason_when_should_stop")

    directives = {"stop_conditions": {"after_cycle": 3}}

    reason = get_stop_reason_from_directives(
        directives, cycle_number=3, last_program="phenix.refine"
    )

    assert reason is not None
    assert "cycle" in reason.lower()

    print("  PASSED")


def test_returns_none_when_should_not_stop():
    """Should return None when no stop condition met."""
    print("Test: returns_none_when_should_not_stop")

    directives = {"stop_conditions": {"after_cycle": 5}}

    reason = get_stop_reason_from_directives(
        directives, cycle_number=2, last_program="phenix.refine"
    )

    assert reason is None

    print("  PASSED")


# =============================================================================
# Tests for format_validation_result
# =============================================================================

def test_formats_modifications():
    """Should format modifications list."""
    print("Test: formats_modifications")

    result = {
        "modifications": ["Added resolution=2.5"],
        "warnings": [],
        "should_stop": False
    }

    output = format_validation_result(result)

    assert "Modifications" in output
    assert "resolution" in output

    print("  PASSED")


def test_formats_warnings():
    """Should format warnings list."""
    print("Test: formats_warnings")

    result = {
        "modifications": [],
        "warnings": ["LLM contradicted directive"],
        "should_stop": False
    }

    output = format_validation_result(result)

    assert "Warnings" in output
    assert "contradicted" in output

    print("  PASSED")


def test_formats_stop():
    """Should format stop reason."""
    print("Test: formats_stop")

    result = {
        "modifications": [],
        "warnings": [],
        "should_stop": True,
        "stop_reason": "Max cycles reached"
    }

    output = format_validation_result(result)

    assert "Stop" in output
    assert "Max cycles" in output

    print("  PASSED")


def test_empty_result():
    """Should handle empty result."""
    print("Test: empty_result")

    result = {
        "modifications": [],
        "warnings": [],
        "should_stop": False
    }

    output = format_validation_result(result)

    assert "No changes" in output

    print("  PASSED")


# =============================================================================
# Test runner
# =============================================================================

def run_all_tests():
    """Run all directive_validator tests."""
    print("=" * 60)
    print("DIRECTIVE VALIDATOR TESTS")
    print("=" * 60)

    # validate_intent tests
    test_empty_directives()
    test_none_directives()
    test_adds_missing_strategy_from_default()
    test_adds_missing_strategy_from_specific()
    test_specific_overrides_default()
    test_warns_and_overrides_llm_contradiction()

    # _apply_program_settings tests
    test_empty_settings()
    test_merge_multiple_settings()

    # _check_stop_conditions tests
    test_no_conditions()
    test_r_free_target_reached()
    test_r_free_target_not_reached()
    test_map_cc_target()

    # _check_max_program_cycles tests
    test_no_max_set()
    test_max_refine_cycles_reached()
    test_max_refine_cycles_not_reached()

    # augment_intent_with_directives tests
    test_augments_missing_only()
    test_handles_stop_program()

    # get_stop_reason_from_directives tests
    test_returns_reason_when_should_stop()
    test_returns_none_when_should_not_stop()

    # format_validation_result tests
    test_formats_modifications()
    test_formats_warnings()
    test_formats_stop()
    test_empty_result()

    print()
    print("=" * 60)
    print("ALL DIRECTIVE VALIDATOR TESTS PASSED!")
    print("=" * 60)


if __name__ == "__main__":
    run_all_tests()
