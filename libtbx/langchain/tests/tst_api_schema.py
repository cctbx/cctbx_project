"""
Tests for API Schema and Client Adapter.

Run with: python tests/tst_api_schema.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from knowledge.api_schema import (
    API_VERSION,
    validate_request,
    validate_response,
    apply_request_defaults,
    create_request,
    create_response,
    create_error_response,
    create_stop_response,
    parse_version,
    is_version_supported,
)

from agent.api_client import (
    build_request_v2,
    parse_response_v2,
    serialize_request,
    deserialize_request,
    serialize_response,
    deserialize_response,
    detect_api_version,
    is_v2_request,
)


# =============================================================================
# TEST HELPERS
# =============================================================================

def assert_equal(actual, expected, message=""):
    if actual != expected:
        raise AssertionError(f"{message}: expected {expected!r}, got {actual!r}")


def assert_true(value, message=""):
    if not value:
        raise AssertionError(f"{message}: expected True, got {value!r}")


def assert_false(value, message=""):
    if value:
        raise AssertionError(f"{message}: expected False, got {value!r}")


def assert_in(item, container, message=""):
    if item not in container:
        raise AssertionError(f"{message}: {item!r} not in {container!r}")


def assert_not_none(value, message=""):
    if value is None:
        raise AssertionError(f"{message}: expected not None")


def assert_is_none(value, message=""):
    if value is not None:
        raise AssertionError(f"{message}: expected None, got {value!r}")


# =============================================================================
# SCHEMA TESTS
# =============================================================================

def test_api_version():
    """Test API version is defined."""
    print("Test: api_version")

    assert_not_none(API_VERSION)
    assert_true(API_VERSION.startswith("2."), "API version should be 2.x")

    print("  PASSED")


def test_parse_version():
    """Test version string parsing."""
    print("Test: parse_version")

    assert_equal(parse_version("2.0"), (2, 0))
    assert_equal(parse_version("2.1"), (2, 1))
    assert_equal(parse_version("2.1.3"), (2, 1, 3))
    assert_equal(parse_version("1.0"), (1, 0))
    assert_equal(parse_version(None), (1, 0), "None should default to 1.0")
    assert_equal(parse_version(""), (1, 0), "Empty should default to 1.0")
    assert_equal(parse_version("invalid"), (1, 0), "Invalid should default to 1.0")

    print("  PASSED")


def test_is_version_supported():
    """Test version support checking."""
    print("Test: is_version_supported")

    assert_true(is_version_supported("1.0"), "1.0 should be supported")
    assert_true(is_version_supported("2.0"), "2.0 should be supported")
    assert_true(is_version_supported(API_VERSION), "Current version should be supported")

    print("  PASSED")


def test_validate_request_minimal():
    """Test validation of minimal valid request."""
    print("Test: validate_request_minimal")

    request = {
        "api_version": "2.0",
        "files": ["/path/to/file.mtz"],
        "cycle_number": 1,
    }

    is_valid, errors = validate_request(request)
    assert_true(is_valid, f"Minimal request should be valid: {errors}")
    assert_equal(len(errors), 0)

    print("  PASSED")


def test_validate_request_missing_required():
    """Test validation catches missing required fields."""
    print("Test: validate_request_missing_required")

    request = {
        "api_version": "2.0",
        # Missing: files, cycle_number
    }

    is_valid, errors = validate_request(request)
    assert_false(is_valid, "Request missing required fields should be invalid")
    assert_true(len(errors) > 0, "Should have error messages")

    print("  PASSED")


def test_validate_request_wrong_types():
    """Test validation catches wrong types."""
    print("Test: validate_request_wrong_types")

    request = {
        "api_version": "2.0",
        "files": "not a list",  # Should be list
        "cycle_number": "1",  # Should be int
    }

    is_valid, errors = validate_request(request)
    assert_false(is_valid, "Request with wrong types should be invalid")
    assert_true(len(errors) >= 2, f"Should have multiple errors: {errors}")

    print("  PASSED")


def test_validate_request_ignores_unknown():
    """Test validation ignores unknown fields by default."""
    print("Test: validate_request_ignores_unknown")

    request = {
        "api_version": "2.0",
        "files": [],
        "cycle_number": 1,
        "unknown_field": "should be ignored",
        "another_unknown": 123,
    }

    is_valid, errors = validate_request(request)
    assert_true(is_valid, f"Unknown fields should be ignored: {errors}")

    print("  PASSED")


def test_validate_request_strict_mode():
    """Test strict validation catches unknown fields."""
    print("Test: validate_request_strict_mode")

    request = {
        "api_version": "2.0",
        "files": [],
        "cycle_number": 1,
        "unknown_field": "should fail in strict mode",
    }

    is_valid, errors = validate_request(request, strict=True)
    assert_false(is_valid, "Unknown fields should fail in strict mode")
    assert_in("unknown_field", str(errors), "Error should mention unknown field")

    print("  PASSED")


def test_validate_response_minimal():
    """Test validation of minimal valid response."""
    print("Test: validate_response_minimal")

    response = {
        "api_version": "2.0",
        "decision": {
            "program": "phenix.refine",
            "command": "phenix.refine model.pdb data.mtz",
        },
    }

    is_valid, errors = validate_response(response)
    assert_true(is_valid, f"Minimal response should be valid: {errors}")

    print("  PASSED")


def test_apply_request_defaults():
    """Test that defaults are applied correctly."""
    print("Test: apply_request_defaults")

    request = {
        "api_version": "2.0",
        "files": ["/path/to/file.mtz"],
        "cycle_number": 1,
    }

    result = apply_request_defaults(request)

    # Check defaults are applied
    assert_equal(result["log_content"], "")
    assert_equal(result["history"], [])
    assert_equal(result["user_advice"], "")
    assert_equal(result["session_state"], {})
    assert_equal(result["settings"], {})

    # Check original values preserved
    assert_equal(result["files"], ["/path/to/file.mtz"])
    assert_equal(result["cycle_number"], 1)

    print("  PASSED")


def test_apply_request_defaults_nested():
    """Test that nested defaults are applied."""
    print("Test: apply_request_defaults_nested")

    request = {
        "api_version": "2.0",
        "files": [],
        "cycle_number": 1,
        "session_state": {
            "resolution": 2.5,
            # Missing: experiment_type, rfree_mtz, best_files
        },
        "settings": {
            "provider": "anthropic",
            # Missing: abort_on_red_flags, abort_on_warnings
        },
    }

    result = apply_request_defaults(request)

    # Check nested defaults
    assert_equal(result["session_state"]["resolution"], 2.5)
    # Note: experiment_type and rfree_mtz have None defaults, so they're NOT added
    # (the function deliberately skips None defaults to avoid adding unnecessary keys)
    assert "experiment_type" not in result["session_state"]
    assert "rfree_mtz" not in result["session_state"]
    # best_files has {} default, so it IS added
    assert_equal(result["session_state"]["best_files"], {})

    assert_equal(result["settings"]["provider"], "anthropic")
    assert_equal(result["settings"]["abort_on_red_flags"], True)
    assert_equal(result["settings"]["abort_on_warnings"], False)

    print("  PASSED")


def test_create_request():
    """Test request creation helper."""
    print("Test: create_request")

    request = create_request(
        files=["/path/to/data.mtz", "/path/to/model.pdb"],
        cycle_number=3,
        log_content="Previous log...",
        history=[{"cycle": 1, "program": "phenix.xtriage"}],
        session_state={"resolution": 2.0},
        user_advice="Use high resolution",
    )

    assert_equal(request["api_version"], API_VERSION)
    assert_equal(len(request["files"]), 2)
    assert_equal(request["cycle_number"], 3)
    assert_equal(request["log_content"], "Previous log...")
    assert_equal(len(request["history"]), 1)
    assert_equal(request["session_state"]["resolution"], 2.0)
    assert_equal(request["user_advice"], "Use high resolution")

    # Validate the created request
    is_valid, errors = validate_request(request)
    assert_true(is_valid, f"Created request should be valid: {errors}")

    print("  PASSED")


def test_create_response():
    """Test response creation helper."""
    print("Test: create_response")

    response = create_response(
        program="phenix.refine",
        command="phenix.refine model.pdb data.mtz",
        reasoning="Model needs refinement",
        strategy={"output_prefix": "refine_001"},
        experiment_type="xray",
    )

    assert_equal(response["api_version"], API_VERSION)
    assert_equal(response["decision"]["program"], "phenix.refine")
    assert_equal(response["decision"]["command"], "phenix.refine model.pdb data.mtz")
    assert_equal(response["decision"]["reasoning"], "Model needs refinement")
    assert_equal(response["decision"]["strategy"]["output_prefix"], "refine_001")
    assert_equal(response["metadata"]["experiment_type"], "xray")
    assert_false(response["stop"])

    # Validate the created response
    is_valid, errors = validate_response(response)
    assert_true(is_valid, f"Created response should be valid: {errors}")

    print("  PASSED")


def test_create_stop_response():
    """Test stop response creation."""
    print("Test: create_stop_response")

    response = create_stop_response(
        stop_reason="converged",
        reasoning="R-free has converged at 0.22",
        final_metrics={"r_free": 0.22, "clashscore": 2.1},
    )

    assert_true(response["stop"])
    assert_equal(response["stop_reason"], "converged")
    assert_equal(response["decision"]["program"], "STOP")
    assert_equal(response["decision"]["command"], "STOP")
    assert_equal(response["metadata"]["final_metrics"]["r_free"], 0.22)

    print("  PASSED")


def test_create_error_response():
    """Test error response creation."""
    print("Test: create_error_response")

    response = create_error_response(
        error_message="Something went wrong",
        debug_log=["Step 1 OK", "Step 2 FAILED"],
    )

    assert_true(response["stop"])
    assert_equal(response["stop_reason"], "error")
    assert_equal(response["error"], "Something went wrong")
    assert_equal(len(response["debug"]["log"]), 2)

    print("  PASSED")


# =============================================================================
# CLIENT ADAPTER TESTS
# =============================================================================

def test_build_request_v2():
    """Test building v2 request from client data."""
    print("Test: build_request_v2")

    request = build_request_v2(
        files=["/path/to/data.mtz"],
        cycle_number=5,
        log_content="Refinement completed...",
        history=[
            {"cycle_number": 1, "program": "phenix.xtriage", "result": "SUCCESS"},
        ],
        session_state={
            "resolution": 2.1,
            "experiment_type": "xray",
            "rfree_mtz": "/path/to/refine_001_data.mtz",
            "best_files": {"model": "/path/to/best.pdb"},
        },
        user_advice="Solve the structure",
        provider="google",
        abort_on_red_flags=True,
    )

    assert_equal(request["api_version"], API_VERSION)
    assert_equal(request["cycle_number"], 5)
    assert_equal(request["session_state"]["resolution"], 2.1)
    assert_equal(request["session_state"]["experiment_type"], "xray")
    assert_equal(request["settings"]["provider"], "google")

    # Validate
    is_valid, errors = validate_request(request)
    assert_true(is_valid, f"Built request should be valid: {errors}")

    print("  PASSED")


def test_parse_response_v2():
    """Test parsing v2 response to legacy format."""
    print("Test: parse_response_v2")

    response = {
        "api_version": "2.0",
        "decision": {
            "program": "phenix.refine",
            "command": "phenix.refine model.pdb data.mtz",
            "reasoning": "Model needs refinement",
            "strategy": {"output_prefix": "refine_001"},
            "confidence": "high",
        },
        "stop": False,
        "metadata": {
            "experiment_type": "xray",
            "workflow_state": "xray_refined",
        },
        "debug": {
            "log": ["PERCEIVE: OK", "PLAN: OK"],
        },
    }

    parsed = parse_response_v2(response)

    # Check legacy format
    assert_not_none(parsed["next_move"])
    assert_equal(parsed["next_move"]["program"], "phenix.refine")
    assert_equal(parsed["next_move"]["command"], "phenix.refine model.pdb data.mtz")
    assert_equal(parsed["program"], "phenix.refine")
    assert_equal(parsed["command"], "phenix.refine model.pdb data.mtz")

    # Check history_record
    assert_not_none(parsed["history_record"])
    assert_equal(parsed["history_record"]["program"], "phenix.refine")
    assert_equal(parsed["history_record"]["experiment_type"], "xray")
    assert_false(parsed["history_record"]["stop"])

    print("  PASSED")


def test_parse_response_v2_stop():
    """Test parsing v2 stop response."""
    print("Test: parse_response_v2_stop")

    response = create_stop_response(
        stop_reason="converged",
        reasoning="Structure is well refined",
    )

    parsed = parse_response_v2(response)

    assert_equal(parsed["program"], "STOP")
    assert_equal(parsed["command"], "STOP")
    assert_true(parsed["history_record"]["stop"])
    assert_equal(parsed["history_record"]["stop_reason"], "converged")

    print("  PASSED")


def test_serialize_deserialize_request():
    """Test request serialization round-trip."""
    print("Test: serialize_deserialize_request")

    original = create_request(
        files=["/path/to/file.mtz"],
        cycle_number=1,
        user_advice="Test advice",
        session_state={"resolution": 2.5},
    )

    # Serialize
    json_str = serialize_request(original)
    assert_true(isinstance(json_str, str))
    assert_in("api_version", json_str)

    # Deserialize
    restored = deserialize_request(json_str)

    # Compare
    assert_equal(restored["api_version"], original["api_version"])
    assert_equal(restored["files"], original["files"])
    assert_equal(restored["cycle_number"], original["cycle_number"])
    assert_equal(restored["user_advice"], original["user_advice"])
    assert_equal(restored["session_state"]["resolution"], 2.5)

    print("  PASSED")


def test_serialize_deserialize_response():
    """Test response serialization round-trip."""
    print("Test: serialize_deserialize_response")

    original = create_response(
        program="phenix.refine",
        command="phenix.refine model.pdb data.mtz",
        reasoning="Needs refinement",
        strategy={"resolution": 3.0},
    )

    # Serialize
    json_str = serialize_response(original)
    assert_true(isinstance(json_str, str))

    # Deserialize
    restored = deserialize_response(json_str)

    # Compare
    assert_equal(restored["decision"]["program"], "phenix.refine")
    assert_equal(restored["decision"]["command"], original["decision"]["command"])
    assert_equal(restored["decision"]["strategy"]["resolution"], 3.0)

    print("  PASSED")


def test_detect_api_version():
    """Test API version detection."""
    print("Test: detect_api_version")

    # v2 request
    assert_equal(detect_api_version({"api_version": "2.0"}), "2.0")
    assert_equal(detect_api_version({"api_version": "2.1"}), "2.1")

    # v1 (no version field)
    assert_equal(detect_api_version({}), "1.0")
    assert_equal(detect_api_version({"some_field": "value"}), "1.0")

    # Invalid
    assert_equal(detect_api_version(None), "1.0")
    assert_equal(detect_api_version("not a dict"), "1.0")

    print("  PASSED")


def test_is_v2_request():
    """Test v2 request detection."""
    print("Test: is_v2_request")

    # v2 requests
    assert_true(is_v2_request({"api_version": "2.0"}))
    assert_true(is_v2_request({"api_version": "2.1"}))
    assert_true(is_v2_request('{"api_version": "2.0"}'))  # JSON string

    # v1 requests
    assert_false(is_v2_request({}))
    assert_false(is_v2_request({"api_version": "1.0"}))
    assert_false(is_v2_request(None))
    assert_false(is_v2_request("invalid json"))

    print("  PASSED")


def test_backwards_compatibility_empty_session_state():
    """Test that empty session_state works (old client scenario)."""
    print("Test: backwards_compatibility_empty_session_state")

    # Old client sends minimal request without session_state
    request = {
        "api_version": "2.0",
        "files": ["/path/to/file.mtz"],
        "cycle_number": 1,
        # No session_state - simulates old client
    }

    # Apply defaults
    result = apply_request_defaults(request)

    # Should have empty but valid session_state
    assert_equal(result["session_state"], {})

    # Should be valid
    is_valid, errors = validate_request(result)
    assert_true(is_valid, f"Request without session_state should be valid: {errors}")

    print("  PASSED")


def test_backwards_compatibility_unknown_fields():
    """Test that new fields from future clients are ignored."""
    print("Test: backwards_compatibility_unknown_fields")

    # Future client sends request with new fields
    request = {
        "api_version": "2.5",  # Future version
        "files": [],
        "cycle_number": 1,
        "new_feature": "something",  # Unknown field
        "session_state": {
            "resolution": 2.0,
            "new_session_field": "ignored",  # Unknown nested field
        },
    }

    # Should still be valid (non-strict mode)
    is_valid, errors = validate_request(request)
    assert_true(is_valid, "Request with unknown fields should be valid")

    print("  PASSED")


def test_session_state_full():
    """Test that full session_state is properly handled."""
    print("Test: session_state_full")

    # Build request with full session state
    request = build_request_v2(
        files=["/path/to/data.mtz", "/path/to/model.pdb"],
        cycle_number=5,
        log_content="Refinement completed...",
        history=[
            {"cycle_number": 1, "program": "phenix.xtriage"},
            {"cycle_number": 2, "program": "phenix.phaser"},
        ],
        session_state={
            "resolution": 2.1,
            "experiment_type": "xray",
            "rfree_mtz": "/path/to/refine_001_data.mtz",
            "best_files": {
                "model": "/path/to/refine_001_001.pdb",
                "data_mtz": "/path/to/refine_001_data.mtz",
            },
        },
        user_advice="Solve the structure",
        provider="google",
    )

    # Validate
    is_valid, errors = validate_request(request)
    assert_true(is_valid, f"Full session_state request should be valid: {errors}")

    # Check session_state fields are present
    assert_equal(request["session_state"]["resolution"], 2.1)
    assert_equal(request["session_state"]["experiment_type"], "xray")
    assert_equal(request["session_state"]["rfree_mtz"], "/path/to/refine_001_data.mtz")
    assert_equal(request["session_state"]["best_files"]["model"], "/path/to/refine_001_001.pdb")

    # Serialize and deserialize - verify nothing is lost
    json_str = serialize_request(request)
    restored = deserialize_request(json_str)

    assert_equal(restored["session_state"]["resolution"], 2.1)
    assert_equal(restored["session_state"]["experiment_type"], "xray")
    assert_equal(restored["session_state"]["rfree_mtz"], "/path/to/refine_001_data.mtz")
    assert_equal(restored["session_state"]["best_files"]["model"], "/path/to/refine_001_001.pdb")

    print("  PASSED")


def test_response_with_metadata():
    """Test that response includes metadata from server."""
    print("Test: response_with_metadata")

    # Create response with all metadata
    response = create_response(
        program="phenix.refine",
        command="phenix.refine model.pdb data.mtz",
        reasoning="Model needs refinement",
        strategy={"output_prefix": "refine_002"},
        experiment_type="xray",
        workflow_state="xray_refined",
        warnings=["Resolution is low"],
        red_flags=[],
        debug_log=["PERCEIVE: OK", "PLAN: OK"],
    )

    # Parse response
    parsed = parse_response_v2(response)

    # Check history_record has metadata
    hr = parsed["history_record"]
    assert_equal(hr["experiment_type"], "xray")
    assert_equal(hr["workflow_state"], "xray_refined")
    assert_equal(hr["warnings"], ["Resolution is low"])
    assert_equal(hr["red_flag_issues"], [])

    print("  PASSED")


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.tst_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
