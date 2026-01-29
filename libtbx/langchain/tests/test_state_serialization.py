"""
Tests for State Serialization and Local/Server Equivalence.

Test 1: State round-trip - verify state can be packaged and unpackaged without changes
Test 2: Local vs Server equivalence - verify identical behavior except for where analysis runs

Run with: python tests/test_state_serialization.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import json
import copy

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Mock libtbx imports
import types
libtbx = types.ModuleType('libtbx')
libtbx.langchain = types.ModuleType('libtbx.langchain')
libtbx.langchain.agent = types.ModuleType('libtbx.langchain.agent')
libtbx.langchain.knowledge = types.ModuleType('libtbx.langchain.knowledge')
sys.modules['libtbx'] = libtbx
sys.modules['libtbx.langchain'] = libtbx.langchain
sys.modules['libtbx.langchain.agent'] = libtbx.langchain.agent
sys.modules['libtbx.langchain.knowledge'] = libtbx.langchain.knowledge

import knowledge.yaml_loader
libtbx.langchain.knowledge.yaml_loader = knowledge.yaml_loader
sys.modules['libtbx.langchain.knowledge.yaml_loader'] = knowledge.yaml_loader

# Import API schema functions
from knowledge.api_schema import (
    apply_request_defaults,
    validate_request,
    validate_response,
)

# Import transport functions
from agent.transport import (
    sanitize_request,
    sanitize_response,
    get_transport_config,
    prepare_request_for_transport,
    process_request_from_transport,
    prepare_response_for_transport,
    process_response_from_transport,
    verify_roundtrip,
)


# =============================================================================
# TEST UTILITIES
# =============================================================================

def assert_equal(a, b, msg=""):
    if a != b:
        raise AssertionError("%s: %r != %r" % (msg, a, b) if msg else "%r != %r" % (a, b))


def assert_true(condition, msg="Assertion failed"):
    if not condition:
        raise AssertionError(msg)


def deep_compare(obj1, obj2, path=""):
    """
    Deep compare two objects, returning list of differences.
    """
    differences = []

    if type(obj1) != type(obj2):
        differences.append("%s: type mismatch %s vs %s" % (path, type(obj1).__name__, type(obj2).__name__))
        return differences

    if isinstance(obj1, dict):
        all_keys = set(obj1.keys()) | set(obj2.keys())
        for key in all_keys:
            new_path = "%s.%s" % (path, key) if path else key
            if key not in obj1:
                differences.append("%s: missing in first" % new_path)
            elif key not in obj2:
                differences.append("%s: missing in second" % new_path)
            else:
                differences.extend(deep_compare(obj1[key], obj2[key], new_path))
    elif isinstance(obj1, list):
        if len(obj1) != len(obj2):
            differences.append("%s: length mismatch %d vs %d" % (path, len(obj1), len(obj2)))
        else:
            for i, (item1, item2) in enumerate(zip(obj1, obj2)):
                differences.extend(deep_compare(item1, item2, "%s[%d]" % (path, i)))
    elif isinstance(obj1, float):
        # Float comparison with tolerance
        if abs(obj1 - obj2) > 1e-9:
            differences.append("%s: %r != %r" % (path, obj1, obj2))
    else:
        if obj1 != obj2:
            differences.append("%s: %r != %r" % (path, obj1, obj2))

    return differences


# =============================================================================
# TEST 1: STATE ROUND-TRIP SERIALIZATION
# =============================================================================

def create_sample_history():
    """Create a realistic history for testing."""
    return [
        {
            "cycle_number": 1,
            "job_id": "1",
            "program": "phenix.xtriage",
            "command": "phenix.xtriage data.mtz",
            "result": "SUCCESS: Analysis complete",
            "summary": "Decision: Analyze data\nCommand: phenix.xtriage data.mtz\nResult: SUCCESS",
            "analysis": {
                "resolution": 2.5,
                "completeness": 0.98,
                "i_over_sigma": 15.2,
            },
            "output_files": ["/path/to/xtriage_output.log"],
        },
        {
            "cycle_number": 2,
            "job_id": "2",
            "program": "phenix.refine",
            "command": "phenix.refine model.pdb data.mtz output.prefix=refine_001",
            "result": "SUCCESS: Refinement complete",
            "summary": "Decision: Refine model\nCommand: phenix.refine ...\nResult: SUCCESS",
            "analysis": {
                "r_work": 0.22,
                "r_free": 0.26,
                "resolution": 2.5,
            },
            "output_files": [
                "/path/to/refine_001_001.pdb",
                "/path/to/refine_001_001.mtz",
            ],
        },
        {
            "cycle_number": 3,
            "job_id": "3",
            "program": "phenix.refine",
            "command": "phenix.refine refine_001_001.pdb data.mtz output.prefix=refine_002",
            "result": "SUCCESS: Refinement complete",
            "summary": "Decision: Continue refinement\nCommand: phenix.refine ...\nResult: SUCCESS",
            "analysis": {
                "r_work": 0.20,
                "r_free": 0.24,
                "resolution": 2.5,
            },
            "output_files": [
                "/path/to/refine_002_001.pdb",
                "/path/to/refine_002_001.mtz",
            ],
        },
    ]


def create_sample_session_state():
    """Create a realistic session state for testing."""
    return {
        "resolution": 2.5,
        "experiment_type": "xray",
        "rfree_mtz": "/path/to/rfree_data.mtz",
        "best_files": {
            "model": "/path/to/refine_002_001.pdb",
            "data_mtz": "/path/to/rfree_data.mtz",
        },
    }


def create_sample_request():
    """Create a complete sample request."""
    return {
        "api_version": "2.0",
        "files": [
            "/path/to/model.pdb",
            "/path/to/data.mtz",
            "/path/to/sequence.fa",
        ],
        "history": create_sample_history(),
        "cycle_number": 4,
        "session_state": create_sample_session_state(),
        "user_advice": "Focus on R-free improvement",
        "settings": {
            "provider": "google",
            "abort_on_red_flags": True,
            "abort_on_warnings": False,
            "max_cycles": 20,
        },
    }


def create_sample_response():
    """Create a complete sample response."""
    return {
        "api_version": "2.0",
        "command": "phenix.refine refine_002_001.pdb data.mtz output.prefix=refine_003",
        "reasoning": "Continuing refinement to improve R-free. Current R-free is 0.24, targeting below 0.22.",
        "stop": False,
        "stop_reason": None,
        "decision": {
            "program": "phenix.refine",
            "command": "phenix.refine refine_002_001.pdb data.mtz output.prefix=refine_003",
            "reasoning": "Continuing refinement to improve R-free",
            "strategy": {"output_prefix": "refine_003"},
            "confidence": "high",
        },
        "history_record": {
            "cycle_number": 4,
            "program": "phenix.refine",
            "command": "phenix.refine refine_002_001.pdb data.mtz output.prefix=refine_003",
            "decision": "Continue refinement",
            "reasoning": "R-free improving, continue optimization",
            "experiment_type": "xray",
        },
        "metadata": {
            "cycle_number": 4,
            "experiment_type": "xray",
            "workflow_state": "xray_refined",
            "valid_programs": ["phenix.refine", "phenix.molprobity", "STOP"],
        },
        "session_state": {
            "resolution": 2.5,
            "experiment_type": "xray",
            "rfree_mtz": "/path/to/rfree_data.mtz",
            "best_files": {
                "model": "/path/to/refine_002_001.pdb",
                "data_mtz": "/path/to/rfree_data.mtz",
            },
        },
        "debug_log": [
            "PERCEIVE: Analyzing cycle 3 results",
            "PLAN: R-free=0.24, continuing refinement",
            "BUILD: output_prefix=refine_003",
        ],
    }


def test_request_round_trip():
    """Test that request can be serialized and deserialized without changes."""
    print("Test: request_round_trip")

    original = create_sample_request()

    # Step 1: Apply defaults (simulates what server does on receive)
    with_defaults = apply_request_defaults(copy.deepcopy(original))

    # Step 2: Serialize to JSON
    json_str = json.dumps(with_defaults, sort_keys=True)

    # Step 3: Deserialize
    restored = json.loads(json_str)

    # Step 4: Apply defaults again (should be idempotent)
    restored_with_defaults = apply_request_defaults(restored)

    # Compare
    differences = deep_compare(with_defaults, restored_with_defaults)

    if differences:
        print("  DIFFERENCES FOUND:")
        for diff in differences[:10]:
            print("    - %s" % diff)
        raise AssertionError("Request round-trip failed: %d differences" % len(differences))

    print("  PASSED")


def test_response_round_trip():
    """Test that response can be serialized and deserialized without changes."""
    print("Test: response_round_trip")

    original = create_sample_response()

    # Step 1: Serialize to JSON (without apply_response_defaults to avoid adding optional fields)
    json_str = json.dumps(original, sort_keys=True)

    # Step 2: Deserialize
    restored = json.loads(json_str)

    # Step 3: Compare key fields that must be preserved
    key_fields = ["api_version", "command", "reasoning", "stop", "stop_reason",
                  "history_record", "metadata", "session_state"]

    for field in key_fields:
        if field in original:
            differences = deep_compare(original.get(field), restored.get(field), field)
            if differences:
                print("  DIFFERENCES in %s:" % field)
                for diff in differences[:5]:
                    print("    - %s" % diff)
                raise AssertionError("Response field '%s' changed" % field)

    print("  PASSED")


def test_history_round_trip():
    """Test that history survives multiple round trips."""
    print("Test: history_round_trip")

    original_history = create_sample_history()

    # Simulate multiple cycles of request/response
    history = copy.deepcopy(original_history)

    for cycle in range(3):
        # Create request with history
        request = {
            "api_version": "2.0",
            "files": ["/path/to/file.pdb"],
            "history": history,
            "cycle_number": len(history) + 1,
            "session_state": {},
            "settings": {},
        }

        # Serialize/deserialize (simulates network transport)
        request_json = json.dumps(request)
        request_restored = json.loads(request_json)

        # Apply defaults
        request_with_defaults = apply_request_defaults(request_restored)

        # Verify history preserved
        assert_equal(len(request_with_defaults["history"]), len(history),
                    "History length changed in cycle %d" % cycle)

        # Add new history entry (simulates agent response)
        new_entry = {
            "cycle_number": len(history) + 1,
            "job_id": str(len(history) + 1),
            "program": "phenix.refine",
            "command": "phenix.refine ...",
            "result": "SUCCESS",
            "analysis": {"r_free": 0.24 - cycle * 0.01},
        }
        history = request_with_defaults["history"] + [new_entry]

    # Final check - all original entries still present
    for i, original_entry in enumerate(original_history):
        restored_entry = history[i]
        for key in original_entry:
            if key in restored_entry:
                differences = deep_compare(original_entry[key], restored_entry[key], key)
                if differences:
                    raise AssertionError("History entry %d.%s changed: %s" % (i, key, differences))

    print("  PASSED")


def test_transport_sanitization_round_trip():
    """Test that transport sanitization preserves data."""
    print("Test: transport_sanitization_round_trip")

    config = get_transport_config()

    # Test request
    original_request = create_sample_request()
    sanitized_request = sanitize_request(copy.deepcopy(original_request), config)

    # Key fields should be preserved
    assert_equal(sanitized_request["api_version"], original_request["api_version"])
    assert_equal(sanitized_request["cycle_number"], original_request["cycle_number"])
    assert_equal(len(sanitized_request["history"]), len(original_request["history"]))

    # Test response
    original_response = create_sample_response()
    sanitized_response = sanitize_response(copy.deepcopy(original_response), config)

    # Key fields should be preserved
    assert_equal(sanitized_response["api_version"], original_response["api_version"])
    assert_equal(sanitized_response["command"], original_response["command"])
    assert_equal(sanitized_response["stop"], original_response["stop"])

    print("  PASSED")


def test_session_state_preservation():
    """Test that session state is preserved across cycles."""
    print("Test: session_state_preservation")

    # Initial session state
    session_state = {
        "resolution": 2.5,
        "experiment_type": "xray",
        "rfree_mtz": "/path/to/rfree.mtz",
        "best_files": {
            "model": "/path/to/model.pdb",
            "data_mtz": "/path/to/data.mtz",
        },
    }

    # Simulate 5 cycles
    for cycle in range(5):
        # Create request
        request = {
            "api_version": "2.0",
            "files": [],
            "history": [],
            "cycle_number": cycle + 1,
            "session_state": copy.deepcopy(session_state),
            "settings": {},
        }

        # Round-trip through JSON
        restored = json.loads(json.dumps(request))
        restored = apply_request_defaults(restored)

        # Verify session state preserved
        for key in session_state:
            if key in restored["session_state"]:
                assert_equal(
                    restored["session_state"][key],
                    session_state[key],
                    "session_state.%s changed in cycle %d" % (key, cycle)
                )

        # Update session state (simulates what happens after each cycle)
        session_state["best_files"]["model"] = "/path/to/refine_%03d.pdb" % (cycle + 1)

    print("  PASSED")


# =============================================================================
# TEST 2: LOCAL VS SERVER EQUIVALENCE
# =============================================================================

def test_request_format_equivalence():
    """Test that local and server use the same request format."""
    print("Test: request_format_equivalence")

    # The request format should be identical whether running locally or via server
    # This tests that build_request_v2 produces valid requests

    from agent.api_client import build_request_v2

    # Simulate what local agent does
    files = ["/path/to/model.pdb", "/path/to/data.mtz"]
    history = create_sample_history()
    cycle_number = 4

    # Build request (as local agent would)
    request = build_request_v2(
        files=files,
        history=history,
        cycle_number=cycle_number,
        session_state=create_sample_session_state(),
        user_advice="Test advice",
        provider="google",
    )

    # Validate request format
    is_valid, errors = validate_request(request)
    assert_true(is_valid, "Request validation failed: %s" % errors)

    # Verify all required fields present
    assert_true("api_version" in request, "Missing api_version")
    assert_true("files" in request, "Missing files")
    assert_true("history" in request, "Missing history")
    assert_true("cycle_number" in request, "Missing cycle_number")

    print("  PASSED")


def test_response_format_equivalence():
    """Test that local and server produce the same response format."""
    print("Test: response_format_equivalence")

    # Both local and server should produce responses matching RESPONSE_V2_SCHEMA

    response = create_sample_response()

    # Validate response format
    is_valid, errors = validate_response(response)
    assert_true(is_valid, "Response validation failed: %s" % errors)

    # Verify all required fields present
    required_fields = ["api_version", "command", "stop"]
    for field in required_fields:
        assert_true(field in response, "Missing required field: %s" % field)

    print("  PASSED")


def test_history_record_format():
    """Test that history_record format is consistent."""
    print("Test: history_record_format")

    # History record should have consistent format whether from local or server

    history_record = {
        "cycle_number": 1,
        "program": "phenix.refine",
        "command": "phenix.refine model.pdb data.mtz",
        "decision": "Refine model",
        "reasoning": "Model needs refinement",
        "experiment_type": "xray",
    }

    # Should serialize cleanly
    json_str = json.dumps(history_record)
    restored = json.loads(json_str)

    # All fields preserved
    for key in history_record:
        assert_equal(restored[key], history_record[key], "history_record.%s" % key)

    print("  PASSED")


def test_graph_state_to_request_conversion():
    """Test conversion from graph state to API request."""
    print("Test: graph_state_to_request_conversion")

    # This simulates what happens when LocalAgent converts internal state to request

    # Internal graph state (what LocalAgent works with)
    graph_state = {
        "available_files": ["/path/to/model.pdb", "/path/to/data.mtz"],
        "history": create_sample_history(),
        "cycle_number": 4,
        "session_info": {
            "experiment_type": "xray",
            "best_files": {"model": "/path/to/best.pdb"},
            "rfree_mtz": "/path/to/rfree.mtz",
        },
        "session_resolution": 2.5,
        "workflow_state": {
            "state": "xray_refined",
            "valid_programs": ["phenix.refine", "STOP"],
        },
    }

    # Convert to API request format
    request = {
        "api_version": "2.0",
        "files": graph_state["available_files"],
        "history": graph_state["history"],
        "cycle_number": graph_state["cycle_number"],
        "session_state": {
            "resolution": graph_state.get("session_resolution"),
            "experiment_type": graph_state["session_info"].get("experiment_type"),
            "rfree_mtz": graph_state["session_info"].get("rfree_mtz"),
            "best_files": graph_state["session_info"].get("best_files", {}),
        },
        "settings": {},
    }

    # Should be valid
    is_valid, errors = validate_request(request)
    assert_true(is_valid, "Converted request invalid: %s" % errors)

    # Should round-trip
    restored = json.loads(json.dumps(request))
    assert_equal(restored["cycle_number"], graph_state["cycle_number"])
    assert_equal(len(restored["history"]), len(graph_state["history"]))

    print("  PASSED")


def test_response_to_graph_state_conversion():
    """Test conversion from API response back to graph state."""
    print("Test: response_to_graph_state_conversion")

    # This simulates what happens when LocalAgent receives a response

    response = create_sample_response()

    # Convert response to graph state updates
    graph_updates = {
        "command": response["command"],
        "reasoning": response["reasoning"],
        "stop": response["stop"],
        "stop_reason": response.get("stop_reason"),
    }

    # History record should be usable
    history_record = response.get("history_record", {})
    if history_record:
        graph_updates["program"] = history_record.get("program")
        graph_updates["decision"] = history_record.get("decision")

    # Session state should be extractable
    session_state = response.get("session_state", {})
    if session_state:
        graph_updates["session_resolution"] = session_state.get("resolution")
        graph_updates["session_info"] = {
            "experiment_type": session_state.get("experiment_type"),
            "rfree_mtz": session_state.get("rfree_mtz"),
            "best_files": session_state.get("best_files", {}),
        }

    # Verify conversion worked
    assert_equal(graph_updates["command"], response["command"])
    assert_equal(graph_updates["program"], "phenix.refine")
    assert_equal(graph_updates["session_resolution"], 2.5)

    print("  PASSED")


def test_full_transport_round_trip():
    """Test complete transport round-trip (simulates local <-> server)."""
    print("Test: full_transport_round_trip")

    # Create a realistic request
    original_request = create_sample_request()

    # Step 1: Prepare for transport (what client does before sending)
    encoded_request, _ = prepare_request_for_transport(original_request, do_encode=True)

    # Step 2: Process from transport (what server does on receive)
    restored_request = process_request_from_transport(encoded_request, was_encoded=True)

    # Step 3: Verify key fields preserved
    assert_equal(restored_request["api_version"], original_request["api_version"])
    assert_equal(restored_request["cycle_number"], original_request["cycle_number"])
    assert_equal(len(restored_request["history"]), len(original_request["history"]))
    assert_equal(len(restored_request["files"]), len(original_request["files"]))

    # Verify history content
    for i, (orig, rest) in enumerate(zip(original_request["history"], restored_request["history"])):
        assert_equal(rest.get("program"), orig.get("program"), "history[%d].program" % i)
        assert_equal(rest.get("cycle_number"), orig.get("cycle_number"), "history[%d].cycle_number" % i)

    # Step 4: Now test response round-trip
    original_response = create_sample_response()

    encoded_response, _ = prepare_response_for_transport(original_response, do_encode=True)
    restored_response = process_response_from_transport(encoded_response, was_encoded=True)

    # Verify key fields preserved
    assert_equal(restored_response["api_version"], original_response["api_version"])
    assert_equal(restored_response["command"], original_response["command"])
    assert_equal(restored_response["stop"], original_response["stop"])

    print("  PASSED")


def test_verify_roundtrip_function():
    """Test the built-in verify_roundtrip function."""
    print("Test: verify_roundtrip_function")

    request = create_sample_request()

    # verify_roundtrip returns (success, message, details)
    result = verify_roundtrip(request)
    success = result[0]
    message = result[1] if len(result) > 1 else ""
    details = result[2] if len(result) > 2 else {}

    if not success:
        print("  Round-trip verification failed:")
        print("    Message: %s" % message)
        print("    Details: %s" % details)
        raise AssertionError("verify_roundtrip failed: %s" % message)

    print("  PASSED")


# =============================================================================
# RUN ALL TESTS
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    from tests.test_utils import run_tests_with_fail_fast
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
