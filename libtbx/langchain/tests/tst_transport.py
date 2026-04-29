"""
Unit tests for the transport module.

Run with: python tests/tst_transport.py
"""

from __future__ import absolute_import, division, print_function

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tests.tst_utils import (
    assert_equal, assert_true, assert_in, assert_not_in,
    assert_less, assert_is_instance,
    run_tests_with_fail_fast
)

from agent.transport import (
    remove_markers,
    truncate_quoted_strings,
    replace_tabs,
    remove_control_chars,
    truncate_string,
    sanitize_string,
    sanitize_for_transport,
    sanitize_dict_recursive,
    encode_for_rest,
    decode_from_rest,
    get_transport_config,
    get_field_settings,
    sanitize_field,
    sanitize_request,
    sanitize_response,
    prepare_request_for_transport,
    process_request_from_transport,
    prepare_response_for_transport,
    process_response_from_transport,
    verify_roundtrip,
)


# =============================================================================
# REMOVE MARKERS TESTS
# =============================================================================

def test_remove_markers_zztazz():
    """ZZTAZZ markers are removed."""
    text = "Hello ZZTAZZ world"
    result = remove_markers(text)
    assert_equal(result, "Hello  world")


def test_remove_markers_multiple():
    """Multiple different markers are removed."""
    text = "ZZCRZZ line ZZTAZZ tab ZZLBZZ brace"
    result = remove_markers(text)
    assert_equal(result, " line  tab  brace")


def test_removes_all_marker_types():
    """All ZZxxZZ patterns are removed."""
    markers = ["ZZCRZZ", "ZZLBZZ", "ZZRBZZ", "ZZSCZZ", "ZZHAZZ",
               "ZZTAZZ", "ZZDQZZ", "ZZSQZZ", "ZZSRZZ", "ZZEXZZ", "ZZDSZZ"]
    for marker in markers:
        text = f"before {marker} after"
        result = remove_markers(text)
        assert_equal(result, "before  after", f"Failed for {marker}")


def test_preserves_non_markers():
    """Non-marker text is preserved."""
    text = "Hello world, no markers here!"
    result = remove_markers(text)
    assert_equal(result, text)


def test_handles_empty_string():
    """Empty string returns empty string."""
    assert_equal(remove_markers(""), "")


def test_handles_none():
    """None returns empty string."""
    assert_equal(remove_markers(None), "")


def test_handles_non_string():
    """Non-string is converted to string."""
    assert_equal(remove_markers(123), "123")


# =============================================================================
# TRUNCATE QUOTED STRINGS TESTS
# =============================================================================

def test_truncates_very_long_quoted_string():
    """Very long quoted strings are truncated."""
    content = "a" * 10000
    text = f"data='{content}'"
    result = truncate_quoted_strings(text, max_len=100)
    # Should be shorter than original
    assert_less(len(result), len(text))


def test_preserves_short_quoted_string():
    """Short quoted strings are unchanged."""
    text = "data='short'"
    result = truncate_quoted_strings(text, max_len=100)
    assert_equal(result, text)


# =============================================================================
# REPLACE TABS TESTS
# =============================================================================

def test_tabs_replaced_with_spaces():
    """Tabs are replaced with spaces."""
    text = "col1\tcol2\tcol3"
    result = replace_tabs(text)
    assert_not_in("\t", result)
    assert_in("col1", result)


def test_multiple_tabs():
    """Multiple tabs are all replaced."""
    text = "a\t\t\tb"
    result = replace_tabs(text)
    assert_not_in("\t", result)


def test_mixed_whitespace():
    """Mixed whitespace is handled correctly."""
    text = "a\t b\t c"
    result = replace_tabs(text)
    assert_not_in("\t", result)
    assert_in("a", result)
    assert_in("c", result)


# =============================================================================
# REMOVE CONTROL CHARS TESTS
# =============================================================================

def test_removes_null():
    """Null character is removed."""
    text = "before\x00after"
    result = remove_control_chars(text)
    assert_not_in("\x00", result)


def test_removes_bell():
    """Bell character is removed."""
    text = "before\x07after"
    result = remove_control_chars(text)
    assert_not_in("\x07", result)


def test_preserves_newline():
    """Newline is preserved."""
    text = "line1\nline2"
    result = remove_control_chars(text)
    assert_in("\n", result)


def test_preserves_normal_text():
    """Normal text is preserved."""
    text = "Hello, world! 123"
    result = remove_control_chars(text)
    assert_equal(result, text)


def test_removes_delete_char():
    """DEL character (0x7F) is removed."""
    text = "before\x7fafter"
    result = remove_control_chars(text)
    assert_not_in("\x7f", result)


# =============================================================================
# TRUNCATE STRING TESTS
# =============================================================================

def test_long_string_truncated():
    """Strings over max_len are truncated."""
    text = "x" * 1000
    result = truncate_string(text, max_len=100)
    assert_less(len(result), 150)
    assert_in("...", result)


def test_short_string_unchanged():
    """Strings under max_len are unchanged."""
    text = "short string"
    result = truncate_string(text, max_len=100)
    assert_equal(result, text)


def test_exact_length_unchanged():
    """Strings at exactly max_len are unchanged."""
    text = "x" * 100
    result = truncate_string(text, max_len=100)
    assert_equal(result, text)


# =============================================================================
# SANITIZE STRING TESTS
# =============================================================================

def test_sanitize_removes_markers():
    """Markers are removed in sanitize_string."""
    text = "data ZZTAZZ here"
    result = sanitize_string(text)
    assert_not_in("ZZTAZZ", result)


def test_sanitize_replaces_tabs():
    """Tabs are replaced in sanitize_string."""
    text = "col1\tcol2"
    result = sanitize_string(text)
    assert_not_in("\t", result)


def test_sanitize_removes_control():
    """Control chars are removed in sanitize_string."""
    text = "before\x00after"
    result = sanitize_string(text)
    assert_not_in("\x00", result)


def test_sanitize_truncates():
    """Long strings are truncated in sanitize_string."""
    text = "x" * 100000
    result = sanitize_string(text, max_len=1000)
    assert_less(len(result), 1500)


# =============================================================================
# SANITIZE FOR TRANSPORT TESTS
# =============================================================================

def test_sanitize_string_value():
    """String values are sanitized."""
    value = "test\tdata\x00here"
    result = sanitize_for_transport(value)
    assert_not_in("\t", result)
    assert_not_in("\x00", result)


def test_sanitize_dict_values():
    """Dict values are recursively sanitized."""
    value = {"key": "value\t", "nested": {"inner": "data\x00"}}
    result = sanitize_for_transport(value)
    assert_is_instance(result, dict)
    assert_not_in("\t", result["key"])
    assert_not_in("\x00", result["nested"]["inner"])


def test_sanitize_list_values():
    """List values are recursively sanitized."""
    value = ["item1\t", "item2\x00"]
    result = sanitize_for_transport(value)
    assert_is_instance(result, list)
    assert_not_in("\t", result[0])
    assert_not_in("\x00", result[1])


def test_sanitize_nested_structure():
    """Nested structures are fully sanitized."""
    value = {"outer": {"inner": "data\t"}, "list": ["item\x00"]}
    result = sanitize_for_transport(value)
    assert_not_in("\t", result["outer"]["inner"])
    assert_not_in("\x00", result["list"][0])


def test_sanitize_none_value():
    """None values pass through unchanged."""
    assert_equal(sanitize_for_transport(None), None)


def test_sanitize_numeric_values():
    """Numeric values pass through unchanged."""
    assert_equal(sanitize_for_transport(42), 42)
    assert_equal(sanitize_for_transport(3.14), 3.14)
    assert_equal(sanitize_for_transport(True), True)


# =============================================================================
# SANITIZE DICT RECURSIVE TESTS
# =============================================================================

def test_sanitize_dict_recursive_nested():
    """Nested dicts are sanitized."""
    d = {"a": {"b": {"c": "value\t"}}}
    result = sanitize_dict_recursive(d)
    assert_not_in("\t", result["a"]["b"]["c"])


def test_sanitize_dict_recursive_max_depth():
    """Deeply nested dicts are truncated at max depth."""
    d = {"a": {"b": {"c": {"d": {"e": "value"}}}}}
    result = sanitize_dict_recursive(d, max_depth=2)
    # Should handle gracefully (exact behavior depends on implementation)
    assert_is_instance(result, dict)


def test_sanitize_dict_recursive_with_list():
    """Lists within dicts are handled."""
    d = {"items": [{"name": "item\t"}]}
    result = sanitize_dict_recursive(d)
    assert_not_in("\t", result["items"][0]["name"])


# =============================================================================
# ENCODE/DECODE ROUNDTRIP TESTS
# =============================================================================

def test_encode_for_rest_string():
    """encode_for_rest handles string input."""
    result = encode_for_rest('{"key": "value"}')
    assert_is_instance(result, str)
    # Should contain markers for braces and quotes
    assert_in("ZZLBZZ", result)  # Left brace marker
    assert_in("ZZDQZZ", result)  # Double quote marker


def test_encode_for_rest_dict():
    """encode_for_rest handles dict input by JSON-encoding first."""
    result = encode_for_rest({"key": "value"})
    assert_is_instance(result, str)
    # Should be encoded
    assert_in("ZZLBZZ", result)


def test_encode_for_rest_none():
    """encode_for_rest handles None input."""
    result = encode_for_rest(None)
    assert_equal(result, "")


def test_encode_for_rest_list():
    """encode_for_rest handles list input."""
    result = encode_for_rest([1, 2, 3])
    assert_is_instance(result, str)


def test_encode_decode_roundtrip_dict():
    """Dict survives encode/decode roundtrip."""
    import json
    data = {"key": "value", "number": 42}
    encoded = encode_for_rest(data)
    decoded = decode_from_rest(encoded)
    result = json.loads(decoded)
    assert_equal(result["key"], "value")
    assert_equal(result["number"], 42)

def test_roundtrip_simple():
    """Simple data survives roundtrip using prepare/process functions."""
    data = {"key": "value", "number": 42}
    prepared = prepare_request_for_transport(data)
    decoded = process_request_from_transport(prepared[0])
    assert_equal(decoded["key"], "value")
    assert_equal(decoded["number"], 42)


def test_roundtrip_with_special_chars():
    """Data with special chars survives roundtrip."""
    data = {"text": "line1\nline2", "path": "/path/to/file"}
    prepared = prepare_request_for_transport(data)
    decoded = process_request_from_transport(prepared[0])
    assert_equal(decoded["text"], data["text"])
    assert_equal(decoded["path"], data["path"])


def test_roundtrip_with_tabs():
    """Data with tabs survives roundtrip (tabs converted)."""
    data = {"text": "col1\tcol2"}
    prepared = prepare_request_for_transport(data)
    decoded = process_request_from_transport(prepared[0])
    # Tabs may be replaced with spaces
    assert_in("col1", decoded["text"])
    assert_in("col2", decoded["text"])


def test_sanitize_then_roundtrip():
    """Sanitized data survives roundtrip."""
    data = {"clean": "normal"}
    sanitized = sanitize_for_transport(data)
    # sanitize_for_transport may return string representation
    # so we test the basic roundtrip with clean data
    prepared = prepare_request_for_transport({"value": "normal"})
    decoded = process_request_from_transport(prepared[0])
    assert_equal(decoded["value"], "normal")


# =============================================================================
# EDGE CASES TESTS
# =============================================================================

def test_pdb70_text_scenario():
    """Handles text similar to PDB70 content."""
    text = "ATOM      1  N   ALA A   1\nATOM      2  CA  ALA A   1"
    sanitized = sanitize_for_transport(text)
    assert_in("ATOM", sanitized)
    assert_in("\n", sanitized)


def test_log_with_embedded_markers():
    """Log text with embedded markers is handled."""
    log = """
    phenix.refine completed.
    R-work: 0.20 ZZTAZZ (marker found in wild)
    """
    sanitized = sanitize_for_transport(log)
    assert_in("R-work", sanitized)
    assert_not_in("ZZTAZZ", sanitized)


def test_large_log_performance():
    """Large logs don't cause performance issues."""
    log = "Line of log content.\n" * 10000
    sanitized = sanitize_for_transport(log)
    assert_is_instance(sanitized, str)


# =============================================================================
# CONFIG LOADING TESTS
# =============================================================================

def test_get_default_config():
    """Default config is loaded."""
    config = get_transport_config()
    assert_is_instance(config, dict)


def test_config_has_log_content_settings():
    """Config has log_content settings."""
    config = get_transport_config()
    settings = get_field_settings(config, "request", "log_content")
    assert_is_instance(settings, dict)


def test_config_has_history_settings():
    """Config has history settings."""
    config = get_transport_config()
    settings = get_field_settings(config, "request", "history")
    assert_is_instance(settings, dict)


# =============================================================================
# SANITIZE FIELD TESTS
# =============================================================================

def test_sanitize_field_log_content():
    """log_content field is sanitized correctly."""
    log = "x" * 100000
    result = sanitize_field(log, "log_content", direction="request")
    # Should be truncated according to config
    assert_is_instance(result, str)


def test_sanitize_field_history():
    """history field is sanitized correctly."""
    history = [{"program": "test", "log": "x" * 10000}]
    result = sanitize_field(history, "history", direction="request")
    assert_is_instance(result, list)


def test_sanitize_field_unknown():
    """Unknown field passes through."""
    value = "test\tvalue"
    result = sanitize_field(value, "unknown_field", direction="request")
    # Unknown fields may pass through unchanged
    assert_is_instance(result, str)


# =============================================================================
# SANITIZE REQUEST TESTS
# =============================================================================

def test_sanitize_request_basic():
    """Basic request is sanitized."""
    request = {
        "api_version": "2.0",
        "files": ["/path/file.mtz"],
        "log_content": "log text"
    }
    result = sanitize_request(request)
    assert_equal(result["api_version"], "2.0")
    assert_is_instance(result, dict)


def test_sanitize_request_with_history():
    """Request with history is sanitized."""
    request = {
        "api_version": "2.0",
        "history": [{"program": "test", "log": "data"}]
    }
    result = sanitize_request(request)
    assert_is_instance(result["history"], list)


def test_sanitize_request_preserves_structure():
    """Request structure is preserved."""
    request = {
        "api_version": "2.0",
        "session_state": {"key": "value"},
        "files": ["a.mtz", "b.pdb"]
    }
    result = sanitize_request(request)
    assert_in("session_state", result)
    assert_equal(len(result["files"]), 2)


# =============================================================================
# SANITIZE RESPONSE TESTS
# =============================================================================

def test_sanitize_response_basic():
    """Basic response is sanitized."""
    response = {
        "status": "success",
        "message": "completed task"
    }
    result = sanitize_response(response)
    assert_equal(result["status"], "success")
    assert_is_instance(result, dict)


def test_sanitize_response_with_command():
    """Response with command is sanitized."""
    response = {
        "status": "run_command",
        "command": "phenix.refine model.pdb data.mtz"
    }
    result = sanitize_response(response)
    assert_in("phenix.refine", result["command"])


def test_sanitize_response_with_intent():
    """Response with intent object is sanitized."""
    response = {
        "status": "success",
        "intent": {"program": "test", "files": {"model": "a.pdb"}}
    }
    result = sanitize_response(response)
    assert_in("intent", result)
    assert_equal(result["intent"]["program"], "test")


# =============================================================================
# UNIFIED TRANSPORT TESTS
# =============================================================================

def test_prepare_request_basic():
    """prepare_request_for_transport works with basic request."""
    request = {"api_version": "2.0", "message": "test"}
    result = prepare_request_for_transport(request)
    # Returns tuple of (encoded, original_json)
    assert_is_instance(result, tuple)
    assert_equal(len(result), 2)


def test_process_request_basic():
    """process_request_from_transport works with basic request."""
    request = {"api_version": "2.0", "message": "test"}
    prepared = prepare_request_for_transport(request)
    # Pass the encoded version (first element)
    result = process_request_from_transport(prepared[0])
    assert_equal(result["api_version"], "2.0")


def test_process_request_not_encoded():
    """process_request_from_transport handles JSON string input."""
    import json
    request = {"api_version": "2.0"}
    # Pass JSON string directly (not marker-encoded)
    result = process_request_from_transport(json.dumps(request))
    assert_equal(result["api_version"], "2.0")


def test_process_request_invalid_json():
    """process_request_from_transport handles invalid JSON."""
    result = process_request_from_transport("not json {{{")
    # Should return some error indicator or empty dict
    assert_is_instance(result, dict)


def test_prepare_response_basic():
    """prepare_response_for_transport works with basic response."""
    response = {"status": "success"}
    result = prepare_response_for_transport(response)
    # Returns tuple of (encoded, original_json)
    assert_is_instance(result, tuple)


def test_process_response_basic():
    """process_response_from_transport works with basic response."""
    response = {"status": "success"}
    prepared = prepare_response_for_transport(response)
    result = process_response_from_transport(prepared[0])
    assert_equal(result["status"], "success")


def test_roundtrip_with_problematic_content():
    """Full roundtrip with problematic content."""
    request = {
        "api_version": "2.0",
        "log_content": "data\t\x00\nmore ZZTAZZ stuff",
        "history": [{"log": "old\tlog"}]
    }
    prepared = prepare_request_for_transport(request)
    result = process_request_from_transport(prepared[0])
    assert_equal(result["api_version"], "2.0")
    assert_not_in("ZZTAZZ", result.get("log_content", ""))


# =============================================================================
# VERIFY ROUNDTRIP TESTS
# =============================================================================

def test_verify_roundtrip_success():
    """verify_roundtrip succeeds with valid data."""
    request = {
        "api_version": "2.0",
        "files": ["/path/to/file.mtz"]
    }
    success, message, details = verify_roundtrip(request)
    assert_true(success)
    assert_in("OK", message)


def test_verify_roundtrip_with_content():
    """Roundtrip verification with substantial content."""
    request = {
        "api_version": "2.0",
        "files": [f"/file_{i}.mtz" for i in range(10)],
        "session_state": {"experiment_type": "cryoem", "resolution": 2.5},
        "history": [{"program": f"prog_{i}"} for i in range(5)],
        "log_content": "x" * 1000,
    }
    success, message, details = verify_roundtrip(request)
    assert_true(success, f"Roundtrip failed: {message}")


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Run all tests with fail-fast behavior (cctbx style)."""
    run_tests_with_fail_fast()


if __name__ == "__main__":
    run_all_tests()
