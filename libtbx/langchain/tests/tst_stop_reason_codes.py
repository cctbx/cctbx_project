"""
STOP_REASON_CODES: parse_assessment validation, think_stop_override, think_file_overrides (P1B).

Tests:
  A. STOP_REASON_CODES constant
  B. parse_assessment: stop_reason_code validation
  C. parse_assessment: file_overrides validation
  D. thinking_agent stop path: think_stop_override populated
  E. thinking_agent stop path: freeform fallback when no code
  F. thinking_agent non-stop path: think_file_overrides populated
  G. thinking_agent non-stop path: empty overrides not set
  H. BUILD consumer: think_file_overrides merged into categorized_files
  I. BUILD consumer: think_file_overrides cleared after use
  J. BUILD consumer: empty overrides — no change to categorized_files
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in [_ROOT, os.path.join(_ROOT, "agent"),
           os.path.join(_ROOT, "knowledge")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

from graph_state import STOP_REASON_CODES
from knowledge.thinking_prompts import parse_assessment


# ============================================================================
# A. STOP_REASON_CODES
# ============================================================================
def test_stop_reason_codes_is_frozenset():
    assert isinstance(STOP_REASON_CODES, frozenset)


def test_stop_reason_codes_contains_expected():
    for code in ("WRONG_MTZ", "WRONG_SPACE_GROUP", "MISMATCHED_SEQUENCE",
                 "NO_SOLUTION_FOUND", "REASONLESS_DIVERGENCE"):
        assert code in STOP_REASON_CODES, "%s missing" % code


def test_stop_reason_codes_no_unknown():
    # Sanity: all codes are upper-case strings
    for code in STOP_REASON_CODES:
        assert isinstance(code, str)
        assert code == code.upper()


# ============================================================================
# B. parse_assessment: stop_reason_code validation
# ============================================================================
def test_parse_assessment_valid_code_kept():
    import json
    raw = json.dumps({
        "action": "stop",
        "analysis": "Wrong MTZ",
        "confidence": "high",
        "stop_reason_code": "WRONG_MTZ",
    })
    result = parse_assessment(raw)
    assert result["stop_reason_code"] == "WRONG_MTZ"


def test_parse_assessment_unknown_code_cleared():
    import json
    raw = json.dumps({
        "action": "stop",
        "analysis": "bad",
        "stop_reason_code": "INVENTED_CODE_XYZ",
    })
    result = parse_assessment(raw)
    assert result["stop_reason_code"] is None, \
        "Unknown code should be cleared, got %r" % result["stop_reason_code"]


def test_parse_assessment_no_code_defaults_none():
    import json
    raw = json.dumps({"action": "stop", "analysis": "expert stop"})
    result = parse_assessment(raw)
    assert result["stop_reason_code"] is None


def test_parse_assessment_all_valid_codes_accepted():
    import json
    for code in STOP_REASON_CODES:
        raw = json.dumps({"action": "stop", "stop_reason_code": code})
        result = parse_assessment(raw)
        assert result["stop_reason_code"] == code, \
            "Code %s was not kept" % code


# ============================================================================
# C. parse_assessment: file_overrides validation
# ============================================================================
def test_parse_assessment_file_overrides_dict_kept():
    import json
    raw = json.dumps({
        "action": "let_run",
        "file_overrides": {"data_mtz": "/work/good.mtz"},
    })
    result = parse_assessment(raw)
    assert result["file_overrides"] == {"data_mtz": "/work/good.mtz"}


def test_parse_assessment_file_overrides_null_normalised():
    import json
    raw = json.dumps({"action": "let_run", "file_overrides": None})
    result = parse_assessment(raw)
    assert result["file_overrides"] == {}


def test_parse_assessment_file_overrides_list_normalised():
    import json
    raw = json.dumps({"action": "let_run", "file_overrides": ["bad"]})
    result = parse_assessment(raw)
    assert result["file_overrides"] == {}


def test_parse_assessment_file_overrides_default_empty():
    import json
    raw = json.dumps({"action": "let_run"})
    result = parse_assessment(raw)
    assert result["file_overrides"] == {}


# ============================================================================
# D-G. thinking_agent stop / non-stop paths  (source-file inspection)
# ============================================================================
_TA_PATH = os.path.join(_ROOT, "agent", "thinking_agent.py")

def _ta_src():
    with open(_TA_PATH) as f:
        return f.read()


def test_thinking_agent_stop_path_sets_think_stop_override():
    src = _ta_src()
    assert "think_stop_override" in src, \
        "think_stop_override not referenced in thinking_agent.py"
    assert "_think_stop_override" in src


def test_thinking_agent_stop_path_uses_code_as_stop_reason():
    src = _ta_src()
    assert "stop_reason = _src" in src, \
        "stop_reason should be assigned _src in the classified stop path"


def test_thinking_agent_stop_path_freeform_fallback():
    src = _ta_src()
    assert "expert: %s" in src, \
        "Freeform 'expert: ...' fallback string missing"


def test_thinking_agent_nonstp_sets_think_file_overrides():
    src = _ta_src()
    assert "think_file_overrides" in src, \
        "think_file_overrides not set in non-stop return path"


def test_thinking_agent_nonstp_uses_file_overrides_var():
    src = _ta_src()
    assert "_file_overrides" in src


# ============================================================================
# H-J. BUILD consumer (_build_with_new_builder) — source-file inspection
# ============================================================================
_GN_PATH = os.path.join(_ROOT, "agent", "graph_nodes.py")

def _gn_build_src():
    """Extract just the _build_with_new_builder function text."""
    with open(_GN_PATH) as f:
        content = f.read()
    start = content.find("def _build_with_new_builder(")
    # next top-level def
    end = content.find("\ndef ", start + 1)
    return content[start:end] if end > start else content[start:]


def test_build_consumer_references_think_file_overrides():
    src = _gn_build_src()
    assert "think_file_overrides" in src, \
        "think_file_overrides not referenced in _build_with_new_builder"


def test_build_consumer_uses_think_overrides_var():
    src = _gn_build_src()
    assert "_think_overrides" in src
    assert "_categorized_files" in src


def test_build_consumer_clears_after_use():
    src = _gn_build_src()
    assert '"think_file_overrides": {}' in src, \
        "think_file_overrides should be cleared (set to {}) after use"


def test_build_consumer_logs_override():
    src = _gn_build_src()
    assert "think_file_override" in src, \
        "Consumer should log the override it applies"


def test_build_consumer_uses_patched_categorized_files():
    src = _gn_build_src()
    assert "categorized_files=_categorized_files" in src, \
        "CommandContext should use _categorized_files (patched), not original"


# ============================================================================
# Runner
# ============================================================================
_TESTS = [
    # A
    test_stop_reason_codes_is_frozenset,
    test_stop_reason_codes_contains_expected,
    test_stop_reason_codes_no_unknown,
    # B
    test_parse_assessment_valid_code_kept,
    test_parse_assessment_unknown_code_cleared,
    test_parse_assessment_no_code_defaults_none,
    test_parse_assessment_all_valid_codes_accepted,
    # C
    test_parse_assessment_file_overrides_dict_kept,
    test_parse_assessment_file_overrides_null_normalised,
    test_parse_assessment_file_overrides_list_normalised,
    test_parse_assessment_file_overrides_default_empty,
    # D-G (source inspection)
    test_thinking_agent_stop_path_sets_think_stop_override,
    test_thinking_agent_stop_path_uses_code_as_stop_reason,
    test_thinking_agent_stop_path_freeform_fallback,
    test_thinking_agent_nonstp_sets_think_file_overrides,
    test_thinking_agent_nonstp_uses_file_overrides_var,
    # H-J
    test_build_consumer_references_think_file_overrides,
    test_build_consumer_uses_think_overrides_var,
    test_build_consumer_clears_after_use,
    test_build_consumer_logs_override,
    test_build_consumer_uses_patched_categorized_files,
]


def run_all_tests():
    for test_fn in _TESTS:
        test_fn()
    print("All %d tests passed." % len(_TESTS))

if __name__ == "__main__":
    passed = failed = 0
    for test_fn in _TESTS:
        try:
            test_fn()
            print("  PASS: %s" % test_fn.__name__)
            passed += 1
        except Exception:
            import traceback
            print("  FAIL: %s" % test_fn.__name__)
            traceback.print_exc()
            failed += 1
    print()
    if failed:
        print("%d/%d tests FAILED." % (failed, passed + failed))
        sys.exit(1)
    else:
        print("All %d tests passed." % passed)
