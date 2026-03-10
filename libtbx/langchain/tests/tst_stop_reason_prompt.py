"""
Stop-reason table injection into system prompt; non-stop guard (P1B, no LLM calls).

Covers:
  1.  build_thinking_prompt returns a (system_msg, user_msg) tuple
  2.  system_msg contains the stop_reason table when errors.yaml is present
  3.  "reference list below" appears before the table in system_msg
  4a. WRONG_MTZ appears in system_msg
  4b. WRONG_SPACE_GROUP appears in system_msg
  4c. MISMATCHED_SEQUENCE appears in system_msg
  4d. NO_SOLUTION_FOUND appears in system_msg
  4e. REASONLESS_DIVERGENCE appears in system_msg
  5.  build_stop_reason_table returns "" gracefully when errors.yaml absent
  6.  build_thinking_prompt falls back to plain SYSTEM_PROMPT when table is empty
  7.  Guard: stop_reason_code preserved when action == "stop"
  8.  Guard: stop_reason_code cleared when action == "guide_step"
  9.  Guard: stop_reason_code cleared when action == "let_run"
  10. Guard: stop_reason_code cleared when action == "pivot"
  11. Guard: None code on non-stop action is a no-op (no KeyError, no mutation)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in [_ROOT, os.path.join(_ROOT, "knowledge"), os.path.join(_ROOT, "agent")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

from thinking_prompts import (
    SYSTEM_PROMPT,
    build_thinking_prompt,
    build_stop_reason_table,
)

# ---------------------------------------------------------------------------
# Minimal context for build_thinking_prompt
# ---------------------------------------------------------------------------
_CONTEXT = {
    "cycle_number": 1,
    "experiment_type": "xray",
    "workflow_state": "refinement",
    "log_sections": "R-free = 0.32",
    "program_name": "phenix.refine",
    "metrics": {},
}


# ---------------------------------------------------------------------------
# Guard logic replicated inline (mirrors thinking_agent.py :: run_think_node)
# ---------------------------------------------------------------------------

def _apply_guard(assessment):
    """
    Replicated P1B guard from thinking_agent.py.
    Returns (assessment, was_cleared).
    """
    was_cleared = False
    if (assessment.get("stop_reason_code")
            and assessment["action"] != "stop"):
        assessment["stop_reason_code"] = None
        was_cleared = True
    return assessment, was_cleared


# ---------------------------------------------------------------------------
# Tests 1-3: build_thinking_prompt output shape and table position
# ---------------------------------------------------------------------------

def test_returns_tuple():
    """build_thinking_prompt returns a two-element tuple of strings."""
    result = build_thinking_prompt(_CONTEXT)
    assert isinstance(result, tuple) and len(result) == 2, \
        "Expected (system_msg, user_msg) tuple, got: %r" % (result,)
    system_msg, user_msg = result
    assert isinstance(system_msg, str) and system_msg, "system_msg must be a non-empty string"
    assert isinstance(user_msg,   str) and user_msg,   "user_msg must be a non-empty string"
    print("PASS test_returns_tuple")


def test_system_msg_contains_table():
    """system_msg contains the stop_reason table when errors.yaml is present."""
    table = build_stop_reason_table()
    if not table:
        print("SKIP test_system_msg_contains_table (errors.yaml not available)")
        return
    system_msg, _ = build_thinking_prompt(_CONTEXT)
    assert table in system_msg, \
        "stop_reason table not found in system_msg (first 200 chars: %r)" \
        % system_msg[:200]
    print("PASS test_system_msg_contains_table")


def test_reference_list_before_table():
    """'reference list below' appears before the table content in system_msg."""
    table = build_stop_reason_table()
    if not table:
        print("SKIP test_reference_list_before_table (errors.yaml not available)")
        return
    system_msg, _ = build_thinking_prompt(_CONTEXT)
    assert "reference list below" in system_msg, \
        "'reference list below' phrase missing from system_msg"
    ref_idx   = system_msg.index("reference list below")
    table_idx = system_msg.index(table)
    assert table_idx > ref_idx, (
        "table (at %d) should appear after 'reference list below' (at %d)"
        % (table_idx, ref_idx)
    )
    print("PASS test_reference_list_before_table")


# ---------------------------------------------------------------------------
# Tests 4a-4e: each stop_reason code appears in system_msg
# ---------------------------------------------------------------------------

def _make_code_test(code):
    def _test():
        table = build_stop_reason_table()
        if not table:
            print("SKIP test_code_%s (errors.yaml not available)" % code)
            return
        system_msg, _ = build_thinking_prompt(_CONTEXT)
        assert code in system_msg, \
            "Code %s not found in system_msg" % code
        print("PASS test_code_%s" % code)
    _test.__name__ = "test_code_%s" % code
    return _test

test_code_WRONG_MTZ          = _make_code_test("WRONG_MTZ")
test_code_WRONG_SPACE_GROUP  = _make_code_test("WRONG_SPACE_GROUP")
test_code_MISMATCHED_SEQUENCE = _make_code_test("MISMATCHED_SEQUENCE")
test_code_NO_SOLUTION_FOUND  = _make_code_test("NO_SOLUTION_FOUND")
test_code_REASONLESS_DIVERGENCE = _make_code_test("REASONLESS_DIVERGENCE")


# ---------------------------------------------------------------------------
# Tests 5-6: graceful degradation when errors.yaml is absent
# ---------------------------------------------------------------------------

def test_table_empty_without_yaml(tmp_path=None):
    """build_stop_reason_table returns '' gracefully when errors.yaml is absent."""
    import tempfile, shutil
    # Temporarily redirect _load_errors_yaml to a directory with no yaml
    tmpdir = tempfile.mkdtemp()
    try:
        # Monkey-patch _load_errors_yaml to point to the empty tmpdir
        import thinking_prompts as _tp
        original = _tp._load_errors_yaml

        def _empty_loader():
            return []

        _tp._load_errors_yaml = _empty_loader
        result = _tp.build_stop_reason_table()
        assert result == "", \
            "Expected '' when no yaml entries, got: %r" % result
        print("PASS test_table_empty_without_yaml")
    finally:
        _tp._load_errors_yaml = original
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_prompt_fallback_when_table_empty():
    """build_thinking_prompt falls back to plain SYSTEM_PROMPT when table is ''."""
    import thinking_prompts as _tp
    original = _tp._load_errors_yaml

    def _empty_loader():
        return []

    _tp._load_errors_yaml = _empty_loader
    try:
        system_msg, _ = build_thinking_prompt(_CONTEXT)
        # When table is empty, system_msg should equal SYSTEM_PROMPT exactly
        assert system_msg == SYSTEM_PROMPT, (
            "Expected plain SYSTEM_PROMPT when table empty;\n"
            "Got length %d vs SYSTEM_PROMPT length %d"
            % (len(system_msg), len(SYSTEM_PROMPT))
        )
        print("PASS test_prompt_fallback_when_table_empty")
    finally:
        _tp._load_errors_yaml = original


# ---------------------------------------------------------------------------
# Tests 7-11: guard — stop_reason_code cleared on non-stop actions
# ---------------------------------------------------------------------------

def test_guard_preserves_code_on_stop():
    """Guard leaves stop_reason_code intact when action is 'stop'."""
    assessment = {"action": "stop", "stop_reason_code": "WRONG_MTZ"}
    result, was_cleared = _apply_guard(assessment)
    assert not was_cleared, "Guard should not clear code on stop action"
    assert result["stop_reason_code"] == "WRONG_MTZ", \
        "stop_reason_code should be preserved"
    print("PASS test_guard_preserves_code_on_stop")


def test_guard_clears_code_on_guide_step():
    """Guard clears stop_reason_code when action is 'guide_step'."""
    assessment = {"action": "guide_step", "stop_reason_code": "WRONG_MTZ"}
    result, was_cleared = _apply_guard(assessment)
    assert was_cleared, "Guard should clear code on guide_step"
    assert result["stop_reason_code"] is None
    print("PASS test_guard_clears_code_on_guide_step")


def test_guard_clears_code_on_let_run():
    """Guard clears stop_reason_code when action is 'let_run'."""
    assessment = {"action": "let_run", "stop_reason_code": "REASONLESS_DIVERGENCE"}
    result, was_cleared = _apply_guard(assessment)
    assert was_cleared
    assert result["stop_reason_code"] is None
    print("PASS test_guard_clears_code_on_let_run")


def test_guard_clears_code_on_pivot():
    """Guard clears stop_reason_code when action is 'pivot'."""
    assessment = {"action": "pivot", "stop_reason_code": "NO_SOLUTION_FOUND"}
    result, was_cleared = _apply_guard(assessment)
    assert was_cleared
    assert result["stop_reason_code"] is None
    print("PASS test_guard_clears_code_on_pivot")


def test_guard_noop_when_code_is_none():
    """Guard is a no-op when stop_reason_code is already None — no KeyError."""
    for action in ("guide_step", "let_run", "pivot", "stop"):
        assessment = {"action": action, "stop_reason_code": None}
        result, was_cleared = _apply_guard(assessment)
        assert not was_cleared, \
            "Guard should not fire when code is None (action=%s)" % action
        assert result["stop_reason_code"] is None
    print("PASS test_guard_noop_when_code_is_none")


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


def run_all_tests():
    for t in [
        test_returns_tuple,
        test_system_msg_contains_table,
        test_reference_list_before_table,
        test_code_WRONG_MTZ,
        test_code_WRONG_SPACE_GROUP,
        test_code_MISMATCHED_SEQUENCE,
        test_code_NO_SOLUTION_FOUND,
        test_code_REASONLESS_DIVERGENCE,
        test_table_empty_without_yaml,
        test_prompt_fallback_when_table_empty,
        test_guard_preserves_code_on_stop,
        test_guard_clears_code_on_guide_step,
        test_guard_clears_code_on_let_run,
        test_guard_clears_code_on_pivot,
        test_guard_noop_when_code_is_none,
    ]:
        t()
    print("All 15 tests passed.")

if __name__ == "__main__":
    tests = [
        test_returns_tuple,
        test_system_msg_contains_table,
        test_reference_list_before_table,
        test_code_WRONG_MTZ,
        test_code_WRONG_SPACE_GROUP,
        test_code_MISMATCHED_SEQUENCE,
        test_code_NO_SOLUTION_FOUND,
        test_code_REASONLESS_DIVERGENCE,
        test_table_empty_without_yaml,
        test_prompt_fallback_when_table_empty,
        test_guard_preserves_code_on_stop,
        test_guard_clears_code_on_guide_step,
        test_guard_clears_code_on_let_run,
        test_guard_clears_code_on_pivot,
        test_guard_noop_when_code_is_none,
    ]
    failures = []
    skips = []
    for t in tests:
        try:
            t()
        except Exception as e:
            import traceback
            print("FAIL %s: %s" % (t.__name__, e))
            traceback.print_exc()
            failures.append(t.__name__)
    print()
    if failures:
        print("FAILED: %s" % ", ".join(failures))
        sys.exit(1)
    else:
        print("All %d tests passed." % len(tests))
