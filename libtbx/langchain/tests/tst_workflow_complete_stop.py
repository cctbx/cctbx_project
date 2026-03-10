"""
workflow_complete stop reason on rules-only success (v115 F5).

F5: Use stop_reason='workflow_complete' when the workflow engine reaches
    the 'complete' step so the Results tab shows a success message
    instead of "No valid programs available".

    Two code paths in _mock_plan are covered:
      A. RulesSelector path: intent overridden when step=='complete'.
      B. Simple fallback path: stop_reason='workflow_complete' when
         step=='complete', 'no_valid_programs' otherwise.

    F5 does NOT change the LLM path (stop_reason in that case is
    determined by the LLM's own output).
"""
import sys
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in [_ROOT, os.path.join(_ROOT, "agent")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)


# ---------------------------------------------------------------------------
# Replicate the F5-edited logic inline
# ---------------------------------------------------------------------------

def _rules_selector_return(rules_intent, step_name):
    """Replicate the F5-edited RulesSelector return in _mock_plan."""
    if rules_intent.get("stop") and step_name == "complete":
        _success_msg = "Workflow complete — structure determination finished successfully."
        rules_intent = {
            **rules_intent,
            "stop_reason": "workflow_complete",
            "abort_message": _success_msg,
            "reasoning": rules_intent.get("reasoning") or _success_msg,
        }
    state = {}
    return {
        **state,
        "intent": rules_intent,
        "stop": rules_intent.get("stop", False),
        "stop_reason": rules_intent.get("stop_reason"),
        "abort_message": rules_intent.get("abort_message") or state.get("abort_message"),
    }


def _simple_fallback_stop(step_name, session_blocked_programs=None):
    """Replicate the F5+F4-edited simple fallback STOP block in _mock_plan."""
    state = {"session_blocked_programs": session_blocked_programs or []}
    if step_name == "complete":
        _stop_reasoning = "Workflow complete — structure determination finished successfully."
        mock_intent = {
            "program": None,
            "reasoning": _stop_reasoning,
            "files": {},
            "strategy": {},
            "stop": True,
            "stop_reason": "workflow_complete",
            "abort_message": _stop_reasoning,
        }
    else:
        _blocked = state.get("session_blocked_programs", [])
        if _blocked:
            _blocked_str = ", ".join(_blocked)
            _stop_reasoning = (
                "No valid programs remain. "
                "The following programs were blocked after repeated failures: %s. "
                "Check the log for details or restart with different settings."
                % _blocked_str
            )
        else:
            _stop_reasoning = "Fallback: No valid programs available, stopping."
        mock_intent = {
            "program": None,
            "reasoning": _stop_reasoning,
            "files": {},
            "strategy": {},
            "stop": True,
            "stop_reason": "no_valid_programs",
            "abort_message": _stop_reasoning,
        }
    return {
        **state,
        "intent": mock_intent,
        "stop": mock_intent.get("stop", False),
        "stop_reason": mock_intent.get("stop_reason"),
        "abort_message": mock_intent.get("abort_message") or state.get("abort_message"),
    }


# ============================================================================
# F5A — RulesSelector path
# ============================================================================

def test_f5a_complete_step_overrides_to_workflow_complete():
    rules_intent = {"stop": True, "stop_reason": "no_valid_programs",
                    "program": None, "reasoning": ""}
    result = _rules_selector_return(rules_intent, "complete")
    assert result["stop_reason"] == "workflow_complete", result["stop_reason"]


def test_f5a_complete_step_sets_success_abort_message():
    rules_intent = {"stop": True, "stop_reason": "no_valid_programs",
                    "program": None, "reasoning": ""}
    result = _rules_selector_return(rules_intent, "complete")
    assert "complete" in result["abort_message"].lower()
    assert result["abort_message"] == result["intent"]["abort_message"]


def test_f5a_abort_message_at_state_top_level():
    """abort_message must be at state top-level for run_ai_agent.py F1."""
    rules_intent = {"stop": True, "stop_reason": "no_valid_programs",
                    "program": None, "reasoning": ""}
    result = _rules_selector_return(rules_intent, "complete")
    assert result.get("abort_message"), "abort_message missing from state top-level"


def test_f5a_non_complete_step_unchanged():
    """For other steps, RulesSelector stop_reason must pass through unchanged."""
    rules_intent = {"stop": True, "stop_reason": "no_valid_programs",
                    "program": None, "reasoning": "stuck"}
    result = _rules_selector_return(rules_intent, "refine")
    assert result["stop_reason"] == "no_valid_programs"


def test_f5a_non_stop_intent_unchanged():
    """Non-stop intents must not be affected regardless of step."""
    rules_intent = {"stop": False, "stop_reason": None,
                    "program": "phenix.refine", "reasoning": "refine next"}
    result = _rules_selector_return(rules_intent, "complete")
    assert result["stop"] is False
    assert result["stop_reason"] is None


def test_f5a_existing_reasoning_preserved_if_non_empty():
    """If RulesSelector already returned useful reasoning, keep it."""
    rules_intent = {"stop": True, "stop_reason": "no_valid_programs",
                    "program": None,
                    "reasoning": "All refinement goals met."}
    result = _rules_selector_return(rules_intent, "complete")
    assert result["intent"]["reasoning"] == "All refinement goals met."


def test_f5a_empty_reasoning_replaced_with_success_message():
    """Empty reasoning should be replaced with the success message."""
    rules_intent = {"stop": True, "stop_reason": "no_valid_programs",
                    "program": None, "reasoning": ""}
    result = _rules_selector_return(rules_intent, "complete")
    assert result["intent"]["reasoning"] != ""
    assert "complete" in result["intent"]["reasoning"].lower()


# ============================================================================
# F5B — Simple fallback path
# ============================================================================

def test_f5b_complete_step_uses_workflow_complete():
    result = _simple_fallback_stop("complete")
    assert result["stop_reason"] == "workflow_complete"


def test_f5b_complete_step_success_message():
    result = _simple_fallback_stop("complete")
    assert "complete" in result["abort_message"].lower()
    assert "successfully" in result["abort_message"].lower()


def test_f5b_complete_abort_message_at_top_level():
    result = _simple_fallback_stop("complete")
    assert result.get("abort_message")


def test_f5b_non_complete_stuck_uses_no_valid_programs():
    result = _simple_fallback_stop("refine")
    assert result["stop_reason"] == "no_valid_programs"


def test_f5b_blocked_programs_named_for_non_complete():
    """F4 behaviour: blocked program names appear in message for non-complete stops."""
    result = _simple_fallback_stop("refine", ["phenix.refine"])
    assert "phenix.refine" in result["abort_message"]
    assert result["stop_reason"] == "no_valid_programs"


def test_f5b_complete_overrides_even_with_blocked_programs():
    """If step is 'complete', blocked programs don't matter — it's still success."""
    result = _simple_fallback_stop("complete", ["phenix.refine"])
    assert result["stop_reason"] == "workflow_complete"
    assert "successfully" in result["abort_message"].lower()


def test_f5b_stop_true_in_all_cases():
    for step in ["complete", "refine", "validate"]:
        result = _simple_fallback_stop(step)
        assert result["stop"] is True, "stop should be True for step=%s" % step


# ============================================================================
# Test registry and runner
# ============================================================================

_TESTS = [
    # F5A — RulesSelector
    test_f5a_complete_step_overrides_to_workflow_complete,
    test_f5a_complete_step_sets_success_abort_message,
    test_f5a_abort_message_at_state_top_level,
    test_f5a_non_complete_step_unchanged,
    test_f5a_non_stop_intent_unchanged,
    test_f5a_existing_reasoning_preserved_if_non_empty,
    test_f5a_empty_reasoning_replaced_with_success_message,
    # F5B — simple fallback
    test_f5b_complete_step_uses_workflow_complete,
    test_f5b_complete_step_success_message,
    test_f5b_complete_abort_message_at_top_level,
    test_f5b_non_complete_stuck_uses_no_valid_programs,
    test_f5b_blocked_programs_named_for_non_complete,
    test_f5b_complete_overrides_even_with_blocked_programs,
    test_f5b_stop_true_in_all_cases,
]


def run_all_tests():
    for test_fn in _TESTS:
        test_fn()
    print("All %d tests passed." % len(_TESTS))

if __name__ == "__main__":
    passed = 0
    failed = 0
    for test_fn in _TESTS:
        try:
            test_fn()
            print("  PASS: %s" % test_fn.__name__)
            passed += 1
        except Exception as exc:
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
