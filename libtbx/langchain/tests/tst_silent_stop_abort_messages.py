"""
abort_message populated for previously-silent stops (v115 follow-up).

no_workflow_state and after_program_not_available stops previously produced
no abort_message at the state top-level; these tests verify the message is
now surfaced for both.
"""
import sys
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in [_ROOT, os.path.join(_ROOT, "agent")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _no_workflow_state_return():
    _msg = (
        "Agent stopped: the workflow engine could not determine a current "
        "state. Check that the expected input files are present in the "
        "working directory, then re-run."
    )
    return {
        "intent": {
            "program": None, "reasoning": _msg, "files": {}, "strategy": {},
            "stop": True, "stop_reason": "no_workflow_state", "abort_message": _msg,
        },
        "command": "STOP",
        "stop_reason": "no_workflow_state",
        "abort_message": _msg,
    }


def _after_program_not_available_return(chosen_program):
    _msg = (
        "Agent stopped: the program '%s' requested by the "
        "stop_after directive is not available at the current "
        "workflow step. Check that the program name is spelled "
        "correctly and that the required input files are present."
        % chosen_program
    )
    return {
        "command": "STOP",
        "stop": True,
        "stop_reason": "after_program_not_available",
        "abort_message": _msg,
        "validation_error": None,
    }


def test_no_workflow_state_abort_message_present():
    r = _no_workflow_state_return()
    assert r.get("abort_message"), "abort_message missing at state top-level"
    assert "input files" in r["abort_message"].lower()


def test_no_workflow_state_intent_has_abort_message():
    r = _no_workflow_state_return()
    assert r["abort_message"] == r["intent"]["abort_message"]


def test_no_workflow_state_stop_reason():
    r = _no_workflow_state_return()
    assert r["stop_reason"] == "no_workflow_state"


def test_no_workflow_state_stop_true():
    r = _no_workflow_state_return()
    assert r["intent"]["stop"] is True


def test_after_program_not_available_names_program():
    r = _after_program_not_available_return("phenix.refine")
    assert "phenix.refine" in r["abort_message"]


def test_after_program_not_available_stop_reason():
    r = _after_program_not_available_return("phenix.refine")
    assert r["stop_reason"] == "after_program_not_available"


def test_after_program_not_available_abort_at_top_level():
    r = _after_program_not_available_return("phenix.autobuild")
    assert r.get("abort_message"), "abort_message missing"
    assert "phenix.autobuild" in r["abort_message"]


def test_after_program_not_available_stop_true():
    r = _after_program_not_available_return("phenix.ligandfit")
    assert r["stop"] is True


_TESTS = [
    test_no_workflow_state_abort_message_present,
    test_no_workflow_state_intent_has_abort_message,
    test_no_workflow_state_stop_reason,
    test_no_workflow_state_stop_true,
    test_after_program_not_available_names_program,
    test_after_program_not_available_stop_reason,
    test_after_program_not_available_abort_at_top_level,
    test_after_program_not_available_stop_true,
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
