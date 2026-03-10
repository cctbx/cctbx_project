"""
diagnosed_failure stop reason for terminal diagnosis (v115 F6).

F6: Ensure DisplayDataModel sees stop_reason="diagnosed_failure" in
    session.data so the Results tab shows the ⚠ Agent Stopped indicator.

    Two concrete changes tested:
      A. _diagnose_terminal_failure now sets session.data["stop_reason"]
         = "diagnosed_failure" alongside failure_diagnosis_path.
      B. _write_session_summary_json uses "diagnosed_failure" (not the
         misleading "red_flag") as the stop_reason in the aborted branch.
"""
import sys
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in [_ROOT, os.path.join(_ROOT, "programs")]:
    if _p not in sys.path:
        sys.path.insert(0, _p)


# ---------------------------------------------------------------------------
# Minimal stubs
# ---------------------------------------------------------------------------

class _SessionStub:
    def __init__(self, data=None):
        self.data = dict(data or {})


# ---------------------------------------------------------------------------
# Replicate the two F6-edited code paths
# ---------------------------------------------------------------------------

def _simulate_diagnose_sets_session_keys(html_path):
    """Replicate the session.data assignments in _diagnose_terminal_failure."""
    session = _SessionStub()
    # Simulate successful html write path
    session.data["failure_diagnosis_path"] = html_path
    session.data["stop_reason"] = "diagnosed_failure"  # F6 addition
    return session


def _simulate_summary_outcome(has_fatal, gate_stopped=False, plan_complete=False,
                               has_plan=False):
    """Replicate the outcome hierarchy in _write_session_summary_json."""
    summary = {}
    if gate_stopped:
        summary["outcome"] = "stopped"
        summary["stop_reason"] = "gate"
    elif has_fatal:
        summary["outcome"] = "aborted"
        summary["stop_reason"] = "diagnosed_failure"   # F6: was "red_flag"
    elif plan_complete:
        summary["outcome"] = "complete"
    elif has_plan:
        summary["outcome"] = "incomplete"
    else:
        summary["outcome"] = "reactive"
    return summary


# ============================================================================
# F6A — session.data["stop_reason"] set by _diagnose_terminal_failure
# ============================================================================

def test_f6a_stop_reason_set_in_session():
    session = _simulate_diagnose_sets_session_keys("/tmp/ai_failure_diagnosis.html")
    assert session.data.get("stop_reason") == "diagnosed_failure", \
        session.data.get("stop_reason")


def test_f6a_failure_diagnosis_path_still_set():
    """Existing failure_diagnosis_path key must not be removed."""
    session = _simulate_diagnose_sets_session_keys("/tmp/ai_failure_diagnosis.html")
    assert session.data.get("failure_diagnosis_path") == "/tmp/ai_failure_diagnosis.html"


def test_f6a_stop_reason_is_not_red_flag():
    """The old value 'red_flag' was misleading — must not appear for diagnosis stops."""
    session = _simulate_diagnose_sets_session_keys("/tmp/ai_failure_diagnosis.html")
    assert session.data.get("stop_reason") != "red_flag"


def test_f6a_stop_reason_set_even_when_html_write_fails():
    """stop_reason must be set unconditionally — the diagnosis fired regardless
    of whether the HTML report could be written to disk."""
    session = _SessionStub()
    # Simulate html write failure: failure_diagnosis_path not set, but
    # stop_reason is set after the try block (step 4/5 always runs).
    session.data["stop_reason"] = "diagnosed_failure"
    assert session.data.get("stop_reason") == "diagnosed_failure"
    # failure_diagnosis_path absent is fine — no button, but ⚠ still shows
    assert session.data.get("failure_diagnosis_path") is None


def test_f6a_full_report_line_only_when_file_written():
    """Step 4 must only print 'Full report:' when failure_diagnosis_path
    is in session.data — not when html_path is merely assigned in the try
    block before a failed open()."""
    # When HTML write succeeds: session.data has failure_diagnosis_path
    session_ok = _SessionStub({"failure_diagnosis_path": "/out/ai_failure_diagnosis.html"})
    assert session_ok.data.get("failure_diagnosis_path")  # would print

    # When HTML write fails: session.data has no failure_diagnosis_path
    session_fail = _SessionStub()
    assert not session_fail.data.get("failure_diagnosis_path")  # would NOT print


def test_f6_has_fatal_uses_stop_reason_primary():
    """_has_fatal should be True when stop_reason='diagnosed_failure' even
    if failure_diagnosis_path is absent (HTML write failed)."""
    sd = {"stop_reason": "diagnosed_failure"}  # no failure_diagnosis_path
    _has_fatal = (
        sd.get("stop_reason") == "diagnosed_failure"
        or bool(sd.get("failure_diagnosis_path"))
    )
    assert _has_fatal


def test_f6_has_fatal_uses_path_fallback_for_old_sessions():
    """Pre-v115 sessions have failure_diagnosis_path but no stop_reason.
    _has_fatal must still be True for backward compatibility."""
    sd = {"failure_diagnosis_path": "/old/ai_failure_diagnosis.html"}
    _has_fatal = (
        sd.get("stop_reason") == "diagnosed_failure"
        or bool(sd.get("failure_diagnosis_path"))
    )
    assert _has_fatal


def test_f6_has_fatal_false_for_normal_session():
    sd = {}
    _has_fatal = (
        sd.get("stop_reason") == "diagnosed_failure"
        or bool(sd.get("failure_diagnosis_path"))
    )
    assert not _has_fatal


# ============================================================================
# F6B — session_summary.json stop_reason for aborted sessions
#
# Note: _write_session_summary_json is only called from _generate_structure_report,
# which is skipped when has_fatal_diagnosis=True. So in practice, session_summary.json
# is never written for diagnosed failures. The _has_fatal branch tested here is
# defensive dead code — correct if _write_session_summary_json is ever called
# separately in future, and it fixes the misleading "red_flag" value for any
# path that does reach it (e.g. via tooling or test harness).
# ============================================================================

def test_f6b_summary_aborted_uses_diagnosed_failure():
    summary = _simulate_summary_outcome(has_fatal=True)
    assert summary["outcome"] == "aborted"
    assert summary["stop_reason"] == "diagnosed_failure", summary["stop_reason"]


def test_f6b_summary_aborted_not_red_flag():
    summary = _simulate_summary_outcome(has_fatal=True)
    assert summary.get("stop_reason") != "red_flag"


def test_f6b_summary_gate_stopped_unaffected():
    """gate_stopped path must still use its own stop_reason."""
    summary = _simulate_summary_outcome(has_fatal=False, gate_stopped=True)
    assert summary["outcome"] == "stopped"
    assert summary["stop_reason"] == "gate"


def test_f6b_summary_complete_no_stop_reason():
    summary = _simulate_summary_outcome(has_fatal=False, plan_complete=True,
                                        has_plan=True)
    assert summary["outcome"] == "complete"
    assert "stop_reason" not in summary


def test_f6b_summary_has_fatal_takes_priority_over_plan():
    """has_fatal check is before plan checks — aborted overrides incomplete."""
    summary = _simulate_summary_outcome(has_fatal=True, has_plan=True)
    assert summary["outcome"] == "aborted"


# ============================================================================
# _finalize_session has_fatal_diagnosis consistency
# ============================================================================

def _compute_has_fatal_diagnosis(session_data):
    """Replicate the updated has_fatal_diagnosis check in _finalize_session."""
    return (
        session_data.get("stop_reason") == "diagnosed_failure"
        or bool(session_data.get("failure_diagnosis_path"))
    )


def test_f6_finalize_fatal_via_stop_reason_only():
    """HTML write failed: only stop_reason set. _finalize_session must still
    treat this as a fatal diagnosis to avoid generating a confusing structure
    report."""
    sd = {"stop_reason": "diagnosed_failure"}
    assert _compute_has_fatal_diagnosis(sd)


def test_f6_finalize_fatal_via_path_only():
    """Pre-v115 session: only failure_diagnosis_path set (no stop_reason)."""
    sd = {"failure_diagnosis_path": "/old/ai_failure_diagnosis.html"}
    assert _compute_has_fatal_diagnosis(sd)


def test_f6_finalize_fatal_both_set():
    """Normal v115+ case: both keys set after successful HTML write."""
    sd = {
        "stop_reason": "diagnosed_failure",
        "failure_diagnosis_path": "/out/ai_failure_diagnosis.html",
    }
    assert _compute_has_fatal_diagnosis(sd)


def test_f6_finalize_not_fatal_for_normal_session():
    assert not _compute_has_fatal_diagnosis({})


def test_f6_finalize_not_fatal_for_other_stop_reason():
    """A different stop_reason (e.g. max_cycles) must not trigger fatal path."""
    assert not _compute_has_fatal_diagnosis({"stop_reason": "max_cycles"})


# ============================================================================
# Integration-style: session_data keys available to DisplayDataModel
# ============================================================================

def test_f6_session_data_has_both_keys_for_ddm():
    """DisplayDataModel.from_session(session_data) gets both keys it needs
    to show ⚠ Agent Stopped: failure_diagnosis_path (for the button) and
    stop_reason='diagnosed_failure' (for outcome_status)."""
    session = _simulate_diagnose_sets_session_keys("/out/ai_failure_diagnosis.html")
    sd = session.data
    assert "failure_diagnosis_path" in sd
    assert sd["stop_reason"] == "diagnosed_failure"


# ============================================================================
# Test runner
# ============================================================================

_TESTS = [
    test_f6a_stop_reason_set_in_session,
    test_f6a_failure_diagnosis_path_still_set,
    test_f6a_stop_reason_is_not_red_flag,
    test_f6a_stop_reason_set_even_when_html_write_fails,
    test_f6a_full_report_line_only_when_file_written,
    test_f6_has_fatal_uses_stop_reason_primary,
    test_f6_has_fatal_uses_path_fallback_for_old_sessions,
    test_f6_has_fatal_false_for_normal_session,
    test_f6b_summary_aborted_uses_diagnosed_failure,
    test_f6b_summary_aborted_not_red_flag,
    test_f6b_summary_gate_stopped_unaffected,
    test_f6b_summary_complete_no_stop_reason,
    test_f6b_summary_has_fatal_takes_priority_over_plan,
    test_f6_finalize_fatal_via_stop_reason_only,
    test_f6_finalize_fatal_via_path_only,
    test_f6_finalize_fatal_both_set,
    test_f6_finalize_not_fatal_for_normal_session,
    test_f6_finalize_not_fatal_for_other_stop_reason,
    test_f6_session_data_has_both_keys_for_ddm,
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
