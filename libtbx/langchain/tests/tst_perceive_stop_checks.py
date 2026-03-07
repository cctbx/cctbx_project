"""Tests for PERCEIVE stop check helpers.

Tests check_directive_stop and check_consecutive_program_cap, which are
standalone functions in perceive_checks.py with no heavy imports.
"""
from __future__ import division, print_function
import sys
import os

# Allow running from the tests/ directory or the project root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from agent.perceive_checks import check_directive_stop, check_consecutive_program_cap


# ---------------------------------------------------------------------------
# check_directive_stop
# ---------------------------------------------------------------------------

def test_no_directives():
    """No directives → no stop."""
    stop, reason = check_directive_stop({}, [], 2)
    assert not stop
    assert reason is None

def test_no_history():
    """Directives present but no history → no stop."""
    d = {"stop_conditions": {"after_program": "phenix.refine"}}
    stop, reason = check_directive_stop(d, [], 2)
    assert not stop

def test_after_program_match():
    """after_program matches last program → NOT a hard stop.

    Since v112.78 (Bug 7), after_program is a minimum-run
    guarantee handled by the PLAN node, not a PERCEIVE
    hard stop.  The LLM decides when to actually stop.
    """
    d = {"stop_conditions": {"after_program": "phenix.refine"}}
    h = [{"program": "phenix.refine", "command": "phenix.refine foo.mtz",
          "result": "SUCCESS"}]
    stop, reason = check_directive_stop(d, h, cycle_number=2)
    assert not stop, (
        "after_program must NOT hard-stop in PERCEIVE "
        "(v112.78: minimum-run guarantee, not hard stop)")

def test_after_program_normalized():
    """after_program='refine' matches phenix.refine → still no hard stop."""
    d = {"stop_conditions": {"after_program": "refine"}}
    h = [{"program": "phenix.refine", "command": "phenix.refine foo.mtz",
          "result": "SUCCESS"}]
    stop, reason = check_directive_stop(d, h, cycle_number=2)
    assert not stop, (
        "after_program must NOT hard-stop even with "
        "normalized name match")

def test_after_program_no_match():
    """after_program doesn't match last program → no stop."""
    d = {"stop_conditions": {"after_program": "phenix.molprobity"}}
    h = [{"program": "phenix.refine", "command": "phenix.refine foo.mtz",
          "result": "SUCCESS"}]
    stop, reason = check_directive_stop(d, h, cycle_number=2)
    assert not stop

def test_after_cycle_met():
    """after_cycle=2, cycle_number=3 (meaning 2 completed) → stop."""
    d = {"stop_conditions": {"after_cycle": 2}}
    h = [{"program": "phenix.refine", "command": "c", "result": "SUCCESS"},
         {"program": "phenix.refine", "command": "c", "result": "SUCCESS"}]
    stop, reason = check_directive_stop(d, h, cycle_number=3)
    assert stop
    assert "after_cycle=2" in reason

def test_after_cycle_not_met():
    """after_cycle=3, cycle_number=2 (meaning 1 completed) → no stop."""
    d = {"stop_conditions": {"after_cycle": 3}}
    h = [{"program": "phenix.refine", "command": "c", "result": "SUCCESS"}]
    stop, reason = check_directive_stop(d, h, cycle_number=2)
    assert not stop

def test_predict_and_build_guard():
    """after_program=predict_and_build with stop_after_predict → suppress stop."""
    d = {"stop_conditions": {"after_program": "phenix.predict_and_build"}}
    h = [{"program": "phenix.predict_and_build",
          "command": "phenix.predict_and_build stop_after_predict=True model.pdb",
          "result": "SUCCESS"}]
    stop, reason = check_directive_stop(d, h, cycle_number=2)
    assert not stop, "Should suppress stop when stop_after_predict used"

def test_predict_and_build_no_guard():
    """after_program=predict_and_build without stop_after_predict → still no hard stop.

    Since v112.78, after_program never hard-stops in
    PERCEIVE regardless of stop_after_predict.
    """
    d = {"stop_conditions": {"after_program": "phenix.predict_and_build"}}
    h = [{"program": "phenix.predict_and_build",
          "command": "phenix.predict_and_build model.pdb data.mtz",
          "result": "SUCCESS"}]
    stop, reason = check_directive_stop(d, h, cycle_number=2)
    assert not stop, (
        "after_program must NOT hard-stop in PERCEIVE "
        "(v112.78: minimum-run guarantee)")

def test_r_free_target_met():
    """r_free below target → stop."""
    d = {"stop_conditions": {"r_free_target": 0.25}}
    h = [{"program": "phenix.refine", "command": "c", "result": "SUCCESS"}]
    metrics = {"r_free": 0.22}
    stop, reason = check_directive_stop(d, h, cycle_number=2,
                                         current_metrics=metrics)
    assert stop
    assert "0.220" in reason
    assert "0.250" in reason

def test_r_free_target_not_met():
    """r_free above target → no stop."""
    d = {"stop_conditions": {"r_free_target": 0.25}}
    h = [{"program": "phenix.refine", "command": "c", "result": "SUCCESS"}]
    metrics = {"r_free": 0.30}
    stop, reason = check_directive_stop(d, h, cycle_number=2,
                                         current_metrics=metrics)
    assert not stop

def test_map_cc_target_met():
    """map_cc above target → stop."""
    d = {"stop_conditions": {"map_cc_target": 0.80}}
    h = [{"program": "phenix.refine", "command": "c", "result": "SUCCESS"}]
    metrics = {"map_cc": 0.85}
    stop, reason = check_directive_stop(d, h, cycle_number=2,
                                         current_metrics=metrics)
    assert stop
    assert "Map CC" in reason

def test_after_cycle_takes_priority():
    """after_cycle fires before after_program check."""
    d = {"stop_conditions": {"after_cycle": 1, "after_program": "phenix.molprobity"}}
    h = [{"program": "phenix.refine", "command": "c", "result": "SUCCESS"}]
    stop, reason = check_directive_stop(d, h, cycle_number=2)
    assert stop
    assert "after_cycle" in reason


# ---------------------------------------------------------------------------
# check_consecutive_program_cap
# ---------------------------------------------------------------------------

def test_no_history_cap():
    """No history → no stop."""
    stop, reason = check_consecutive_program_cap([])
    assert not stop

def test_below_cap():
    """2 consecutive → below default threshold of 3 → no stop."""
    h = [
        {"program": "phenix.refine", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"},
    ]
    stop, reason = check_consecutive_program_cap(h)
    assert not stop

def test_at_cap():
    """3 consecutive → at threshold → stop."""
    h = [
        {"program": "phenix.refine", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"},
    ]
    stop, reason = check_consecutive_program_cap(h)
    assert stop
    assert "phenix.refine" in reason
    assert "3" in reason

def test_above_cap():
    """5 consecutive → above threshold → stop."""
    h = [{"program": "phenix.refine", "result": "SUCCESS"}] * 5
    stop, reason = check_consecutive_program_cap(h)
    assert stop

def test_mixed_programs():
    """Different programs → no stop."""
    h = [
        {"program": "phenix.refine", "result": "SUCCESS"},
        {"program": "phenix.molprobity", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"},
    ]
    stop, reason = check_consecutive_program_cap(h)
    assert not stop

def test_failed_runs_ignored():
    """Only SUCCESS runs count — failures don't break the streak."""
    h = [
        {"program": "phenix.refine", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "FAILED: something"},
        {"program": "phenix.refine", "result": "SUCCESS"},
    ]
    stop, reason = check_consecutive_program_cap(h)
    assert stop, "3 SUCCESS runs of refine, failure is skipped"

def test_different_before_streak():
    """Different program before 3 consecutive → stop (only tail matters)."""
    h = [
        {"program": "phenix.molprobity", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"},
        {"program": "phenix.refine", "result": "SUCCESS"},
    ]
    stop, reason = check_consecutive_program_cap(h)
    assert stop

def test_custom_threshold():
    """Custom max_consecutive=5 → 3 doesn't trigger."""
    h = [{"program": "phenix.refine", "result": "SUCCESS"}] * 3
    stop, reason = check_consecutive_program_cap(h, max_consecutive=5)
    assert not stop

    h5 = [{"program": "phenix.refine", "result": "SUCCESS"}] * 5
    stop, reason = check_consecutive_program_cap(h5, max_consecutive=5)
    assert stop


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

def run():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    failed = 0
    first_failure = None
    for t in tests:
        try:
            t()
            passed += 1
            print("  PASS: %s" % t.__name__)
        except Exception as e:
            failed += 1
            print("  FAIL: %s — %s" % (t.__name__, e))
            if first_failure is None:
                first_failure = "%s — %s" % (t.__name__, e)
    print("\n%d passed, %d failed" % (passed, failed))
    if failed > 0:
        raise AssertionError(
            "First failure: %s" % first_failure)


if __name__ == "__main__":
    try:
        run()
    except AssertionError:
        sys.exit(1)
