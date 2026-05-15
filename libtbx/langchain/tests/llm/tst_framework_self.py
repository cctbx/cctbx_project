"""
Deterministic tests for the LLM-test framework's early-termination
logic.  These tests use a mock run_one_fn and don't make any LLM
calls; they verify the framework's verdict mathematics.

This file IS appropriate to run as part of the deterministic test
suite (unlike the rest of tests/llm/ which need API keys).  However,
to keep the convention clean we still ship it under tests/llm/ — the
user can opt in to running it via the run_llm_tests.py runner with
the special suite name "framework_self_test".

Run standalone:
    phenix.python tests/llm/tst_framework_self.py
"""

from __future__ import absolute_import, division, print_function

import os
import sys

# Set a fake env var so the framework's provider detection returns
# something for the test setup.  The mock run_one_fn never actually
# touches the LLM.
os.environ.setdefault("GOOGLE_API_KEY", "test-key-for-framework-tests")

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from framework import (Scenario, RunOutcome, run_with_early_termination,
                       _required_passes)


def _mock_runner(verdicts):
    """Build a mock run_one_fn that returns canned verdicts in order.

    verdicts: list of items, each is True (pass), False (semantic fail),
              or an Exception instance (error).
    """
    state = {"i": 0}

    def run_one(scenario, provider, run_index):
        item = verdicts[state["i"]]
        state["i"] += 1
        if isinstance(item, Exception):
            return RunOutcome(
                run_index=run_index, elapsed_s=0.0,
                passed=False, why=str(item),
                raw_output="", parsed=None,
                error="MockError: %s" % item)
        return RunOutcome(
            run_index=run_index, elapsed_s=0.0,
            passed=bool(item), why="mock",
            raw_output="", parsed={"mock": item}, error=None)
    return run_one


# =====================================================================
# Capability tests
# =====================================================================

def test_capability_first_pass_early_stops():
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="capability", max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([True, True, True, True, True]))
    assert v.result == "PASS"
    assert v.n_runs_executed == 1
    print("PASS test_capability_first_pass_early_stops")


def test_capability_all_fail_returns_fail():
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="capability", max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google", _mock_runner([False]*5))
    assert v.result == "FAIL"
    assert v.n_runs_executed == 5
    print("PASS test_capability_all_fail_returns_fail")


def test_capability_late_pass_works():
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="capability", max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([False, False, False, False, True]))
    assert v.result == "PASS"
    assert v.n_runs_executed == 5
    print("PASS test_capability_late_pass_works")


def test_capability_all_error_returns_error():
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="capability", max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([Exception("api")]*5))
    assert v.result == "ERROR"
    assert v.n_errored == 5
    print("PASS test_capability_all_error_returns_error")


def test_capability_errors_then_pass():
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="capability", max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([Exception("e"),
                                                 Exception("e"),
                                                 True, True, True]))
    assert v.result == "PASS"
    assert v.n_passed == 1
    assert v.n_errored == 2
    print("PASS test_capability_errors_then_pass")


# =====================================================================
# Reliability tests
# =====================================================================

def test_reliability_required_passes_math():
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.8, max_runs=5,
                  input="x", expected_fn=None)
    assert _required_passes(sc) == 4

    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.9, max_runs=5,
                  input="x", expected_fn=None)
    assert _required_passes(sc) == 5

    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.6, max_runs=10,
                  input="x", expected_fn=None)
    assert _required_passes(sc) == 6
    print("PASS test_reliability_required_passes_math")


def test_reliability_08_early_pass_at_4():
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.8, max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google", _mock_runner([True]*5))
    assert v.result == "PASS"
    assert v.n_runs_executed == 4
    print("PASS test_reliability_08_early_pass_at_4")


def test_reliability_08_early_fail_at_2_fails():
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.8, max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([False, False, True, True, True]))
    assert v.result == "FAIL"
    assert v.n_runs_executed == 2
    print("PASS test_reliability_08_early_fail_at_2_fails")


def test_reliability_08_1fail_4pass_succeeds():
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.8, max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([False, True, True, True, True]))
    assert v.result == "PASS"
    assert v.n_passed == 4
    print("PASS test_reliability_08_1fail_4pass_succeeds")


def test_reliability_09_one_fail_means_fail():
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.9, max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([False, True, True, True, True]))
    assert v.result == "FAIL"
    assert v.n_runs_executed == 1
    print("PASS test_reliability_09_one_fail_means_fail")


def test_reliability_09_all_pass():
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.9, max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google", _mock_runner([True]*5))
    assert v.result == "PASS"
    assert v.n_runs_executed == 5
    print("PASS test_reliability_09_all_pass")


# =====================================================================
# Error-handling fix (P1) — errors should NOT count toward early-fail
# =====================================================================

def test_reliability_09_error_then_passes_returns_error():
    """If run 1 errors then runs 2-5 all pass: 4 pass + 1 error +
    0 semantic fail.  Per the corrected verdict logic, this is ERROR
    not FAIL: no run was ever observed making a wrong decision; we
    just couldn't gather the 5th data point due to a transient error.

    This is Tom's "3/5 with [2 err]" production case generalized."""
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.9, max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([Exception("api"),
                                                 True, True, True, True]))
    assert v.result == "ERROR"
    assert v.n_runs_executed == 5  # critical: did NOT bail at run 1
    assert v.n_passed == 4
    assert v.n_errored == 1
    assert v.n_failed == 0
    print("PASS test_reliability_09_error_then_passes_returns_error")


def test_reliability_09_error_then_fail_early_stops():
    """Run 1 errors, run 2 semantically fails.  Best-case after run 2:
    0 pass + 3 remaining + 1 err = 4 < 5 → early-fail at run 2."""
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.9, max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([Exception("api"),
                                                 False, False, False, False]))
    assert v.result == "FAIL"
    assert v.n_runs_executed == 2
    print("PASS test_reliability_09_error_then_fail_early_stops")


def test_reliability_08_two_errors_three_passes_returns_error():
    """0.8 needs 4/5.  3 pass + 2 errors + 0 semantic fail.  Per the
    corrected verdict logic, this is ERROR not FAIL: we never
    observed the LLM making a wrong choice — the 2 missing data
    points were transient API issues, not LLM misbehavior.

    This is EXACTLY the case Tom hit in his 'dont_retry_same_failure
    3/5 [2 err]' production run.  Old logic would FAIL this; new
    logic correctly returns ERROR."""
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.8, max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([Exception("e"),
                                                 Exception("e"),
                                                 True, True, True]))
    assert v.result == "ERROR"
    assert v.n_runs_executed == 5
    assert v.n_passed == 3
    assert v.n_errored == 2
    assert v.n_failed == 0
    print("PASS test_reliability_08_two_errors_three_passes_returns_error")


def test_reliability_08_fail_with_errors_still_fail_if_unreachable():
    """0.8 needs 4/5.  1 pass + 2 fail + 2 err.  Even if both errors
    had been passes, we'd have 3 passes — still under 4 required.
    Combined with 2 confirmed semantic failures, this is confirmed
    unreliable: FAIL."""
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.8, max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([True, False, False,
                                                 Exception("e"),
                                                 Exception("e")]))
    # Note: this scenario early-fails at run 3 (1 pass + 2 fail + 2 remaining
    # = best-case 3 < 4 required, treating errors as could-have-passed).
    # The 2 errors don't actually get to execute.
    assert v.result == "FAIL"
    assert v.n_failed >= 2
    print("PASS test_reliability_08_fail_with_errors_still_fail_if_unreachable")


def test_reliability_08_one_fail_some_errors_returns_error():
    """0.8 needs 4/5.  2 pass + 1 fail + 2 err.  Best case if errors
    were passes: 4 passes — could have reached threshold.  We have one
    confirmed failure but can't say with confidence that the LLM is
    unreliable.  ERROR."""
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.8, max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([True, True,
                                                 Exception("e"),
                                                 False,
                                                 Exception("e")]))
    assert v.result == "ERROR"
    assert v.n_passed == 2
    assert v.n_failed == 1
    assert v.n_errored == 2
    print("PASS test_reliability_08_one_fail_some_errors_returns_error")


def test_reliability_all_errors_returns_error_not_fail():
    """All 5 runs error → verdict is ERROR (not FAIL), because we have
    no semantic information to distinguish.  Also: does not early-stop
    on first error (would only stop if best_case < required, which
    would require remaining + errors < required, which is impossible
    when we treat errors as indeterminate)."""
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="reliability", threshold=0.9, max_runs=5,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google",
                                   _mock_runner([Exception("e")]*5))
    assert v.result == "ERROR"
    assert v.n_runs_executed == 5
    assert v.n_errored == 5
    print("PASS test_reliability_all_errors_returns_error_not_fail")


# =====================================================================
# Edge cases
# =====================================================================

def test_max_runs_zero_degenerate():
    """max_runs=0 is degenerate but must not crash."""
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="capability", max_runs=0,
                  input="x", expected_fn=None)
    v = run_with_early_termination(sc, "google", _mock_runner([]))
    assert v.result == "FAIL"
    assert v.n_runs_executed == 0
    print("PASS test_max_runs_zero_degenerate")


def test_unknown_test_type_raises():
    """An unrecognized test_type should fail loudly."""
    sc = Scenario(name="t", description="d", decision_point="dp",
                  test_type="not_a_real_type", max_runs=5,
                  input="x", expected_fn=None)
    try:
        run_with_early_termination(sc, "google", _mock_runner([True]))
    except ValueError:
        print("PASS test_unknown_test_type_raises")
        return
    raise AssertionError("expected ValueError")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    """Run all framework self-tests; raise on first failure."""
    test_fns = sorted(
        ((k, v) for k, v in globals().items()
         if k.startswith("test_") and callable(v)),
        key=lambda kv: kv[0])
    passed = 0
    for name, fn in test_fns:
        try:
            fn()
            passed += 1
        except Exception as e:
            print("FAIL %s: %s" % (name, e))
            raise
    print()
    print("Framework self-test: %d / %d passed" % (passed, len(test_fns)))


if __name__ == "__main__":
    run_all_tests()
