"""
Tests for v116.12 Fix #1: Cryo-EM stop logic checks validation_done.

The X-ray path in analyze_metrics_trend correctly defers auto-stop until
validation (phenix.molprobity / phenix.model_vs_data) has been run.  Pre-
v116.12, the cryo-EM path immediately set should_stop=True when CC > 0.70
regardless of whether validation had been run.  This caused AF_7mjs to
auto-stop after refinement (cycle 4, CC=0.78) and skip the planned
validation stage (Stage 5, molprobity).

These tests verify:
  - When CC > 0.70 and validation NOT done: should_stop=False, suggest_validation=True
  - When CC > 0.70 and phenix.molprobity in history: should_stop=True
  - When CC > 0.70 and phenix.validation_cryoem in history: should_stop=True
  - When CC <= 0.70: behaviour unchanged (no validation check needed yet)
  - X-ray path is not affected by this change

Run with:
  python tst_cryoem_stop_validation.py
  or
  phenix.python tst_cryoem_stop_validation.py
"""
import sys
import os

# Allow running standalone without libtbx.langchain on the path.
HERE = os.path.dirname(os.path.abspath(__file__))
_ROOTS = (
    os.path.abspath(os.path.join(HERE, "..")),
    os.path.abspath(HERE),
)
for _root in _ROOTS:
    if _root not in sys.path:
        sys.path.insert(0, _root)

# Stub libtbx so the import chain inside metrics_analyzer works without PHENIX.
if "libtbx" not in sys.modules:
    import types as _types
    _libtbx = _types.ModuleType("libtbx")
    _langchain = _types.ModuleType("libtbx.langchain")
    _agent = _types.ModuleType("libtbx.langchain.agent")
    sys.modules["libtbx"] = _libtbx
    sys.modules["libtbx.langchain"] = _langchain
    sys.modules["libtbx.langchain.agent"] = _agent

try:
    from libtbx.langchain.agent.metrics_analyzer import analyze_metrics_trend
except ImportError:
    from agent.metrics_analyzer import analyze_metrics_trend


# =============================================================================
# Test infrastructure
# =============================================================================

PASS = 0
FAIL = 0
FAILURES = []


def assert_eq(actual, expected, msg):
    global PASS, FAIL
    if actual == expected:
        PASS += 1
    else:
        FAIL += 1
        FAILURES.append("%s: expected %r, got %r" % (msg, expected, actual))


def assert_true(condition, msg):
    global PASS, FAIL
    if condition:
        PASS += 1
    else:
        FAIL += 1
        FAILURES.append("%s: condition was False" % msg)


def assert_false(condition, msg):
    global PASS, FAIL
    if not condition:
        PASS += 1
    else:
        FAIL += 1
        FAILURES.append("%s: condition was True" % msg)


def _rsr(cc, program="phenix.real_space_refine"):
    """Build a real_space_refine metrics history entry."""
    return {"program": program, "map_cc": cc}


def _molprobity():
    return {"program": "phenix.molprobity", "clashscore": 5.0}


def _validation_cryoem():
    return {"program": "phenix.validation_cryoem", "map_cc": 0.78}


# =============================================================================
# Tests
# =============================================================================

def test_primary_af7mjs_regression():
    """The exact AF_7mjs scenario: refinement reached CC > 0.70,
    validation NOT done. Should NOT auto-stop.

    NOTE: The uploaded metrics_analyzer.py has an early-return for
    `len(cc_values) < 2` ("first refinement" branch).  Tom's actually
    deployed version differs (it auto-stopped on AF_7mjs which had
    only 1 RSR cycle — see session JSON).  This test uses 2 RSR cycles
    so it hits the SUCCESS branch in the uploaded file; the SEMANTIC
    fix (validation_done check) applies identically to both versions.
    """
    history = [
        {"program": "phenix.mtriage"},
        {"program": "phenix.resolve_cryo_em"},
        {"program": "phenix.predict_and_build", "map_cc": 0.61},
        _rsr(0.72),
        _rsr(0.7832),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    assert_false(r.get("should_stop"),
                 "AF_7mjs: should_stop must be False when CC>0.70 and validation not done")
    assert_eq(r.get("recommendation"), "validate",
              "AF_7mjs: recommendation should be 'validate'")
    assert_true(r.get("suggest_validation"),
                "AF_7mjs: suggest_validation should be True")


def test_cc_above_target_with_molprobity_done():
    """When CC > 0.70 AND phenix.molprobity has been run, auto-stop is OK."""
    history = [
        _rsr(0.62),
        _rsr(0.78),
        _molprobity(),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    assert_true(r.get("should_stop"),
                "CC>0.70 + molprobity done: should_stop must be True")
    assert_eq(r.get("recommendation"), "stop",
              "CC>0.70 + molprobity done: recommendation should be 'stop'")


def test_cc_above_target_with_validation_cryoem_done():
    """When CC > 0.70 AND phenix.validation_cryoem has been run, auto-stop is OK."""
    history = [
        _rsr(0.62),
        _rsr(0.78),
        _validation_cryoem(),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    assert_true(r.get("should_stop"),
                "CC>0.70 + validation_cryoem done: should_stop must be True")
    assert_eq(r.get("recommendation"), "stop",
              "CC>0.70 + validation_cryoem done: recommendation should be 'stop'")


def test_cc_above_target_validation_not_done():
    """When CC > 0.70 AND no validation in history: defer to validation."""
    history = [
        _rsr(0.62),
        _rsr(0.78),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    assert_false(r.get("should_stop"),
                 "CC>0.70 + no validation: should_stop must be False")
    assert_eq(r.get("recommendation"), "validate",
              "CC>0.70 + no validation: recommendation should be 'validate'")
    assert_true(r.get("suggest_validation"),
                "CC>0.70 + no validation: suggest_validation should be True")
    # Reason should mention validation
    reason = r.get("reason", "")
    assert_true("validation" in reason.lower() or "validate" in reason.lower(),
                "Reason should mention validation: got %r" % reason)


def test_cc_above_target_no_validation_trend_summary():
    """The trend_summary should signal validation is needed (for UI)."""
    history = [
        _rsr(0.62),
        _rsr(0.78),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    summary = r.get("trend_summary", "")
    assert_true("TARGET REACHED" in summary,
                "trend_summary should contain 'TARGET REACHED': got %r" % summary)
    assert_true("VALIDATE" in summary,
                "trend_summary should contain 'VALIDATE': got %r" % summary)


def test_cc_above_target_with_validation_trend_summary():
    """When validation IS done, trend_summary should NOT say VALIDATE BEFORE STOPPING."""
    history = [
        _rsr(0.62),
        _rsr(0.78),
        _molprobity(),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    summary = r.get("trend_summary", "")
    assert_false("VALIDATE BEFORE STOPPING" in summary,
                 "With validation done, trend_summary should NOT urge validation: got %r" % summary)
    assert_true("ABOVE TARGET" in summary,
                "With validation done, trend_summary should say ABOVE TARGET: got %r" % summary)


def test_cc_just_above_threshold_validation_not_done():
    """CC just above 0.70 (e.g., 0.71): still need validation."""
    history = [
        _rsr(0.65),
        _rsr(0.71),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    assert_false(r.get("should_stop"),
                 "CC=0.71 (just above) + no validation: should_stop must be False")


def test_cc_exactly_at_threshold():
    """CC = 0.70 exactly (not strictly above): falls through to other checks.
    The condition is `latest_cc > 0.70` (strict), so CC=0.70 doesn't hit
    the SUCCESS branch."""
    history = [
        _rsr(0.65),
        _rsr(0.70),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    # CC=0.70 (not strictly >), so the SUCCESS branch doesn't fire.
    # Result depends on PLATEAU/EXCESSIVE/WORSENING — none of those should
    # match with 2 cycles. So no should_stop, no suggest_validation.
    assert_false(r.get("should_stop"),
                 "CC=0.70 (not strictly >): SUCCESS branch should not fire")
    assert_false(r.get("suggest_validation"),
                 "CC=0.70: suggest_validation should not be set")


def test_cc_below_target_one_cycle():
    """First RSR cycle with CC < 0.70: no stop conditions apply."""
    history = [
        _rsr(0.55),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    assert_false(r.get("should_stop"),
                 "First RSR with CC<0.70: should_stop must be False")
    # First refinement special case
    assert_true("first refinement" in r.get("trend_summary", "").lower(),
                "First RSR: trend_summary should mention 'first refinement'")


def test_cc_below_target_multiple_cycles_no_stop():
    """Multiple RSR cycles with CC still below 0.70: no auto-stop."""
    history = [
        _rsr(0.55),
        _rsr(0.60),
        _rsr(0.63),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    assert_false(r.get("should_stop"),
                 "CC<0.70 through multiple cycles: should_stop must be False")


def test_excessive_rsr_still_stops_regardless_of_validation():
    """If we hit 5+ consecutive RSR cycles, EXCESSIVE stop still fires.
    This is a separate stop condition (loop guard), independent of the
    validation_done check."""
    history = [_rsr(0.60 + 0.005 * i) for i in range(5)]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    assert_true(r.get("should_stop"),
                "5 consecutive RSR cycles: EXCESSIVE stop must fire even without validation")
    assert_true("EXCESSIVE" in r.get("reason", ""),
                "EXCESSIVE: reason must mention 'EXCESSIVE'")


def test_xray_path_unchanged_validation_done():
    """The X-ray path was already correct; verify our cryo-EM change didn't break it.
    R-free well below target AND molprobity done: should_stop=True.

    For resolution=2.0Å: dynamic_target = 0.20, success_threshold = 0.18.
    Use r_free=0.15 to be clearly below the success threshold.
    """
    history = [
        {"program": "phenix.refine", "r_free": 0.28},
        {"program": "phenix.refine", "r_free": 0.15},
        {"program": "phenix.molprobity", "clashscore": 5.0},
    ]
    r = analyze_metrics_trend(history, experiment_type="xray", resolution=2.0)
    assert_true(r.get("should_stop"),
                "X-ray: R-free<target + molprobity done: should_stop must be True")


def test_xray_path_unchanged_validation_not_done():
    """X-ray: R-free below target but NO validation -> should_stop=False.

    Use r_free=0.15 (well below success_threshold=0.18 for 2.0Å).
    """
    history = [
        {"program": "phenix.refine", "r_free": 0.28},
        {"program": "phenix.refine", "r_free": 0.15},
    ]
    r = analyze_metrics_trend(history, experiment_type="xray", resolution=2.0)
    assert_false(r.get("should_stop"),
                 "X-ray: R-free<target + no validation: should_stop must be False")
    assert_true(r.get("suggest_validation"),
                "X-ray: suggest_validation should be True")


def test_phenix_model_vs_data_NOT_a_cryoem_validation():
    """phenix.model_vs_data is the X-ray validation program; it should NOT
    satisfy the cryo-EM validation_done check.  Cryo-EM uses phenix.molprobity
    and phenix.validation_cryoem."""
    history = [
        _rsr(0.62),
        _rsr(0.78),
        {"program": "phenix.model_vs_data", "r_free": 0.20},  # X-ray validation
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    assert_false(r.get("should_stop"),
                 "Cryo-EM with only model_vs_data (X-ray validation): should_stop must be False")
    assert_true(r.get("suggest_validation"),
                "Cryo-EM: model_vs_data does not count as cryo-EM validation")


def test_validation_at_any_history_position_counts():
    """Validation done EARLIER in history (not the latest cycle) still counts.
    This matches the X-ray behaviour (any/all in history)."""
    history = [
        _rsr(0.62),
        _molprobity(),       # Validation in the middle
        _rsr(0.75),
        _rsr(0.78),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    assert_true(r.get("should_stop"),
                "Validation earlier in history still counts: should_stop must be True")


def test_consecutive_count_includes_non_rsr_separator():
    """consecutive_rsr counts consecutive real_space_refine at the END of
    history. Validation (molprobity) is not RSR, so it breaks consecutive."""
    history = [
        _rsr(0.60),
        _rsr(0.65),
        _rsr(0.72),
        _molprobity(),  # breaks consecutive streak
        _rsr(0.78),
    ]
    r = analyze_metrics_trend(history, experiment_type="cryoem")
    # Validation IS done, so should_stop=True for the SUCCESS branch
    assert_true(r.get("should_stop"),
                "After validation, single RSR with CC>0.70: should_stop must be True (SUCCESS)")
    # And consecutive_rsr should be 1 (only the last entry)
    assert_eq(r.get("consecutive_rsr"), 1,
              "consecutive_rsr should be 1 (broken by molprobity)")


# =============================================================================
# Test runner (matches the project's tst_*.py convention)
# =============================================================================

def run_all_tests():
    try:
        from libtbx.langchain.tests.tst_utils import (
            run_tests_with_fail_fast)
    except ImportError:
        try:
            from tests.tst_utils import run_tests_with_fail_fast
        except ImportError:
            _standalone_runner()
            return
    run_tests_with_fail_fast()


def _standalone_runner():
    test_fns = [v for k, v in sorted(globals().items())
                if k.startswith("test_") and callable(v)]
    for fn in test_fns:
        try:
            fn()
        except Exception as e:
            global FAIL, FAILURES
            FAIL += 1
            FAILURES.append("%s: raised %s: %s" % (
                fn.__name__, type(e).__name__, e))
    print("%d passed, %d failed" % (PASS, FAIL))
    if FAILURES:
        for f in FAILURES:
            print("  FAIL: %s" % f)
    if FAIL:
        sys.exit(1)


if __name__ == "__main__":
    _standalone_runner()
