"""Tests for v116.13 cryo-EM validation_done check in metric_evaluator.

Bug: v116.12 patched `agent/metrics_analyzer.py` thinking it was
the active code path.  It wasn't.  graph_nodes.py sets
`USE_YAML_METRICS = True` (default).  With this flag set,
`analyze_metrics_trend()` in metrics_analyzer.py delegates immediately
to `analyze_refinement_trend()` in `agent/metric_evaluator.py` and
returns:

    if use_yaml_evaluator and analyze_refinement_trend is not None:
        try:
            return analyze_refinement_trend(...)
        except Exception:
            ... # fallback to hardcoded (where v116.12 fix was applied)

The hardcoded fallback path almost never runs.  The actual production
path is `MetricEvaluator.analyze_trend()`, which had the same
asymmetry as the hardcoded path: the X-ray branch (line 403)
checked validation_done before allowing auto-stop, but the cryo-EM
branch (line ~500) did not.

This caused AF_7mjs to auto-stop after refinement (CC=0.827)
without running Stage 5 (phenix.molprobity), even with v116.12
"deployed" — the deployment was correct but the code was dead.

Root cause: missing validation_done check in
MetricEvaluator.analyze_trend() cryo-EM SUCCESS branch.

Fix: add the validation_done check, mirroring the X-ray pattern at
line 403.  Accepts both phenix.molprobity (general validation,
shared with X-ray) and phenix.validation_cryoem (cryo-EM-specific).
phenix.model_vs_data is X-ray-only; does NOT count for cryo-EM.

These tests verify:
  - AF_7mjs regression: CC>target + no validation → must NOT stop
  - validation programs recognized for cryo-EM (molprobity, validation_cryoem)
  - X-ray validation program (model_vs_data) does NOT count for cryo-EM
  - X-ray path unchanged
  - trend_summary "TARGET REACHED" format preserved
  - Edge cases: CC exactly at threshold (strict >), CC below target,
    validation earlier in history, MetricEvaluator class direct usage
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

# Stub libtbx so the import chain works without PHENIX.
# Both libtbx.langchain.knowledge.yaml_loader and the agent module
# are stubbed; the yaml_loader stub returns production-matching
# thresholds for r_free and map_cc.
if "libtbx" not in sys.modules:
    import types as _types
    for _mod in (
        "libtbx", "libtbx.langchain",
        "libtbx.langchain.agent",
        "libtbx.langchain.knowledge",
        "libtbx.langchain.knowledge.yaml_loader",
    ):
        sys.modules[_mod] = _types.ModuleType(_mod)

    yl = sys.modules["libtbx.langchain.knowledge.yaml_loader"]

    def _get_metric_threshold(name, level=None, resolution=None):
        # Production thresholds from metrics.yaml.  Cryo-EM
        # map_cc.acceptable = 0.70 is the success threshold;
        # get_target("map_cc") returns this via
        # get_threshold(metric, "acceptable").
        thresholds = {
            "r_free": {"acceptable": 0.20, "good": 0.18, "excellent": 0.15,
                       "target": 0.20},
            "map_cc": {"acceptable": 0.70, "good": 0.80, "excellent": 0.90,
                       "target": 0.70},
        }
        if name in thresholds:
            if level in thresholds[name]:
                return thresholds[name][level]
            return thresholds[name].get("target")
        return None

    yl.get_metric_threshold = _get_metric_threshold
    yl.get_metric_direction = lambda name: (
        "lower_is_better" if name == "r_free" else "higher_is_better")
    yl.is_metric_good = lambda name, val, resolution=None: False
    yl.is_metric_acceptable = lambda name, val, resolution=None: True
    yl.load_metrics = lambda: {}

try:
    from libtbx.langchain.agent.metric_evaluator import (
        analyze_refinement_trend, MetricEvaluator)
except ImportError:
    from agent.metric_evaluator import (
        analyze_refinement_trend, MetricEvaluator)


# =============================================================================
# Test helpers
# =============================================================================

def _rsr(cc):
    """Build a phenix.real_space_refine metrics_history entry."""
    return {"program": "phenix.real_space_refine", "map_cc": cc}


def _molprobity():
    """Build a phenix.molprobity metrics_history entry."""
    return {"program": "phenix.molprobity", "clashscore": 5.0}


def _validation_cryoem():
    """Build a phenix.validation_cryoem metrics_history entry."""
    return {"program": "phenix.validation_cryoem", "map_cc": 0.78}


# =============================================================================
# Primary AF_7mjs regression test
# =============================================================================

def test_af7mjs_regression():
    """The AF_7mjs scenario: cycle 4 RSR produced CC=0.827,
    validation NOT done.  Pre-v116.13 this auto-stopped before
    running Stage 5 (phenix.molprobity).  Post-fix should defer."""
    print("Test: af7mjs_regression")
    history = [
        {"program": "phenix.mtriage"},
        {"program": "phenix.resolve_cryo_em"},
        {"program": "phenix.predict_and_build", "map_cc": 0.61},
        _rsr(0.827),
    ]
    r = analyze_refinement_trend(history, experiment_type="cryoem")

    assert r.get("should_stop") is False, (
        "AF_7mjs regression: should_stop must be False when CC>target "
        "and validation not done.  Got: %r" % r.get("should_stop"))
    assert r.get("recommendation") == "validate", (
        "AF_7mjs: recommendation should be 'validate'.  Got: %r"
        % r.get("recommendation"))
    assert r.get("suggest_validation") is True, (
        "AF_7mjs: suggest_validation should be True.  Got: %r"
        % r.get("suggest_validation"))
    print("  PASS")


# =============================================================================
# Validation program recognition
# =============================================================================

def test_cc_above_target_with_molprobity_done():
    """CC>target + phenix.molprobity in history → should_stop=True.
    Matches the pre-v116.13 SUCCESS behavior for the happy path."""
    print("Test: cc_above_target_with_molprobity_done")
    history = [
        _rsr(0.65),
        _rsr(0.78),
        _molprobity(),
    ]
    r = analyze_refinement_trend(history, experiment_type="cryoem")

    assert r.get("should_stop") is True, (
        "CC>target + molprobity done: should_stop must be True.  Got: %r"
        % r.get("should_stop"))
    assert r.get("recommendation") == "stop", (
        "CC>target + molprobity done: recommendation should be 'stop'.  "
        "Got: %r" % r.get("recommendation"))
    print("  PASS")


def test_cc_above_target_with_validation_cryoem_done():
    """CC>target + phenix.validation_cryoem in history → should_stop=True.
    validation_cryoem is the cryo-EM-specific validation program."""
    print("Test: cc_above_target_with_validation_cryoem_done")
    history = [
        _rsr(0.65),
        _rsr(0.78),
        _validation_cryoem(),
    ]
    r = analyze_refinement_trend(history, experiment_type="cryoem")

    assert r.get("should_stop") is True, (
        "CC>target + validation_cryoem done: should_stop must be True.  "
        "Got: %r" % r.get("should_stop"))
    print("  PASS")


def test_phenix_model_vs_data_not_cryoem_validation():
    """phenix.model_vs_data is X-ray-only; should NOT count for cryo-EM
    validation.  Mirrors the X-ray branch which lists model_vs_data
    but not validation_cryoem."""
    print("Test: phenix_model_vs_data_not_cryoem_validation")
    history = [
        _rsr(0.65),
        _rsr(0.78),
        {"program": "phenix.model_vs_data", "r_free": 0.20},
    ]
    r = analyze_refinement_trend(history, experiment_type="cryoem")

    assert r.get("should_stop") is False, (
        "Cryo-EM with only model_vs_data (X-ray val): should NOT stop.  "
        "Got: %r" % r.get("should_stop"))
    assert r.get("suggest_validation") is True, (
        "Cryo-EM: model_vs_data does not count as cryo-EM validation; "
        "suggest_validation should be True.  Got: %r"
        % r.get("suggest_validation"))
    print("  PASS")


# =============================================================================
# Validation history positioning
# =============================================================================

def test_validation_earlier_in_history_counts():
    """Validation done earlier in history still counts.  Uses
    any() over full metrics_history, not just the latest cycle."""
    print("Test: validation_earlier_in_history_counts")
    history = [
        _rsr(0.62),
        _molprobity(),       # Validation in the middle
        _rsr(0.78),
    ]
    r = analyze_refinement_trend(history, experiment_type="cryoem")

    assert r.get("should_stop") is True, (
        "Validation earlier in history still counts: should_stop must "
        "be True.  Got: %r" % r.get("should_stop"))
    print("  PASS")


# =============================================================================
# Defer-to-validation behavior
# =============================================================================

def test_cc_above_target_validation_not_done():
    """CC>target + no validation → should_stop=False with validation hint."""
    print("Test: cc_above_target_validation_not_done")
    history = [_rsr(0.65), _rsr(0.78)]
    r = analyze_refinement_trend(history, experiment_type="cryoem")

    assert r.get("should_stop") is False, (
        "CC>target + no validation: should_stop must be False.  Got: %r"
        % r.get("should_stop"))
    assert r.get("recommendation") == "validate", (
        "CC>target + no validation: recommendation should be 'validate'.  "
        "Got: %r" % r.get("recommendation"))
    assert r.get("suggest_validation") is True, (
        "CC>target + no validation: suggest_validation should be True.  "
        "Got: %r" % r.get("suggest_validation"))
    print("  PASS")


def test_reason_mentions_validation():
    """When validation not done, the reason should hint why we're not stopping."""
    print("Test: reason_mentions_validation")
    history = [_rsr(0.65), _rsr(0.78)]
    r = analyze_refinement_trend(history, experiment_type="cryoem")
    reason = r.get("reason", "")

    assert "validation" in reason.lower() or "validate" in reason.lower(), (
        "reason should mention validation when deferring auto-stop.  "
        "Got: %r" % reason)
    print("  PASS")


def test_trend_summary_format_preserved():
    """trend_summary 'TARGET REACHED' format is preserved exactly,
    so the cycle box display doesn't change.  The validation advisory
    lives in reason/recommendation, not trend_summary."""
    print("Test: trend_summary_format_preserved")
    history = [_rsr(0.65), _rsr(0.78)]
    r = analyze_refinement_trend(history, experiment_type="cryoem")
    summary = r.get("trend_summary", "")

    assert "TARGET REACHED" in summary, (
        "trend_summary must contain 'TARGET REACHED' (preserved format).  "
        "Got: %r" % summary)
    print("  PASS")


# =============================================================================
# Boundary / edge cases
# =============================================================================

def test_cc_just_above_threshold():
    """CC just above target (0.71) still defers to validation —
    same as CC=0.78.  The threshold is binary, not graded."""
    print("Test: cc_just_above_threshold")
    history = [_rsr(0.65), _rsr(0.71)]
    r = analyze_refinement_trend(history, experiment_type="cryoem")

    assert r.get("should_stop") is False, (
        "CC=0.71 (just above target) + no validation: should_stop must "
        "be False.  Got: %r" % r.get("should_stop"))
    print("  PASS")


def test_cc_exactly_at_threshold_no_success():
    """CC = 0.70 exactly does NOT strictly exceed target (uses `>`).
    SUCCESS branch doesn't fire; falls through to PLATEAU/EXCESSIVE/continue."""
    print("Test: cc_exactly_at_threshold_no_success")
    history = [_rsr(0.65), _rsr(0.70)]
    r = analyze_refinement_trend(history, experiment_type="cryoem")

    # SUCCESS branch should NOT fire for CC=0.70 (strict >).
    # If it doesn't fire, suggest_validation is not set by it.
    assert r.get("suggest_validation", False) is False, (
        "CC=0.70 (not strictly >): SUCCESS branch should not fire; "
        "suggest_validation should not be set.  Got: %r"
        % r.get("suggest_validation"))
    print("  PASS")


def test_cc_below_target_no_success_suggestion():
    """CC < target: should_stop=False, no validation suggestion
    from SUCCESS branch (it doesn't fire)."""
    print("Test: cc_below_target_no_success_suggestion")
    history = [_rsr(0.55), _rsr(0.65)]
    r = analyze_refinement_trend(history, experiment_type="cryoem")

    assert r.get("should_stop") is False, (
        "CC<target: should_stop must be False.  Got: %r"
        % r.get("should_stop"))
    assert r.get("suggest_validation", False) is False, (
        "CC<target: suggest_validation must not be set by SUCCESS branch "
        "(it doesn't fire below target).  Got: %r"
        % r.get("suggest_validation"))
    print("  PASS")


# =============================================================================
# X-ray path unchanged (regression guard)
# =============================================================================

def test_xray_path_unchanged_validation_done():
    """X-ray path was already correct; verify cryo-EM change didn't
    break it.  R-free=0.15 with resolution=2.0Å (target=0.20, threshold=0.18)
    + molprobity done → should_stop=True."""
    print("Test: xray_path_unchanged_validation_done")
    history = [
        {"program": "phenix.refine", "r_free": 0.28},
        {"program": "phenix.refine", "r_free": 0.15},
        {"program": "phenix.molprobity", "clashscore": 5.0},
    ]
    r = analyze_refinement_trend(history, experiment_type="xray",
                                   resolution=2.0)

    assert r.get("should_stop") is True, (
        "X-ray: R-free<threshold + molprobity done: should_stop must be "
        "True.  Got: %r" % r.get("should_stop"))
    print("  PASS")


def test_xray_path_unchanged_validation_not_done():
    """X-ray: R-free<threshold + no validation → should_stop=False
    (X-ray path's existing validation_done check, unchanged)."""
    print("Test: xray_path_unchanged_validation_not_done")
    history = [
        {"program": "phenix.refine", "r_free": 0.28},
        {"program": "phenix.refine", "r_free": 0.15},
    ]
    r = analyze_refinement_trend(history, experiment_type="xray",
                                   resolution=2.0)

    assert r.get("should_stop") is False, (
        "X-ray: R-free<threshold + no validation: should_stop must be "
        "False (existing X-ray behavior preserved).  Got: %r"
        % r.get("should_stop"))
    print("  PASS")


# =============================================================================
# Class-level direct usage (the fix is in MetricEvaluator.analyze_trend, not
# just the analyze_refinement_trend wrapper)
# =============================================================================

def test_using_evaluator_class_directly():
    """Verify the change is in MetricEvaluator.analyze_trend, not
    just the module-level wrapper.  Callers that use the class
    directly (e.g. phenix_ai internals) must see the same behavior."""
    print("Test: using_evaluator_class_directly")
    ev = MetricEvaluator()
    history = [_rsr(0.65), _rsr(0.78)]
    r = ev.analyze_trend(history, experiment_type="cryoem")

    assert r.get("should_stop") is False, (
        "MetricEvaluator.analyze_trend direct call: CC>target + no "
        "validation → should_stop must be False.  Got: %r"
        % r.get("should_stop"))
    print("  PASS")


# =============================================================================
# Existing fields preserved
# =============================================================================

def test_consecutive_field_still_set():
    """The fix shouldn't remove existing result fields like
    consecutive_rsr.  Regression guard against accidental dict
    rebuilding."""
    print("Test: consecutive_field_still_set")
    history = [_rsr(0.62), _rsr(0.66), _rsr(0.78)]
    r = analyze_refinement_trend(history, experiment_type="cryoem")

    assert r.get("consecutive_rsr") is not None, (
        "consecutive_rsr field must still exist in the result.  "
        "Got keys: %r" % list(r.keys()))
    print("  PASS")


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
    passed = 0
    failed = 0
    for fn in test_fns:
        try:
            fn()
            passed += 1
        except Exception as e:
            print("  FAIL: %s" % e)
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    _standalone_runner()
