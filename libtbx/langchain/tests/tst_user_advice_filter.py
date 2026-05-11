"""Tests for v116.10 Phase 1 fix to _apply_user_advice in rules_selector.

The bug: `_apply_user_advice` was over-aggressively treating any
"stop" mention as an immediate-stop request, truncating
`valid_programs` to `[STOP]`.  Phrasings like "predict and stop"
or "refine then stop" are stop-TARGET requests (do X, then stop),
not immediate-stop requests.

The fix: extend `stop_condition_patterns` with "and stop",
"then stop", and ", stop".

These tests run by stubbing the libtbx import path so the module
can be imported and exercised standalone.  All tests are pure
functional checks of the filter logic.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import types


# --- Path setup + libtbx stubs ---------------------------------------------
# rules_selector imports from libtbx.langchain.agent.*.  We stub those
# packages so the import succeeds outside a PHENIX environment.

_HERE = os.path.dirname(os.path.abspath(__file__))
_AGENT_DIR = os.path.normpath(os.path.join(_HERE, "..", "agent"))
if _AGENT_DIR not in sys.path:
    sys.path.insert(0, _AGENT_DIR)


def _install_stubs():
    """Install stub modules for libtbx.langchain.agent.* dependencies
    of rules_selector.py.  Only the symbols actually imported are
    provided; the constructors return objects with attributes the
    test path doesn't touch."""

    class _StubWorkflowEngine:
        def __init__(self):
            pass

    class _StubMetricEvaluator:
        def __init__(self):
            pass

    class _StubProgramRegistry:
        def __init__(self):
            pass

        def get_user_advice_keywords(self, prog):
            return []  # Tests use generic programs with no keywords

    # Build the nested package structure
    def _ensure(modname):
        if modname not in sys.modules:
            sys.modules[modname] = types.ModuleType(modname)
        return sys.modules[modname]

    _ensure("libtbx")
    _ensure("libtbx.langchain")
    _ensure("libtbx.langchain.agent")

    we = _ensure("libtbx.langchain.agent.workflow_engine")
    we.WorkflowEngine = _StubWorkflowEngine

    me = _ensure("libtbx.langchain.agent.metric_evaluator")
    me.MetricEvaluator = _StubMetricEvaluator

    pr = _ensure("libtbx.langchain.agent.program_registry")
    pr.ProgramRegistry = _StubProgramRegistry

    # rules_selector also references program_registry in select_next_action;
    # we don't call that here.


# Only install stubs if libtbx isn't already importable (i.e. when
# running outside a PHENIX environment)
try:
    import libtbx.langchain.agent.rules_selector as _real
    _real_imported = True
except ImportError:
    _install_stubs()
    _real_imported = False


try:
    from libtbx.langchain.agent.rules_selector import RulesSelector
except ImportError:
    from rules_selector import RulesSelector


# =====================================================================
# Test fixtures
# =====================================================================
# A typical xray_initial valid_programs after the workflow engine has
# computed it but before _apply_user_advice runs.
_VALID_INITIAL = ["phenix.xtriage", "STOP"]

# A multi-program valid set (post-analyze)
_VALID_ANALYZED = ["phenix.predict_and_build", "phenix.phaser",
                   "phenix.autosol", "STOP"]


def _make_selector():
    """Construct a RulesSelector for testing.

    The selector's __init__ calls WorkflowEngine(), MetricEvaluator(),
    and ProgramRegistry().  In the standalone test environment these
    are stubs; in PHENIX they're real.  Either way, _apply_user_advice
    only uses self.program_registry.get_user_advice_keywords().
    """
    return RulesSelector()


# =====================================================================
# SECTION A: The bug — "X and stop" was being truncated to [STOP]
# =====================================================================

def test_predict_and_stop_preserves_valid_programs():
    """'predict and stop' must NOT truncate valid_programs to [STOP].

    Regression test for the v116.10 bug: this advice phrasing was
    treated as immediate stop because "and" alone isn't in
    sequencing_words and no stop-condition pattern matched.
    """
    print("Test: predict_and_stop_preserves_valid_programs")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_INITIAL, "predict and stop")
    assert result == _VALID_INITIAL, (
        "'predict and stop' should preserve valid_programs intact, "
        "got %s (expected %s)" % (result, _VALID_INITIAL))
    print("  PASS")


def test_refine_then_stop_preserves_valid_programs():
    """'refine then stop' is a sequence with a stop target."""
    print("Test: refine_then_stop_preserves_valid_programs")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_ANALYZED, "refine then stop")
    # "refine" matches phenix.refine — wait, refine isn't in
    # _VALID_ANALYZED, so no single-program return.  Goes to
    # stop branch.  With fix: "then stop" matches → not immediate.
    assert "STOP" in result, "STOP should remain valid"
    assert len(result) > 1, (
        "'refine then stop' should NOT truncate to [STOP], got %s"
        % result)
    print("  PASS")


def test_comma_stop_preserves_valid_programs():
    """', stop' is a comma-separated stop target.

    Tests two sub-cases:
    (a) 'xtriage, stop' — 'xtriage' matches phenix.xtriage in the
        single-step branch, which returns early before reaching the
        stop branch.  The result is [phenix.xtriage] (no STOP).
    (b) 'predict, stop' — 'predict' doesn't match any valid program,
        so we fall through to the stop branch.  This is the case
        that exercises the new ', stop' pattern.
    """
    print("Test: comma_stop_preserves_valid_programs")
    sel = _make_selector()

    # Sub-case (a): program-name match in single-step branch
    result_a = sel._apply_user_advice(_VALID_INITIAL, "xtriage, stop")
    assert "phenix.xtriage" in result_a, (
        "'xtriage, stop' should keep phenix.xtriage available, got %s"
        % result_a)

    # Sub-case (b): the ', stop' pattern path
    result_b = sel._apply_user_advice(_VALID_INITIAL, "predict, stop")
    assert "STOP" in result_b, "STOP should remain available"
    assert len(result_b) > 1, (
        "'predict, stop' should NOT truncate to [STOP], got %s"
        % result_b)
    print("  PASS")


def test_xtriage_and_stop_preserves_valid_programs():
    """Variant: a real program name + 'and stop'.

    Note: the single-step branch may match 'xtriage' to
    phenix.xtriage and return early.  In that case the stop-
    keyword branch never runs.  This test verifies that either
    path produces a non-truncated result.
    """
    print("Test: xtriage_and_stop_preserves_valid_programs")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_INITIAL, "xtriage and stop")
    # Either: single-step picks xtriage → returns [phenix.xtriage]
    # Or: falls through, hits stop branch, "and stop" matches → preserved
    # Both outcomes preserve xtriage as a runnable option:
    assert "phenix.xtriage" in result, (
        "'xtriage and stop' should keep phenix.xtriage available, got %s"
        % result)
    print("  PASS")


def test_predict_then_stop_preserves_valid_programs():
    """Variant of 'predict and stop' using 'then'."""
    print("Test: predict_then_stop_preserves_valid_programs")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_INITIAL, "predict then stop")
    assert result == _VALID_INITIAL, (
        "'predict then stop' should preserve valid_programs, got %s"
        % result)
    print("  PASS")


# =====================================================================
# SECTION B: Regression tests — immediate-stop requests still work
# =====================================================================

def test_stop_now_returns_stop_only():
    """'stop now' is an immediate-stop request.  Must still truncate."""
    print("Test: stop_now_returns_stop_only")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_INITIAL, "stop now")
    assert result == ["STOP"], (
        "'stop now' should return [STOP], got %s" % result)
    print("  PASS")


def test_please_stop_returns_stop_only():
    """'please stop' is an immediate-stop request."""
    print("Test: please_stop_returns_stop_only")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_INITIAL, "please stop")
    assert result == ["STOP"], (
        "'please stop' should return [STOP], got %s" % result)
    print("  PASS")


def test_bare_stop_returns_stop_only():
    """A bare 'stop' is an immediate-stop request."""
    print("Test: bare_stop_returns_stop_only")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_INITIAL, "stop")
    assert result == ["STOP"], (
        "'stop' should return [STOP], got %s" % result)
    print("  PASS")


# =====================================================================
# SECTION C: Existing stop-condition patterns still work
# =====================================================================

def test_stop_after_X_preserves_valid_programs():
    """'stop after refine' — the original supported pattern."""
    print("Test: stop_after_X_preserves_valid_programs")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_ANALYZED,
                                    "stop after refine completes")
    assert "STOP" in result
    assert len(result) > 1, (
        "Existing 'stop after X' should preserve valid_programs, got %s"
        % result)
    print("  PASS")


def test_stop_when_X_preserves_valid_programs():
    """'stop when R-free converges' — existing pattern."""
    print("Test: stop_when_X_preserves_valid_programs")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_ANALYZED,
                                    "stop when R-free converges")
    assert "STOP" in result
    assert len(result) > 1
    print("  PASS")


def test_stop_at_X_preserves_valid_programs():
    """'stop at R=0.20' — existing pattern."""
    print("Test: stop_at_X_preserves_valid_programs")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_ANALYZED,
                                    "stop at R-free 0.20")
    assert "STOP" in result
    assert len(result) > 1
    print("  PASS")


# =====================================================================
# SECTION D: No-stop advice (control)
# =====================================================================

def test_no_stop_mention_no_change():
    """Advice without 'stop' shouldn't touch the stop branch at all."""
    print("Test: no_stop_mention_no_change")
    sel = _make_selector()
    result = sel._apply_user_advice(
        _VALID_ANALYZED, "use phaser for molecular replacement")
    # "phaser" matches phenix.phaser in single-step branch
    # → returns [phenix.phaser]
    assert result == ["phenix.phaser"], (
        "'use phaser' should select phaser, got %s" % result)
    print("  PASS")


def test_empty_advice_no_filter():
    """Empty advice → no filtering."""
    print("Test: empty_advice_no_filter")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_ANALYZED, "")
    assert result == _VALID_ANALYZED
    print("  PASS")


# =====================================================================
# SECTION E: Edge cases
# =====================================================================

def test_stop_and_predict_still_treated_as_immediate():
    """'stop and predict' has 'and' but not in 'X and stop' form.

    The 'and stop' substring requires 'and' to be immediately before
    'stop'.  'stop and predict' does NOT contain 'and stop'.
    """
    print("Test: stop_and_predict_still_treated_as_immediate")
    sel = _make_selector()
    advice = "stop and predict"
    # No program name match — none of valid programs in advice.
    # Falls to stop branch.  "and stop" NOT in advice (the order
    # is "stop and", not "and stop").
    # is_stop_condition = False → returns [STOP]
    result = sel._apply_user_advice(_VALID_INITIAL, advice)
    assert result == ["STOP"], (
        "'stop and predict' should be treated as immediate stop, got %s"
        % result)
    print("  PASS")


def test_capitalized_stop_phrasing():
    """The function lowercases advice before matching — uppercase still works."""
    print("Test: capitalized_stop_phrasing")
    sel = _make_selector()
    result = sel._apply_user_advice(_VALID_INITIAL, "Predict And Stop")
    assert result == _VALID_INITIAL, (
        "Mixed-case 'Predict And Stop' should preserve valid_programs, got %s"
        % result)
    print("  PASS")


def test_multiword_advice_with_stop_target():
    """Realistic phrasing: 'do alphafold prediction and stop'."""
    print("Test: multiword_advice_with_stop_target")
    sel = _make_selector()
    advice = "do alphafold prediction and stop"
    result = sel._apply_user_advice(_VALID_INITIAL, advice)
    assert "STOP" in result
    assert len(result) > 1, (
        "Realistic prediction-and-stop should preserve valid_programs, got %s"
        % result)
    print("  PASS")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    """Run all tests using cctbx fail-fast convention."""
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
    """Simple runner when tst_utils is not available."""
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
