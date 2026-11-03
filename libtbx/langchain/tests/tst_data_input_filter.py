"""Tests for v116.10 Phase 4b — xtriage/mtriage filter when no data.

The bug: xtriage and mtriage are added unconditionally by the analyze
step in workflows.yaml, but they cannot run without their respective
data inputs.  For sequence-only sessions (X-ray with .seq but no .mtz),
valid_programs would include phenix.xtriage, the LLM might pick it,
and the command builder would generate a command that fails at runtime.

The fix: a new helper `_filter_programs_missing_data_inputs` removes
xtriage when no .mtz data is present, and removes mtriage when no map
is present.  Called from `get_valid_programs` after the run_once filter
and before the directive-aware path.

These tests verify the helper in isolation.  Integration tests against
the tutorial corpus (run_29_openai, run_29_ollama) verify the end-to-end
behavior.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import types


# --- Path setup + libtbx stubs ---------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
_AGENT_DIR = os.path.normpath(os.path.join(_HERE, "..", "agent"))
if _AGENT_DIR not in sys.path:
    sys.path.insert(0, _AGENT_DIR)


def _install_stubs():
    """Install stub modules for libtbx imports.

    workflow_engine.py imports from libtbx.langchain.knowledge.yaml_loader
    and libtbx.langchain.agent.program_registry.  The filter we're testing
    does not call into either, but Python still needs the imports to
    resolve at module load.
    """

    class _Stub:
        def __init__(self):
            pass

    def _ensure(modname):
        if modname not in sys.modules:
            sys.modules[modname] = types.ModuleType(modname)
        return sys.modules[modname]

    _ensure("libtbx")
    _ensure("libtbx.langchain")
    _ensure("libtbx.langchain.agent")
    _ensure("libtbx.langchain.knowledge")

    yl = _ensure("libtbx.langchain.knowledge.yaml_loader")
    yl.get_workflow_steps = lambda et: {}
    yl.get_workflow_targets = lambda et, m: None
    yl.get_metric_threshold = lambda et, m: None
    yl.get_program = lambda p: None

    pr = _ensure("libtbx.langchain.agent.program_registry")
    pr.ProgramRegistry = _Stub


try:
    from libtbx.langchain.agent.workflow_engine import WorkflowEngine
except ImportError:
    _install_stubs()
    try:
        from libtbx.langchain.agent.workflow_engine import WorkflowEngine
    except ImportError:
        from workflow_engine import WorkflowEngine


# Construct the engine.  Heavy dependencies are stubbed but ProgramRegistry
# is built-in shape, so this just works.
def _make_engine():
    return WorkflowEngine()


# =====================================================================
# SECTION A: X-ray xtriage filtering
# =====================================================================

def test_xtriage_removed_when_no_xray_data():
    """xtriage at xray_initial with no .mtz → removed."""
    print("Test: xtriage_removed_when_no_xray_data")
    engine = _make_engine()
    valid = ["phenix.xtriage", "STOP"]
    context = {
        "has_sequence": True,
        "has_data_mtz": False,
        "has_phased_data_mtz": False,
    }
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert "phenix.xtriage" not in result, (
        "xtriage should be removed when no .mtz data present\n"
        "  Result: %s" % result)
    assert "STOP" in result, (
        "STOP should remain after xtriage removal\n  Result: %s" % result)
    print("  PASS")


def test_xtriage_kept_when_data_mtz_present():
    """xtriage at xray_initial with has_data_mtz → stays."""
    print("Test: xtriage_kept_when_data_mtz_present")
    engine = _make_engine()
    valid = ["phenix.xtriage", "STOP"]
    context = {"has_data_mtz": True}
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert "phenix.xtriage" in result, (
        "xtriage should be preserved when .mtz is present\n"
        "  Result: %s" % result)
    print("  PASS")


def test_xtriage_kept_when_phased_data_mtz_present():
    """xtriage with phased data .mtz (intensities embedded) → stays.

    Phased MTZ files include the original intensity data alongside the
    phase columns.  xtriage can analyze these for resolution, twinning,
    etc.  The prediction-only allowance at _check_program_prerequisites
    treats has_phased_data_mtz as equivalent to has_data_mtz, and this
    filter must match that contract.
    """
    print("Test: xtriage_kept_when_phased_data_mtz_present")
    engine = _make_engine()
    valid = ["phenix.xtriage", "STOP"]
    context = {
        "has_data_mtz": False,
        "has_phased_data_mtz": True,
    }
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert "phenix.xtriage" in result, (
        "xtriage should be preserved when phased .mtz is present\n"
        "  Result: %s" % result)
    print("  PASS")


def test_xtriage_kept_when_both_data_types_present():
    """xtriage with both has_data_mtz and has_phased_data_mtz → stays."""
    print("Test: xtriage_kept_when_both_data_types_present")
    engine = _make_engine()
    valid = ["phenix.xtriage", "STOP"]
    context = {
        "has_data_mtz": True,
        "has_phased_data_mtz": True,
    }
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert "phenix.xtriage" in result, (
        "xtriage should be preserved when both data types present\n"
        "  Result: %s" % result)
    print("  PASS")


# =====================================================================
# SECTION B: Cryo-EM mtriage filtering
# =====================================================================

def test_mtriage_removed_when_no_map():
    """mtriage at cryoem_initial with no map → removed."""
    print("Test: mtriage_removed_when_no_map")
    engine = _make_engine()
    valid = ["phenix.mtriage", "STOP"]
    context = {
        "has_sequence": True,
        "has_full_map": False,
        "has_half_map": False,
    }
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert "phenix.mtriage" not in result, (
        "mtriage should be removed when no map present\n"
        "  Result: %s" % result)
    assert "STOP" in result
    print("  PASS")


def test_mtriage_kept_when_full_map_present():
    """mtriage with full_map → stays."""
    print("Test: mtriage_kept_when_full_map_present")
    engine = _make_engine()
    valid = ["phenix.mtriage", "STOP"]
    context = {"has_full_map": True}
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert "phenix.mtriage" in result, (
        "mtriage should be preserved when full_map present\n"
        "  Result: %s" % result)
    print("  PASS")


def test_mtriage_kept_when_half_map_present():
    """mtriage with half_map → stays.

    Half-maps alone are enough for mtriage analysis even without a
    processed full map.  The map_correlations program needs a full
    map, but mtriage handles half-maps.
    """
    print("Test: mtriage_kept_when_half_map_present")
    engine = _make_engine()
    valid = ["phenix.mtriage", "STOP"]
    context = {
        "has_full_map": False,
        "has_half_map": True,
    }
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert "phenix.mtriage" in result, (
        "mtriage should be preserved when half_map present\n"
        "  Result: %s" % result)
    print("  PASS")


# =====================================================================
# SECTION C: No-op cases (programs not present, other programs preserved)
# =====================================================================

def test_no_xtriage_in_valid_no_change():
    """If xtriage isn't in valid_programs, filter is a no-op."""
    print("Test: no_xtriage_in_valid_no_change")
    engine = _make_engine()
    valid = ["phenix.refine", "phenix.molprobity", "STOP"]
    context = {"has_data_mtz": False}  # no data, but no xtriage to remove
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert result == valid, (
        "Filter should be a no-op when xtriage absent\n"
        "  Expected: %s\n  Got: %s" % (valid, result))
    print("  PASS")


def test_no_mtriage_in_valid_no_change():
    """If mtriage isn't in valid_programs, filter is a no-op."""
    print("Test: no_mtriage_in_valid_no_change")
    engine = _make_engine()
    valid = ["phenix.real_space_refine", "phenix.molprobity", "STOP"]
    context = {"has_full_map": False}
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert result == valid, (
        "Filter should be a no-op when mtriage absent\n"
        "  Expected: %s\n  Got: %s" % (valid, result))
    print("  PASS")


def test_other_programs_preserved_when_xtriage_removed():
    """Removing xtriage must not affect other programs in valid."""
    print("Test: other_programs_preserved_when_xtriage_removed")
    engine = _make_engine()
    # Hypothetical: xtriage gets included with prediction-only candidates
    valid = ["phenix.xtriage", "phenix.predict_and_build", "STOP"]
    context = {
        "has_sequence": True,
        "has_data_mtz": False,
        "has_phased_data_mtz": False,
    }
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert "phenix.xtriage" not in result, (
        "xtriage should be removed\n  Result: %s" % result)
    assert "phenix.predict_and_build" in result, (
        "predict_and_build should remain\n  Result: %s" % result)
    assert "STOP" in result
    print("  PASS")


def test_other_programs_preserved_when_mtriage_removed():
    """Removing mtriage must not affect other programs in valid."""
    print("Test: other_programs_preserved_when_mtriage_removed")
    engine = _make_engine()
    valid = ["phenix.mtriage", "phenix.predict_and_build", "STOP"]
    context = {
        "has_sequence": True,
        "has_full_map": False,
        "has_half_map": False,
    }
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert "phenix.mtriage" not in result
    assert "phenix.predict_and_build" in result
    assert "STOP" in result
    print("  PASS")


# =====================================================================
# SECTION D: Edge cases
# =====================================================================

def test_empty_valid_programs():
    """Empty valid_programs returns empty list."""
    print("Test: empty_valid_programs")
    engine = _make_engine()
    result = engine._filter_programs_missing_data_inputs([], {})
    assert result == [], "Empty input should give empty output, got %s" % result
    print("  PASS")


def test_input_list_not_mutated():
    """The filter returns a NEW list; the input is unchanged.

    Important: callers may have references to the original list and
    shouldn't see surprise mutations.
    """
    print("Test: input_list_not_mutated")
    engine = _make_engine()
    original = ["phenix.xtriage", "STOP"]
    snapshot = list(original)
    context = {"has_data_mtz": False}
    result = engine._filter_programs_missing_data_inputs(original, context)
    assert original == snapshot, (
        "Input list was mutated.\n  Before: %s\n  After: %s"
        % (snapshot, original))
    assert "phenix.xtriage" not in result, "xtriage should be removed in result"
    print("  PASS")


def test_context_with_none_values():
    """Context dict where data flags are explicitly None (not missing).

    .get() returns None for missing keys; if the flag exists as None,
    bool() correctly treats it as falsy.
    """
    print("Test: context_with_none_values")
    engine = _make_engine()
    valid = ["phenix.xtriage", "STOP"]
    context = {
        "has_data_mtz": None,
        "has_phased_data_mtz": None,
    }
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert "phenix.xtriage" not in result, (
        "None-valued data flags should be treated as no-data\n"
        "  Result: %s" % result)
    print("  PASS")


def test_both_xtriage_and_mtriage_filtered_independently():
    """Hypothetical: both xtriage and mtriage present (cross-modality).

    This shouldn't happen in practice (one analyze step per experiment
    type) but the filter handles each independently.
    """
    print("Test: both_xtriage_and_mtriage_filtered_independently")
    engine = _make_engine()
    valid = ["phenix.xtriage", "phenix.mtriage", "STOP"]
    # X-ray data present, no map
    context = {
        "has_data_mtz": True,
        "has_full_map": False,
        "has_half_map": False,
    }
    result = engine._filter_programs_missing_data_inputs(valid, context)
    assert "phenix.xtriage" in result, (
        "xtriage should be kept (data present)\n  Result: %s" % result)
    assert "phenix.mtriage" not in result, (
        "mtriage should be removed (no map)\n  Result: %s" % result)
    print("  PASS")


def test_diag_flag_does_not_change_result():
    """The _diag flag only adds printing; results are identical."""
    print("Test: diag_flag_does_not_change_result")
    engine = _make_engine()
    valid = ["phenix.xtriage", "STOP"]
    context = {"has_data_mtz": False}
    r1 = engine._filter_programs_missing_data_inputs(valid, context, _diag=False)
    r2 = engine._filter_programs_missing_data_inputs(valid, context, _diag=True)
    assert r1 == r2, (
        "_diag should not affect output: %s vs %s" % (r1, r2))
    print("  PASS")


# =====================================================================
# SECTION E: Directive-side defense via _check_program_prerequisites
#
# These tests guard against `_apply_directives` re-adding xtriage/mtriage
# after `_filter_programs_missing_data_inputs` has removed them.  The
# prerequisite check is the directive-side complement to the filter:
# both paths must enforce the same requirement.
# =====================================================================

def test_prereq_check_rejects_xtriage_when_no_data():
    """_check_program_prerequisites returns False for xtriage with no data.

    This ensures `_apply_directives` won't re-add xtriage via the
    after_program or program_settings paths when data is absent.
    """
    print("Test: prereq_check_rejects_xtriage_when_no_data")
    engine = _make_engine()
    context = {
        "has_sequence": True,
        "has_data_mtz": False,
        "has_phased_data_mtz": False,
    }
    result = engine._check_program_prerequisites(
        "phenix.xtriage", context, "analyze")
    assert result is False, (
        "xtriage prereq should fail with no data, got %r" % result)
    print("  PASS")


def test_prereq_check_allows_xtriage_when_data_present():
    """_check_program_prerequisites returns True for xtriage with data."""
    print("Test: prereq_check_allows_xtriage_when_data_present")
    engine = _make_engine()
    context = {"has_data_mtz": True}
    result = engine._check_program_prerequisites(
        "phenix.xtriage", context, "analyze")
    assert result is True, (
        "xtriage prereq should pass with data, got %r" % result)
    print("  PASS")


def test_prereq_check_allows_xtriage_when_phased_data_present():
    """_check_program_prerequisites accepts phased data for xtriage."""
    print("Test: prereq_check_allows_xtriage_when_phased_data_present")
    engine = _make_engine()
    context = {
        "has_data_mtz": False,
        "has_phased_data_mtz": True,
    }
    result = engine._check_program_prerequisites(
        "phenix.xtriage", context, "analyze")
    assert result is True, (
        "xtriage prereq should pass with phased data, got %r" % result)
    print("  PASS")


def test_prereq_check_rejects_mtriage_when_no_map():
    """_check_program_prerequisites returns False for mtriage with no map."""
    print("Test: prereq_check_rejects_mtriage_when_no_map")
    engine = _make_engine()
    context = {
        "has_full_map": False,
        "has_half_map": False,
    }
    result = engine._check_program_prerequisites(
        "phenix.mtriage", context, "analyze")
    assert result is False, (
        "mtriage prereq should fail with no map, got %r" % result)
    print("  PASS")


def test_prereq_check_allows_mtriage_when_full_map_present():
    """_check_program_prerequisites returns True for mtriage with full map."""
    print("Test: prereq_check_allows_mtriage_when_full_map_present")
    engine = _make_engine()
    context = {"has_full_map": True}
    result = engine._check_program_prerequisites(
        "phenix.mtriage", context, "analyze")
    assert result is True, (
        "mtriage prereq should pass with full map, got %r" % result)
    print("  PASS")


def test_prereq_check_allows_mtriage_when_half_map_present():
    """_check_program_prerequisites accepts half map for mtriage."""
    print("Test: prereq_check_allows_mtriage_when_half_map_present")
    engine = _make_engine()
    context = {
        "has_full_map": False,
        "has_half_map": True,
    }
    result = engine._check_program_prerequisites(
        "phenix.mtriage", context, "analyze")
    assert result is True, (
        "mtriage prereq should pass with half map, got %r" % result)
    print("  PASS")


# =====================================================================
# Runner
# =====================================================================

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
