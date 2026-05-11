"""Tests for v116.10 Phase 6b — sequence-only state routing.

The bug: When a user uploads only a sequence (no diffraction data),
`_detect_xray_step` returns step="analyze" because `xtriage_done` is
False.  The state becomes "xray_initial" and valid_programs is
[phenix.xtriage] from YAML.  But xtriage can't run without .mtz data,
so the workflow is stuck.

The v116.10 prediction-only allowance in _check_program_prerequisites
was a workaround: it lets predict_and_build join valid_programs even at
the analyze step.  But that's a state-shape mismatch — the workflow is
nominally at "analyze" but actually running a model-building step.

The fix: _detect_xray_step routes sequence-only sessions directly to
obtain_model, where predict_and_build is the natural program (has:
sequence condition in YAML).  State becomes "xray_analyzed", which
correctly reflects what's about to happen.

These tests verify the routing logic in isolation.
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
    # Minimal step definitions covering everything _detect_xray_step
    # references.  Only the step names matter for the routing tests;
    # the contents are returned via _make_step_result.
    yl.get_workflow_steps = lambda et: {
        "analyze": {"description": "Analyze data"},
        "probe_placement": {"description": "Probe placement"},
        "obtain_model": {"description": "Obtain model"},
        "molecular_replacement": {"description": "MR"},
        "experimental_phasing": {"description": "Experimental phasing"},
        "build_from_phases": {"description": "Build from phases"},
        "refine": {"description": "Refine"},
        "combine_ligand": {"description": "Combine ligand"},
        "validate": {"description": "Validate"},
        "complete": {"description": "Complete", "stop": True},
    }
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


def _make_engine():
    return WorkflowEngine()


def _steps():
    """Return the same step dict the stubbed yaml_loader returns."""
    from libtbx.langchain.knowledge import yaml_loader
    return yaml_loader.get_workflow_steps("xray")


# =====================================================================
# SECTION A: The Phase 6b path — sequence-only sessions route to obtain_model
# =====================================================================

def test_sequence_only_routes_to_obtain_model():
    """has_sequence + no .mtz → obtain_model (not analyze)."""
    print("Test: sequence_only_routes_to_obtain_model")
    engine = _make_engine()
    context = {
        "has_sequence": True,
        "has_data_mtz": False,
        "has_phased_data_mtz": False,
        # xtriage_done absent (False by default)
    }
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] == "obtain_model", (
        "Expected step='obtain_model', got %r\n  Full result: %r"
        % (result.get("step"), result))
    print("  PASS")


def test_sequence_only_state_maps_to_xray_analyzed():
    """Phase 6b's routing maps to the xray_analyzed state name."""
    print("Test: sequence_only_state_maps_to_xray_analyzed")
    engine = _make_engine()
    context = {
        "has_sequence": True,
        "has_data_mtz": False,
    }
    result = engine._detect_xray_step(_steps(), context)
    state = engine._map_step_to_state(result["step"], "xray")
    assert state == "xray_analyzed", (
        "Expected state='xray_analyzed', got %r" % state)
    print("  PASS")


# =====================================================================
# SECTION B: Phase 6b does NOT fire when data is present
# =====================================================================

def test_sequence_with_data_mtz_still_analyzes():
    """has_sequence + has_data_mtz → analyze (existing behavior)."""
    print("Test: sequence_with_data_mtz_still_analyzes")
    engine = _make_engine()
    context = {
        "has_sequence": True,
        "has_data_mtz": True,
        "has_phased_data_mtz": False,
    }
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] == "analyze", (
        "Sequence + data must still go to analyze, got step=%r"
        % result.get("step"))
    print("  PASS")


def test_sequence_with_phased_data_mtz_still_analyzes():
    """has_sequence + has_phased_data_mtz → analyze.

    Phased MTZ files have intensity columns xtriage can analyze, so
    we must NOT skip analyze for them.
    """
    print("Test: sequence_with_phased_data_mtz_still_analyzes")
    engine = _make_engine()
    context = {
        "has_sequence": True,
        "has_data_mtz": False,
        "has_phased_data_mtz": True,
    }
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] == "analyze", (
        "Sequence + phased data must still go to analyze, got step=%r"
        % result.get("step"))
    print("  PASS")


# =====================================================================
# SECTION C: Phase 6b does NOT fire without a sequence
# =====================================================================

def test_no_sequence_no_data_still_analyzes():
    """no sequence + no data → analyze (xtriage filter will then STOP).

    Without a sequence we have nothing to predict, so the appropriate
    behavior is to enter the analyze step and let the Phase 4b filter
    + the empty-valid-programs fallback produce [STOP].
    """
    print("Test: no_sequence_no_data_still_analyzes")
    engine = _make_engine()
    context = {
        "has_sequence": False,
        "has_data_mtz": False,
        "has_phased_data_mtz": False,
    }
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] == "analyze", (
        "No sequence + no data must go to analyze (so the downstream "
        "STOP path fires), got step=%r" % result.get("step"))
    print("  PASS")


def test_sequence_explicitly_false():
    """has_sequence=False (explicitly) routes the same as missing."""
    print("Test: sequence_explicitly_false")
    engine = _make_engine()
    context = {
        "has_sequence": False,
        "has_data_mtz": False,
    }
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] == "analyze"
    print("  PASS")


# =====================================================================
# SECTION D: Phase 6b does NOT fire when xtriage has already run
# =====================================================================

def test_xtriage_done_proceeds_past_phase_6b():
    """xtriage_done=True → Phase 6b doesn't fire (analysis is complete).

    After xtriage runs, the workflow should follow normal post-analyze
    logic, not be re-routed by Phase 6b.
    """
    print("Test: xtriage_done_proceeds_past_phase_6b")
    engine = _make_engine()
    context = {
        "xtriage_done": True,
        "has_sequence": True,
        "has_data_mtz": False,
        # No model, no prediction, no other state → falls through to
        # "Step 2: Need model" → obtain_model
    }
    result = engine._detect_xray_step(_steps(), context)
    # The expected outcome here is obtain_model (Step 2 fall-through),
    # but reached via different reasoning than the Phase 6b path.
    # The "reason" string distinguishes them.
    assert result["step"] == "obtain_model", (
        "After xtriage_done, expected obtain_model (Step 2 path), "
        "got step=%r" % result.get("step"))
    # Make sure it's the Step 2 reason, not the Phase 6b reason
    reason = result.get("reason", "")
    assert "Sequence-only" not in reason, (
        "Should NOT route via Phase 6b after xtriage_done; reason=%r"
        % reason)
    print("  PASS")


# =====================================================================
# SECTION E: Edge cases
# =====================================================================

def test_sequence_only_with_predicted_model_still_routes():
    """Phase 6b fires even when has_predicted_model=True.

    This is deliberate.  Without diffraction data, the user's
    predicted model can't be placed (MR requires .mtz).  Routing to
    obtain_model means the workflow either re-predicts (if
    predict_full_done=False — re-prediction with the same sequence
    yields essentially the same model) or stops gracefully (if
    predict_full_done=True — the YAML `not_done` condition blocks
    re-prediction, valid_programs becomes empty, falls to STOP).
    Both outcomes are better than silently stalling at the analyze
    step the way the pre-Phase-6b code did.
    """
    print("Test: sequence_only_with_predicted_model_still_routes")
    engine = _make_engine()
    context = {
        "has_sequence": True,
        "has_data_mtz": False,
        "has_phased_data_mtz": False,
        "has_predicted_model": True,
        "has_placed_model": False,
    }
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] == "obtain_model", (
        "Phase 6b should fire even with has_predicted_model=True; "
        "got step=%r" % result.get("step"))
    print("  PASS")


def test_context_with_none_for_data_flags():
    """None-valued data flags are correctly treated as missing."""
    print("Test: context_with_none_for_data_flags")
    engine = _make_engine()
    context = {
        "has_sequence": True,
        "has_data_mtz": None,
        "has_phased_data_mtz": None,
    }
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] == "obtain_model", (
        "None-valued data flags should route via Phase 6b, got step=%r"
        % result.get("step"))
    print("  PASS")


def test_result_structure_is_well_formed():
    """The returned step result has the expected keys."""
    print("Test: result_structure_is_well_formed")
    engine = _make_engine()
    context = {"has_sequence": True, "has_data_mtz": False}
    result = engine._detect_xray_step(_steps(), context)
    # Whatever _make_step_result returns, it must have at least "step"
    assert "step" in result, (
        "Result must have 'step' key, got %r" % list(result.keys()))
    # And a reason for diagnostics
    assert "reason" in result, (
        "Result must have 'reason' key, got %r" % list(result.keys()))
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
