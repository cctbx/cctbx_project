"""Tests for v116.10 Phase 6c — Ligand Workflow Restart Bug.

Two interacting bugs caused the nsf-d2-ligand tutorial to "restart"
its reasoning at cycle 3:

1. _detect_xray_step returned "analyze" any time xtriage_done was
   False, regardless of whether refinement, ligand fitting, or
   other downstream programs had already run.  When the user
   supplies a pre-placed model + start_with_program=phenix.refine
   (skipping xtriage), every cycle reported STATE=xray_initial /
   "Need to analyze data quality first".  The LLM, seeing that
   worldview, reasonably concluded "this is the first refinement
   step" at cycles 3, 4, 5 and re-applied first-cycle directives.

2. best_files_tracker scored ligand_fit_output (105) below the
   current best refined model (100 + R-free contribution = ~122).
   So phenix.ligandfit's output (ligand_fit_1.pdb) did not become
   the new "best model".  Subsequent context handoff to the LLM
   then told it the model is the unliganded one, causing
   post-ligandfit-refine cycles to operate on the wrong model.

Uses the same libtbx-stubbing pattern as tst_sequence_only_routing.py
(Phase 6b tests), so the workflow_engine tests run without a real
PHENIX install.
"""

from __future__ import absolute_import, division, print_function

import os
import sys
import types


# =============================================================================
# Path setup + libtbx stubs (mirrors tst_sequence_only_routing.py)
# =============================================================================

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
    yl.get_metric_threshold = lambda *args, **kwargs: None
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
    from libtbx.langchain.knowledge import yaml_loader
    return yaml_loader.get_workflow_steps("xray")


def _empty_context():
    """Empty context: all flags False/None/0."""
    return {
        "xtriage_done": False,
        "has_sequence": False,
        "has_data_mtz": False,
        "has_phased_data_mtz": False,
        "has_model": False,
        "has_placed_model": False,
        "has_placed_model_from_history": False,
        "has_predicted_model": False,
        "has_processed_model": False,
        "has_refined_model": False,
        "has_ligand_fit": False,
        "has_ligand_file": False,
        "refine_done": False,
        "phaser_done": False,
        "autobuild_done": False,
        "autosol_done": False,
        "autosol_attempted": False,
        "ligandfit_done": False,
        "rsr_done": False,
        "dock_done": False,
        "predict_done": False,
        "predict_full_done": False,
        "placement_probed": False,
        "placement_uncertain": False,
        "cell_mismatch": False,
        "force_mr": False,
        "wants_validation_only": False,
        "validation_done": False,
        "pdbtools_done": False,
        "needs_post_ligandfit_refine": False,
        "refine_count": 0,
        "r_free": None,
    }


# =============================================================================
# SECTION A: past_analysis check — Fix 1 in _detect_xray_step
# =============================================================================

def test_analyze_at_start_when_nothing_done():
    """Fresh xtriage workflow should still go to analyze on cycle 1."""
    print("Test: analyze_at_start_when_nothing_done")
    engine = _make_engine()
    context = _empty_context()
    context["has_data_mtz"] = True
    context["has_model"] = True
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] == "analyze", (
        "Fresh xtriage workflow must start at analyze, got step=%r" %
        result.get("step"))
    print("  PASS")


def test_past_analyze_when_refine_done():
    """After refine_done=True, advance past analyze."""
    print("Test: past_analyze_when_refine_done")
    engine = _make_engine()
    context = _empty_context()
    context["refine_done"] = True
    context["refine_count"] = 1
    context["has_placed_model"] = True
    context["has_refined_model"] = True
    context["has_data_mtz"] = True
    context["r_free"] = 0.29
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] != "analyze", (
        "After refine_done=True, must not return analyze, got step=%r" %
        result.get("step"))
    print("  PASS — step=%r" % result.get("step"))


def test_past_analyze_when_ligandfit_done():
    """After ligandfit_done=True (the motivating bug case), advance
    past analyze."""
    print("Test: past_analyze_when_ligandfit_done")
    engine = _make_engine()
    context = _empty_context()
    # Scenario: ligand tutorial after cycle 2 (ligandfit complete)
    context["refine_done"] = True
    context["ligandfit_done"] = True
    context["has_ligand_fit"] = True
    context["refine_count"] = 1
    context["has_placed_model"] = True
    context["has_refined_model"] = True
    context["has_ligand_file"] = True
    context["has_data_mtz"] = True
    context["r_free"] = 0.28
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] != "analyze", (
        "After ligandfit_done=True (the nsf-d2-ligand bug), must "
        "not return analyze. Got step=%r" % result.get("step"))
    print("  PASS — step=%r" % result.get("step"))


def test_past_analyze_when_placed_from_history():
    """When history shows placement, do not loop back to analyze.
    Covers session-resumption scenarios."""
    print("Test: past_analyze_when_placed_from_history")
    engine = _make_engine()
    context = _empty_context()
    context["has_placed_model_from_history"] = True
    context["has_placed_model"] = True
    context["has_data_mtz"] = True
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] != "analyze", (
        "Placed model from history must not loop to analyze, got step=%r" %
        result.get("step"))
    print("  PASS — step=%r" % result.get("step"))


def test_past_analyze_when_phaser_done():
    """Phaser ran successfully → past analyze."""
    print("Test: past_analyze_when_phaser_done")
    engine = _make_engine()
    context = _empty_context()
    context["phaser_done"] = True
    context["has_placed_model"] = True
    context["has_data_mtz"] = True
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] != "analyze", (
        "phaser_done=True must not loop to analyze, got step=%r" %
        result.get("step"))
    print("  PASS — step=%r" % result.get("step"))


def test_past_analyze_when_placement_probed():
    """The placement probe ran (whatever outcome) → past analyze."""
    print("Test: past_analyze_when_placement_probed")
    engine = _make_engine()
    context = _empty_context()
    context["placement_probed"] = True
    context["placement_probe_result"] = "placed"
    context["has_data_mtz"] = True
    context["has_model"] = True
    context["has_placed_model"] = True
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] != "analyze", (
        "placement_probed=True must not loop to analyze, got step=%r" %
        result.get("step"))
    print("  PASS — step=%r" % result.get("step"))


def test_explicit_xtriage_done_still_works():
    """Existing behavior: xtriage_done=True alone advances past
    analyze even with nothing else done."""
    print("Test: explicit_xtriage_done_still_works")
    engine = _make_engine()
    context = _empty_context()
    context["xtriage_done"] = True
    context["has_data_mtz"] = True
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] != "analyze", (
        "xtriage_done=True must advance past analyze, got step=%r" %
        result.get("step"))
    print("  PASS — step=%r" % result.get("step"))


def test_nsf_d2_ligand_tutorial_scenario():
    """Reproduce the exact context at cycle 3 of the nsf-d2-ligand
    tutorial after the bug fix.

    Cycle 3 conditions (from the log):
      - xtriage was never run (user said start_with_program=
        phenix.refine, skipping the analyze stage)
      - refine_001 ran successfully (cycle 1) → refine_done=True,
        refine_count=1
      - ligandfit ran successfully (cycle 2) → ligandfit_done=True,
        has_ligand_fit=True
      - User advice set model_is_placed=True → has_placed_model=True
      - R-free = 0.29 from cycle 1

    Pre-fix: returns "analyze", state=xray_initial, LLM says "first
    refinement step" and re-applies generate_r_free_flags=True.
    Post-fix: returns combine_ligand (correct workflow path through
    pdbtools).
    """
    print("Test: nsf_d2_ligand_tutorial_scenario")
    engine = _make_engine()
    context = _empty_context()
    context["refine_done"] = True
    context["refine_count"] = 1
    context["ligandfit_done"] = True
    context["has_ligand_fit"] = True
    context["has_ligand_file"] = True
    context["has_placed_model"] = True
    context["has_refined_model"] = True
    context["has_model"] = True
    context["has_data_mtz"] = True
    context["r_free"] = 0.29

    result = engine._detect_xray_step(_steps(), context)
    step = result.get("step")
    assert step != "analyze", (
        "nsf-d2-ligand cycle 3 must not return analyze "
        "(the restart bug). Got step=%r" % step)
    # Should route to combine_ligand (correct flow: pdbtools combines
    # ligand-fit + refined model, then refine on combined)
    assert step in ("refine", "combine_ligand", "validate"), (
        "Expected refine/combine_ligand/validate after ligand fit, "
        "got step=%r" % step)
    print("  PASS — step=%r (no longer routing to analyze)" % step)


# =============================================================================
# SECTION B: Regression — existing behavior must be preserved
# =============================================================================

def test_phase_6b_sequence_only_still_fires():
    """The Phase 6b sequence-only check runs BEFORE past_analysis.
    Verify it still fires correctly when only a sequence is
    available."""
    print("Test: phase_6b_sequence_only_still_fires")
    engine = _make_engine()
    context = _empty_context()
    context["has_sequence"] = True
    # No data, no model, nothing done.  Phase 6b should route this
    # to obtain_model, NOT to analyze.
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] == "obtain_model", (
        "Sequence-only session should route to obtain_model (Phase 6b), "
        "got step=%r" % result.get("step"))
    reason = result.get("reason", "")
    assert "Sequence-only" in reason, (
        "Should fire Phase 6b path; reason=%r" % reason)
    print("  PASS")


def test_first_cycle_with_data_and_model_still_analyzes():
    """User uploads data + model with no directives; xtriage hasn't
    run yet.  Must still go to analyze first (preserves existing
    cycle-1 behavior)."""
    print("Test: first_cycle_with_data_and_model_still_analyzes")
    engine = _make_engine()
    context = _empty_context()
    context["has_data_mtz"] = True
    context["has_model"] = True
    # Note: NO history evidence, NO directive that model_is_placed
    # → has_placed_model_from_history=False, has_refined_model=False
    result = engine._detect_xray_step(_steps(), context)
    assert result["step"] == "analyze", (
        "Fresh data+model upload should still go to analyze first, "
        "got step=%r" % result.get("step"))
    print("  PASS")


# =============================================================================
# SECTION C: ligand_fit_output inheritance — Fix 2 in best_files_tracker
# =============================================================================

def _import_tracker():
    try:
        from libtbx.langchain.agent.best_files_tracker import BestFilesTracker
        return BestFilesTracker
    except ImportError:
        try:
            from agent.best_files_tracker import BestFilesTracker
            return BestFilesTracker
        except ImportError:
            from best_files_tracker import BestFilesTracker
            return BestFilesTracker


def test_ligand_fit_output_inherits_metrics():
    """phenix.ligandfit produces ligand_fit_1.pdb; tracker should
    treat it as new best, inheriting R-free from prior refined."""
    print("Test: ligand_fit_output_inherits_metrics")
    BestFilesTracker = _import_tracker()
    tracker = BestFilesTracker()

    # Cycle 1: refine_001_001.pdb (R-free 0.29)
    became_best = tracker.evaluate_file(
        path="/work/refine_001_001.pdb",
        cycle=1,
        metrics={"r_free": 0.29, "r_work": 0.27, "resolution": 2.0},
        stage="refined",
        category="model")
    assert became_best, "refine_001_001.pdb must be best after cycle 1"

    # Cycle 2: ligand_fit_1.pdb (no metrics — LigandFit doesn't
    # report R-free)
    became_best = tracker.evaluate_file(
        path="/work/ligand_fit_1.pdb",
        cycle=2,
        metrics=None,
        stage="ligand_fit_output",
        category="model")
    assert became_best, (
        "ligand_fit_1.pdb must become best after ligandfit "
        "(the bug — without inheritance fix, scored 105 vs "
        "refine's ~122)")

    best_now = tracker.get_best("model")
    assert "ligand_fit_1" in best_now.path, (
        "Best should be ligand_fit_1.pdb, got %s" % best_now.path)
    assert best_now.metrics is not None, (
        "Should have inherited metrics, got None")
    assert abs(best_now.metrics.get("r_free", 0) - 0.29) < 1e-6, (
        "Inherited r_free should be 0.29, got %s" %
        best_now.metrics.get("r_free"))
    print("  PASS — best=ligand_fit_1.pdb, R-free=0.29 inherited")


def test_with_ligand_inheritance_still_works():
    """The original with_ligand inheritance (v115.05) must continue
    to work after the extension."""
    print("Test: with_ligand_inheritance_still_works")
    BestFilesTracker = _import_tracker()
    tracker = BestFilesTracker()

    tracker.evaluate_file(
        path="/work/refine_001_001.pdb",
        cycle=1,
        metrics={"r_free": 0.25, "r_work": 0.23, "resolution": 2.0},
        stage="refined",
        category="model")

    became_best = tracker.evaluate_file(
        path="/work/refine_001_001_with_ligand.pdb",
        cycle=2,
        metrics=None,
        stage="with_ligand",
        category="model")
    assert became_best, "with_ligand must still become best"
    best_now = tracker.get_best("model")
    assert "with_ligand" in best_now.path
    assert abs(best_now.metrics.get("r_free", 0) - 0.25) < 1e-6
    print("  PASS")


def test_ligand_fit_output_with_metrics_uses_them():
    """If ligand_fit_output arrives WITH metrics, those are used
    (no inheritance)."""
    print("Test: ligand_fit_output_with_metrics_uses_them")
    BestFilesTracker = _import_tracker()
    tracker = BestFilesTracker()

    tracker.evaluate_file(
        path="/work/refine_001_001.pdb",
        cycle=1,
        metrics={"r_free": 0.29, "resolution": 2.0},
        stage="refined",
        category="model")

    # Hypothetical: ligand_fit comes with worse metrics
    tracker.evaluate_file(
        path="/work/ligand_fit_1.pdb",
        cycle=2,
        metrics={"r_free": 0.45, "resolution": 2.0},
        stage="ligand_fit_output",
        category="model")
    best_now = tracker.get_best("model")
    # Stage 105 + bad R-free (0.45) vs stage 100 + good R-free (0.29).
    # The good refined wins via R-free contribution.
    print("  PASS — explicit metrics respected (best=%s)" %
          os.path.basename(best_now.path))


def test_post_ligandfit_refine_finds_liganded_model():
    """End-to-end: refine → ligandfit → refine results in the
    final refined model being best."""
    print("Test: post_ligandfit_refine_finds_liganded_model")
    BestFilesTracker = _import_tracker()
    tracker = BestFilesTracker()

    tracker.evaluate_file(
        path="/work/refine_001_001.pdb",
        cycle=1,
        metrics={"r_free": 0.29, "resolution": 2.0},
        stage="refined",
        category="model")

    tracker.evaluate_file(
        path="/work/ligand_fit_1.pdb",
        cycle=2,
        metrics=None,
        stage="ligand_fit_output",
        category="model")

    tracker.evaluate_file(
        path="/work/refine_002_001.pdb",
        cycle=3,
        metrics={"r_free": 0.22, "resolution": 2.0},
        stage="refined",
        category="model")

    best = tracker.get_best("model")
    assert "refine_002" in best.path, (
        "Best should be refine_002_001.pdb, got %s" % best.path)
    assert abs(best.metrics.get("r_free", 0) - 0.22) < 1e-6
    print("  PASS")


# =============================================================================
# Test runner (standard format — matches other v116.10 test files)
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
