"""Tests for `_apply_directives` after_program_done branch under the
v117 stop_refactor architecture (stop_after_requested-gated).

ARCHITECTURE (v117):
  The `elif after_program_done:` branch in `_apply_directives` is now
  gated by `stop_conditions.stop_after_requested`:

    stop_after_requested=True  + after_program_done → wipe to [STOP]
    stop_after_requested=False + after_program_done → no modification
                                                      (plan progresses)

  This supersedes:
    - v112.78 "append STOP" everywhere (replaced by gated wipe)
    - v116.17 validate-step guard (no longer special-cased; the wipe
      only fires when user-explicit stop is requested, and at that
      point wiping at the validate step is what the user asked for)
    - skip_validation carve-out (the flag is sufficient signal now)

  Reference: `plans/stop_refactor_plan.md`, `info2.dat` (May 16 design
  discussion).

WHAT REMAINS UNCHANGED:
  - The v116.17 validate-step guard CODE is still present in the file
    for explicit logging when at validate step + after_program_done
    + validation pending.  The guard's WIPE behavior is replaced by
    the gated logic; the diagnostic branch lives on.
  - Bug 3 retry path (failed refine still retries).
  - Post-ligandfit exemption.
  - Min-run guarantee for not-yet-done after_program.

TEST GROUPS (M1-M7):
  M1: stop_after_requested=True + done → wipe (user-explicit stop fires)
  M2: stop_after_requested=False + done → no wipe (plan progression)
  M3: AF_7mjs preservation — validate step + after_program_done +
      no user-explicit stop → validation programs preserved (plan
      progression interpretation)
  M4: same for X-ray parallel path
  M5: Bug-3-style retry — failed refine still retries even with
      stop_after_requested=True
  M6: Not-yet-done after_program — min-run guarantee applies regardless
      of flag
  M7: context=None doesn't crash
"""

import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
_ROOTS = (
    os.path.abspath(os.path.join(HERE, "..")),
    os.path.abspath(HERE),
)
for _root in _ROOTS:
    if _root not in sys.path:
        sys.path.insert(0, _root)

# Stub libtbx so imports work standalone.
if "libtbx" not in sys.modules:
    import types as _types
    for _mod in (
        "libtbx", "libtbx.langchain",
        "libtbx.langchain.agent",
        "libtbx.langchain.knowledge",
        "libtbx.langchain.knowledge.yaml_loader",
        "libtbx.langchain.knowledge.program_registration",
        "libtbx.langchain.agent.program_registry",
    ):
        sys.modules[_mod] = _types.ModuleType(_mod)

    # Minimal yaml_loader stub
    import yaml as _yaml

    # Find programs.yaml and workflows.yaml relative to this test
    # (typically tests/ is sibling to agent/ and knowledge/)
    _knowledge_dir = None
    for _candidate in (
        os.path.join(HERE, "..", "knowledge"),
        os.path.join(HERE, "..", "..", "knowledge"),
    ):
        if os.path.isdir(_candidate):
            _knowledge_dir = os.path.abspath(_candidate)
            break

    # If knowledge YAMLs aren't found, provide hard-coded minimums
    # so the test still runs.  The guard logic does not require
    # workflow_engine to consult YAML for these tests.
    if _knowledge_dir and os.path.isfile(
            os.path.join(_knowledge_dir, "programs.yaml")):
        with open(os.path.join(_knowledge_dir, "programs.yaml")) as _f:
            _PROGRAMS = _yaml.safe_load(_f)
    else:
        _PROGRAMS = {}
    if _knowledge_dir and os.path.isfile(
            os.path.join(_knowledge_dir, "workflows.yaml")):
        with open(os.path.join(_knowledge_dir, "workflows.yaml")) as _f:
            _WORKFLOWS = _yaml.safe_load(_f)
    else:
        _WORKFLOWS = {}

    yl = sys.modules["libtbx.langchain.knowledge.yaml_loader"]
    yl.get_program = lambda n: _PROGRAMS.get(n)
    yl.get_all_programs = lambda: _PROGRAMS

    def _get_metric_threshold(name, level=None, resolution=None):
        thresholds = {
            "r_free": {"acceptable": 0.20},
            "map_cc": {"acceptable": 0.70},
        }
        if name in thresholds and level in thresholds[name]:
            return thresholds[name][level]
        return None
    yl.get_metric_threshold = _get_metric_threshold
    yl.get_metric_direction = lambda n: "lower" if n == "r_free" else "higher"
    yl.is_metric_good = lambda *a, **k: False
    yl.is_metric_acceptable = lambda *a, **k: True
    yl.load_metrics = lambda: {}
    yl.get_workflow_steps = lambda et: _WORKFLOWS.get(et, {}).get("steps", {})
    yl.get_workflow_targets = lambda et: _WORKFLOWS.get(et, {}).get("targets", {})
    yl.load_workflows = lambda: _WORKFLOWS

    pr = sys.modules["libtbx.langchain.knowledge.program_registration"]
    pr.get_program_done_flag_map = lambda: {
        name: (info.get("done_tracking") or {}).get("flag")
        for name, info in _PROGRAMS.items()
        if isinstance(info, dict)
        and (info.get("done_tracking") or {}).get("flag")
    }

    preg = sys.modules["libtbx.langchain.agent.program_registry"]
    class _PR:
        def __init__(self): pass
        def get_user_advice_keywords(self, prog): return []
    preg.ProgramRegistry = _PR


from agent.workflow_engine import WorkflowEngine


def _af7mjs_context():
    """AF_7mjs cryo-EM scenario: RSR has run once successfully,
    validation pending."""
    return {
        "rsr_count": 1,
        "refine_count": 0,
        "validation_done": False,
        "map_cc": 0.786,
        "has_refined_model": True,
        "has_placed_model": True,
        "experiment_type": "cryoem",
        "last_program": "phenix.real_space_refine",
        # Bug 3 override needs this to recognize successful completion.
        # Without it, after_program_done flips to False for retry.
        "successful_programs": {"phenix.real_space_refine"},
    }


def _xray_context():
    """X-ray scenario: refine has run successfully, validation pending."""
    return {
        "refine_count": 1,
        "rsr_count": 0,
        "validation_done": False,
        "r_free": 0.24,
        "has_refined_model": True,
        "has_placed_model": True,
        "experiment_type": "xray",
        "last_program": "phenix.refine",
        "successful_programs": {"phenix.refine"},
    }


# =============================================================================
# M1: stop_after_requested=True + after_program_done → wipe
# =============================================================================

def test_M1a_user_explicit_stop_validate_step_wipes_to_stop():
    """When user explicitly requested stop (stop_after_requested=True)
    and the target program has run, valid_programs is wiped to [STOP]
    even at the validate step.  (Option (a) per info2.dat: hard stop,
    no validation.)"""
    print("[test] M1a: user-explicit stop + validate step -> wipe to [STOP]")
    engine = WorkflowEngine()
    directives = {"stop_conditions": {
        "after_program": "phenix.real_space_refine",
        "stop_after_requested": True,
    }}
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem",
                "phenix.map_correlations"]
    result = engine._apply_directives(
        valid_in, directives, "validate", _af7mjs_context(),
        experiment_type="cryoem")
    assert result == ["STOP"], (
        "User-explicit stop should wipe to [STOP] even at validate "
        "step, got: %r" % result)
    print("  PASS")


def test_M1b_user_explicit_stop_refine_step_wipes_to_stop():
    """Same: user-explicit stop wipes at refine step too."""
    print("[test] M1b: user-explicit stop + refine step -> wipe to [STOP]")
    engine = WorkflowEngine()
    directives = {"stop_conditions": {
        "after_program": "phenix.real_space_refine",
        "stop_after_requested": True,
    }}
    valid_in = ["phenix.real_space_refine", "phenix.dock_in_map"]
    result = engine._apply_directives(
        valid_in, directives, "refine", _af7mjs_context(),
        experiment_type="cryoem")
    assert result == ["STOP"], (
        "User-explicit stop should wipe to [STOP] at any step, got: %r"
        % result)
    print("  PASS")


def test_M1c_user_explicit_stop_complete_step_wipes_to_stop():
    """Same at the complete step."""
    print("[test] M1c: user-explicit stop + complete step -> wipe to [STOP]")
    engine = WorkflowEngine()
    directives = {"stop_conditions": {
        "after_program": "phenix.real_space_refine",
        "stop_after_requested": True,
    }}
    valid_in = ["phenix.real_space_refine"]
    result = engine._apply_directives(
        valid_in, directives, "complete", _af7mjs_context(),
        experiment_type="cryoem")
    assert result == ["STOP"], (
        "User-explicit stop at complete step -> [STOP], got: %r" % result)
    print("  PASS")


# =============================================================================
# M2: stop_after_requested=False/absent + after_program_done -> no wipe
# =============================================================================

def test_M2a_plan_progression_validate_step_no_wipe():
    """When stop_after_requested is False/absent, after_program is a
    plan-progression hint.  At the validate step with after_program
    done, the validation programs must remain available."""
    print("[test] M2a: plan progression + validate step -> no wipe")
    engine = WorkflowEngine()
    directives = {"stop_conditions": {
        "after_program": "phenix.real_space_refine",
    }}
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem",
                "phenix.map_correlations"]
    result = engine._apply_directives(
        valid_in, directives, "validate", _af7mjs_context(),
        experiment_type="cryoem")
    for p in ("phenix.molprobity", "phenix.validation_cryoem",
              "phenix.map_correlations"):
        assert p in result, (
            "%s should remain in result for plan-progression case, "
            "got: %r" % (p, result))
    assert result != ["STOP"], (
        "Plan progression should not wipe to [STOP], got: %r" % result)
    print("  PASS")


def test_M2b_plan_progression_flag_false_no_wipe():
    """Same with stop_after_requested=False explicitly set."""
    print("[test] M2b: stop_after_requested=False + validate -> no wipe")
    engine = WorkflowEngine()
    directives = {"stop_conditions": {
        "after_program": "phenix.real_space_refine",
        "stop_after_requested": False,
    }}
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem"]
    result = engine._apply_directives(
        valid_in, directives, "validate", _af7mjs_context(),
        experiment_type="cryoem")
    assert "phenix.molprobity" in result
    print("  PASS")


# =============================================================================
# M3: AF_7mjs - the original cryo-EM regression
# =============================================================================

def test_M3_af7mjs_cryoem_validate_preserves_validation():
    """AF_7mjs scenario: cryo-EM tutorial, the extractor (mis)set
    after_program=phenix.real_space_refine from preprocessed prose.
    No user-explicit stop was requested.  At the validate step with
    RSR done, validation programs must remain available."""
    print("[test] M3: AF_7mjs cryo-EM validate preserves validation")
    engine = WorkflowEngine()
    directives = {"stop_conditions": {
        "after_program": "phenix.real_space_refine",
    }}
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem",
                "phenix.map_correlations"]
    result = engine._apply_directives(
        valid_in, directives, "validate", _af7mjs_context(),
        experiment_type="cryoem")
    assert "phenix.molprobity" in result
    assert "phenix.validation_cryoem" in result
    assert "phenix.map_correlations" in result
    assert result != ["STOP"]
    print("  PASS - result: %r" % result)


# =============================================================================
# M4: X-ray parallel
# =============================================================================

def test_M4_xray_validate_preserves_validation():
    """X-ray parallel to M3."""
    print("[test] M4: X-ray validate step preserves validation")
    engine = WorkflowEngine()
    directives = {"stop_conditions": {
        "after_program": "phenix.refine",
    }}
    valid_in = ["phenix.molprobity", "phenix.model_vs_data",
                "phenix.map_correlations"]
    result = engine._apply_directives(
        valid_in, directives, "validate", _xray_context(),
        experiment_type="xray")
    assert "phenix.molprobity" in result
    assert "phenix.model_vs_data" in result
    assert result != ["STOP"]
    print("  PASS - result: %r" % result)


# =============================================================================
# M5: Bug-3 retry path under stop_after_requested=True
# =============================================================================

def test_M5_bug3_retry_still_works_with_stop_after_requested():
    """When the target program failed and is set up to retry, the
    Bug 3 override fires BEFORE the wipe - the failed program must
    still be available for retry even with stop_after_requested=True.

    Property called out in info2.dat: failed refine still retries
    because the Bug 3 override fires inside after_program_done
    calculation, so after_program_done is False and the wipe doesn't
    apply.
    """
    print("[test] M5: Bug-3 retry works with stop_after_requested=True")
    engine = WorkflowEngine()
    directives = {"stop_conditions": {
        "after_program": "phenix.real_space_refine",
        "stop_after_requested": True,
    }}
    context = {
        "rsr_count": 0,
        "refine_count": 0,
        "validation_done": False,
        "map_cc": 0.0,
        "has_refined_model": False,
        "has_placed_model": True,
        "experiment_type": "cryoem",
        "last_program": None,
    }
    valid_in = ["phenix.real_space_refine", "phenix.molprobity"]
    result = engine._apply_directives(
        valid_in, directives, "refine", context, experiment_type="cryoem")
    assert "phenix.real_space_refine" in result, (
        "Target program must remain available when not yet completed, "
        "got: %r" % result)
    assert result != ["STOP"], (
        "Should not wipe when after_program not yet done, got: %r" % result)
    print("  PASS - result: %r" % result)


# =============================================================================
# M6: Min-run guarantee unchanged
# =============================================================================

def test_M6a_min_run_target_preserved_when_not_done():
    """When after_program is set but not yet done, the target program
    is preserved at the front of valid_programs regardless of the
    stop_after_requested flag."""
    print("[test] M6a: min-run guarantee preserves target when not done")
    engine = WorkflowEngine()
    for flag_state in [True, False, None]:
        directives = {"stop_conditions": {
            "after_program": "phenix.refine",
        }}
        if flag_state is not None:
            directives["stop_conditions"]["stop_after_requested"] = flag_state
        context = {
            "refine_count": 0,
            "rsr_count": 0,
            "validation_done": False,
            "r_free": None,
            "has_refined_model": False,
            "has_placed_model": True,
            "experiment_type": "xray",
            "last_program": None,
        }
        valid_in = ["phenix.refine", "phenix.molprobity"]
        result = engine._apply_directives(
            valid_in, directives, "refine", context,
            experiment_type="xray")
        assert "phenix.refine" in result, (
            "Target program missing when not yet done (flag=%r): %r"
            % (flag_state, result))
        assert result != ["STOP"], (
            "Wipe should not fire when target not yet done (flag=%r): %r"
            % (flag_state, result))
    print("  PASS")


def test_M6b_no_after_program_no_change():
    """When after_program is not set at all, _apply_directives makes
    no after_program-related modifications."""
    print("[test] M6b: no after_program -> no after_program-related change")
    engine = WorkflowEngine()
    directives = {"stop_conditions": {}}
    context = {
        "refine_count": 1,
        "rsr_count": 0,
        "validation_done": False,
        "r_free": 0.24,
        "has_refined_model": True,
        "has_placed_model": True,
        "experiment_type": "xray",
        "last_program": "phenix.refine",
    }
    valid_in = ["phenix.molprobity", "phenix.model_vs_data"]
    result = engine._apply_directives(
        valid_in, directives, "validate", context, experiment_type="xray")
    assert result != ["STOP"]
    assert "phenix.molprobity" in result
    print("  PASS")


# =============================================================================
# M7: context=None safety
# =============================================================================

def test_M7_context_none_does_not_crash():
    """_apply_directives must handle context=None gracefully."""
    print("[test] M7: context=None does not crash")
    engine = WorkflowEngine()
    directives = {"stop_conditions": {
        "after_program": "phenix.refine",
        "stop_after_requested": True,
    }}
    valid_in = ["phenix.refine", "phenix.molprobity"]
    try:
        result = engine._apply_directives(
            valid_in, directives, "refine", None, experiment_type="xray")
    except Exception as e:
        raise AssertionError(
            "_apply_directives raised with context=None: %s" % e)
    assert isinstance(result, list)
    print("  PASS - no crash, result: %r" % result)


def run_all_tests():
    tests = [
        test_M1a_user_explicit_stop_validate_step_wipes_to_stop,
        test_M1b_user_explicit_stop_refine_step_wipes_to_stop,
        test_M1c_user_explicit_stop_complete_step_wipes_to_stop,
        test_M2a_plan_progression_validate_step_no_wipe,
        test_M2b_plan_progression_flag_false_no_wipe,
        test_M3_af7mjs_cryoem_validate_preserves_validation,
        test_M4_xray_validate_preserves_validation,
        test_M5_bug3_retry_still_works_with_stop_after_requested,
        test_M6a_min_run_target_preserved_when_not_done,
        test_M6b_no_after_program_no_change,
        test_M7_context_none_does_not_crash,
    ]
    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            passed += 1
        except Exception as e:
            failed += 1
            print("  FAIL: %s" % e)
    print()
    print("%d passed, %d failed" % (passed, failed))
    return failed == 0


if __name__ == "__main__":
    sys.exit(0 if run_all_tests() else 1)
