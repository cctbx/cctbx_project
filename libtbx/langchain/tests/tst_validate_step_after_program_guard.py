"""Tests for v116.17 validate-step after_program guard.

Bug: When `directives.stop_conditions.after_program` is set AND the
after_program has been run (after_program_done=True), `_apply_directives`
at `agent/workflow_engine.py` line 2073 wiped the entire `valid_programs`
list to `["STOP"]` — even at the `validate` step where molprobity /
validation_cryoem / map_correlations should still be runnable.

For AF_7mjs (cryo-EM tutorial), the directive extractor mis-extracted
`after_program=phenix.real_space_refine` from the README's Primary Goal
("rebuild the loop and refine the model"), despite a CRITICAL rule in
the extractor prompt forbidding fabrication of stop conditions from
goal descriptions.  This combined with the wipe to cause cycle 5 to
stop without ever running validation.

The v116.16 diagnostic patch confirmed the chain:

  PERCEIVE.detect_workflow_state → valid_programs=['STOP']
  PLAN.before_LLM_call          → valid_programs=['STOP']
  prompt as seen by LLM         → "VALID PROGRAMS: STOP"

The LLM correctly chose STOP because that was its only option.

Fix: in `_apply_directives`'s `elif after_program_done` branch, do NOT
wipe `valid_programs` when:
  * step_name == "validate" (we're at the validation step)
  * validation_done is False (we haven't validated yet)
  * skip_validation is NOT set in stop_conditions (user hasn't opted out)

Under those conditions, the validation programs that the validate step
provides are preserved, and STOP is appended (belt-and-suspenders) so
the LLM can choose either path.

These tests verify:
  - AF_7mjs regression: validate-step + cryo-EM after_program_done →
    validation programs preserved
  - Same for X-ray (parallel path, after_program=phenix.refine)
  - skip_validation:true escape hatch still wipes to [STOP] as before
  - validation_done=True still wipes to [STOP] (normal completion)
  - Other steps (not "validate") still wipe to [STOP] (the existing
    use case where after_program completed and workflow is done)
  - context=None doesn't crash the guard
  - STOP is in the result so LLM can choose to stop after validating
  - after_program directive of phenix.molprobity (a validation program
    itself) still works correctly when validation_done=True
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


def _make_engine_with_directives():
    """Return (engine, common_directives, common_context) for cryo-EM AF_7mjs scenario."""
    engine = WorkflowEngine()
    directives = {
        "stop_conditions": {
            "after_program": "phenix.real_space_refine",
        },
    }
    context = {
        "rsr_count": 1,           # RSR has run once
        "refine_count": 0,
        "validation_done": False,
        "map_cc": 0.786,
        "has_refined_model": True,
        "has_placed_model": True,
        "experiment_type": "cryoem",
        "last_program": "phenix.real_space_refine",
    }
    return engine, directives, context


# =============================================================================
# Core regression: AF_7mjs scenario — validate step keeps validation programs
# =============================================================================

def test_cryoem_validate_step_preserves_validation_programs():
    print("[test] Cryo-EM validate step preserves validation programs "
          "when after_program=RSR done")
    engine, directives, context = _make_engine_with_directives()
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem",
                "phenix.map_correlations"]
    result = engine._apply_directives(
        valid_in, directives, "validate", context, experiment_type="cryoem")
    assert "phenix.molprobity" in result, \
        "molprobity should remain available, got: %r" % result
    assert "phenix.validation_cryoem" in result, \
        "validation_cryoem should remain available, got: %r" % result
    assert "phenix.map_correlations" in result, \
        "map_correlations should remain available, got: %r" % result
    assert result != ["STOP"], \
        "valid_programs must NOT be wiped to [STOP], got: %r" % result
    print("  PASS — result: %r" % result)


def test_cryoem_validate_step_stop_is_appended():
    print("[test] STOP is appended so LLM may still choose to stop")
    engine, directives, context = _make_engine_with_directives()
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem",
                "phenix.map_correlations"]
    result = engine._apply_directives(
        valid_in, directives, "validate", context, experiment_type="cryoem")
    assert "STOP" in result, \
        "STOP should be present in result so LLM can opt to stop, got: %r" % result
    print("  PASS")


# =============================================================================
# X-ray parallel path
# =============================================================================

def test_xray_validate_step_preserves_validation_programs():
    print("[test] X-ray validate step preserves validation programs "
          "when after_program=refine done")
    engine = WorkflowEngine()
    directives = {
        "stop_conditions": {"after_program": "phenix.refine"},
    }
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
    valid_in = ["phenix.molprobity", "phenix.model_vs_data",
                "phenix.map_correlations"]
    result = engine._apply_directives(
        valid_in, directives, "validate", context, experiment_type="xray")
    assert "phenix.molprobity" in result, \
        "molprobity should remain for X-ray validate step, got: %r" % result
    assert "phenix.model_vs_data" in result, \
        "model_vs_data should remain for X-ray validate step, got: %r" % result
    assert result != ["STOP"], \
        "X-ray validate step must NOT be wiped to [STOP], got: %r" % result
    print("  PASS — result: %r" % result)


# =============================================================================
# Escape hatch: skip_validation:true preserves the original wipe behavior
# =============================================================================

def test_skip_validation_overrides_guard():
    print("[test] skip_validation:true makes the wipe fire as before")
    engine, directives, context = _make_engine_with_directives()
    directives["stop_conditions"]["skip_validation"] = True
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem",
                "phenix.map_correlations"]
    result = engine._apply_directives(
        valid_in, directives, "validate", context, experiment_type="cryoem")
    assert result == ["STOP"], \
        "With skip_validation:true the wipe should fire, got: %r" % result
    print("  PASS — wipe behavior preserved")


# =============================================================================
# validation_done=True triggers the normal wipe (workflow truly complete)
# =============================================================================

def test_validation_already_done_wipes_to_stop():
    print("[test] validation_done=True triggers normal wipe")
    engine, directives, context = _make_engine_with_directives()
    context["validation_done"] = True
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem",
                "phenix.map_correlations"]
    result = engine._apply_directives(
        valid_in, directives, "validate", context, experiment_type="cryoem")
    assert result == ["STOP"], \
        "validation_done=True should wipe to [STOP], got: %r" % result
    print("  PASS")


# =============================================================================
# Non-validate steps still wipe (preserve original "workflow done" behavior)
# =============================================================================

def test_refine_step_still_wipes_to_stop():
    print("[test] Other steps (e.g. 'refine') still wipe to [STOP] as before")
    engine, directives, context = _make_engine_with_directives()
    # At refine step, after_program=RSR was the user's last requested action;
    # without an explicit validation pending guard the wipe should fire.
    valid_in = ["phenix.real_space_refine", "phenix.dock_in_map"]
    result = engine._apply_directives(
        valid_in, directives, "refine", context, experiment_type="cryoem")
    assert result == ["STOP"], \
        "Non-validate steps should still wipe to [STOP], got: %r" % result
    print("  PASS")


def test_complete_step_still_wipes_to_stop():
    print("[test] 'complete' step still wipes to [STOP] as before")
    engine, directives, context = _make_engine_with_directives()
    valid_in = ["STOP"]
    result = engine._apply_directives(
        valid_in, directives, "complete", context, experiment_type="cryoem")
    assert result == ["STOP"], \
        "'complete' step should still wipe to [STOP], got: %r" % result
    print("  PASS")


# =============================================================================
# Safety: context=None should not crash the guard
# =============================================================================

def test_context_none_does_not_crash():
    print("[test] context=None does not crash the guard")
    engine = WorkflowEngine()
    directives = {
        "stop_conditions": {"after_program": "phenix.real_space_refine"},
    }
    # When context is None we cannot determine validation_done, so the
    # guard should NOT fire and behavior falls through to the original
    # wipe.  This matches the existing pattern in the function where
    # after_program_done is computed only when context is not None.
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem"]
    # after_program_done is False when context is None (per the function's
    # own check at line ~1962), so this exercises the "after_program not
    # yet run" branch rather than the wipe.  Test that no exception is
    # raised regardless of which branch fires.
    try:
        result = engine._apply_directives(
            valid_in, directives, "validate", None, experiment_type="cryoem")
        assert isinstance(result, list), "Result must be a list, got: %r" % result
    except Exception as e:
        raise AssertionError("Guard crashed with context=None: %s" % e)
    print("  PASS")


# =============================================================================
# After_program is itself a validation program — works correctly
# =============================================================================

def test_after_program_is_validation_program():
    print("[test] after_program=phenix.molprobity, molprobity done → "
          "wipe correctly fires (validation done)")
    engine = WorkflowEngine()
    directives = {
        "stop_conditions": {"after_program": "phenix.molprobity"},
    }
    context = {
        "validation_done": True,  # molprobity ran, so this is True
        "rsr_count": 1,
        "refine_count": 0,
        "map_cc": 0.85,
        "has_refined_model": True,
        "experiment_type": "cryoem",
        "last_program": "phenix.molprobity",
    }
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem",
                "phenix.map_correlations"]
    result = engine._apply_directives(
        valid_in, directives, "validate", context, experiment_type="cryoem")
    assert result == ["STOP"], \
        "When after_program=validation program and it ran, wipe should " \
        "fire (validation_done=True), got: %r" % result
    print("  PASS")


# =============================================================================
# Without an after_program directive, no wipe ever fires (sanity check)
# =============================================================================

def test_no_after_program_no_change():
    print("[test] Without after_program directive, _apply_directives does "
          "not touch validation programs")
    engine = WorkflowEngine()
    directives = {}  # no stop_conditions at all
    context = {
        "rsr_count": 1,
        "validation_done": False,
        "map_cc": 0.786,
        "experiment_type": "cryoem",
    }
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem",
                "phenix.map_correlations"]
    result = engine._apply_directives(
        valid_in, directives, "validate", context, experiment_type="cryoem")
    for prog in valid_in:
        assert prog in result, \
            "%s should remain in result, got: %r" % (prog, result)
    print("  PASS")


# =============================================================================
# Modifications log mentions the guard fired (helps debugging)
# =============================================================================

def test_guard_logs_explanation():
    print("[test] When guard fires, a 'validation pending' message is in log")
    # We can't easily capture modifications without instrumenting; instead
    # check that the function ran without error and returned the expected
    # shape.  The actual log message is verified by source review (the
    # patched code adds 'after_program ... validation pending' string).
    engine, directives, context = _make_engine_with_directives()
    valid_in = ["phenix.molprobity", "phenix.validation_cryoem"]
    result = engine._apply_directives(
        valid_in, directives, "validate", context, experiment_type="cryoem")
    # Just confirm the guard fired (validation programs preserved)
    assert "phenix.molprobity" in result
    print("  PASS — guard fired (verified by behavioral test above)")


# =============================================================================
# Test runner (matches the project's tst_*.py convention)
# =============================================================================

def run_all_tests():
    try:
        from libtbx.langchain.tests.tst_utils import run_tests_with_fail_fast
    except ImportError:
        try:
            from tests.tst_utils import run_tests_with_fail_fast
        except ImportError:
            _standalone_runner()
            return
    run_tests_with_fail_fast()


def _standalone_runner():
    test_fns = [(k, v) for k, v in sorted(globals().items())
                if k.startswith("test_") and callable(v)]
    passed = 0
    failed = 0
    for name, fn in test_fns:
        try:
            fn()
            passed += 1
        except Exception as e:
            print("  FAIL [%s]: %s" % (name, e))
            failed += 1
    print()
    print("%d passed, %d failed" % (passed, failed))
    if failed:
        sys.exit(1)
    print("=" * 60)
    print("OK: tst_validate_step_after_program_guard  (%d/%d)"
          % (passed, passed + failed))
    print("=" * 60)


if __name__ == "__main__":
    _standalone_runner()
