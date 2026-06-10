"""
Tests for v116.12 Fix #2: PLAN AUTO-STOP defense-in-depth.

Adds two things to graph_nodes.py PLAN node:
  (a) An additional elif that suppresses AUTO-STOP when workflow_engine
      has identified step="validate" with validation_done=False.  This
      provides defense-in-depth when plan_has_pending_stages doesn't
      fire (e.g. no expert plan, or the flag didn't propagate).
  (b) A diagnostic context dump in the AUTO-STOP path that captures
      plan_has_pending_stages, step, validation_done, after_program,
      and experiment_type.  Helps debug future AUTO-STOP regressions
      by recording exactly what state led to the stop.

These tests directly exercise the plan() function's elif chain using
synthetic state dicts.  They isolate the suppression logic from the
LLM / rules selection that runs after the chain.

Run with:
  python tst_plan_autostop_validation_suppression.py
  or
  phenix.python tst_plan_autostop_validation_suppression.py
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

# Comprehensive libtbx stub for graph_nodes import chain
if "libtbx" not in sys.modules:
    import types as _types
    for _mod in (
        "libtbx", "libtbx.langchain",
        "libtbx.langchain.agent",
        "libtbx.langchain.agent.pattern_manager",
        "libtbx.langchain.agent.event_log",
        "libtbx.langchain.agent.template_builder",
        "libtbx.langchain.agent.metrics_analyzer",
        "libtbx.langchain.agent.workflow_state",
        "libtbx.langchain.agent.rules_selector",
        "libtbx.langchain.agent.perceive_checks",
        "libtbx.langchain.knowledge",
        "libtbx.langchain.knowledge.prompts_hybrid",
    ):
        sys.modules[_mod] = _types.ModuleType(_mod)

    # Stub out symbols needed by graph_nodes imports
    pm = sys.modules["libtbx.langchain.agent.pattern_manager"]
    pm.patterns = {}
    pm.extract_cycle_number = lambda path, default=0: default
    pm.extract_all_numbers = lambda path: ()

    ev = sys.modules["libtbx.langchain.agent.event_log"]
    class _ET:
        DEBUG = "debug"
        STOP_DECISION = "stop_decision"
        PROGRAM_SELECTED = "program_selected"
        WORKFLOW_STATE_DETECTED = "workflow_state_detected"
        VALID_PROGRAMS = "valid_programs"
        METRICS_TREND = "metrics_trend"
    ev.EventType = _ET

    tb = sys.modules["libtbx.langchain.agent.template_builder"]
    class _TB: pass
    tb.TemplateBuilder = _TB

    ma = sys.modules["libtbx.langchain.agent.metrics_analyzer"]
    ma.derive_metrics_from_history = lambda *a, **k: {}
    ma.analyze_metrics_trend = lambda *a, **k: {}
    ma.get_latest_resolution = lambda *a, **k: None

    ws = sys.modules["libtbx.langchain.agent.workflow_state"]
    ws.detect_workflow_state = lambda **k: {"state": "test", "step_info": {"step": ""}, "context": {}, "valid_programs": []}
    ws.validate_program_choice = lambda *a, **k: (True, "")

    ph = sys.modules["libtbx.langchain.knowledge.prompts_hybrid"]
    ph.get_planning_prompt = lambda *a, **k: ("", "")

    rs = sys.modules["libtbx.langchain.agent.rules_selector"]
    rs.select_action_by_rules = lambda *a, **k: {"program": None}

    pc = sys.modules["libtbx.langchain.agent.perceive_checks"]
    pc.check_directive_stop = lambda *a, **k: (False, None, None)
    pc.check_consecutive_program_cap = lambda *a, **k: (False, None, None)


# Import the plan function from graph_nodes
try:
    from libtbx.langchain.agent.graph_nodes import plan as plan_node
except ImportError:
    from agent.graph_nodes import plan as plan_node


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


def _make_state(should_stop=True,
                reason="SUCCESS: Map-model CC (0.78) above target",
                trend_summary="Map CC: 0.78 - TARGET REACHED",
                plan_has_pending_stages=False,
                step="validate",
                validation_done=False,
                after_program=None,
                user_wants_ligandfit=False,
                experiment_type="cryoem",
                advice_changed=False,
                extra_ws_context=None):
    """Build a synthetic state dict for the PLAN node tests."""
    _ws_context = {
        "validation_done": validation_done,
        "user_wants_ligandfit": user_wants_ligandfit,
    }
    if extra_ws_context:
        _ws_context.update(extra_ws_context)
    return {
        "stop": False,
        "session_info": {
            "plan_has_pending_stages": plan_has_pending_stages,
            "advice_changed": advice_changed,
            "force_retry_program": None,
        },
        "metrics_trend": {
            "should_stop": should_stop,
            "reason": reason,
            "trend_summary": trend_summary,
        },
        "workflow_state": {
            "experiment_type": experiment_type,
            "step_info": {"step": step},
            "context": _ws_context,
        },
        "directives": {
            "stop_conditions": {"after_program": after_program} if after_program else {},
        },
        "history": [],
        "debug_log": [],
        "events": [],
    }


def _has_log_substring(result_state, substring):
    """True if any debug_log entry contains the given substring."""
    for line in result_state.get("debug_log", []):
        if substring in str(line):
            return True
    return False


# =============================================================================
# Tests for the new elif: workflow_state-based validation suppression
# =============================================================================

def test_validate_step_suppresses_autostop():
    """When step='validate' and validation_done=False, AUTO-STOP should be suppressed."""
    state = _make_state(
        should_stop=True,
        plan_has_pending_stages=False,  # plan_has_pending_stages OFF
        step="validate",
        validation_done=False,
    )
    result = plan_node(state)
    # Should NOT have stopped
    assert_false(result.get("stop", False),
                 "step=validate + validation_done=False: should NOT stop")
    # metrics_trend.should_stop should be cleared
    assert_false(result.get("metrics_trend", {}).get("should_stop", True),
                 "metrics_trend.should_stop should be False after suppression")
    # Debug log should mention the suppression
    assert_true(_has_log_substring(result, "Suppressing AUTO-STOP"),
                "debug_log should mention 'Suppressing AUTO-STOP'")
    assert_true(_has_log_substring(result, "validate"),
                "debug_log should mention 'validate'")


def test_validate_step_with_validation_already_done_does_not_suppress():
    """If validation_done=True, even with step='validate', AUTO-STOP fires."""
    state = _make_state(
        should_stop=True,
        plan_has_pending_stages=False,
        step="validate",
        validation_done=True,   # ← validation already done
    )
    result = plan_node(state)
    # AUTO-STOP should fire (no suppression)
    assert_true(result.get("stop", False),
                "step=validate but validation_done=True: SHOULD stop")


def test_non_validate_step_does_not_suppress():
    """When step is something other than 'validate', the new elif doesn't fire."""
    state = _make_state(
        should_stop=True,
        plan_has_pending_stages=False,
        step="complete",   # ← not 'validate'
        validation_done=False,
    )
    result = plan_node(state)
    # AUTO-STOP should fire
    assert_true(result.get("stop", False),
                "step=complete (not 'validate'): SHOULD stop")


def test_plan_pending_takes_priority_over_validate_step():
    """plan_has_pending_stages is checked FIRST (earlier elif).  Even with
    step='validate', the plan_has_pending_stages suppression should fire
    if set."""
    state = _make_state(
        should_stop=True,
        plan_has_pending_stages=True,   # ← takes priority
        step="validate",
        validation_done=False,
    )
    result = plan_node(state)
    assert_false(result.get("stop", False),
                 "plan_has_pending_stages=True: should NOT stop (matched first)")
    # The suppression log should mention "plan has pending stages",
    # not the validate-step variant.
    assert_true(_has_log_substring(result, "plan has pending stages"),
                "Should match the plan_has_pending_stages branch, not the new validate-step branch")


def test_after_program_takes_priority_over_validate_step():
    """after_program elif comes BEFORE the new validate-step elif."""
    state = _make_state(
        should_stop=True,
        plan_has_pending_stages=False,
        step="validate",
        validation_done=False,
        after_program="phenix.refine",
    )
    # Make sure after_program is not "done" — make history have a different program
    state["history"] = [{"program": "phenix.mtriage"}]
    result = plan_node(state)
    # after_program suppression should fire (not the validate-step one)
    assert_false(result.get("stop", False),
                 "after_program suppression takes priority: should NOT stop")


def test_no_should_stop_means_no_chain_runs():
    """When should_stop=False, the entire elif chain is skipped.
    PLAN proceeds to normal LLM/rules planning."""
    state = _make_state(
        should_stop=False,
        plan_has_pending_stages=False,
        step="validate",
        validation_done=False,
    )
    result = plan_node(state)
    # Doesn't matter if it stops or not — chain isn't entered.
    # The key assertion: NO "Suppressing AUTO-STOP" log entry from our branches.
    assert_false(_has_log_substring(result, "Suppressing AUTO-STOP"),
                 "When should_stop=False, no suppression logs should appear")


# =============================================================================
# Tests for the diagnostic context dump
# =============================================================================

def test_autostop_context_logged_on_stop():
    """When AUTO-STOP fires (else branch), a diagnostic context dump
    should be in debug_log."""
    state = _make_state(
        should_stop=True,
        plan_has_pending_stages=False,
        step="complete",   # nothing suppresses
        validation_done=True,
        after_program=None,
        experiment_type="cryoem",
    )
    result = plan_node(state)
    assert_true(result.get("stop", False),
                "Should have stopped (precondition for diagnostic dump)")
    # Diagnostic dump should be in the log
    assert_true(_has_log_substring(result, "AUTO-STOP context:"),
                "AUTO-STOP context dump should appear in debug_log")
    # Should include key fields
    assert_true(_has_log_substring(result, "plan_has_pending_stages="),
                "Context dump should include plan_has_pending_stages")
    assert_true(_has_log_substring(result, "validation_done="),
                "Context dump should include validation_done")
    assert_true(_has_log_substring(result, "step="),
                "Context dump should include step")
    assert_true(_has_log_substring(result, "experiment_type="),
                "Context dump should include experiment_type")
    assert_true(_has_log_substring(result, "after_program="),
                "Context dump should include after_program")


def test_autostop_context_NOT_logged_on_suppression():
    """The diagnostic context is only logged when AUTO-STOP actually fires.
    If suppression triggers, no diagnostic dump."""
    state = _make_state(
        should_stop=True,
        plan_has_pending_stages=True,  # Will suppress
        step="refine",
        validation_done=False,
    )
    result = plan_node(state)
    assert_false(result.get("stop", False),
                 "Should not have stopped (suppression)")
    assert_false(_has_log_substring(result, "AUTO-STOP context:"),
                 "Diagnostic dump should NOT appear when suppression fired")


def test_autostop_context_records_correct_values():
    """The dump should accurately reflect the values that gated the chain."""
    state = _make_state(
        should_stop=True,
        plan_has_pending_stages=False,
        step="complete",
        validation_done=True,
        after_program=None,
        experiment_type="xray",
    )
    result = plan_node(state)
    # Find the dump line
    dump_line = None
    for line in result.get("debug_log", []):
        if "AUTO-STOP context:" in str(line):
            dump_line = str(line)
            break
    assert_true(dump_line is not None, "Found the dump line")
    if dump_line:
        # Verify each value (use repr to be exact)
        assert_true("plan_has_pending_stages=False" in dump_line,
                    "Dump should show plan_has_pending_stages=False")
        assert_true("validation_done=True" in dump_line,
                    "Dump should show validation_done=True")
        assert_true("step='complete'" in dump_line,
                    "Dump should show step='complete'")
        assert_true("after_program=None" in dump_line,
                    "Dump should show after_program=None")
        assert_true("experiment_type='xray'" in dump_line,
                    "Dump should show experiment_type='xray'")


# =============================================================================
# End-to-end scenario tests
# =============================================================================

def test_af7mjs_scenario_suppressed_by_new_branch():
    """The full AF_7mjs scenario: cryoem, refined, validation pending.
    Even if plan_has_pending_stages didn't fire (the original bug),
    the new validate-step branch catches it."""
    state = _make_state(
        should_stop=True,
        reason="SUCCESS: Map-model CC (0.783) above target (0.70)",
        trend_summary="Map CC: 0.783 - TARGET REACHED",
        plan_has_pending_stages=False,  # Simulate the original bug
        step="validate",                # workflow_engine knows validation needed
        validation_done=False,
        experiment_type="cryoem",
    )
    result = plan_node(state)
    assert_false(result.get("stop", False),
                 "AF_7mjs scenario: new validate-step branch should suppress AUTO-STOP")
    assert_true(_has_log_substring(result, "validate"),
                "Log should mention validate-step suppression")


def test_af7mjs_scenario_double_protection():
    """The full AF_7mjs scenario with BOTH suppressors active.
    plan_has_pending_stages fires first.  This is the post-fix happy path."""
    state = _make_state(
        should_stop=True,
        reason="SUCCESS: Map-model CC (0.783) above target (0.70)",
        plan_has_pending_stages=True,
        step="validate",
        validation_done=False,
        experiment_type="cryoem",
    )
    result = plan_node(state)
    assert_false(result.get("stop", False),
                 "Should not stop (plan_has_pending_stages fires first)")


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
            import traceback
            FAILURES.append("%s: raised %s: %s\n%s" % (
                fn.__name__, type(e).__name__, e,
                traceback.format_exc()))
    print("%d passed, %d failed" % (PASS, FAIL))
    if FAILURES:
        for f in FAILURES:
            print("  FAIL: %s" % f)
    if FAIL:
        sys.exit(1)


if __name__ == "__main__":
    _standalone_runner()
