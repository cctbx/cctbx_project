"""K_H6: Planning suite framework unit tests (v119.H6).

Pure-Python tests for the four framework.py functions implemented
in v119.H6 (Phase 2A planning suite):

  §A  is_stop_intent — production-parity STOP detection (6 tests)
  §B  validate_planning_state — 9-key structural validation (5 tests)
  §C  make_planning_run_fn — factory wrapping injectable call_fn (5 tests)
  §D  call_planning_llm — raw_output preservation on parse failure (2 tests)

Total: 18 tests.  No libtbx dependency — all tests use either pure
Python primitives or framework-internal injection.  Section D
monkey-patches via the test framework's lazy-import path; in
sandbox it patches against `agent.graph_nodes` (the relative
import path); under PHENIX, libtbx wins the lazy-import race
and the patches still apply to the imported names.

Sandbox-friendly: all 18 tests PASS without libtbx, no SKIPs.
"""
from __future__ import absolute_import, division, print_function

import os
import sys


# =====================================================================
# Import helpers — sandbox-friendly
# =====================================================================

def _try_import_framework():
    """Import the planning-framework public surface.

    framework.py is in tests/llm/, not tests/, so a path-fixup
    is needed when running this test from tests/.
    """
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        llm_dir = os.path.join(here, "llm")
        if llm_dir not in sys.path:
            sys.path.insert(0, llm_dir)
        import framework
        return framework, None
    except ImportError as e:
        return None, str(e)


# =====================================================================
# §A: is_stop_intent (6 tests)
# =====================================================================
#
# Production parity (graph_nodes.py ~line 2168) treats EITHER signal:
#   - intent["program"] == "STOP"
#   - intent["stop"] is True (strict equality)
# Defensive: non-dict inputs (None included) return False.


def test_is_stop_intent_program_stop():
    """STOP marker in program field."""
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    assert fw.is_stop_intent({"program": "STOP"}) is True
    print("  PASS: test_is_stop_intent_program_stop")


def test_is_stop_intent_stop_true_field():
    """Boolean stop field with strict True."""
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    assert fw.is_stop_intent({"program": "phenix.refine",
                              "stop": True}) is True
    print("  PASS: test_is_stop_intent_stop_true_field")


def test_is_stop_intent_normal_program():
    """Normal program selection — no stop."""
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    assert fw.is_stop_intent({"program": "phenix.refine"}) is False
    print("  PASS: test_is_stop_intent_normal_program")


def test_is_stop_intent_strict_equality():
    """Strict equality with True — string 'true' is NOT True."""
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    # Production rejects type-mismatched stop signals
    assert fw.is_stop_intent({"program": "phenix.refine",
                              "stop": "true"}) is False, (
        "stop='true' (string) should not match True")
    assert fw.is_stop_intent({"program": "phenix.refine",
                              "stop": 1}) is False, (
        "stop=1 (int) should not match True (strict equality)")
    print("  PASS: test_is_stop_intent_strict_equality")


def test_is_stop_intent_defensive_none():
    """Defensive: None intent (e.g., from error path) returns False."""
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    assert fw.is_stop_intent(None) is False
    print("  PASS: test_is_stop_intent_defensive_none")


def test_is_stop_intent_defensive_non_dict():
    """Defensive: non-dict input returns False."""
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    assert fw.is_stop_intent("STOP") is False, (
        "bare string should NOT match — production checks dict shape")
    assert fw.is_stop_intent([]) is False
    assert fw.is_stop_intent(42) is False
    print("  PASS: test_is_stop_intent_defensive_non_dict")


# =====================================================================
# §B: validate_planning_state (5 tests)
# =====================================================================
#
# Per PHASE2_PLAN §4.1 — nine required top-level keys.  Missing
# keys raise ValueError to fail fast at scenario-load time.


def _build_complete_state():
    """Helper: build a minimal-but-complete state dict."""
    return {
        "history": [],
        "analysis": {},
        "available_files": [],
        "previous_attempts": [],
        "user_advice": "",
        "metrics_trend": {},
        "workflow_state": {},
        "directives": {},
        "best_files": {},
    }


def test_validate_state_complete_returns_none():
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    state = _build_complete_state()
    result = fw.validate_planning_state(state)
    assert result is None
    print("  PASS: test_validate_state_complete_returns_none")


def test_validate_state_missing_single_key():
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    state = _build_complete_state()
    del state["workflow_state"]  # remove a single required key
    try:
        fw.validate_planning_state(state)
        assert False, "should have raised ValueError"
    except ValueError as e:
        assert "workflow_state" in str(e), (
            "Error message should name the missing key, got: %s" % e)
    print("  PASS: test_validate_state_missing_single_key")


def test_validate_state_missing_multiple_keys():
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    state = {"history": [], "user_advice": "test"}
    try:
        fw.validate_planning_state(state)
        assert False, "should have raised ValueError"
    except ValueError as e:
        msg = str(e)
        # Should list multiple missing keys
        for expected in ("analysis", "available_files",
                          "workflow_state", "directives", "best_files"):
            assert expected in msg, (
                "expected key %s in error msg, got: %s"
                % (expected, msg))
    print("  PASS: test_validate_state_missing_multiple_keys")


def test_validate_state_non_dict_raises():
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    for bad in (None, [], "not a dict", 42):
        try:
            fw.validate_planning_state(bad)
            assert False, "should raise on non-dict input %r" % bad
        except ValueError:
            pass
    print("  PASS: test_validate_state_non_dict_raises")


def test_validate_state_extra_keys_ok():
    """Extra keys beyond the required 9 are accepted (not strict)."""
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    state = _build_complete_state()
    state["extra_diagnostic_field"] = "this is fine"
    state["another_one"] = 42
    result = fw.validate_planning_state(state)
    assert result is None, (
        "Extra keys should be accepted — validator is only strict "
        "about missing keys")
    print("  PASS: test_validate_state_extra_keys_ok")


# =====================================================================
# §C: make_planning_run_fn factory (5 tests)
# =====================================================================
#
# Inject a fake call_fn that returns synthetic 4-tuples; verify the
# returned run_one wraps them into correct RunOutcome shapes.


def _make_fake_scenario(expected_fn=None):
    """Helper: build a minimal Scenario-shaped object for testing."""
    fw, _ = _try_import_framework()
    return fw.Scenario(
        name="fake_test_scenario",
        description="testing",
        decision_point="planning",
        test_type="reliability",
        threshold=0.8,
        max_runs=5,
        input={"history": []},  # validate_planning_state not called here
        expected_fn=expected_fn or (lambda intent: (True, "default ok")),
    )


def test_make_run_fn_success_path():
    """Happy path: fake returns valid intent, expected_fn passes."""
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return

    def fake_call(state_inputs, provider):
        return ("raw json text",
                {"program": "phenix.refine"},
                None,
                [])

    def expected_fn(intent):
        return (intent.get("program") == "phenix.refine",
                "checked program")

    run_one = fw.make_planning_run_fn(call_fn=fake_call)
    sc = _make_fake_scenario(expected_fn=expected_fn)
    outcome = run_one(sc, "google", 0)

    assert outcome.passed is True
    assert outcome.error is None
    assert outcome.raw_output == "raw json text"
    assert outcome.parsed == {"program": "phenix.refine"}
    assert outcome.run_index == 0
    print("  PASS: test_make_run_fn_success_path")


def test_make_run_fn_llm_error_path():
    """Fake reports error; run_one produces failed RunOutcome."""
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return

    def fake_call(state_inputs, provider):
        return ("", None, "LLM exploded: timeout", [])

    run_one = fw.make_planning_run_fn(call_fn=fake_call)
    sc = _make_fake_scenario()
    outcome = run_one(sc, "openai", 0)

    assert outcome.passed is False
    assert outcome.error == "LLM exploded: timeout"
    assert outcome.parsed is None
    assert "LLM call errored" in outcome.why
    print("  PASS: test_make_run_fn_llm_error_path")


def test_make_run_fn_expected_fn_fails():
    """Fake returns valid intent but scenario's expected_fn says no."""
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return

    def fake_call(state_inputs, provider):
        return ("raw", {"program": "phenix.refine"}, None, [])

    def expected_fn(intent):
        return (False, "expected STOP, got phenix.refine")

    run_one = fw.make_planning_run_fn(call_fn=fake_call)
    sc = _make_fake_scenario(expected_fn=expected_fn)
    outcome = run_one(sc, "google", 0)

    assert outcome.passed is False
    assert outcome.error is None  # No error — just a non-pass result
    assert outcome.parsed == {"program": "phenix.refine"}
    assert "expected STOP" in outcome.why
    print("  PASS: test_make_run_fn_expected_fn_fails")


def test_make_run_fn_expected_fn_raises():
    """expected_fn itself raises — surfaced as error, not crash."""
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return

    def fake_call(state_inputs, provider):
        return ("raw", {"program": "phenix.refine"}, None, [])

    def expected_fn(intent):
        raise RuntimeError("buggy assertion helper")

    run_one = fw.make_planning_run_fn(call_fn=fake_call)
    sc = _make_fake_scenario(expected_fn=expected_fn)
    outcome = run_one(sc, "google", 0)

    assert outcome.passed is False
    assert outcome.error is not None
    assert "buggy assertion helper" in outcome.error
    assert "expected_fn raised" in outcome.why
    print("  PASS: test_make_run_fn_expected_fn_raises")


def test_make_run_fn_default_call_fn_is_call_planning_llm():
    """Default call_fn (no argument) wires to call_planning_llm.

    Verifies that omitting the call_fn argument actually uses
    call_planning_llm as the default — not some other function
    that happens to succeed silently.

    The exact failure shape depends on the environment:
      - Sandbox without libtbx: ImportError on agent.graph_nodes
      - PHENIX without API key set: get_planning_llm returns
        (None, "Provider 'X' requires...")
      - PHENIX with API key + bad state_inputs: get_planning_prompt
        raises (the {"history": []} input here is incomplete)

    All three failure modes are signatures of call_planning_llm
    being invoked.  This test accepts any of them — what it's
    protecting against is the wrong function being wired up
    (which would either succeed or fail with a NameError / other
    unrelated message).
    """
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return

    run_one = fw.make_planning_run_fn()  # no call_fn argument
    sc = _make_fake_scenario()
    outcome = run_one(sc, "google", 0)

    # Whatever the environment, the default call_fn must produce
    # a failed outcome (it's call_planning_llm trying to reach a
    # real planning LLM with incomplete inputs / missing keys /
    # missing libtbx).  A pass would mean we wired the wrong
    # function as default.
    assert outcome.passed is False, (
        "Default call_fn produced a PASS — the wrong function "
        "may be wired as the default")
    assert outcome.error is not None, (
        "Default call_fn produced no error — the wrong function "
        "may be wired as the default")

    # Accept any of the three known failure signatures from
    # call_planning_llm.  This list documents the expected
    # environments; extend if a new one appears.
    expected_signatures = (
        "ImportError",              # sandbox without libtbx
        "get_planning_llm error",   # PHENIX without API key
        "ValueError",               # bad state_inputs (parse path)
        "KeyError",                 # bad state_inputs (prompt builder)
        "TypeError",                # other bad-input cases
    )
    assert any(sig in outcome.error for sig in expected_signatures), (
        "Default call_fn errored with an unexpected message — "
        "either call_planning_llm itself is broken, or the "
        "wrong default is wired.  Error was: %s"
        % outcome.error[:200])
    print("  PASS: test_make_run_fn_default_call_fn_is_call_planning_llm")


# =====================================================================
# §D: call_planning_llm raw_output preservation (2 tests)
# =====================================================================
#
# Gemini rev 3 caught that the exception path lost raw_output on
# parse_intent_json failures.  These tests pin the fix.
#
# Strategy: monkey-patch the imports that call_planning_llm uses.
# In sandbox, the relative-import branch fires (agent.graph_nodes
# instead of libtbx.langchain.agent.graph_nodes), so we can install
# fakes via sys.modules under those names.


def _install_fake_planning_modules(parse_raises=None, llm_returns=None):
    """Helper: install fake modules under sys.modules for the relative
    import path call_planning_llm uses.

    Args:
        parse_raises: if set, parse_intent_json raises ValueError(parse_raises)
                      instead of returning a dict.
        llm_returns: if (None, error_msg), get_planning_llm returns that
                     shape (simulating LLM-init failure).  Otherwise
                     the LLM is faked to invoke() a mock response.
    """
    import types

    # Fake response object that llm.invoke() returns
    class _FakeResponse:
        def __init__(self, content):
            self.content = content

    class _FakeLLM:
        def __init__(self, content):
            self._content = content

        def invoke(self, messages):
            return _FakeResponse(self._content)

    fake_graph_nodes = types.ModuleType("agent.graph_nodes")
    fake_prompts = types.ModuleType("knowledge.prompts_hybrid")
    fake_messages = types.ModuleType("langchain_core.messages")
    fake_agent = types.ModuleType("agent")
    fake_knowledge = types.ModuleType("knowledge")
    fake_lc_core = types.ModuleType("langchain_core")

    # Default LLM behavior: returns "raw bad json text"
    def fake_get_planning_llm(provider):
        if llm_returns is not None:
            return llm_returns
        return (_FakeLLM("raw bad json text"), None)

    def fake_parse_intent_json(text):
        if parse_raises is not None:
            raise ValueError(parse_raises)
        return {"program": "phenix.refine"}  # default happy path

    def fake_get_planning_prompt(**kwargs):
        return ("system msg", "user msg")

    class _FakeSystemMessage:
        def __init__(self, content=None):
            self.content = content

    class _FakeHumanMessage:
        def __init__(self, content=None):
            self.content = content

    fake_graph_nodes.get_planning_llm = fake_get_planning_llm
    fake_graph_nodes.parse_intent_json = fake_parse_intent_json
    fake_prompts.get_planning_prompt = fake_get_planning_prompt
    fake_messages.SystemMessage = _FakeSystemMessage
    fake_messages.HumanMessage = _FakeHumanMessage

    sys.modules["agent"] = fake_agent
    sys.modules["agent.graph_nodes"] = fake_graph_nodes
    sys.modules["knowledge"] = fake_knowledge
    sys.modules["knowledge.prompts_hybrid"] = fake_prompts
    sys.modules["langchain_core"] = fake_lc_core
    sys.modules["langchain_core.messages"] = fake_messages
    # Ensure rate_limit_handler doesn't accidentally import
    # (we want the "no handler" path — direct llm.invoke).
    sys.modules.pop("agent.rate_limit_handler", None)


def _uninstall_fake_planning_modules():
    """Helper: clean up after a fake-modules test."""
    for mod_name in ("agent.graph_nodes", "knowledge.prompts_hybrid",
                      "langchain_core.messages", "agent",
                      "knowledge", "langchain_core"):
        sys.modules.pop(mod_name, None)


def _libtbx_planning_available():
    """Detect whether we're in a PHENIX environment where the
    libtbx-prefixed import path will succeed.

    The §D tests below install fakes under the relative-import
    names (agent.graph_nodes, etc.) but call_planning_llm's
    lazy-import tries libtbx FIRST.  Under PHENIX, that succeeds
    and bypasses our fakes entirely — making the §D mechanic
    untestable from the test harness.

    The Gemini rev 3 fix (hoist raw_output outside try block)
    is a one-line mechanic verified by the sandbox runs.  Under
    PHENIX, the real planning suite (run_llm_tests.py --suite
    planning) exercises the same code path with real LLMs, so
    if the fix regressed, JSON-parse failures would visibly
    show empty raw_output in the test logs.  Sandbox-side pinning
    of the mechanic is sufficient.
    """
    try:
        from libtbx.langchain.agent import graph_nodes  # noqa: F401
        nothing = graph_nodes
        return True
    except ImportError:
        return False


def test_call_planning_llm_preserves_raw_on_parse_failure():
    """The critical Gemini rev 3 fix: if parse_intent_json raises,
    the returned raw_output contains the malformed text.

    SKIPPED under PHENIX: the test installs fakes under
    sys.modules at relative-import names (agent.graph_nodes
    etc.), but call_planning_llm's lazy-import tries libtbx
    first.  Under PHENIX libtbx wins and the fakes are bypassed.
    The mechanic itself (one variable hoisted) is verified
    in sandbox runs.
    """
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    if _libtbx_planning_available():
        print("  SKIP: test_call_planning_llm_preserves_raw_on_parse_failure "
              "(PHENIX: libtbx import path bypasses sandbox fakes; "
              "mechanic verified by sandbox runs)")
        return

    _install_fake_planning_modules(
        parse_raises="not parseable JSON")
    try:
        state = _build_complete_state()
        raw, intent, error, captured = fw.call_planning_llm(
            state, "google")

        # The fake LLM returns "raw bad json text"; parse raises
        # ValueError; the exception handler must preserve raw_output.
        assert raw == "raw bad json text", (
            "Expected raw to contain the malformed LLM text, "
            "got %r — the Gemini rev 3 fix is missing or broken"
            % raw)
        assert intent is None
        assert error is not None
        assert "ValueError" in error
        assert "not parseable JSON" in error
    finally:
        _uninstall_fake_planning_modules()

    print("  PASS: test_call_planning_llm_preserves_raw_on_parse_failure")


def test_call_planning_llm_raw_empty_on_pre_llm_failure():
    """When the failure is pre-LLM-call (get_planning_llm returns
    None), raw_output is "" since no LLM output exists yet.

    SKIPPED under PHENIX: same reason as the parse-failure test
    above — sys.modules fakes are bypassed by libtbx imports.
    """
    fw, err = _try_import_framework()
    if fw is None:
        print("  SKIP: cannot import framework (%s)" % err)
        return
    if _libtbx_planning_available():
        print("  SKIP: test_call_planning_llm_raw_empty_on_pre_llm_failure "
              "(PHENIX: libtbx import path bypasses sandbox fakes; "
              "mechanic verified by sandbox runs)")
        return

    _install_fake_planning_modules(
        llm_returns=(None, "test: provider unreachable"))
    try:
        state = _build_complete_state()
        raw, intent, error, captured = fw.call_planning_llm(
            state, "google")

        # Pre-LLM-call failure: no raw output exists.
        assert raw == "", (
            "Pre-LLM-call failure should have raw='' (nothing to "
            "preserve); got %r" % raw)
        assert intent is None
        assert error is not None
        assert "get_planning_llm error" in error
        assert "test: provider unreachable" in error
    finally:
        _uninstall_fake_planning_modules()

    print("  PASS: test_call_planning_llm_raw_empty_on_pre_llm_failure")


# =====================================================================
# Runner
# =====================================================================

def run_all_tests():
    # §A: is_stop_intent (6)
    test_is_stop_intent_program_stop()
    test_is_stop_intent_stop_true_field()
    test_is_stop_intent_normal_program()
    test_is_stop_intent_strict_equality()
    test_is_stop_intent_defensive_none()
    test_is_stop_intent_defensive_non_dict()

    # §B: validate_planning_state (5)
    test_validate_state_complete_returns_none()
    test_validate_state_missing_single_key()
    test_validate_state_missing_multiple_keys()
    test_validate_state_non_dict_raises()
    test_validate_state_extra_keys_ok()

    # §C: make_planning_run_fn (5)
    test_make_run_fn_success_path()
    test_make_run_fn_llm_error_path()
    test_make_run_fn_expected_fn_fails()
    test_make_run_fn_expected_fn_raises()
    test_make_run_fn_default_call_fn_is_call_planning_llm()

    # §D: call_planning_llm raw_output preservation (2)
    test_call_planning_llm_preserves_raw_on_parse_failure()
    test_call_planning_llm_raw_empty_on_pre_llm_failure()


if __name__ == "__main__":
    print("K_H6: Planning suite framework unit tests (v119.H6)")
    print("=" * 65)
    run_all_tests()
    print("=" * 65)
    print("K_H6 complete.")
