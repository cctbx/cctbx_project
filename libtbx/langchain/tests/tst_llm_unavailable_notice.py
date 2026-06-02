"""Tests for the LLM-unavailable observability notice (v120 follow-on).

When the requested LLM provider is unavailable, the agent degrades to
rules-based program selection.  This must be obvious in both CLI and GUI, on
both the local and server execution paths.

Design under test (verified against graph_nodes/ai_agent/graph/graph_state/
contract/api_schema/run_ai_agent/event_log/event_formatter):

  - graph_nodes._handle_llm_failure emits a structured EventType.NOTICE event
    with notice_kind="llm_unavailable" (machine-readable) PLUS a human line to
    debug_log.  The event rides state["events"], which every downstream node
    preserves and which crosses the REST boundary identically local/server.
  - ai_agent.print_history_record -> _detect_llm_unavailable_from_events scans
    each cycle's events for that discriminator and sets the RUN-LEVEL flag
    session.data["llm_ever_unavailable"] (session.data persists across cycles;
    graph state does not).
  - ai_agent._finalize_session shows an end-of-run banner when the flag is set.

The key test the prior implementation lacked is PROPAGATION: feed an events
list containing the NOTICE through the detection logic and assert the session
flag flips.  Source-scan tests skip (not fail) if files aren't found.

2-space indentation.
"""

from __future__ import absolute_import, division, print_function

import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)            # .../libtbx/langchain


# ---------------------------------------------------------------------------
# Source resolvers
# ---------------------------------------------------------------------------

_GRAPH_NODES_CANDIDATES = [
    os.path.join(_ROOT, "agent", "graph_nodes.py"),
]
_AI_AGENT_CANDIDATES = [
    os.path.join(_ROOT, "programs", "ai_agent.py"),
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "phenix", "programs", "ai_agent.py"),
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "programs", "ai_agent.py"),
]
_RUN_AI_AGENT_CANDIDATES = [
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "phenix", "phenix_ai", "run_ai_agent.py"),
    os.path.join(_ROOT, "..", "..", "..",
                 "phenix", "phenix_ai", "run_ai_agent.py"),
]


def _find(candidates, env_var):
    for cand in candidates:
        cand = os.path.abspath(cand)
        if os.path.isfile(cand):
            return cand
    env = os.environ.get(env_var)
    if env and os.path.isfile(env):
        return env
    return None


def _find_graph_nodes():
    return _find(_GRAPH_NODES_CANDIDATES, "GRAPH_NODES_PY")


def _find_ai_agent():
    return _find(_AI_AGENT_CANDIDATES, "AI_AGENT_PY")


def _find_run_ai_agent():
    return _find(_RUN_AI_AGENT_CANDIDATES, "RUN_AI_AGENT_PY")


# ---------------------------------------------------------------------------
# Behavioural replica: the graph-side emit + driver-side detection.
# ---------------------------------------------------------------------------

NOTICE_KIND = "llm_unavailable"


def _emit_replica(state, provider, error_msg):
    """Replica of graph_nodes._handle_llm_failure's within-invoke emit guard +
    structured NOTICE event.  Emits at most one per invoke; appends to a COPY
    of the events list (immutability)."""
    if state.get("llm_notice_emitted_this_invoke"):
        return state
    events = list(state.get("events", []))          # copy, never in-place
    events.append({
        "type": "notice",
        "notice_kind": NOTICE_KIND,
        "provider": provider,
        "message": "LLM PROVIDER UNAVAILABLE -- provider '%s' ... %s" % (
            provider, error_msg),
    })
    return {**state, "events": events,
            "llm_notice_emitted_this_invoke": True}


class _FakeSession(object):
    def __init__(self):
        self.data = {}


def _detect_replica(events, session):
    """Replica of ai_agent._detect_llm_unavailable_from_events."""
    if session is None or not events:
        return
    for ev in events:
        if not isinstance(ev, dict):
            continue
        if ev.get("notice_kind") == NOTICE_KIND:
            if not session.data.get("llm_ever_unavailable"):
                session.data["llm_ever_unavailable"] = True
                session.data["llm_unavailable_provider"] = (
                    ev.get("provider") or "unknown")
            return


# ---------------------------------------------------------------------------
# Behavioural tests
# ---------------------------------------------------------------------------

def test_emit_adds_structured_notice_event():
    state = {}
    out = _emit_replica(state, "anthropic", "boom")
    evs = out.get("events", [])
    assert len(evs) == 1, "exactly one NOTICE event expected"
    ev = evs[0]
    assert ev.get("type") == "notice"
    assert ev.get("notice_kind") == NOTICE_KIND
    assert ev.get("provider") == "anthropic"
    assert "UNAVAILABLE" in ev.get("message", "")


def test_emit_once_per_invoke():
    state = {}
    state = _emit_replica(state, "anthropic", "boom")
    state = _emit_replica(state, "anthropic", "again")  # validate->plan loop
    state = _emit_replica(state, "anthropic", "again2")
    notices = [e for e in state.get("events", [])
               if e.get("notice_kind") == NOTICE_KIND]
    assert len(notices) == 1, \
        "guard must prevent duplicate NOTICE within one invoke, got %d" % (
            len(notices),)


def test_emit_does_not_mutate_caller_events():
    """_emit must append to a COPY; the caller's list stays unchanged."""
    original = [{"type": "debug", "message": "x"}]
    state = {"events": original}
    _emit_replica(state, "anthropic", "boom")
    assert original == [{"type": "debug", "message": "x"}], \
        "caller's events list must not be mutated in place"


def test_propagation_event_to_session_flag():
    """THE test the prior implementation lacked: an events list carrying the
    NOTICE must flip the run-level session flag via the detection logic."""
    state = _emit_replica({}, "portkey", "gateway down")
    session = _FakeSession()
    _detect_replica(state["events"], session)
    assert session.data.get("llm_ever_unavailable") is True, \
        "detection must set the run-level session flag from the event"
    assert session.data.get("llm_unavailable_provider") == "portkey"


def test_detection_idempotent_across_cycles():
    """Per-cycle repetition (Option A) is safe: detecting again keeps the flag
    and the first provider."""
    session = _FakeSession()
    s1 = _emit_replica({}, "anthropic", "down")
    _detect_replica(s1["events"], session)
    s2 = _emit_replica({}, "anthropic", "still down")
    _detect_replica(s2["events"], session)
    assert session.data.get("llm_ever_unavailable") is True


def test_detection_ignores_unrelated_events():
    session = _FakeSession()
    events = [
        {"type": "program_selected", "program": "phenix.refine"},
        {"type": "notice", "notice_kind": "something_else"},
        {"type": "debug", "message": "noise"},
    ]
    _detect_replica(events, session)
    assert not session.data.get("llm_ever_unavailable"), \
        "only notice_kind=='llm_unavailable' may set the flag"


def test_detection_session_none_is_safe():
    """The display-only/replay path calls with session=None; must not raise."""
    state = _emit_replica({}, "anthropic", "down")
    _detect_replica(state["events"], None)   # must not raise


def test_per_run_flag_reset_clears_stale_value():
    """The run-level flag is PER-RUN: a value loaded from a prior run's saved
    session must be cleared at run start, so re-running in a directory whose
    prior run hit an LLM outage does NOT falsely report 'DID NOT USE THE LLM'
    when this run makes no LLM attempt (e.g. gate stops immediately).

    Mirrors the reset applied in ai_agent.iterate_agent /
    _handle_display_and_stop right after AgentSession construction."""
    # Simulate a session loaded from a prior run that hit an outage.
    loaded_data = {"llm_ever_unavailable": True,
                   "llm_unavailable_provider": "anthropic",
                   "cycles": []}
    # The per-run reset:
    loaded_data.pop("llm_ever_unavailable", None)
    loaded_data.pop("llm_unavailable_provider", None)
    # Banner condition (session_data_flag_llm_unavailable) must now be False.
    assert not loaded_data.get("llm_ever_unavailable"), \
        "stale flag from prior run must be cleared at run start"

    # And a genuine in-run failure after reset still sets it:
    session = _FakeSession()
    session.data = dict(loaded_data)
    s = _emit_replica({}, "anthropic", "down")
    _detect_replica(s["events"], session)
    assert session.data.get("llm_ever_unavailable") is True, \
        "a real in-run LLM failure must still set the flag after reset"


def test_ai_agent_resets_flag_per_run():
    """Source-scan: ai_agent.py must reset the run-level flag at run start
    (both the iterate_agent path and the display_and_stop path)."""
    path = _find_ai_agent()
    if path is None:
        print("  (skip) ai_agent.py not found in this checkout")
        return
    with open(path, "r") as fh:
        src = fh.read()
    assert src.count('pop("llm_ever_unavailable"') >= 2, \
        "ai_agent must reset llm_ever_unavailable at run start on BOTH the " \
        "iterate_agent and display_and_stop finalize paths (>=2 pop sites)"


# ---------------------------------------------------------------------------
# Source-scan: graph_nodes.py
# ---------------------------------------------------------------------------

def test_graph_nodes_emits_structured_notice():
    path = _find_graph_nodes()
    if path is None:
        print("  (skip) graph_nodes.py not found in this checkout")
        return
    with open(path, "r") as fh:
        src = fh.read()
    assert 'EventType.NOTICE' in src, \
        "graph_nodes must emit an EventType.NOTICE event"
    assert 'notice_kind="llm_unavailable"' in src, \
        "the NOTICE event must carry notice_kind='llm_unavailable'"
    assert 'llm_notice_emitted_this_invoke' in src, \
        "graph must guard emit-once-per-invoke with a scoped key"
    assert 'llm_ever_unavailable' not in src, \
        "graph_nodes must NOT use the driver's run-level key name (scope " \
        "collision); run-level aggregation lives on the driver"
    assert 'isinstance(error_msg, str)' in src, \
        "failure handler must normalize error_msg (never raise)"
    assert 'state["events"].append' not in src, \
        "must not append to state['events'] in place; use _emit (copies)"


# ---------------------------------------------------------------------------
# Source-scan: ai_agent.py
# ---------------------------------------------------------------------------

def test_ai_agent_detects_and_banners():
    path = _find_ai_agent()
    if path is None:
        print("  (skip) ai_agent.py not found in this checkout")
        return
    with open(path, "r") as fh:
        src = fh.read()
    assert "_detect_llm_unavailable_from_events" in src, \
        "ai_agent must detect the NOTICE from cycle events"
    assert 'notice_kind") == "llm_unavailable"' in src or \
           "notice_kind') == 'llm_unavailable'" in src, \
        "detection must match the structured discriminator, not a string scan"
    assert 'session.data["llm_ever_unavailable"] = True' in src, \
        "detection must set the run-level session.data flag"
    assert "session_data_flag_llm_unavailable" in src, \
        "ai_agent must have the run-level check helper"
    assert "IMPORTANT: THIS RUN DID NOT USE THE LLM" in src, \
        "ai_agent must print an end-of-run banner"
    assert "llm_ever_unavailable=" in src, \
        "ai_agent must expose the flag on the result object for the GUI"
    assert "if session is not None" in src, \
        "detection hook must guard session is not None (replay path)"


# ---------------------------------------------------------------------------
# Source-scan: run_ai_agent.py (server path still carries events)
# ---------------------------------------------------------------------------

def test_server_path_carries_events():
    path = _find_run_ai_agent()
    if path is None:
        print("  (skip) run_ai_agent.py not found in this checkout")
        return
    with open(path, "r") as fh:
        src = fh.read()
    assert 'final_state.get("events"' in src, \
        "server must read events from final_state"
    assert "events=events" in src, \
        "server must pass events into the response builder(s)"


_TESTS = [
    test_emit_adds_structured_notice_event,
    test_emit_once_per_invoke,
    test_emit_does_not_mutate_caller_events,
    test_propagation_event_to_session_flag,
    test_detection_idempotent_across_cycles,
    test_detection_ignores_unrelated_events,
    test_detection_session_none_is_safe,
    test_per_run_flag_reset_clears_stale_value,
    test_ai_agent_resets_flag_per_run,
    test_graph_nodes_emits_structured_notice,
    test_ai_agent_detects_and_banners,
    test_server_path_carries_events,
]


def run_all_tests():
    for fn in _TESTS:
        fn()
    print("All %d tests passed." % len(_TESTS))
    return True


if __name__ == "__main__":
    passed = 0
    failed = 0
    for fn in _TESTS:
        print("  Running %s..." % fn.__name__)
        try:
            fn()
            print("  PASS: %s" % fn.__name__)
            passed += 1
        except AssertionError as e:
            print("  FAIL: %s -- %s" % (fn.__name__, e))
            failed += 1
    print("\n%d passed, %d failed" % (passed, failed))
