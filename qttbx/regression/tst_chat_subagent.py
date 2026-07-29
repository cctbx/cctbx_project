"""run_subagent runs a child conversation and returns its output.

The child is a full AgentSession at depth+1. Because depth>0 disables
autosave, the runner must persist the SubagentRecord itself; because
context_tokens is a per-agent monotonic peak, the child's usage must land on
the record and never on the parent; and because a child has no user, every
approval it hits must resolve rather than park the parent's worker.
"""

import shutil
import sys
import tempfile
import threading
from pathlib import Path

from libtbx import group_args
from libtbx.utils import Sorry, format_cpu_times, null_out

sys.path.insert(0, str(Path(__file__).resolve().parent))

from qttbx.widgets.chat.agent.base import ToolSpec
from qttbx.widgets.chat.agent.conversation import (
  Conversation, ContentBlock, Message, TokenUsage, now)
from qttbx.widgets.chat.agent.errors import CancelToken, TurnCancelled
# Usage reaches a message as its OWN event -- TurnDone carries only
# stop_reason/finish -- so a scripted turn emits TokenUsage separately, exactly
# as the backends do. Aliased as in tst_chat_session, because the stored
# conversation.TokenUsage shares the name.
from qttbx.widgets.chat.agent.events import (
  CANCELLED, CAP, CLEAN, ERROR, TRUNCATED,
  AgentError, ImageEmitted, TextDelta, ToolResultObserved, ToolUseRequested,
  TokenUsage as TokenUsageEvent, TurnDone)
from qttbx.widgets.chat.agent.profile import Profile
from qttbx.widgets.chat.agent.session import AgentSession
from qttbx.widgets.chat.agent.storage import ConversationStorage
from qttbx.widgets.chat.agent.subagent import (
  RUN_SUBAGENT_DESCRIPTION, register_run_subagent, run_child)
from qttbx.widgets.chat.agent.tools import ToolPolicy, ToolRegistry

from tst_chat_session import FakeAgent


# A child that cannot answer an approval used to park the parent's worker
# forever, so every run_child here goes through _run_child, which fails loudly
# instead of hanging the whole suite.
_DEADLINE_S = 20.0


def _parent(project_dir, **profile_kwargs):
  """Build a depth-0 parent session with a real storage."""
  storage = ConversationStorage(Path(project_dir), log=null_out())
  conv = Conversation.new(profile_name="parent", model="parent-model")
  storage.save(conv)
  kwargs = dict(subagents_max_depth=1)
  kwargs.update(profile_kwargs)
  return AgentSession(
    agent=FakeAgent([]), conversation=conv, storage=storage,
    tools=ToolRegistry(log=null_out()), policy=ToolPolicy(default="allow"),
    profile=Profile(name="parent", model="parent-model", **kwargs),
    depth=0, log=null_out()), storage


class _CancellingAgent(FakeAgent):
  """Replays *script*, then raises ``TurnCancelled`` instead of a further turn.

  ``TurnCancelled`` is the one SEMANTIC CANCEL shape that escapes ``run_turn``
  -- not the only exception that can (a crashing agent's own error escapes
  too, and must keep escaping; see
  :func:`exercise_a_child_crash_propagates`). ``AgentSession`` catches
  ``TurnCancelled`` only inside ``_dispatch_and_build_results``, so an unwind
  reaching ``run_turn`` from anywhere else -- an agent that aborts its stream
  on Stop -- leaves the caller holding it. That caller is the parent's own
  tool loop, where an escape costs the whole SubagentRecord: the persisted
  transcript and the GUI's only handle on the run, not merely the text.

  Raises BETWEEN turns: the previous turn is already committed and nothing is
  in flight. For the partway abort see :class:`_MidStreamCancellingAgent`.
  """

  def stream_turn(self, conversation, tools, cancel):
    if self._cursor >= len(self.turn_scripts):
      raise TurnCancelled()
    return FakeAgent.stream_turn(self, conversation, tools, cancel)


class _MidStreamCancellingAgent(FakeAgent):
  """Replays *script*, aborting mid-stream during the LAST scripted turn.

  The shape a backend makes when Stop lands mid-response: events are already
  delivered when the unwind starts, and it escapes ``_collect_one_response``
  before ``run_turn`` can ever append the in-progress assistant message -- so
  everything streamed in the aborted iteration is lost and only the ALREADY
  COMMITTED iterations remain. Script the aborting turn WITHOUT a ``TurnDone``:
  the abort is what ends it.
  """

  def stream_turn(self, conversation, tools, cancel):
    for event in FakeAgent.stream_turn(self, conversation, tools, cancel):
      yield event
    if self._cursor >= len(self.turn_scripts):
      raise TurnCancelled()


class _CrashingAgent(FakeAgent):
  """Fails with a genuine (non-cancel) error, the way a broken backend does."""

  def stream_turn(self, conversation, tools, cancel):
    raise RuntimeError("child agent exploded")


class _LateCrashingAgent(FakeAgent):
  """Replays *script*, then crashes instead of a further turn.

  The shape that costs the most: a child that ran, measured and reported, and
  only THEN hit a genuine (non-cancel) failure. ``_CrashingAgent`` dies before
  producing anything, so it cannot tell a runner that persists a crashed
  child's transcript from one that drops it.
  """

  def stream_turn(self, conversation, tools, cancel):
    if self._cursor >= len(self.turn_scripts):
      raise RuntimeError("child agent exploded")
    return FakeAgent.stream_turn(self, conversation, tools, cancel)


class _ClosableAgent(FakeAgent):
  """Records ``close()``, as every real backend agent implements it.

  The API agents close an HTTP connection pool; ``ClaudeCodeAgent.close()``
  tears down a ``claude`` CLI subprocess (which itself hosts the child's MCP
  servers) and a daemon asyncio loop thread. Neither has a finalizer, so
  nothing collects them if this is never called.
  """

  def __init__(self, script):
    FakeAgent.__init__(self, script)
    self.closed = 0

  def close(self):
    self.closed += 1


class _FakeConnection:
  """A STARTED ``McpServerConnection`` stand-in that records ``stop()``."""

  def __init__(self, name="srv"):
    self.name = name
    self.stopped = 0

  def stop(self):
    self.stopped += 1


class _UnstoppableConnection(_FakeConnection):
  """Its ``stop()`` raises, as a server whose subprocess already died can."""

  def stop(self):
    _FakeConnection.stop(self)
    raise RuntimeError("server already gone")


def _bundle(script, model="child-model", policy=None, tools=None,
            profile=None, agent=None, connections=None, notice=None,
            measurement_servers=()):
  """A ChildBundle whose agent replays *script* (or is *agent*, if given).

  *measurement_servers* is the phenix-side factory's inventory of the MCP
  servers a child was really wired up with. It is REQUIRED on a real bundle --
  the gate refuses one without it, because guessing means choosing between
  never flagging and always flagging -- so the default here is an empty list,
  the honest reading for a child wired to no measurement server at all. Pass
  ``False`` to omit the field and exercise that refusal.
  """
  bundle = group_args(
    agent=agent if agent is not None else FakeAgent(script),
    tools=tools if tools is not None else ToolRegistry(log=null_out()),
    policy=policy if policy is not None else ToolPolicy(default="allow"),
    profile=(profile if profile is not None else
             Profile(name="child", model=model, subagents_max_depth=0)),
    profile_name="child", model=model, backend="fake",
    connections=list(connections or []), notice=notice or "")
  if measurement_servers is not False:
    bundle.measurement_servers = list(measurement_servers)
  return bundle


def _run_child(*args, **kwargs):
  """``run_child`` on a daemon thread, with a deadline.

  A regression that re-parks the child on an unanswerable approval must FAIL
  the test, not wedge it: run_child holds the parent's worker thread and polls
  nothing while parked, so calling it inline would hang forever.
  """
  deadline = kwargs.pop("deadline_s", _DEADLINE_S)
  out = {}

  def _go():
    try:
      out["record"] = run_child(*args, **kwargs)
    except BaseException as exc:                   # re-raised on the caller
      out["error"] = exc

  thread = threading.Thread(target=_go, daemon=True)
  thread.start()
  thread.join(deadline)
  assert not thread.is_alive(), (
    "run_child did not return within %.0fs -- a child that cannot answer an "
    "approval parked the parent's worker" % deadline)
  if "error" in out:
    raise out["error"]
  return out["record"]


def _tool_result(messages, tool_use_id):
  """``(is_error, text)`` of the block answering *tool_use_id*, else None."""
  for message in messages:
    for block in message.content:
      if (block.type == "tool_result"
          and block.data.get("tool_use_id") == tool_use_id):
        text = "\n".join(b.data.get("text", "")
                         for b in block.data.get("content", []))
        return bool(block.data.get("is_error")), text
  return None


def exercise_child_runs_and_returns_final_text():
  """The child's final text comes back on the record."""
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    bundle = _bundle([[TextDelta(text="FINDINGS: none"),
                       TurnDone(stop_reason="end_turn")]])
    record = _run_child(parent, "review this", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    assert record.final_text == "FINDINGS: none", record.final_text
    assert record.parent_conversation_id == parent.conv.meta.id
    assert record.parent_tool_use_id == "t1"
    assert record.model == "child-model", record.model
    # A child that ran to its own end is complete, and says so in the field
    # code reads -- the empty string, never a stand-in for "don't know".
    assert record.incomplete_reason == "", record.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_record_is_persisted_despite_depth_gating():
  """depth>0 disables autosave, so the runner must store the record."""
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    bundle = _bundle([[TextDelta(text="done"),
                       TurnDone(stop_reason="end_turn")]])
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    loaded = storage.load_subagent(parent.conv.meta.id, record.sub_id)
    assert loaded.final_text == "done", loaded.final_text
    assert len(loaded.messages) >= 2, len(loaded.messages)
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_child_usage_never_touches_the_parent():
  """context_tokens is a per-agent peak; the parent's must stay untouched."""
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    # Seed a real usage-bearing parent turn: with an EMPTY parent transcript
    # the before/after comparison below would only ever match [] == [].
    parent.conv.append(Message(
      role="assistant", timestamp=now(),
      content=[ContentBlock(type="text", data={"text": "parent turn"})],
      usage=TokenUsage(input=5, output=1, context_tokens=1234)))
    before = [m.usage.context_tokens for m in parent.conv.messages
              if m.usage is not None]
    assert before == [1234], before
    bundle = _bundle([[TextDelta(text="x"),
                       TokenUsageEvent(input=10, context_tokens=900000),
                       TurnDone(stop_reason="end_turn")]])
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    after = [m.usage.context_tokens for m in parent.conv.messages
             if m.usage is not None]
    assert after == before, (before, after)
    assert record.token_usage.context_tokens == 900000, record.token_usage
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_costs_sum_but_context_peaks():
  """The four cost fields sum across iterations; context_tokens takes the max.

  Needs a child with TWO usage-bearing turns: with one, sum and max coincide
  and every mix-up of the two rules survives.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    # 'noop' is unregistered, so its dispatch fails with a Sorry that the
    # session turns into an error tool_result -- a second iteration without
    # needing a real tool.
    bundle = _bundle([
      [ToolUseRequested(id="tu1", name="noop", input={}),
       TokenUsageEvent(input=100, output=10, cache_read=4, cache_creation=1,
                       context_tokens=500000),
       TurnDone(stop_reason="tool_use")],
      [TextDelta(text="done"),
       TokenUsageEvent(input=200, output=20, cache_read=8, cache_creation=2,
                       context_tokens=600000),
       TurnDone(stop_reason="end_turn")]])
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    usage = record.token_usage
    assert usage.input == 300, usage
    assert usage.output == 30, usage
    assert usage.cache_read == 12, usage
    assert usage.cache_creation == 3, usage
    # 600000, NOT the 1100000 a sum would give: context_tokens is a peak.
    assert usage.context_tokens == 600000, usage
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_child_depth_is_derived_from_the_parent():
  """The child runs at parent.depth+1, so its own spawn attempt is refused.

  Derived, never assumed: a child left at the parent's depth would pass its
  own max_depth gate (unbounded recursion) and re-enable autosave on its
  throwaway conversation.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    spawned = []

    def _build(profile_name, model, backend):
      spawned.append((profile_name, model, backend))
      return _bundle([[TextDelta(text="GRANDCHILD RAN"),
                       TurnDone(stop_reason="end_turn")]])

    child_tools = ToolRegistry(log=null_out())
    register_run_subagent(child_tools, _build)
    bundle = _bundle(
      [[ToolUseRequested(id="tu1", name="run_subagent",
                         input={"task": "spawn a grandchild"}),
        TurnDone(stop_reason="tool_use")],
       [TextDelta(text="refused"), TurnDone(stop_reason="end_turn")]],
      tools=child_tools,
      # max_depth 1 (not 0) so the gate distinguishes depth 1 from depth 0.
      profile=Profile(name="child", model="child-model",
                      subagents_max_depth=1))
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    result = _tool_result(record.messages, "tu1")
    assert result is not None, "the grandchild tool_use was never answered"
    is_error, text = result
    assert is_error, text
    assert "depth" in text.lower(), text
    assert spawned == [], spawned
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_child_gets_its_own_event_sink():
  """No child event may reach the parent's on_event."""
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    seen = []
    parent.on_event = seen.append
    bundle = _bundle([[TextDelta(text="child text"),
                       TurnDone(stop_reason="end_turn")]])
    _run_child(parent, "task", bundle, CancelToken(), tool_use_id="t1",
               max_turns=5)
    assert seen == [], seen
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_an_ask_policy_is_denied_not_parked():
  """A child has no user, so an 'ask' must resolve rather than block.

  ToolPolicy's default IS 'ask', so the most natural child bundle walks
  straight into this. Parking would wedge the parent's worker past even the
  parent's Stop (the shared CancelToken is not polled while parked).
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    bundle = _bundle(
      [[ToolUseRequested(id="tu1", name="noop", input={}),
        TurnDone(stop_reason="tool_use")],
       [TextDelta(text="gave up"), TurnDone(stop_reason="end_turn")]],
      policy=ToolPolicy())                         # the DEFAULT: 'ask'
    assert bundle.policy.resolve("noop") == "ask"
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    is_error, text = _tool_result(record.messages, "tu1")
    assert is_error, text
    assert "denied" in text.lower(), text
    # The child's own words survive verbatim...
    assert record.final_text.startswith("gave up"), record.final_text
    # ...and are not the whole story: it was refused every tool it tried, so
    # "gave up" alone would read as a considered conclusion.
    assert "denied every tool" in record.final_text, record.final_text
    assert "noop" in record.final_text, record.final_text
    assert record.incomplete_reason == "tools_denied", record.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_the_destructive_floor_is_denied_not_parked():
  """The allow->ask rewrite for an allow_remember=False tool must resolve too.

  This path never consults the policy's decision, so coercing the child's
  policy away from 'ask' would leave it parking; the denial has to sit on the
  child's approval coordinator, which every path into _await_approval crosses.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    ran = []
    child_tools = ToolRegistry(log=null_out())
    child_tools.register_builtin(
      spec=ToolSpec(name="force_close", description="d",
                    input_schema={"type": "object"}),
      handler=lambda **kw: ran.append(True) or "closed",
      risk="destructive", allow_remember=False)
    bundle = _bundle(
      [[ToolUseRequested(id="tu1", name="force_close", input={}),
        TurnDone(stop_reason="tool_use")],
       [TextDelta(text="gave up"), TurnDone(stop_reason="end_turn")]],
      tools=child_tools,
      policy=ToolPolicy(default="allow"))          # allow, yet still asks
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    is_error, text = _tool_result(record.messages, "tu1")
    assert is_error, text
    assert "denied" in text.lower(), text
    assert ran == [], "a destructive child tool ran without any approval"
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_tool_returns_the_child_final_text():
  """The builtin hands build_child its args and returns the child's text."""
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    calls = []

    def _build(profile_name, model, backend):
      calls.append((profile_name, model, backend))
      return _bundle([[TextDelta(text="CHILD SAYS OK"),
                       TurnDone(stop_reason="end_turn")]])

    register_run_subagent(parent.tools, _build)
    result = parent.tools.invoke_builtin(
      "run_subagent",
      {"task": "measure the model", "profile": "reviewer",
       "model": "m", "backend": "b"},
      cancel=CancelToken(), session=parent, tool_use_id="t1")
    assert result == "CHILD SAYS OK", result     # the child's text, not the task
    assert calls == [("reviewer", "m", "b")], calls
    assert parent.tools.risk_of("run_subagent") == "write"
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_tool_passes_the_parent_turn_cap():
  """The PARENT profile's subagents_default_max_turns caps the child.

  Asserted at the tool boundary -- the string the parent MODEL reads -- so the
  cap is pinned where it is acted on, not merely on the record.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp, subagents_default_max_turns=1)
    script = [[TextDelta(text="FIRST"),
               ToolUseRequested(id="tu1", name="noop", input={}),
               TurnDone(stop_reason="tool_use")],
              [TextDelta(text="SECOND"), TurnDone(stop_reason="end_turn")]]
    register_run_subagent(
      parent.tools, lambda profile_name, model, backend: _bundle(script))
    result = parent.tools.invoke_builtin(
      "run_subagent", {"task": "measure"}, cancel=CancelToken(),
      session=parent, tool_use_id="t1")
    # Capped after one iteration; without the cap the second turn would run.
    assert result.startswith("FIRST"), result
    assert "SECOND" not in result, result
    # ...and the parent must be told it was cut off, or a child that stopped
    # early reads exactly like one that finished with nothing to report.
    assert "turn cap" in result.lower(), result
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_tool_refuses_when_disabled_or_taskless():
  """subagents_enabled=False and an empty task are both refused.

  build_child hands back a RUNNABLE child, so a dropped gate spawns one
  successfully and trips this test's own assertion -- rather than being killed
  incidentally by a FakeAgent that ran out of scripted turns.
  """
  tmp = tempfile.mkdtemp()

  def _build(profile_name, model, backend):
    return _bundle([[TextDelta(text="SPAWNED"),
                     TurnDone(stop_reason="end_turn")]])

  try:
    parent, storage = _parent(tmp, subagents_enabled=False)
    register_run_subagent(parent.tools, _build)
    try:
      parent.tools.invoke_builtin(
        "run_subagent", {"task": "anything"}, cancel=CancelToken(),
        session=parent, tool_use_id="t1")
    except Sorry as exc:
      assert "disabled" in str(exc).lower(), str(exc)
    else:
      raise AssertionError("a disabled profile still spawned a subagent")

    enabled, _ = _parent(tmp)
    register_run_subagent(enabled.tools, _build)
    for args in ({}, {"task": ""}):
      try:
        enabled.tools.invoke_builtin(
          "run_subagent", args, cancel=CancelToken(), session=enabled,
          tool_use_id="t1")
      except Sorry as exc:
        assert "task" in str(exc).lower(), str(exc)
      else:
        raise AssertionError("an empty task still spawned a subagent: %r"
                             % (args,))
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_direct_call_falls_back_to_25_turns():
  """max_turns=None with a profile that lacks the field caps at 25, not 1."""
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    bundle = _bundle(
      [[TextDelta(text="FIRST"),
        ToolUseRequested(id="tu1", name="noop", input={}),
        TurnDone(stop_reason="tool_use")],
       [TextDelta(text="SECOND"), TurnDone(stop_reason="end_turn")]],
      profile=group_args())                        # no subagents_* fields
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=None)
    assert record.final_text == "SECOND", record.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_the_childs_own_cap_outranks_the_callers():
  """A cap the child profile DECLARED beats the number the caller passed.

  Both call sites pass the PARENT profile's ``subagents_default_max_turns``,
  and every parent profile has one -- the dataclass fills it in at 25 whether
  or not the author wrote it -- so a child that declares its own cap could
  never reach ``max_turns=None`` and its number was inert on every API
  backend. A reviewer profile asking for 40 ran capped at the parent's 25,
  which truncates the review and blocks the readiness claim it was spawned to
  settle.

  ``bundle.max_turns`` is the child's DECLARED cap only: the factory writes it
  from the child's own JSON, so ``None`` (no declaration) still leaves the
  caller's number in force -- pinned by the second half here, or "the child
  always wins" would silently uncap every child whose profile says nothing.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    script = [[TextDelta(text="FIRST"),
               ToolUseRequested(id="tu1", name="noop", input={}),
               TurnDone(stop_reason="tool_use")],
              [TextDelta(text="SECOND"), TurnDone(stop_reason="end_turn")]]
    bundle = _bundle(list(script))
    bundle.max_turns = 1                       # what the child's JSON declared
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=25)
    assert record.final_text.startswith("FIRST"), record.final_text
    assert "SECOND" not in record.final_text, record.final_text
    assert record.incomplete_reason == "turn_cap", record.incomplete_reason
    # Declared nothing -> the caller's number still governs, so the same
    # script runs to its second turn.
    loose = _bundle(list(script))
    loose.max_turns = None
    record = _run_child(parent, "task", loose, CancelToken(),
                        tool_use_id="t2", max_turns=25)
    assert record.final_text == "SECOND", record.final_text
    assert not record.incomplete_reason, record.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_turn_cap_marks_the_record_incomplete():
  """Hitting the cap must be visible on the record, not silently clean.

  A child cut off mid-task that reads as a clean finish is indistinguishable
  from one that finished with nothing to say -- and the marker must ARRIVE
  ALONGSIDE whatever the child did produce, never in place of it.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)

    # Three tool-calling turns, but a cap of 2.
    def _turn():
      return [TextDelta(text="working"),
              ToolUseRequested(id="tu1", name="noop", input={}),
              TurnDone(stop_reason="tool_use")]

    bundle = _bundle([_turn(), _turn(), _turn()])
    bundle.tools.register_builtin(
      spec=ToolSpec(name="noop", description="noop",
                    input_schema={"type": "object"}),
      handler=lambda name, input, cancel, session, tool_use_id: "ok",
      risk="read")
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=2)
    assert "turn cap" in record.final_text.lower(), record.final_text
    assert "working" in record.final_text, record.final_text
    # The prose above is what the parent MODEL reads; this is what code reads.
    assert record.incomplete_reason == "turn_cap", record.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_quoted_cap_marker_is_not_a_cap():
  """A child that merely QUOTES the cap marker finished; it was not cut off.

  The declared consumer is a review agent, and reviewing the subagent module
  itself is a plausible task for it -- so the marker phrase turns up in the
  transcript of a perfectly clean run. A whole-transcript substring scan calls
  that run capped AND grafts the unrelated earlier text onto its findings, so
  the parent reads a false "cut off" over a corrupted report. Only the
  session's own marker, which is appended as the LAST message, counts.

  All three halves of that rule are load-bearing, so all three are exercised:
  an earlier message quoting the marker (here OPENING with it, so scanning the
  whole transcript misreads it however carefully it matches), a final message
  quoting it mid-sentence (so matching loosely misreads it however narrowly the
  scan is positioned), and a final ASSISTANT message that is the marker and
  nothing else (so a scan that checks position and text but not ROLE misreads
  it -- the session writes its marker as a USER message, and the same words in
  an assistant turn are the child talking about caps, not being stopped by
  one).
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    quote = ("[Subagent stopped at turn cap (N)] is what subagent.py appends "
             "at the cap; I checked that it always lands last.")
    bundle = _bundle(
      [[TextDelta(text=quote),
        ToolUseRequested(id="tu1", name="noop", input={}),
        TurnDone(stop_reason="tool_use")],
       [TextDelta(text="FINDINGS: none. Review complete."),
        TurnDone(stop_reason="end_turn")]])
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    assert record.incomplete_reason == "", record.incomplete_reason
    assert record.final_text == "FINDINGS: none. Review complete.", \
        record.final_text

    # ...and quoted in the child's OWN closing words, which is where a review
    # of this module would most naturally put it. Only a marker the session
    # appended STANDS ALONE as the message; a quote is embedded in a sentence.
    inline = ("Complete. subagent.py appends '%sN)]' at the cap."
              % "[Subagent stopped at turn cap (")
    record = _run_child(
      parent, "task",
      _bundle([[TextDelta(text=inline), TurnDone(stop_reason="end_turn")]]),
      CancelToken(), tool_use_id="t2", max_turns=5)
    assert record.incomplete_reason == "", record.incomplete_reason
    assert record.final_text == inline, record.final_text

    # ...and as the child's ENTIRE closing message, which a review agent asked
    # "what exactly does the session append?" answers with verbatim. Position
    # and text both match the real thing here; only the ROLE differs, and the
    # role is what makes a marker the SESSION's signal rather than the child's
    # words. Read as a cap, this clean run would be reported to the parent as
    # cut off.
    alone = "[Subagent stopped at turn cap (25)]"
    record = _run_child(
      parent, "task",
      _bundle([[TextDelta(text=alone), TurnDone(stop_reason="end_turn")]]),
      CancelToken(), tool_use_id="t3", max_turns=5)
    assert record.incomplete_reason == "", record.incomplete_reason
    assert record.final_text == alone, record.final_text

    # ...and the TASK quoting it, which is the one user-role message that can
    # end up last: a child that crashes before writing anything leaves its
    # transcript at the task alone. Position and role both match the real
    # marker here; only "is the message the marker AND NOTHING ELSE" tells them
    # apart. Matched by prefix instead, the child's own instructions are
    # returned to the parent as the child's findings.
    quoted_task = ("[Subagent stopped at turn cap (5)] -- explain when "
                   "subagent.py appends this and who reads it.")
    solo, solo_storage = _parent(tmp)            # its own, empty, record dir
    try:
      _run_child(solo, quoted_task, _bundle([], agent=_CrashingAgent([])),
                 CancelToken(), tool_use_id="t4", max_turns=5)
    except RuntimeError:
      pass
    else:
      raise AssertionError("the crashing child stopped raising")
    sub_dir = (solo_storage.root / "conversations" / solo.conv.meta.id
               / "subagents")
    crashed = sorted(sub_dir.glob("*.json"))
    assert len(crashed) == 1, [p.name for p in crashed]
    loaded = solo_storage.load_subagent(solo.conv.meta.id, crashed[0].stem)
    assert "turn cap" not in loaded.final_text.lower(), loaded.final_text
    assert "explain when" not in loaded.final_text, loaded.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_capped_child_still_reports_what_it_found():
  """A child cut off AFTER measuring must still deliver its findings.

  The most damaging way this module can fail, because it fails silently and in
  the safe-looking direction: the parent is told the child was cut off, so it
  has no reason to doubt the (empty) findings that came with the notice.

  ``AgentSession`` commits each completed tool iteration and opens a FRESH
  assistant message for the next round. When the cap fires on that next round
  the residual is appended carrying no text -- empty on the claude_code shape
  below, or holding only ``tool_use`` blocks on an API backend -- and the cap
  marker goes on after it. A ``_final_text`` that reads the last assistant
  MESSAGE therefore found the residual, found no text in it, and returned
  nothing: a child that measured a structure, found real clashes and ran out of
  turns reported the cap marker and NOT ONE of its findings.

  Both no-text residual shapes are exercised, since they arrive by different
  routes and a fix for one need not cover the other.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    # claude_code shape: the iteration is committed from ToolResultObserved,
    # and the residual assistant message is appended completely empty.
    bundle = _bundle([[TextDelta(text="FINDINGS: clashes at A:12, A:44"),
                       ToolUseRequested(id="tu1", name="mcp__probe__measure",
                                        input={}),
                       ToolResultObserved(tool_use_id="tu1",
                                          content="2 clashes", is_error=False),
                       TurnDone(stop_reason="tool_use")]])
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=1)
    assert "clashes at A:12, A:44" in record.final_text, record.final_text
    assert "turn cap" in record.final_text.lower(), record.final_text
    assert record.incomplete_reason == "turn_cap", record.incomplete_reason
    # The residual really is text-free, or this exercise proves nothing.
    assert record.messages[-2].role == "assistant", record.messages[-2].role
    assert not [b for b in record.messages[-2].content if b.type == "text"], \
        record.messages[-2].content

    # API shape: the residual carries the next round's tool_use blocks and no
    # text at all, so "has content" is not the same question as "said
    # something".
    bundle = _bundle([[TextDelta(text="FINDINGS: density unmodelled at B:7"),
                       ToolUseRequested(id="tu1", name="noop", input={}),
                       TurnDone(stop_reason="tool_use")],
                      [ToolUseRequested(id="tu2", name="noop", input={}),
                       TurnDone(stop_reason="tool_use")]])
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t2", max_turns=2)
    assert "density unmodelled at B:7" in record.final_text, record.final_text
    assert record.incomplete_reason == "turn_cap", record.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_crashed_child_keeps_its_transcript():
  """A crash must still leave a record: the error escapes, the work does not.

  ``depth > 0`` disables autosave AND checkpointing, so this module's own
  ``store_subagent`` call is the ONLY thing that ever writes a child's messages
  to disk. Letting a non-cancel exception propagate from ``run_turn`` skipped
  every line below it, so a child that ran twenty turns, measured plenty and
  then hit a broken backend lost its entire transcript -- the parent got a bare
  tool error and the GUI got no handle on the run at all.

  The error must still ESCAPE (see
  :func:`exercise_a_child_crash_propagates`): a crash swallowed into a returned
  record reads as a child that finished with nothing to say, which is the
  confusion this module exists to prevent. So both halves are asserted here --
  the raise reaches the caller, and the transcript is on disk when it does.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    script = [[TextDelta(text="FINDINGS: ligand density is absent"),
               ToolUseRequested(id="tu1", name="noop", input={}),
               TurnDone(stop_reason="tool_use")]]
    bundle = _bundle(script, agent=_LateCrashingAgent(list(script)))
    try:
      _run_child(parent, "task", bundle, CancelToken(), tool_use_id="t1",
                 max_turns=5)
    except RuntimeError as exc:
      assert "exploded" in str(exc), str(exc)
    else:
      raise AssertionError("a crashed child no longer raises")
    sub_dir = (storage.root / "conversations" / parent.conv.meta.id
               / "subagents")
    saved = sorted(sub_dir.glob("*.json")) if sub_dir.exists() else []
    assert len(saved) == 1, [p.name for p in saved]
    loaded = storage.load_subagent(parent.conv.meta.id, saved[0].stem)
    # A crash is its own disposition -- never "" (finished) and never
    # "cancelled" (which the parent asked for).
    assert loaded.incomplete_reason == "crashed", loaded.incomplete_reason
    assert "ligand density is absent" in loaded.final_text, loaded.final_text
    assert "crashed" in loaded.final_text.lower(), loaded.final_text
    assert "exploded" in loaded.final_text, loaded.final_text
    assert len(loaded.messages) >= 2, len(loaded.messages)
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_the_child_transcript_names_who_wrote_it():
  """Every assistant message the child produced names its model and backend.

  The parent stamps both on its own assistant messages; a child's were stamped
  from what the SESSION could see -- ``agent.model`` and ``profile.backend`` --
  and neither is authoritative for a child. The factory resolves the pair
  before building the agent and substitutes a different backend (running the
  parent's model) when the requested one has no usable credentials, while the
  child's ``Profile`` goes on naming the backend that was ASKED for. A stored
  transcript that cannot say which backend really answered defeats the whole
  point of commissioning a cross-model review.

  The bundle's resolved pair is the same one ``SubagentRecord.model`` carries,
  so the per-message stamp agrees with the record around it instead of being a
  third answer.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    # The profile deliberately names a DIFFERENT backend from the bundle's, as
    # a credentials substitution leaves it: the resolved pair must win.
    bundle = _bundle([[TextDelta(text="reviewed"),
                       ToolUseRequested(id="tu1", name="noop", input={}),
                       TurnDone(stop_reason="tool_use")],
                      [TextDelta(text="done"),
                       TurnDone(stop_reason="end_turn")]],
                     profile=Profile(name="child", model="asked-for-model",
                                     backend="anthropic",
                                     subagents_max_depth=0))
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    assistants = [m for m in record.messages if m.role == "assistant"]
    assert len(assistants) >= 2, len(assistants)
    for message in assistants:
      assert message.model == bundle.model, (message.model, bundle.model)
      assert message.backend == bundle.backend, (message.backend,
                                                 bundle.backend)
    # Survives the round trip, which is the point: the record on disk is what
    # a later reader consults.
    loaded = storage.load_subagent(parent.conv.meta.id, record.sub_id)
    stamped = [m for m in loaded.messages if m.role == "assistant"]
    assert stamped and all(m.backend == bundle.backend for m in stamped), \
        [(m.role, m.backend) for m in loaded.messages]
    # The user's own message is not an assistant turn and stays unstamped.
    assert all(m.model is None for m in loaded.messages if m.role == "user"), \
        [(m.role, m.model) for m in loaded.messages]
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_child_image_lands_under_the_parent_conversation():
  """A child's images must be resolvable from the conversation that keeps them.

  The child ``Conversation`` is never saved -- its transcript lives inside the
  ``SubagentRecord``, which is stored under the PARENT -- so attachments filed
  under the child's id created a ``conversations/<child_id>/attachments/``
  directory for a conversation that does not exist, and the record's image
  blocks, resolved against the conversation the record is stored under, pointed
  at nothing. A figure the child produced to support its findings was
  unviewable, which for a map or a density picture is the evidence itself.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    bundle = _bundle([[TextDelta(text="see the difference map"),
                       ImageEmitted(data=b"\x89PNG not-really", caption="map",
                                    mime="image/png"),
                       TurnDone(stop_reason="end_turn")]])
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    shas = [b.data.get("attachment_sha256") for m in record.messages
            for b in m.content if b.type == "image"]
    assert len(shas) == 1 and shas[0], shas
    # Resolvable from the parent id, which is the only id the record carries.
    assert storage.load_attachment(record.parent_conversation_id, shas[0]) \
        == b"\x89PNG not-really"
    # ...and no orphan directory was minted for the unsaved child.
    conv_dirs = sorted(p.name for p in
                       (storage.root / "conversations").iterdir() if p.is_dir())
    assert conv_dirs == [parent.conv.meta.id], conv_dirs
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_the_final_text_walk_back_is_bounded():
  """Skipping the residual must not reach back into superseded rounds.

  ``_final_text`` walks back over text-free assistant messages so a child cut
  off after real findings still delivers them. Unbounded, that walk crosses
  however many SILENT tool rounds the child ran -- a model that calls tools
  without narrating them writes no text round after round -- and returns its
  OPENING line as the report. "Let me start by looking at the model", handed to
  the parent as the findings, is worse than the empty string the walk was added
  to avoid, because it reads like an answer and the cap marker beside it makes
  it look like a truncated one.

  One step is the whole of what is needed and the whole of what is safe:
  ``_collect_one_response`` opens exactly ONE fresh assistant message per
  completed iteration, so at most one text-free residual can follow the child's
  last real turn.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    # Round 1 speaks; rounds 2 and 3 call tools in silence; the cap fires and
    # leaves the residual. Round 1's line is three rounds superseded.
    bundle = _bundle([[TextDelta(text="Let me start by looking at the model."),
                       ToolUseRequested(id="tu1", name="noop", input={}),
                       TurnDone(stop_reason="tool_use")],
                      [ToolUseRequested(id="tu2", name="noop", input={}),
                       TurnDone(stop_reason="tool_use")],
                      [ToolUseRequested(id="tu3", name="noop", input={}),
                       TurnDone(stop_reason="tool_use")]])
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=3)
    assert "Let me start by looking" not in record.final_text, (
      "an opening line the child superseded three rounds ago was returned as "
      "its report: %r" % record.final_text)
    assert record.incomplete_reason == "turn_cap", record.incomplete_reason
    assert "turn cap" in record.final_text.lower(), record.final_text
    # The control, and the reason the bound is 1 and not 0: ONE residual is
    # exactly what a cap leaves, and the round before it must still be read.
    bundle = _bundle([[TextDelta(text="FINDINGS: density unmodelled at B:7"),
                       ToolUseRequested(id="tu1", name="noop", input={}),
                       TurnDone(stop_reason="tool_use")],
                      [ToolUseRequested(id="tu2", name="noop", input={}),
                       TurnDone(stop_reason="tool_use")]])
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t2", max_turns=2)
    assert "density unmodelled at B:7" in record.final_text, record.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_provider_executed_tool_is_a_measurement():
  """A child that measured through server-side tools measured SOMETHING.

  Provider-executed tools -- web search and its siblings, opted into through a
  profile's ``server_tools`` -- run at the provider and are recorded by
  ``_accumulate`` as ``server_tool_use`` / ``server_tool_result`` blocks, not as
  the ``tool_result`` the dispatch loop writes. ``_any_tool_succeeded`` scanned
  only the latter, so such a child was reported as having measured NOTHING: a
  false accusation, the mirror of the failure this backstop exists to catch,
  and one that makes a good review look broken. It stayed hidden while the flag
  was gated on a denial record, since a child using provider tools is denied
  nothing.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)

    class _ServerToolAgent(FakeAgent):
      """Emits the two blocks a provider-executed tool leaves behind."""

      def __init__(self, content):
        FakeAgent.__init__(self, [])
        self._content = content

      def stream_turn(self, conversation, tools, cancel):
        from qttbx.widgets.chat.agent.events import (
          ServerToolResult, ServerToolUsed)
        yield ServerToolUsed(id="s1", name="web_search", input={})
        yield ServerToolResult(tool_use_id="s1", content=self._content)
        yield TextDelta(text="FINDINGS: the deposited model differs at A:12.")
        yield TurnDone(stop_reason="end_turn", finish=CLEAN)

    ok = {"type": "web_search_tool_result",
          "content": [{"type": "web_search_result", "title": "PDB 1yjp"}]}
    bundle = _bundle([], agent=_ServerToolAgent(ok),
                     measurement_servers=["probe"])
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    assert not record.incomplete_reason, (
      "a child that measured through provider-executed tools was accused of "
      "measuring nothing: %r" % record.incomplete_reason)
    assert "measured nothing" not in record.final_text, record.final_text
    # ...and the provider's own ERROR shape is still not a measurement, or the
    # widening would have turned every failed search into evidence.
    bad = {"type": "web_search_tool_result_error", "error_code": "max_uses"}
    bundle = _bundle([], agent=_ServerToolAgent(bad),
                     measurement_servers=["probe"])
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t2", max_turns=5)
    assert record.incomplete_reason == "no_measurement", \
        record.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_late_provider_failure_is_not_a_turn_cap():
  """A completed ROUND is not a completed TURN.

  ``AgentSession`` runs the API backends' tool loop over several rounds, and
  every completed round reports ``TurnDone(stop_reason='tool_use')`` before the
  next begins. A ``saw_turn_done`` latch that only ever SET therefore found
  itself already true when the provider failed on round two: the no-terminal
  arm was skipped, ``_outcome_for_finish`` read round one's ``TOOL_USE`` as a
  cap, and the parent was told to raise a turn cap the child never reached --
  for a rate limit or an outage the record never mentioned. Raising the cap
  reproduces the failure exactly.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    bundle = _bundle([[TextDelta(text="Checking. "),
                       ToolUseRequested(id="tu1", name="noop", input={}),
                       TurnDone(stop_reason="tool_use")],
                      [TextDelta(text="FINDINGS: nothing yet"),
                       AgentError(message="rate limited", recoverable=True)]],
                     measurement_servers=["probe"])
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=9)
    assert record.incomplete_reason == "error", (
      "a provider failure after a completed round was reported as %r; the "
      "parent would raise a cap that was never reached"
      % record.incomplete_reason)
    assert "rate limited" in record.final_text, record.final_text
    assert "turn cap" not in record.final_text.lower(), record.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_denial_survives_another_reason_for_stopping():
  """Two things went wrong, and the parent needs both.

  A child refused every tool AND cut off by a provider failure has two
  independent faults, and they call for opposite responses: the failure says
  retry, the denial says the profile grants this child nothing and retrying
  reproduces it exactly. Ordered as one more ``elif`` arm, the zero-measurement
  reading was simply LOST whenever another arm matched first -- so the parent
  saw only the transient half and retried into the permanent one.

  ``incomplete_reason`` keeps the terminal cause, since that is what ended the
  child; the denial rides in ``final_text``, which is the only thing about the
  child the parent MODEL ever sees.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    bundle = _bundle([[TextDelta(text="Checking. "),
                       ToolUseRequested(id="tu1", name="probe_measure",
                                        input={}),
                       TurnDone(stop_reason="tool_use")],
                      [TextDelta(text="FINDINGS: looks fine"),
                       AgentError(message="rate limited", recoverable=True)]],
                     policy=ToolPolicy(default="deny"),
                     measurement_servers=["probe"])
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=9)
    assert record.incomplete_reason, record.incomplete_reason
    low = record.final_text.lower()
    assert "rate limited" in low, record.final_text
    assert "denied every tool" in low, (
      "the denial was lost behind the provider failure, so the parent is told "
      "to retry a spawn whose real fault is its policy: %r" % record.final_text)
    assert "probe_measure" in record.final_text, record.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_one_denial_is_recorded_once():
  """Each denial is recorded by the layer that made it, exactly once.

  An unanswerable ``ask`` was counted TWICE: once by
  ``_HeadlessApprovals.open``, which decides it, and again by the dispatch
  loop, which recorded every ``deny`` decision whatever settled it -- and
  labelled the second copy "tool policy", which is the one thing it was not.
  For a child the difference between "your profile denies this tool" and "this
  tool needs an approval and a child has nobody to give one" is the difference
  between editing a policy and not spawning that way at all.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    from qttbx.widgets.chat.agent.subagent import _HeadlessApprovals
    seen = []
    real_record = _HeadlessApprovals.record_denial

    def _spy(self, tool_name, reason="approval required, no user available"):
      seen.append((tool_name, reason))
      return real_record(self, tool_name, reason)

    _HeadlessApprovals.record_denial = _spy
    try:
      # An 'ask' the child cannot answer: decided inside open().
      bundle = _bundle([[ToolUseRequested(id="tu1", name="probe_measure",
                                          input={}),
                         TurnDone(stop_reason="tool_use")],
                        [TextDelta(text="done"),
                         TurnDone(stop_reason="end_turn", finish=CLEAN)]],
                       policy=ToolPolicy(default="ask"),
                       measurement_servers=["probe"])
      record = _run_child(parent, "review", bundle, CancelToken(),
                          tool_use_id="t1", max_turns=5)
      assert len(seen) == 1, (
        "one unanswerable approval was recorded %d times: %s"
        % (len(seen), seen))
      assert "no user" in seen[0][1], seen
      assert record.incomplete_reason == "tools_denied", \
          record.incomplete_reason
      # A POLICY deny is settled without opening a request, so it is still
      # recorded -- once, by the resolver, with its own reason.
      seen[:] = []
      bundle = _bundle([[ToolUseRequested(id="tu1", name="probe_measure",
                                          input={}),
                         TurnDone(stop_reason="tool_use")],
                        [TextDelta(text="done"),
                         TurnDone(stop_reason="end_turn", finish=CLEAN)]],
                       policy=ToolPolicy(default="deny"),
                       measurement_servers=["probe"])
      record = _run_child(parent, "review", bundle, CancelToken(),
                          tool_use_id="t2", max_turns=5)
      assert len(seen) == 1, (
        "one policy denial was recorded %d times: %s" % (len(seen), seen))
      assert seen[0] == ("probe_measure", "tool policy"), seen
      assert record.incomplete_reason == "tools_denied", \
          record.incomplete_reason
    finally:
      _HeadlessApprovals.record_denial = real_record
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_child_reads_back_the_images_it_was_filed():
  """What the session WRITES an attachment under is what the agent READS.

  Filing a child's images under the PARENT's conversation is right -- the child
  conversation is never saved -- but it left the agents resolving image blocks
  against ``conversation.meta.id``, the CHILD's. Writing to one id and reading
  from another means a child handed an image by a tool cannot send it back on
  its next round: the load raises, mid-turn, on bytes that were on disk the
  whole time. ``AgentSession`` now tells the agent which id it files under.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    seen = {}

    class _ImageReadingAgent(FakeAgent):
      """Emits an image, then reads it back the way a real backend does."""

      def __init__(self):
        FakeAgent.__init__(self, [])
        self._round = 0

      def stream_turn(self, conversation, tools, cancel):
        self._round += 1
        if self._round == 1:
          yield ImageEmitted(data=b"\x89PNG evidence", caption="map",
                             mime="image/png")
          yield ToolUseRequested(id="tu1", name="noop", input={})
          yield TurnDone(stop_reason="tool_use")
          return
        from qttbx.widgets.chat.agent.base import load_image_attachment
        conv_id = self.attachment_conv_id_for(conversation)
        seen["conv_id"] = conv_id
        for message in conversation.messages:
          for block in message.content:
            if block.type == "image":
              # Caught, so the split is REPORTED rather than escaping as the
              # bare "Attachment not found" a user would see mid-turn.
              try:
                seen["bytes"] = load_image_attachment(
                  storage, conv_id, block, Sorry)[1]
              except Sorry as exc:
                seen["error"] = str(exc)
        yield TextDelta(text="FINDINGS: the map supports it.")
        yield TurnDone(stop_reason="end_turn", finish=CLEAN)

    bundle = _bundle([], agent=_ImageReadingAgent(),
                     measurement_servers=["probe"])
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    assert "error" not in seen, (
      "the child could not read back an image the session had just filed for "
      "it (%s): the write used %r and the read used %r"
      % (seen["error"], parent.conv.meta.id, seen.get("conv_id")))
    assert seen.get("conv_id") == parent.conv.meta.id, (
      "the agent resolved the child's images against %r while the session "
      "filed them under %r" % (seen.get("conv_id"), parent.conv.meta.id))
    assert seen.get("bytes") == b"\x89PNG evidence", seen
    assert "the map supports it" in record.final_text, record.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_non_string_task_is_refused_not_crashed():
  """A ``task`` that is not a string is refused with a Sorry, naming its type.

  The schema asks for a string; a schema is a request, not a guarantee, and a
  model emitting an object or an array here is an ordinary mis-generation. The
  old ``(args.get("task") or "").strip()`` met it with a bare
  ``AttributeError`` -- an unhandled crash inside a tool dispatch, where this
  layer contracts to raise ``Sorry`` -- carrying a message ("'dict' object has
  no attribute 'strip'") the model cannot act on. Naming the type it sent makes
  the retry an informed one.
  """
  tmp = tempfile.mkdtemp()

  def _build(profile_name, model, backend):
    raise AssertionError("a non-string task still spawned a child")

  try:
    parent, storage = _parent(tmp)
    register_run_subagent(parent.tools, _build)
    for bad in ({"instructions": "review"}, ["review"], 7, 1.5, True):
      try:
        parent.tools.invoke_builtin(
          "run_subagent", {"task": bad}, cancel=CancelToken(),
          session=parent, tool_use_id="t1")
      except Sorry as exc:
        assert "task" in str(exc), str(exc)
        assert type(bad).__name__ in str(exc), (str(exc), type(bad).__name__)
      else:
        raise AssertionError("a %s task was accepted" % type(bad).__name__)
    # A real string still works, so the guard refuses only what it should.
    # A fresh session, because ToolRegistry.register_builtin does not overwrite
    # an existing name and the refusing _build above is still installed.
    ok, _ = _parent(tmp)
    register_run_subagent(
      ok.tools,
      lambda profile_name, model, backend: _bundle(
        [[TextDelta(text="SPAWNED"), TurnDone(stop_reason="end_turn")]]))
    assert ok.tools.invoke_builtin(
      "run_subagent", {"task": "  review this  "}, cancel=CancelToken(),
      session=ok, tool_use_id="t2") == "SPAWNED"
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_the_two_turn_cap_keys_mean_one_thing_each():
  """``subagents.default_max_turns`` caps CHILDREN; ``max_turns`` caps oneself.

  The two readings used to share one key. A parent's
  ``subagents_default_max_turns`` capped the children it spawned, and the same
  field read off a CHILD profile was taken as that child's own cap -- one key
  carrying opposite meanings depending on who held it, which a profile in the
  middle of a spawn chain would have to carry both of at once.

  The split, pinned in three parts:

  - a child profile's ``subagents_default_max_turns`` no longer caps the child.
    It is that profile's default for the grandchildren it would spawn, and
    nothing here spawns any, so it must be inert;
  - a child profile's own ``max_turns`` DOES cap it, outranking the caller's
    default -- the reviewer case, where a profile that needs 40 turns was
    silently held to the parent's 25;
  - a caller naming nothing and a child declaring nothing land on the module's
    own default, borrowing neither profile's ``subagents`` number.
  """
  tmp = tempfile.mkdtemp()

  def _two_turns():
    return [[TextDelta(text="FIRST"),
             ToolUseRequested(id="tu1", name="noop", input={}),
             TurnDone(stop_reason="tool_use")],
            [TextDelta(text="SECOND"), TurnDone(stop_reason="end_turn")]]

  try:
    parent, storage = _parent(tmp)
    # A child whose subagents block says 1: about ITS children, not itself.
    inert = _bundle(_two_turns(),
                    profile=Profile(name="child", model="m",
                                    subagents_default_max_turns=1,
                                    subagents_max_depth=0))
    record = _run_child(parent, "task", inert, CancelToken(),
                        tool_use_id="t1", max_turns=None)
    assert record.final_text == "SECOND", record.final_text
    assert record.incomplete_reason == "", record.incomplete_reason

    # The child's OWN key, which does cap it -- and beats the caller's number.
    own = _bundle(_two_turns(),
                  profile=Profile(name="child", model="m", max_turns=1,
                                  subagents_max_depth=0))
    record = _run_child(parent, "task", own, CancelToken(),
                        tool_use_id="t2", max_turns=25)
    assert record.final_text.startswith("FIRST"), record.final_text
    assert "SECOND" not in record.final_text, record.final_text
    assert record.incomplete_reason == "turn_cap", record.incomplete_reason

    # ...and the other direction: a child asking for MORE than the caller's
    # default gets it, which is the case the reviewer profile exists for.
    roomy = _bundle(_two_turns(),
                    profile=Profile(name="child", model="m", max_turns=40,
                                    subagents_max_depth=0))
    record = _run_child(parent, "task", roomy, CancelToken(),
                        tool_use_id="t3", max_turns=1)
    assert record.final_text == "SECOND", record.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_an_unusable_turn_cap_is_refused_before_the_child_runs():
  """A cap below 1 is refused, not obeyed and not silently clamped.

  ``AgentSession`` tests ``iterations >= max_turns`` on the FIRST pass, so a
  declared cap of ``0`` or a negative -- inert while nothing read the child's
  declaration -- cuts the child off before it dispatches a single tool. That is
  a child which measured nothing for a configuration reason nobody would think
  to suspect, and its record says only ``turn_cap``.

  Refusal beats clamping: a clamp to 1 runs the child under a number its author
  never wrote, one turn is itself a zero-measurement child on any real task, and
  the misconfiguration stays hidden while still costing a provider call. The
  refusal names the key so the profile's author can fix it.

  Whatever ``build_child`` opened is still RELEASED on the refusal -- a spawn
  that never ran must not leak an agent or an MCP server subprocess.
  """
  tmp = tempfile.mkdtemp()

  def _script():
    return [[TextDelta(text="FINDINGS: real work"),
             ToolUseRequested(id="tu1", name="noop", input={}),
             TurnDone(stop_reason="tool_use")],
            [TextDelta(text="done"), TurnDone(stop_reason="end_turn")]]

  try:
    parent, storage = _parent(tmp)
    # The child's own declaration, and the caller's default: both are capped
    # numbers and both must be validated, or the hole simply moves.
    for declared, caller in ((0, 25), (-5, 25), (None, 0), (None, -1)):
      agent = _ClosableAgent(_script())
      conn = _FakeConnection()
      bundle = _bundle(_script(), agent=agent, connections=[conn])
      if declared is not None:
        bundle.max_turns = declared
      try:
        _run_child(parent, "task", bundle, CancelToken(), tool_use_id="t1",
                   max_turns=caller)
      except Sorry as exc:
        assert "at least 1" in str(exc), str(exc)
        assert "max_turns" in str(exc), str(exc)
      else:
        raise AssertionError("cap declared=%r caller=%r was accepted"
                             % (declared, caller))
      assert agent.closed == 1, (declared, caller, agent.closed)
      assert conn.stopped == 1, (declared, caller, conn.stopped)

    # A non-integer is refused the same way rather than reaching run_turn's
    # comparison, where it would raise a raw TypeError from inside the loop.
    for bad in ("5", 2.5, True):
      bundle = _bundle(_script())
      bundle.max_turns = bad
      try:
        _run_child(parent, "task", bundle, CancelToken(), tool_use_id="t1",
                   max_turns=25)
      except Sorry as exc:
        assert "integer" in str(exc), str(exc)
      else:
        raise AssertionError("a %r cap was accepted" % (bad,))

    # Nothing was recorded for any refused spawn: no half-run child exists.
    sub_dir = (storage.root / "conversations" / parent.conv.meta.id
               / "subagents")
    assert not sub_dir.exists() or not list(sub_dir.glob("*.json")), \
        sorted(p.name for p in sub_dir.glob("*.json"))
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_child_crash_propagates():
  """A crashed child is not a cancelled one: its error must still escape.

  The cancel guard catches ``TurnCancelled`` and nothing else on purpose.
  Widened to ``except Exception`` it would convert a crash into a recorded,
  non-error, empty-text child -- indistinguishable from one that ran fine
  and found nothing -- and every other exercise in this file would still pass.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    bundle = _bundle([], agent=_CrashingAgent([]))
    try:
      _run_child(parent, "task", bundle, CancelToken(), tool_use_id="t1",
                 max_turns=5)
    except RuntimeError as exc:
      assert "exploded" in str(exc), str(exc)
    else:
      raise AssertionError("a crashed child was swallowed as a clean record")
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_cancel_returns_partial_text_not_an_exception():
  """A cancelled child yields what it had; the parent must not see a raise.

  Both cancel shapes: the ordinary one, where the shared token trips mid-run
  and ``_run_turn`` short-circuits, and the unwinding one, where the turn
  raises ``TurnCancelled`` past every handler inside the session.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)

    # 1. The token trips during a tool call -- the shape a user's Stop makes.
    def _stop(name, input, cancel, session, tool_use_id):
      cancel.set()
      return "stopped"

    cancel = CancelToken()
    child_tools = ToolRegistry(log=null_out())
    child_tools.register_builtin(
      spec=ToolSpec(name="stopper", description="d",
                    input_schema={"type": "object"}),
      handler=_stop, risk="read")
    bundle = _bundle(
      [[TextDelta(text="partial"),
        ToolUseRequested(id="tu1", name="stopper", input={}),
        TurnDone(stop_reason="tool_use")],
       [TextDelta(text="NEVER"), TurnDone(stop_reason="end_turn")]],
      tools=child_tools)
    record = _run_child(parent, "task", bundle, cancel, tool_use_id="t1",
                        max_turns=5)
    assert record.final_text.startswith("partial"), record.final_text
    assert record.incomplete_reason == "cancelled", record.incomplete_reason

    # 2. The turn unwinds as TurnCancelled: still a record, still the text.
    script = [[TextDelta(text="partial"),
               ToolUseRequested(id="tu1", name="noop", input={}),
               TurnDone(stop_reason="tool_use")]]
    raising = _bundle([], agent=_CancellingAgent(script))
    record = _run_child(parent, "task", raising, CancelToken(),
                        tool_use_id="t2", max_turns=5)
    assert record.final_text.startswith("partial"), record.final_text
    assert record.incomplete_reason == "cancelled", record.incomplete_reason
    loaded = storage.load_subagent(parent.conv.meta.id, record.sub_id)
    assert loaded.final_text == record.final_text, loaded.final_text
    assert loaded.incomplete_reason == "cancelled", loaded.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_cancel_is_never_silently_a_clean_finish():
  """Every cancel shape must reach the parent MODEL as prose, not just code.

  The parent model reads exactly one thing about a child: the tool-result
  string. A cancelled child whose text is byte-identical to a clean finish --
  or, worse, EMPTY -- reads as "ran fine, nothing to report", which is the one
  confusion this whole feature exists to prevent. (An empty text block is also
  a provider 400.) So the record's prose carries the cancel too.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)

    # A child cancelled with NOTHING committed: its text would otherwise be ''.
    nothing = _bundle([], agent=_MidStreamCancellingAgent(
      [[TextDelta(text="half a thou")]]))
    record = _run_child(parent, "task", nothing, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    assert record.final_text.strip(), "a cancelled child returned empty text"
    assert "cancel" in record.final_text.lower(), record.final_text
    assert record.incomplete_reason == "cancelled", record.incomplete_reason

    # ...and one whose partial answer would read exactly like a clean finish.
    clean_looking = _bundle(
      [[TextDelta(text="FINDINGS: none"), TurnDone(stop_reason="end_turn")]])
    complete = _run_child(parent, "task", clean_looking, CancelToken(),
                          tool_use_id="t2", max_turns=5)
    cancelled_script = [[TextDelta(text="FINDINGS: none"),
                         ToolUseRequested(id="tu1", name="noop", input={}),
                         TurnDone(stop_reason="tool_use")]]
    cancelled = _run_child(
      parent, "task", _bundle([], agent=_CancellingAgent(cancelled_script)),
      CancelToken(), tool_use_id="t3", max_turns=5)
    assert cancelled.final_text != complete.final_text, cancelled.final_text
    assert cancelled.incomplete_reason == "cancelled"
    assert complete.incomplete_reason == ""
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_mid_stream_abort_keeps_its_committed_work():
  """An abort PARTWAY through a stream keeps every committed iteration.

  The fixture the cancel test names: the agent aborts with events already
  delivered, so the unwind escapes ``_collect_one_response`` before ``run_turn``
  can append the in-progress message. Raising between turns (the easier shape)
  never exercises that -- there the previous turn is already committed. What
  must survive is the record itself plus the earlier iterations' transcript.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    # Turn 1 completes and commits; turn 2 aborts after one delta.
    script = [[TextDelta(text="partial finding"),
               ToolUseRequested(id="tu1", name="noop", input={}),
               TurnDone(stop_reason="tool_use")],
              [TextDelta(text="never committed")]]
    bundle = _bundle([], agent=_MidStreamCancellingAgent(script))
    record = _run_child(parent, "task", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    assert record.incomplete_reason == "cancelled", record.incomplete_reason
    assert record.final_text.startswith("partial finding"), record.final_text
    assert "never committed" not in record.final_text, record.final_text
    # The committed transcript survives: user task, assistant turn 1, results.
    assert len(record.messages) >= 3, len(record.messages)
    said = [b.data.get("text", "") for m in record.messages
            if m.role == "assistant" for b in m.content if b.type == "text"]
    assert "partial finding" in said, said
    loaded = storage.load_subagent(parent.conv.meta.id, record.sub_id)
    assert loaded.incomplete_reason == "cancelled", loaded.incomplete_reason
    assert len(loaded.messages) == len(record.messages)
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_capped_claude_code_child_is_not_a_clean_one():
  """The DEFAULT backend's cap fires where AgentSession cannot see it.

  ``claude_code`` runs its agentic loop inside the SDK subprocess, so
  ``AgentSession``'s ``max_turns`` never fires for it and the transcript marker
  the API backends leave is never written. The SDK reports its own cap as
  ``ResultMessage(subtype='error_max_turns')``, which the backend maps to
  ``TurnDone(finish=CAP)`` -- an event ``_accumulate`` records only the
  ``stop_reason`` half of, and which a no-op child sink then discarded. Result,
  measured: a capped child came back byte-identical to a clean one, on the
  backend every shipped profile inherits.

  Asserted as the comparison itself, not merely as a field value: the property
  is that the two runs DIFFER, and both halves (the field code reads and the
  prose the parent model reads) have to carry it.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    same_words = "I checked the R-free flags..."
    capped = _run_child(
      parent, "review",
      # stop_reason "" is what the SDK sends with an error_max_turns result:
      # a keyed-on-stop_reason test would pass while the bug was live.
      _bundle([[TextDelta(text=same_words),
                TurnDone(stop_reason="", finish=CAP)]]),
      CancelToken(), tool_use_id="t1", max_turns=25)
    clean = _run_child(
      parent, "review",
      _bundle([[TextDelta(text=same_words),
                TurnDone(stop_reason="end_turn", finish=CLEAN)]]),
      CancelToken(), tool_use_id="t2", max_turns=25)
    assert clean.incomplete_reason == "", clean.incomplete_reason
    assert clean.final_text == same_words, clean.final_text
    assert capped.incomplete_reason == "turn_cap", capped.incomplete_reason
    assert capped.final_text != clean.final_text, capped.final_text
    assert same_words in capped.final_text, capped.final_text
    assert "cap" in capped.final_text.lower(), capped.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_every_non_clean_finish_reaches_the_parent():
  """Cap is not the only state the discarded ``finish`` was hiding.

  The same dropped event swallowed a backend ERROR result (rate-limited,
  ``error_during_execution``) and, on the API backends, the states
  ``finish_for_stop_reason`` deliberately FAILS CLOSED on: ``content_filter``,
  ``pause_turn``, a future provider value. Each must arrive under its own name
  -- a fail-closed mapper feeding a consumer that only recognizes a fixed set
  would put every unknown back into the "clean" bucket it was built to escape.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    cases = [(ERROR, "error", "error"),
             (TRUNCATED, "truncated", "output-token"),
             ("content_filter", "content_filter", "content_filter"),
             ("pause_turn", "pause_turn", "pause_turn")]
    for i, (finish, reason, phrase) in enumerate(cases):
      record = _run_child(
        parent, "review",
        _bundle([[TextDelta(text="FINDINGS: none"),
                  TurnDone(stop_reason="", finish=finish)]]),
        CancelToken(), tool_use_id="t%d" % i, max_turns=25)
      assert record.incomplete_reason == reason, (finish,
                                                  record.incomplete_reason)
      assert "FINDINGS: none" in record.final_text, record.final_text
      assert phrase in record.final_text.lower(), (finish, record.final_text)
    # ...and CANCELLED reported through the sink alone (no raise, no
    # 'cancelled' stop_reason) is still a cancel.
    record = _run_child(
      parent, "review",
      _bundle([[TextDelta(text="half"), TurnDone(stop_reason="",
                                                 finish=CANCELLED)]]),
      CancelToken(), tool_use_id="t9", max_turns=25)
    assert record.incomplete_reason == "cancelled", record.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_an_empty_child_answer_is_never_an_empty_result():
  """`""` is both a provider 400 and a lie about what happened.

  A ``claude_code`` ERROR result carries no assistant text at all, so the
  child's closing text is genuinely empty (measured). The registry path then
  emits an empty text block and the SDK path a ``_text_result("")`` -- the same
  empty block the cancel branch was already careful to avoid. Only the cancel
  branch was protected.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    # An error with nothing said: the shape actually measured.
    errored = _run_child(
      parent, "review", _bundle([[TurnDone(stop_reason="", finish=ERROR)]]),
      CancelToken(), tool_use_id="t1", max_turns=25)
    assert errored.final_text.strip(), "an errored child returned empty text"
    assert errored.incomplete_reason == "error", errored.incomplete_reason
    # ...and a CLEAN turn that simply said nothing: no reason to report, but
    # still not an empty string on the wire.
    silent = _run_child(
      parent, "review", _bundle([[TurnDone(stop_reason="end_turn",
                                           finish=CLEAN)]]),
      CancelToken(), tool_use_id="t2", max_turns=25)
    assert silent.final_text.strip(), "a silent child returned empty text"
    assert silent.incomplete_reason == "", silent.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_provider_error_is_never_a_clean_finish():
  """A turn that ends on an ``AgentError`` emits NO ``TurnDone`` at all.

  Every API backend surfaces a provider failure that way -- anthropic's
  ``APIError``/``TransportError`` arms, openai's, gemini's ``_map_error``, and
  claude_code's own input-extraction guard all ``yield AgentError`` and return.
  ``AgentSession`` then exits its loop on a ``stop_reason`` of ``''``, so the
  sink latches nothing and ``incomplete_reason`` stays ``''`` -- which
  ``SubagentRecord`` documents as "it ran to its own end".

  Measured as the comparison, not as a field: a review cut off by a rate limit
  halfway through returned a record byte-identical in shape to one that
  reviewed everything and found nothing.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    words = "R-free looks fine so far"
    errored = _run_child(
      parent, "review",
      _bundle([[TextDelta(text=words),
                AgentError(message="rate limited", recoverable=True,
                           kind="rate_limited")]]),
      CancelToken(), tool_use_id="t1", max_turns=5)
    clean = _run_child(
      parent, "review",
      _bundle([[TextDelta(text=words), TurnDone(stop_reason="end_turn",
                                                finish=CLEAN)]]),
      CancelToken(), tool_use_id="t2", max_turns=5)
    assert clean.incomplete_reason == "", clean.incomplete_reason
    assert errored.incomplete_reason == "error", errored.incomplete_reason
    # Both halves: the field code reads, and the prose the parent MODEL reads.
    assert errored.final_text != clean.final_text, errored.final_text
    assert words in errored.final_text, errored.final_text
    assert "rate limited" in errored.final_text, errored.final_text
    # ...and it survives the round trip to disk.
    loaded = storage.load_subagent(parent.conv.meta.id, errored.sub_id)
    assert loaded.incomplete_reason == "error", loaded.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_stream_that_just_stops_is_never_a_clean_finish():
  """No terminal event and no error either: still not a finished child.

  The other shape of the same gap -- a backend generator that returns without
  saying anything about how the turn ended. There is nothing to report a reason
  from, so it is reported under its own name rather than defaulting to clean:
  a table of KNOWN bad finishes cannot catch an absent one.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    record = _run_child(
      parent, "review",
      _bundle([[TextDelta(text="FINDINGS: none")]]),   # no TurnDone at all
      CancelToken(), tool_use_id="t1", max_turns=5)
    assert record.incomplete_reason == "no_terminal_event", \
        record.incomplete_reason
    assert "FINDINGS: none" in record.final_text, record.final_text
    assert "did not finish" in record.final_text.lower(), record.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_the_skill_tools_do_not_count_as_measurement():
  """A child pre-allowed its own harness has still measured nothing.

  Every child is pre-allowed ``load_skill`` / ``read_skill_file`` /
  ``list_skill_files`` so an on_demand skill can deliver its operating
  instructions. Those three read the child's OWN prompt material -- they
  measure nothing about the thing under review -- so a reviewer that loaded its
  harness and was then denied every real tool satisfied a
  "did any tool succeed" test while producing exactly the confident,
  unverified report the backstop exists to catch.

  Skill-ness is read off the registry (``source_of`` -> ``'skill'``), the same
  structural signal the factory uses to decide which names to pre-allow -- not
  a hardcoded name list that a renamed or added skill tool walks past.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    child_tools = ToolRegistry(log=null_out())
    child_tools.register_skill_tool(
      spec=ToolSpec(name="load_skill", description="d",
                    input_schema={"type": "object"}),
      handler=lambda name, input: "SKILL BODY")
    script = [[ToolUseRequested(id="tu1", name="load_skill", input={}),
               ToolUseRequested(id="tu2", name="phenix_get_results", input={}),
               TurnDone(stop_reason="tool_use")],
              [TextDelta(text="FINDINGS: the model is fine"),
               TurnDone(stop_reason="end_turn")]]
    # The production shape exactly: every shipped profile leaves
    # tool_policy_default 'ask' and the factory's _allow_skill_tools pre-allows
    # the three skill names, so the skill tool runs and every measurement tool
    # is denied by _HeadlessApprovals.
    policy = ToolPolicy(default="ask")
    policy.per_tool["load_skill"] = "allow"
    bundle = _bundle(script, tools=child_tools, policy=policy)
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    # The skill tool DID succeed ...
    ok, text = _tool_result(record.messages, "tu1")
    assert not ok and "SKILL BODY" in text, text
    # ... and that must not be mistaken for a measurement.
    assert record.incomplete_reason == "tools_denied", record.incomplete_reason
    assert "measured nothing" in record.final_text, record.final_text
    assert "phenix_get_results" in record.final_text, record.final_text

    # The same child refused by POLICY rather than by an unanswerable 'ask':
    # a 'deny' is settled before any approval is opened, so the denial record
    # was empty and the backstop had nothing to fire on.
    denying = ToolPolicy(default="deny")
    denying.per_tool["load_skill"] = "allow"
    record = _run_child(parent, "review",
                        _bundle(script, tools=child_tools, policy=denying),
                        CancelToken(), tool_use_id="t2", max_turns=5)
    assert record.incomplete_reason == "tools_denied", record.incomplete_reason
    assert "phenix_get_results" in record.final_text, record.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_one_real_tool_result_is_still_a_measurement():
  """The other direction: the flag must not fire on a child that DID measure.

  A reviewer denied one destructive tool while reading five files measured
  plenty; flagging it would make the flag meaningless. So a single non-skill
  tool result clears it -- including one a backend ran in its own subprocess
  and merely reported (``ToolResultObserved``), which never passes through the
  session's dispatch and so is not in any registry the child holds.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    child_tools = ToolRegistry(log=null_out())
    child_tools.register_builtin(
      spec=ToolSpec(name="grep_model", description="d",
                    input_schema={"type": "object"}),
      handler=lambda name, input, cancel, session, tool_use_id: "3 outliers",
      risk="read")
    bundle = _bundle(
      [[ToolUseRequested(id="tu1", name="grep_model", input={}),
        ToolUseRequested(id="tu2", name="wipe_disk", input={}),
        TurnDone(stop_reason="tool_use")],
       [TextDelta(text="FINDINGS: 3 outliers"), TurnDone(stop_reason="end_turn")]],
      tools=child_tools, policy=ToolPolicy(default="allow"))
    bundle.policy.per_tool["wipe_disk"] = "deny"
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    assert record.incomplete_reason == "", record.incomplete_reason
    assert "measured nothing" not in record.final_text, record.final_text

    # ...and the claude_code shape: the tool ran inside the SDK subprocess and
    # was only OBSERVED, so nothing about it is registry-resolvable.
    observed = _bundle(
      [[ToolUseRequested(id="tu1", name="mcp__phenix__phenix_get_results",
                         input={}),
        ToolResultObserved(tool_use_id="tu1", content="R-free 0.21",
                           name="mcp__phenix__phenix_get_results"),
        TextDelta(text="FINDINGS: none"), TurnDone(stop_reason="end_turn")]],
      tools=ToolRegistry(log=null_out()))
    record = _run_child(parent, "review", observed, CancelToken(),
                        tool_use_id="t2", max_turns=5)
    assert record.incomplete_reason == "", record.incomplete_reason

    # The same shape from the phenix factory, which DOES say what the child was
    # equipped with -- so the flag is armed here, and must still not fire. The
    # server the SDK named it under is one the profile asked for, which is the
    # whole difference between this and a bookkeeping result.
    equipped = _bundle(
      [[ToolUseRequested(id="tu1", name="mcp__phenix__phenix_get_results",
                         input={}),
        ToolResultObserved(tool_use_id="tu1", content="R-free 0.21",
                           name="mcp__phenix__phenix_get_results"),
        TextDelta(text="FINDINGS: none"), TurnDone(stop_reason="end_turn")]],
      measurement_servers=["phenix"])
    record = _run_child(parent, "review", equipped, CancelToken(),
                        tool_use_id="t3", max_turns=5)
    assert record.incomplete_reason == "", record.incomplete_reason
    assert "measured nothing" not in record.final_text, record.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_an_equipped_child_that_measured_nothing_is_never_clean():
  """The backstop turns on the OUTCOME, not on how the outcome came about.

  It used to be gated on ``approvals.denied_tools``, which made it unreachable
  for the two commonest shapes -- neither of which records a denial:

  - every call came back as an ERROR (an unknown tool, a server that never
    started, a backing program that failed). Nothing was refused, so nothing
    reached the denial record and the check behind it was never consulted;
  - the child never called a tool AT ALL. It was equipped, it decided not to
    measure, and it wrote its report anyway.

  Both produced ``incomplete_reason=''`` and no marker, i.e. a record whose
  reader concludes the child ran to its own end and found nothing wrong. The
  gate is now "did anything it could measure with come back", and the denial
  record is kept only to say WHICH tools were refused.

  Read on the two routes a child's equipment arrives by, because each is
  invisible from the other: an API child's MCP tools are registry entries
  sourced ``mcp:<server>``, while a ``claude_code`` child's live in its SDK
  subprocess and appear only as the bundle's ``measurement_servers``.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    # 1. Errored, on the claude_code route: allowed by policy, dispatched, and
    #    answered with an error. No denial anywhere.
    errored = _bundle(
      [[ToolUseRequested(id="tu1", name="mcp__phenix__phenix_get_results",
                         input={}),
        ToolResultObserved(tool_use_id="tu1", name="mcp__phenix__"
                           "phenix_get_results",
                           content="No such tool available", is_error=True),
        TextDelta(text="FINDINGS: the model looks fine"),
        TurnDone(stop_reason="end_turn")]],
      measurement_servers=["phenix"])
    record = _run_child(parent, "review", errored, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    assert record.incomplete_reason == "no_measurement", \
        record.incomplete_reason
    assert "measured nothing" in record.final_text, record.final_text
    assert "the model looks fine" in record.final_text, record.final_text

    # 2. Never called a tool at all, on the API route: equipment read off the
    #    registry, so this holds for a bundle that names no servers.
    child_tools = ToolRegistry(log=null_out())
    child_tools.register_mcp_tool(
      spec=ToolSpec(name="phenix_get_results", description="d",
                    input_schema={"type": "object"}),
      server_name="phenix", handler=lambda name, input, cancel: "R-free 0.21")
    silent = _bundle([[TextDelta(text="FINDINGS: the model looks fine"),
                       TurnDone(stop_reason="end_turn")]],
                     tools=child_tools)
    record = _run_child(parent, "review", silent, CancelToken(),
                        tool_use_id="t2", max_turns=5)
    assert record.incomplete_reason == "no_measurement", \
        record.incomplete_reason

    # 3. The same silent child with NO inventory on either route. Nothing told
    #    this module what it held, so it says nothing -- the pre-existing
    #    behavior, and the only honest answer.
    unknown = _bundle([[TextDelta(text="FINDINGS: the model looks fine"),
                        TurnDone(stop_reason="end_turn")]])
    record = _run_child(parent, "review", unknown, CancelToken(),
                        tool_use_id="t3", max_turns=5)
    assert record.incomplete_reason == "", record.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_backends_own_bookkeeping_server_is_not_a_measurement():
  """A tool the child's PROFILE never asked for measures nothing for it.

  Every backend wires in servers of its own on top of the profile's.
  ``claude_code`` builds an in-process ``phenix_chat`` server on EVERY agent,
  carrying ``phenix_get_job_history`` -- auto-approved ahead of every policy
  check, and (with no provider wired, which is every child: only the chat
  window ever wires one) returning ``"(no jobs recorded)"`` as a SUCCESSFUL
  result. That satisfied "did any tool succeed" for a child that had been
  denied every real tool and measured nothing at all.

  Decided STRUCTURALLY -- is this tool's server one the profile asked for? --
  and not by a list of tool names, which lives in the other repo and would have
  to be kept in step with it forever. So a bookkeeping tool added to that
  server tomorrow is excluded the day it is added, and the same rule covers
  every backend's housekeeping server without naming any of them.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    bookkeeping = "mcp__phenix_chat__phenix_get_job_history"
    measure = "mcp__phenix__phenix_get_results"
    bundle = _bundle(
      [[ToolUseRequested(id="tu1", name=bookkeeping, input={}),
        ToolResultObserved(tool_use_id="tu1", name=bookkeeping,
                           content="(no jobs recorded)", is_error=False),
        ToolUseRequested(id="tu2", name=measure, input={}),
        ToolResultObserved(tool_use_id="tu2", name=measure,
                           content="Denied by tool policy.", is_error=True),
        TextDelta(text="FINDINGS: the model looks fine"),
        TurnDone(stop_reason="end_turn")]],
      measurement_servers=["phenix"])
    record = _run_child(parent, "review", bundle, CancelToken(),
                        tool_use_id="t1", max_turns=5)
    assert record.incomplete_reason == "no_measurement", \
        record.incomplete_reason
    assert "measured nothing" in record.final_text, record.final_text
    # ... and the profile's OWN server, named the same way, still counts.
    ok = _bundle(
      [[ToolUseRequested(id="tu1", name=measure, input={}),
        ToolResultObserved(tool_use_id="tu1", name=measure,
                           content="R-free 0.21", is_error=False),
        TextDelta(text="FINDINGS: none"), TurnDone(stop_reason="end_turn")]],
      measurement_servers=["phenix"])
    record = _run_child(parent, "review", ok, CancelToken(),
                        tool_use_id="t2", max_turns=5)
    assert record.incomplete_reason == "", record.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_backend_side_denial_reaches_the_record():
  """A backend that denies tools ITSELF must still feed the denial record.

  ``claude_code`` resolves every permission inside its own SDK callback: a
  child's tool calls never reach ``_HeadlessApprovals``, so the whole
  denied-every-tool backstop was dead on the DEFAULT backend -- the one every
  shipped profile inherits. ``run_child`` therefore hands such an agent a sink
  to report into, and the record is built from the union.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)

    class _SelfDenyingAgent(FakeAgent):
      """Denies in its own permission layer, as ClaudeCodeAgent does."""

      def __init__(self, script):
        FakeAgent.__init__(self, script)
        self.sink = None

      def set_denial_sink(self, sink):
        self.sink = sink

      def stream_turn(self, conversation, tools, cancel):
        if self.sink is not None:
          self.sink("mcp__phenix__phenix_get_results")
        return FakeAgent.stream_turn(self, conversation, tools, cancel)

    agent = _SelfDenyingAgent([[TextDelta(text="FINDINGS: none"),
                               TurnDone(stop_reason="end_turn")]])
    record = _run_child(parent, "review", _bundle([], agent=agent),
                        CancelToken(), tool_use_id="t1", max_turns=5)
    assert agent.sink is not None, "run_child never wired the denial sink"
    assert record.incomplete_reason == "tools_denied", record.incomplete_reason
    assert "phenix_get_results" in record.final_text, record.final_text
    assert "FINDINGS: none" in record.final_text, record.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_the_child_releases_its_agent_and_its_servers():
  """A child that SUCCEEDS must release what was opened for it.

  ``build_child``'s own guard covers a FAILED build only, so a successful spawn
  kept one provider client / ``claude`` CLI subprocess (plus that subprocess's
  own MCP servers and a daemon loop thread) and one server subprocess per
  session-side connection -- for the life of the parent session, once per
  spawn, with no finalizer to collect them. The declared consumer is a review
  agent invoked repeatedly in a long-lived session.

  Release is deterministic and covers the FAILURE paths too: a crashed child
  raises past this point, and a release that only ran on the happy path would
  leak exactly when things are going wrong.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    agent = _ClosableAgent([[TextDelta(text="done"),
                             TurnDone(stop_reason="end_turn")]])
    conns = [_FakeConnection("a"), _FakeConnection("b")]
    record = _run_child(
      parent, "task", _bundle([], agent=agent, connections=conns),
      CancelToken(), tool_use_id="t1", max_turns=5)
    assert record.incomplete_reason == "", record.incomplete_reason
    assert agent.closed == 1, agent.closed
    assert [c.stopped for c in conns] == [1, 1], conns

    # A crash keeps propagating (that contract is exercised elsewhere) AND
    # still releases.
    crash_agent = _CrashingAgent([])
    crash_agent.closed = 0
    crash_agent.close = lambda: setattr(crash_agent, "closed",
                                        crash_agent.closed + 1)
    crash_conns = [_FakeConnection("c")]
    try:
      _run_child(parent, "task",
                 _bundle([], agent=crash_agent, connections=crash_conns),
                 CancelToken(), tool_use_id="t2", max_turns=5)
    except RuntimeError:
      pass
    else:
      raise AssertionError("the crash stopped propagating")
    assert crash_agent.closed == 1, crash_agent.closed
    assert crash_conns[0].stopped == 1, crash_conns[0].stopped

    # One server that cannot be stopped must not strand the ones after it, nor
    # turn a finished child into a failed tool call.
    ok_after = _FakeConnection("after")
    noisy = _ClosableAgent([[TextDelta(text="done"),
                             TurnDone(stop_reason="end_turn")]])
    record = _run_child(
      parent, "task",
      _bundle([], agent=noisy,
              connections=[_UnstoppableConnection("bad"), ok_after]),
      CancelToken(), tool_use_id="t3", max_turns=5)
    assert record.final_text == "done", record.final_text
    assert ok_after.stopped == 1, ok_after.stopped
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_build_notice_leads_the_final_text():
  """How the child was really built has to reach the party that reads it.

  ``build_child`` falls back to the parent's own backend/model when the
  requested one has no usable credentials -- the right call, since a review on
  the parent's model beats no review. But the parent MODEL reads exactly one
  thing about a child, ``final_text``, so a substitution reported only to the
  session log returns a review the parent commissioned from an INDEPENDENT
  second model, unattributed, from the parent's own. That is this feature's own
  argument for the cancel marker, in the one case where independence is the
  entire point.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    notice = "[Subagent note: ran on claude_code/opus, not google/gemini]"
    record = _run_child(
      parent, "review",
      _bundle([[TextDelta(text="FINDINGS: none"),
                TurnDone(stop_reason="end_turn")]], notice=notice),
      CancelToken(), tool_use_id="t1", max_turns=5)
    # FIRST: the caveat has to precede the findings it qualifies.
    assert record.final_text.startswith(notice), record.final_text
    assert "FINDINGS: none" in record.final_text, record.final_text
    # A substitution is not an incomplete run -- the child did finish.
    assert record.incomplete_reason == "", record.incomplete_reason
    # ...and a child that also said nothing still reports the no-text marker,
    # rather than letting the notice stand in for an answer.
    silent = _run_child(
      parent, "review",
      _bundle([[TurnDone(stop_reason="end_turn")]], notice=notice),
      CancelToken(), tool_use_id="t2", max_turns=5)
    assert silent.final_text.startswith(notice), silent.final_text
    assert len(silent.final_text) > len(notice), silent.final_text
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_child_usage_lands_on_the_record_and_not_on_the_parent():
  """The child's spend belongs to the record, and only to the record.

  It was once filed a second time on the parent session
  (``add_subagent_usage``), which nothing read -- two sinks for one number,
  and the redundant one survived because this test pinned it. The record is
  the persisted home, so that is what is asserted.

  The id still matters: ``AgentSession`` mints its own ``sub_id`` at depth>0,
  and a record carrying a SECOND id would name a child nothing can join back
  to the transcript the GUI shows.
  """
  from qttbx.widgets.chat.agent import subagent as subagent_mod
  tmp = tempfile.mkdtemp()
  sessions = []

  class _RecordingSession(AgentSession):
    """Captures the child session so its OWN sub_id is observable."""

    def __init__(self, *args, **kwargs):
      AgentSession.__init__(self, *args, **kwargs)
      sessions.append(self)

  real_session = subagent_mod.AgentSession
  try:
    parent, storage = _parent(tmp)
    subagent_mod.AgentSession = _RecordingSession
    try:
      record = _run_child(
        parent, "task",
        _bundle([[TextDelta(text="x"),
                  TokenUsageEvent(input=10, output=2, context_tokens=900000),
                  TurnDone(stop_reason="end_turn")]]),
        CancelToken(), tool_use_id="t1", max_turns=5)
    finally:
      subagent_mod.AgentSession = real_session
    assert len(sessions) == 1, len(sessions)
    # ONE id for one child. The session mints its own at depth>0 and the
    # rollup is keyed by it, so a record minting a SECOND files the usage
    # under a name nothing can join back to the record the GUI shows.
    assert record.sub_id == sessions[0].sub_id, (record.sub_id,
                                                 sessions[0].sub_id)
    usage = record.token_usage
    assert usage.input == 10 and usage.output == 2, usage
    # Filed UNDER the child, not folded into the parent's own peak: the parent
    # agent's context reading must stay its own (a false handoff nudge is what
    # the child's separate agent exists to prevent).
    assert usage.context_tokens == 900000, usage
    assert [m.usage for m in parent.conv.messages
            if m.usage is not None] == [], "child usage landed on the parent"
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_failed_save_keeps_the_child_result():
  """A disk problem must not destroy a turn's work.

  Every sibling save in this subsystem is try/except + log for this reason.
  Unguarded, a full disk / read-only project converts a child that finished
  into a bare tool error AND throws away its text -- the expensive half of the
  loss, since the record is reconstructible and the child's minutes of work are
  not.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)

    def _explode(conv_id, record):
      raise OSError("[Errno 28] No space left on device")

    storage.store_subagent = _explode
    record = _run_child(
      parent, "task",
      _bundle([[TextDelta(text="FINDINGS: none"),
                TurnDone(stop_reason="end_turn")]]),
      CancelToken(), tool_use_id="t1", max_turns=5)
    assert record.final_text == "FINDINGS: none", record.final_text
    assert record.incomplete_reason == "", record.incomplete_reason
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_a_blank_task_is_refused_on_both_tool_paths():
  """``"   "`` is an empty task, whichever backend family asks.

  claude_code's SDK tool strips before its emptiness check; the registry
  builtin did not, so a whitespace-only task was refused on one backend family
  and spawned a child with no instructions on the other. Same tool, same
  refusal.
  """
  tmp = tempfile.mkdtemp()
  try:
    parent, storage = _parent(tmp)
    spawned = []

    def _build(profile_name, model, backend):
      spawned.append(True)
      return _bundle([[TextDelta(text="SPAWNED"),
                       TurnDone(stop_reason="end_turn")]])

    register_run_subagent(parent.tools, _build)
    for blank in ("   ", "\t", "\n "):
      try:
        parent.tools.invoke_builtin(
          "run_subagent", {"task": blank}, cancel=CancelToken(),
          session=parent, tool_use_id="t1")
      except Sorry as exc:
        assert "task" in str(exc).lower(), str(exc)
      else:
        raise AssertionError("a blank task spawned a child: %r" % (blank,))
    assert spawned == [], spawned
    # ...and the surrounding whitespace of a REAL task is trimmed, not sent:
    # what run_child receives is the child's ENTIRE prompt.
    from qttbx.widgets.chat.agent import subagent as subagent_mod
    seen = []
    real_run_child = subagent_mod.run_child

    def _capturing(parent_session, task, *args, **kwargs):
      seen.append(task)
      return real_run_child(parent_session, task, *args, **kwargs)

    subagent_mod.run_child = _capturing
    try:
      result = parent.tools.invoke_builtin(
        "run_subagent", {"task": "  review this  "}, cancel=CancelToken(),
        session=parent, tool_use_id="t2")
    finally:
      subagent_mod.run_child = real_run_child
    assert result == "SPAWNED", result
    assert seen == ["review this"], seen
  finally:
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_the_tool_description_is_a_shared_constant():
  """One tool, two backend surfaces, one wording.

  ``ClaudeCodeAgent`` registers ``run_subagent`` on the SDK from its own copy
  of the description, so a verbatim duplicate drifts silently -- and this
  wording is what tells the model the child starts blind and cannot ask, i.e.
  what it needs to write a self-contained task. Exported so the phenix side
  reads the same object rather than a matching string.
  """
  registry = ToolRegistry(log=null_out())
  register_run_subagent(registry, lambda *a: None)
  spec = [s for s in registry.specs() if s.name == "run_subagent"][0]
  assert spec.description is RUN_SUBAGENT_DESCRIPTION, spec.description
  assert "no user" in RUN_SUBAGENT_DESCRIPTION, RUN_SUBAGENT_DESCRIPTION


def exercise():
  exercise_child_runs_and_returns_final_text()
  exercise_record_is_persisted_despite_depth_gating()
  exercise_child_usage_never_touches_the_parent()
  exercise_costs_sum_but_context_peaks()
  exercise_child_depth_is_derived_from_the_parent()
  exercise_child_gets_its_own_event_sink()
  exercise_an_ask_policy_is_denied_not_parked()
  exercise_the_destructive_floor_is_denied_not_parked()
  exercise_tool_returns_the_child_final_text()
  exercise_tool_passes_the_parent_turn_cap()
  exercise_tool_refuses_when_disabled_or_taskless()
  exercise_direct_call_falls_back_to_25_turns()
  exercise_the_childs_own_cap_outranks_the_callers()
  exercise_turn_cap_marks_the_record_incomplete()
  exercise_a_quoted_cap_marker_is_not_a_cap()
  exercise_a_capped_child_still_reports_what_it_found()
  exercise_a_crashed_child_keeps_its_transcript()
  exercise_the_child_transcript_names_who_wrote_it()
  exercise_a_child_image_lands_under_the_parent_conversation()
  exercise_the_final_text_walk_back_is_bounded()
  exercise_a_provider_executed_tool_is_a_measurement()
  exercise_a_late_provider_failure_is_not_a_turn_cap()
  exercise_a_denial_survives_another_reason_for_stopping()
  exercise_one_denial_is_recorded_once()
  exercise_a_child_reads_back_the_images_it_was_filed()
  exercise_a_non_string_task_is_refused_not_crashed()
  exercise_the_two_turn_cap_keys_mean_one_thing_each()
  exercise_an_unusable_turn_cap_is_refused_before_the_child_runs()
  exercise_a_child_crash_propagates()
  exercise_cancel_returns_partial_text_not_an_exception()
  exercise_cancel_is_never_silently_a_clean_finish()
  exercise_a_mid_stream_abort_keeps_its_committed_work()
  exercise_a_capped_claude_code_child_is_not_a_clean_one()
  exercise_every_non_clean_finish_reaches_the_parent()
  exercise_an_empty_child_answer_is_never_an_empty_result()
  exercise_a_provider_error_is_never_a_clean_finish()
  exercise_a_stream_that_just_stops_is_never_a_clean_finish()
  exercise_the_skill_tools_do_not_count_as_measurement()
  exercise_one_real_tool_result_is_still_a_measurement()
  exercise_an_equipped_child_that_measured_nothing_is_never_clean()
  exercise_a_backends_own_bookkeeping_server_is_not_a_measurement()
  exercise_a_backend_side_denial_reaches_the_record()
  exercise_the_child_releases_its_agent_and_its_servers()
  exercise_a_build_notice_leads_the_final_text()
  exercise_child_usage_lands_on_the_record_and_not_on_the_parent()
  exercise_a_failed_save_keeps_the_child_result()
  exercise_a_blank_task_is_refused_on_both_tool_paths()
  exercise_the_tool_description_is_a_shared_constant()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
