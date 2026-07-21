import io
import json
import shutil
import tempfile
import threading
from collections import namedtuple
from pathlib import Path

from libtbx.utils import format_cpu_times, null_out
from qttbx.widgets.chat.agent.base import (
  Agent, ToolSpec)
from qttbx.widgets.chat.agent.conversation import (
  ContentBlock, Conversation, Message, TokenUsage, now)
from qttbx.widgets.chat.agent.errors import CancelToken, TurnCancelled
from qttbx.widgets.chat.agent.events import (
  AskUserQuestionRequested, ServerToolResult, ServerToolUsed, TextDelta,
  ToolResultObserved, ToolUseRequested, TurnDone,
  TokenUsage as TokenUsageEvent)
from qttbx.widgets.chat.agent.profile import Profile
from qttbx.widgets.chat.agent.session import AgentSession
from qttbx.widgets.chat.agent.storage import ConversationStorage
from qttbx.widgets.chat.agent.tools import (
  ToolPolicy, ToolRegistry, ToolApprovalResponse, _Cancelled)


# ---- shared FakeAgent ------------------------------------------------------

class FakeAgent(Agent):
  """Replays a scripted sequence of event lists. Each call to stream_turn
  yields events from the next list; raises if exhausted."""
  name = "fake"
  model = "fake"

  def __init__(self, turn_scripts):
    self.turn_scripts = list(turn_scripts)
    self._cursor = 0

  def stream_turn(self, conversation, tools, cancel):
    if self._cursor >= len(self.turn_scripts):
      raise RuntimeError("FakeAgent: ran out of scripted turns")
    events = self.turn_scripts[self._cursor]
    self._cursor += 1
    for ev in events:
      yield ev

  def resolve_credentials(self, cli_override=None):
    return "fake-key"

  def credentials_dialog_class(self):
    return object


def _new_test_session(turn_scripts, project_dir=None, max_depth=1, agent=None,
                      clock=None, autosave_interval_s=5.0, depth=0):
  """Build a session with an in-memory queue and minimal profile."""
  if project_dir is None:
    project_dir = tempfile.mkdtemp()
  storage = ConversationStorage(Path(project_dir), log=null_out())
  conv = Conversation.new(profile_name="fake", model="fake")
  storage.save(conv)
  registry = ToolRegistry(log=null_out())
  policy = ToolPolicy(default="allow")
  profile = Profile(name="fake", model="fake",
                    subagents_max_depth=max_depth)
  session = AgentSession(
    agent=agent if agent is not None else FakeAgent(turn_scripts),
    conversation=conv,
    storage=storage,
    tools=registry,
    policy=policy,
    profile=profile,
    depth=depth,
    log=null_out(),
    clock=clock,
    autosave_interval_s=autosave_interval_s)
  return session, project_dir


def _answer_surfaced_approvals(session, decision="deny"):
  """Answer each ToolApprovalRequest the session surfaces from a SEPARATE
  thread, modelling production: the worker thread registers + emits the card and
  then blocks in _await_approval's fut.result() WITHOUT holding the lock, while
  the GUI thread submits the decision. on_event now runs UNDER the coordinator
  lock, so a test must never submit inline from it (that re-enters the lock on
  the worker thread and deadlocks) -- it hands the request to the answering
  thread instead.

  Returns (captured, on_event, stop): install on_event as session.on_event,
  inspect captured after the calls, and call stop() to join the thread.
  """
  import queue as _queue
  import threading as _threading
  from qttbx.widgets.chat.agent.tools import (
    ToolApprovalRequest, ToolApprovalResponse)
  captured = []
  surfaced = _queue.Queue()
  done = _threading.Event()

  def _on_event(ev):
    captured.append(ev)
    if isinstance(ev, ToolApprovalRequest):
      surfaced.put(ev)

  def _answer():
    while not done.is_set():
      try:
        req = surfaced.get(timeout=0.05)
      except _queue.Empty:
        continue
      # Carry the surfaced id (the real card echoes req.request_id) so
      # _await_approval matches it instead of dropping it as a stale click.
      session.approvals.submit(ToolApprovalResponse(
        request_id=req.request_id, decision=decision))

  th = _threading.Thread(target=_answer, daemon=True)
  th.start()

  def _stop():
    done.set()
    th.join(timeout=2)

  return captured, _on_event, _stop


# ---- tests -----------------------------------------------------------------

def exercise_simple_text_turn():
  """One model turn with text only, stop_reason='end_turn'."""
  session, tmp = _new_test_session([
    [TextDelta(text="Hello"),
     TextDelta(text=", world."),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    cancel = CancelToken()
    user_msg = Message(role="user",
                       content=[ContentBlock(type="text", data={"text": "hi"})],
                       timestamp=now())
    assistant = session.run_turn(user_msg, cancel)
    assert assistant.stop_reason == "end_turn"
    assert assistant.content[0].type == "text"
    assert "Hello, world." in assistant.content[0].data["text"]
  finally:
    shutil.rmtree(tmp)


def exercise_assistant_messages_stamped_with_model_and_backend():
  """Each assistant message records the model (from the agent) and the
  backend (from the profile) that produced it, so a reloaded conversation
  shows what generated each turn."""
  session, tmp = _new_test_session([
    [TextDelta(text="hi"), TurnDone(stop_reason="end_turn")],
  ])
  try:
    session.profile.backend = "anthropic"          # distinct from default
    cancel = CancelToken()
    user_msg = Message(role="user",
                       content=[ContentBlock(type="text", data={"text": "hi"})],
                       timestamp=now())
    session.run_turn(user_msg, cancel)
    assistants = [m for m in session.conv.messages if m.role == "assistant"]
    assert assistants, "expected an assistant message"
    for m in assistants:
      assert m.model == "fake", m.model            # FakeAgent.model
      assert m.backend == "anthropic", m.backend   # profile.backend
  finally:
    shutil.rmtree(tmp)


def exercise_reconciles_meta_model_and_backend_on_continue():
  """Continuing a conversation under a different model/backend updates the
  conversation meta to the active one (the per-message stamp preserves the
  full per-turn history)."""
  session, tmp = _new_test_session([
    [TextDelta(text="hi"), TurnDone(stop_reason="end_turn")],
  ])
  try:
    # A conversation created earlier under a different model/backend.
    session.conv.meta.model = "old-model"
    session.conv.meta.backend = "old-backend"
    session.profile.backend = "anthropic"          # the active backend now
    cancel = CancelToken()
    user_msg = Message(role="user",
                       content=[ContentBlock(type="text", data={"text": "hi"})],
                       timestamp=now())
    session.run_turn(user_msg, cancel)
    assert session.conv.meta.model == "fake", session.conv.meta.model
    assert session.conv.meta.backend == "anthropic", session.conv.meta.backend
  finally:
    shutil.rmtree(tmp)


def exercise_tool_use_loop_completes():
  """Turn 1 emits a tool_use; we register an allow-policy builtin tool;
  session dispatches it; turn 2 emits final text."""
  session, tmp = _new_test_session([
    # Turn 1: text + tool_use, stop_reason='tool_use'
    [TextDelta(text="Calling tool..."),
     ToolUseRequested(id="t1", name="echo", input={"text": "hi back"}),
     TurnDone(stop_reason="tool_use")],
    # Turn 2: final text, stop_reason='end_turn'
    [TextDelta(text="Done."),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    # Register echo tool with allow policy (default already 'allow').
    session.tools.register_builtin(
      ToolSpec(name="echo", description="echo",
               input_schema={"type": "object"}),
      handler=lambda name, input, cancel, session, tool_use_id: input["text"],
      risk="write")
    cancel = CancelToken()
    user_msg = Message(role="user",
                       content=[ContentBlock(type="text", data={"text": "hi"})],
                       timestamp=now())
    final = session.run_turn(user_msg, cancel)
    assert final.stop_reason == "end_turn"
    # Conversation should have: user, assistant(turn1), user(tool_result), assistant(turn2)
    roles = [m.role for m in session.conv.messages]
    assert roles == ["user", "assistant", "user", "assistant"]
    # tool_result message contains the dispatched tool's result
    tr_msg = session.conv.messages[2]
    assert tr_msg.content[0].type == "tool_result"
    assert tr_msg.content[0].data["tool_use_id"] == "t1"
    assert tr_msg.content[0].data["is_error"] is False
  finally:
    shutil.rmtree(tmp)


def exercise_cancel_after_tool_use_does_not_orphan_tool_use():
  """[Major] Stop hit after the model streamed a tool_use but BEFORE dispatch
  must not leave an unmatched tool_use in the saved transcript: an assistant
  tool_use with no following tool_result is a non-recoverable provider 400 on
  the next turn's replay, which permanently wedges the conversation
  (anthropic/openai/portkey/gemini). The cancel early-return must answer the
  pending tool_uses with a tool_result."""
  cancel = CancelToken()
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="echo", input={"text": "hi"}),
     TurnDone(stop_reason="tool_use")]])
  orig_stream = session.agent.stream_turn

  def _stream(conversation, tools, c):
    for ev in orig_stream(conversation, tools, c):
      yield ev
    cancel.set()                 # Stop pressed just as the tool_use finished
  session.agent.stream_turn = _stream
  try:
    user_msg = Message(role="user",
                       content=[ContentBlock(type="text", data={"text": "go"})],
                       timestamp=now())
    assistant = session.run_turn(user_msg, cancel)
    assert assistant.stop_reason == "cancelled", assistant.stop_reason
    # The saved transcript must NOT orphan the tool_use: the final message is a
    # user tool_result answering t1.
    last = session.conv.messages[-1]
    assert last.role == "user", [m.role for m in session.conv.messages]
    assert any(b.type == "tool_result" and b.data.get("tool_use_id") == "t1"
               for b in last.content), [(b.type, b.data) for b in last.content]
  finally:
    shutil.rmtree(tmp)


def exercise_cancel_during_streaming_does_not_orphan_tool_use():
  """[Major] The REALISTIC Stop-during-streaming case: the backend (Anthropic)
  emits a COMPLETED tool_use block and THEN a TurnDone(stop_reason='cancelled')
  -- so the assistant message carries a pending tool_use but its stop_reason is
  'cancelled', not 'tool_use'. run_turn's `stop_reason != 'tool_use'` early
  return must STILL answer the pending tool_use, or it orphans the transcript
  -> non-recoverable replay 400. (A previous fix only covered the narrower
  stop_reason=='tool_use'-then-cancel window.)"""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="echo", input={"text": "hi"}),
     TurnDone(stop_reason="cancelled")]])
  try:
    cancel = CancelToken()
    user_msg = Message(role="user",
                       content=[ContentBlock(type="text", data={"text": "go"})],
                       timestamp=now())
    assistant = session.run_turn(user_msg, cancel)
    assert assistant.stop_reason == "cancelled", assistant.stop_reason
    last = session.conv.messages[-1]
    assert last.role == "user", [m.role for m in session.conv.messages]
    assert any(b.type == "tool_result" and b.data.get("tool_use_id") == "t1"
               for b in last.content), [(b.type, b.data) for b in last.content]
  finally:
    shutil.rmtree(tmp)


def exercise_cancel_during_multi_tool_streaming_answers_all():
  """The window widens for parallel/multi-tool turns: Stop after several
  tool_use blocks streamed (then a cancelled TurnDone) must answer EVERY
  pending tool_use, not just one -- an unmatched tool_use of any of them
  orphans the transcript."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="echo", input={"text": "a"}),
     ToolUseRequested(id="t2", name="echo", input={"text": "b"}),
     ToolUseRequested(id="t3", name="echo", input={"text": "c"}),
     TurnDone(stop_reason="cancelled")]])
  try:
    cancel = CancelToken()
    user_msg = Message(role="user",
                       content=[ContentBlock(type="text", data={"text": "go"})],
                       timestamp=now())
    session.run_turn(user_msg, cancel)
    answered = set()
    for m in session.conv.messages:
      for b in m.content:
        if b.type == "tool_result":
          answered.add(b.data.get("tool_use_id"))
    assert answered == {"t1", "t2", "t3"}, answered
  finally:
    shutil.rmtree(tmp)


def exercise_builtin_handler_receives_tool_use_id():
  """A builtin tool handler is invoked with tool_use_id == the call id,
  passed directly through invoke_builtin (no session-held state needed --
  this is why AgentSession keeps no current_tool_use_id attribute)."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="echo", input={"text": "hi"}),
     TurnDone(stop_reason="tool_use")],
    [TextDelta(text="done"), TurnDone(stop_reason="end_turn")],
  ])
  try:
    seen = []

    def _handler(name, input, cancel, session, tool_use_id):
      seen.append(tool_use_id)
      return input["text"]

    session.tools.register_builtin(
      ToolSpec(name="echo", description="echo",
               input_schema={"type": "object"}),
      handler=_handler, risk="write")
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "hi"})], timestamp=now())
    session.run_turn(user_msg, cancel)
    assert seen == ["t1"], seen
  finally:
    shutil.rmtree(tmp)


def exercise_tool_denied_returns_is_error():
  """Policy='deny' for the tool → tool_result is_error=True; the model's
  next turn sees the denial."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="echo", input={}),
     TurnDone(stop_reason="tool_use")],
    [TextDelta(text="OK, denied."),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    session.policy = ToolPolicy(default="deny")
    session.tools.register_builtin(
      ToolSpec(name="echo", description="echo",
               input_schema={"type": "object"}),
      handler=lambda **kw: "should-not-run",
      risk="write")
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "hi"})], timestamp=now())
    session.run_turn(user_msg, cancel)
    tr_msg = session.conv.messages[2]
    assert tr_msg.content[0].data["is_error"] is True
    assert "denied" in tr_msg.content[0].data["content"][0].data["text"].lower()
  finally:
    shutil.rmtree(tmp)


def exercise_cancel_while_awaiting_approval():
  """The deadlock-fix path: worker parked on the coordinator future; the
  runner cancels the turn via cancel_turn(); the worker raises
  TurnCancelled, the turn ends."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="echo", input={}),
     TurnDone(stop_reason="tool_use")],
  ])
  try:
    session.policy = ToolPolicy(default="ask")
    session.tools.register_builtin(
      ToolSpec(name="echo", description="echo",
               input_schema={"type": "object"}),
      handler=lambda **kw: "x", risk="write")

    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "hi"})], timestamp=now())

    # Track when an approval has been requested
    requested = threading.Event()
    original_await = session._await_approval

    def _await_with_signal(req):
      requested.set()
      return original_await(req)

    session._await_approval = _await_with_signal

    # Run the turn in a worker thread.
    result_box = []
    error_box = []

    def _run():
      try:
        result_box.append(session.run_turn(user_msg, cancel))
      except Exception as e:
        error_box.append(e)

    t = threading.Thread(target=_run)
    t.start()

    # Wait until the worker is parked on the approval queue.
    assert requested.wait(timeout=5.0), "approval was never requested"

    # GUI side: set cancel + push sentinel (Section 10.4 mechanism).
    cancel.set()
    session.approvals.cancel_turn()          # was: approval_queue.put(_Cancelled())

    t.join(timeout=5.0)
    assert not t.is_alive(), "worker did not unblock — deadlock"
    # No exception escaped (TurnCancelled is caught inside dispatch).
    assert not error_box, "unexpected exception: %s" % error_box
    # Turn unwound: tool_result has is_error=True with cancellation message.
    tr_msg = session.conv.messages[2]
    assert tr_msg.content[0].data["is_error"] is True
    assert "cancel" in tr_msg.content[0].data["content"][0].data["text"].lower()
  finally:
    shutil.rmtree(tmp)


def exercise_max_turns_cap():
  """run_turn with max_turns=2 should stop and append the cap marker even
  if the model keeps emitting tool_use."""
  # FakeAgent will yield tool_use repeatedly; we cap at 2 iterations.
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="echo", input={}),
     TurnDone(stop_reason="tool_use")],
    [ToolUseRequested(id="t2", name="echo", input={}),
     TurnDone(stop_reason="tool_use")],
    # Without the cap, the session would request a 3rd turn here.
  ])
  try:
    session.tools.register_builtin(
      ToolSpec(name="echo", description="echo",
               input_schema={"type": "object"}),
      handler=lambda **kw: "ok", risk="write")
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "loop"})], timestamp=now())
    session.run_turn(user_msg, cancel, max_turns=2)
    # Last message in the conversation should be the cap marker.
    last = session.conv.messages[-1]
    assert last.role == "user"
    text = last.content[0].data["text"]
    assert "Subagent stopped at turn cap" in text
  finally:
    shutil.rmtree(tmp)


def exercise_add_subagent_usage_aggregates():
  session, tmp = _new_test_session([[TurnDone(stop_reason="end_turn")]])
  try:
    session.add_subagent_usage("sa_1",
                               TokenUsage(input=100, output=50))
    session.add_subagent_usage("sa_2",
                               TokenUsage(input=200, output=70))
    assert session._subagent_usage_by_id["sa_1"].input == 100
    assert session._subagent_usage_by_id["sa_2"].output == 70
  finally:
    shutil.rmtree(tmp)


def exercise_token_usage_event_carries_all_fields_to_message():
  """A TokenUsage event accumulates into the assistant message's usage as a
  stored TokenUsage with every field preserved (the event -> stored
  conversion drops nothing). Dataclass equality compares all fields, so a
  usage field added later is covered automatically."""
  session, tmp = _new_test_session([
    [TextDelta(text="hi"),
     TokenUsageEvent(input=11, output=22, cache_read=33, cache_creation=44),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "hi"})], timestamp=now())
    assistant = session.run_turn(user_msg, cancel)
    assert isinstance(assistant.usage, TokenUsage), type(assistant.usage)
    assert assistant.usage == TokenUsage(
      input=11, output=22, cache_read=33, cache_creation=44), assistant.usage
  finally:
    shutil.rmtree(tmp)


def exercise_deny_and_stop_ends_turn():
  """User selects 'Deny + Stop' on an approval prompt: the tool gets
  is_error=true, the cancel token is set, run_turn finalizes the in-flight
  assistant message with stop_reason='cancelled', and the model is NOT
  asked for another response (Section 4.4 / Section 10.2)."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="echo", input={}),
     TurnDone(stop_reason="tool_use")],
    # FakeAgent has only one scripted turn — if run_turn incorrectly
    # iterates after deny_and_stop, this would raise RuntimeError("ran out
    # of scripted turns") and the assertion below would fail.
  ])
  try:
    session.policy = ToolPolicy(default="ask")
    session.tools.register_builtin(
      ToolSpec(name="echo", description="echo",
               input_schema={"type": "object"}),
      handler=lambda **kw: "should-not-run",
      risk="write")
    # Deliver Deny+Stop from a separate thread once the card is surfaced -- the
    # production split: the worker blocks in _await_approval's fut.result()
    # while the GUI thread submits (on_event runs under the coordinator lock, so
    # answering inline from it would re-enter and deadlock).
    _captured, on_event, stop = _answer_surfaced_approvals(
      session, decision="deny_and_stop")
    session.on_event = on_event
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "hi"})], timestamp=now())
    assistant = session.run_turn(user_msg, cancel)
    stop()
    # The assistant message gets finalized with stop_reason="cancelled"
    # (turn unwinds without asking the model again).
    assert assistant.stop_reason == "cancelled"
    # tool_result message is appended with is_error=True.
    tr_msg = session.conv.messages[2]
    assert tr_msg.content[0].data["is_error"] is True
    assert "denied" in tr_msg.content[0].data["content"][0].data["text"].lower()
    # cancel token was set by deny_and_stop.
    assert cancel.is_set()
  finally:
    shutil.rmtree(tmp)


def exercise_server_tool_events_accumulate_without_dispatch():
  """ServerToolUsed / ServerToolResult events represent API-executed
  tools (web_search, web_fetch, code_execution). The session must
  append corresponding ContentBlocks to the assistant message so they
  persist + replay, but MUST NOT add them to tool_calls (no client
  dispatch). The turn ends with whatever stop_reason the model
  returned -- the server tool blocks do NOT trigger another loop
  iteration."""
  session, tmp = _new_test_session([
    [TextDelta(text="Let me search. "),
     ServerToolUsed(id="srvtu_1", name="web_search",
                    input={"query": "phenix.refine"}),
     ServerToolResult(tool_use_id="srvtu_1",
                      content={"type": "web_search_tool_result",
                               "results": ["result-1", "result-2"]}),
     TextDelta(text="Found two results."),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "search for phenix"})],
      timestamp=now())
    assistant = session.run_turn(user_msg, cancel)
    assert assistant.stop_reason == "end_turn"
    types_seen = [b.type for b in assistant.content]
    assert "server_tool_use" in types_seen, types_seen
    assert "server_tool_result" in types_seen, types_seen
    use_block = next(b for b in assistant.content
                     if b.type == "server_tool_use")
    res_block = next(b for b in assistant.content
                     if b.type == "server_tool_result")
    assert use_block.data["id"] == "srvtu_1"
    assert use_block.data["name"] == "web_search"
    assert use_block.data["input"] == {"query": "phenix.refine"}
    assert res_block.data["tool_use_id"] == "srvtu_1"
    assert res_block.data["content"]["results"] == ["result-1", "result-2"]
    # Only the user message + this single assistant message exist --
    # no tool_result message (no client dispatch happened).
    assert len(session.conv.messages) == 2, len(session.conv.messages)
  finally:
    shutil.rmtree(tmp)


# ---- claude_code observed-result answering -------------------------------
# The claude_code backend runs its OWN tool loop inside the SDK subprocess and
# surfaces each executed result as a ToolResultObserved -- its tool_use blocks
# still reach the session as ToolUseRequested (populating tool_calls), but the
# turn ends with a normal end_turn, NOT 'tool_use'. The session must answer the
# pending tool_uses with their REAL observed results instead of fabricating
# 'Cancelled by user.' errors, which otherwise surface on resume as bogus
# 'result (error)' entries.


def exercise_observed_results_answer_completed_claude_code_tools():
  """[Major] A claude_code turn (tool_use answered by ToolResultObserved, then a
  normal end_turn) must NOT fabricate a 'Cancelled by user.' result for the
  completed tool, and must commit the iteration as an (assistant tool_use, user
  tool_result) pair whose REAL answer keeps persistable_prefix from trimming the
  committed assistant -- the turn then ends in the final assistant text message
  (the unified, API-matching shape), not a skipped/trailing unanswered
  tool_use."""
  from qttbx.widgets.chat.agent.storage import persistable_prefix
  session, tmp = _new_test_session([
    [TextDelta(text="Checking. "),
     ToolUseRequested(id="t1", name="phenix_get_status", input={"job": "1"}),
     ToolResultObserved(tool_use_id="t1", content="Job 1: done",
                        name="phenix_get_status", input={"job": "1"}),
     TextDelta(text="It finished."),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "status?"})], timestamp=now())
    assistant = session.run_turn(user_msg, cancel)
    assert assistant.stop_reason == "end_turn", assistant.stop_reason
    # The tool_use is answered with the REAL result, is_error=False.
    tr_msgs = [m for m in session.conv.messages if m.role == "user"
               and m.content and m.content[0].type == "tool_result"]
    assert len(tr_msgs) == 1, [m.role for m in session.conv.messages]
    block = tr_msgs[0].content[0]
    assert block.data["tool_use_id"] == "t1", block.data
    assert block.data["is_error"] is False, "completed tool marked as error"
    text = block.data["content"][0].data["text"]
    assert "Job 1: done" in text, text
    # No fabricated cancel anywhere in the conversation.
    for m in session.conv.messages:
      for b in m.content:
        if b.type == "tool_result":
          assert "Cancelled by user" not in \
            b.data["content"][0].data["text"], b.data
    # The whole turn survives persistence in the unified (API-matching) shape:
    # the completed iteration committed as an (assistant tool_use, user
    # tool_result) pair and the trailing text landed as the final assistant, so
    # the turn ends in that final assistant message. Because the committed
    # tool_use is answered by a REAL tool_result (not a skipped/trailing
    # unanswered tool_use), persistable_prefix trims NOTHING -- a skipped answer
    # would leave a trailing unanswered assistant tool_use and drop that
    # committed assistant turn.
    kept = persistable_prefix(session.conv.messages)
    assert len(kept) == len(session.conv.messages), [m.role for m in kept]
    roles = [m.role for m in kept]
    assert roles == ["user", "assistant", "user", "assistant"], roles
    assert kept[-1].role == "assistant", roles
    assert any(m.role == "assistant" for m in kept), [m.role for m in kept]
    session.storage.save(session.conv)
    p = session.storage.conv_dir(session.conv.meta.id) / "messages.json"
    reloaded = json.load(open(p, encoding="utf-8"))["messages"]
    amsg = next(m for m in reloaded if m["role"] == "assistant")
    assert any(b["type"] == "tool_use" for b in amsg["content"]), amsg
    assert "Cancelled by user" not in json.dumps(reloaded), \
      "fabricated cancel persisted to messages.json"
  finally:
    shutil.rmtree(tmp)


def exercise_partial_observed_answers_mix_real_and_cancelled():
  """A claude_code turn cut short after one of two tools ran: the observed tool
  gets its real result, the un-observed one still gets 'Cancelled by user.' so
  the transcript never orphans a tool_use."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="phenix_get_status", input={}),
     ToolResultObserved(tool_use_id="t1", content="ran", name="phenix_get_status",
                        input={}),
     ToolUseRequested(id="t2", name="phenix_get_status", input={}),
     TurnDone(stop_reason="cancelled")],
  ])
  try:
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "go"})], timestamp=now())
    session.run_turn(user_msg, cancel)
    answers = {}
    for m in session.conv.messages:
      for b in m.content:
        if b.type == "tool_result":
          answers[b.data["tool_use_id"]] = (
            b.data["is_error"], b.data["content"][0].data["text"])
    assert set(answers) == {"t1", "t2"}, answers
    assert answers["t1"] == (False, "ran"), answers["t1"]
    assert answers["t2"][0] is True and \
      "Cancelled by user" in answers["t2"][1], answers["t2"]
  finally:
    shutil.rmtree(tmp)


def exercise_observed_error_result_preserves_is_error():
  """A claude_code tool that FAILED in its SDK subprocess -- a ToolResultObserved
  with is_error=True -- must be answered with an is_error tool_result, not
  recorded as a success, so a resumed conversation shows the failure rather than
  a muted 'success' cell."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="Bash", input={"command": "false"}),
     ToolResultObserved(tool_use_id="t1", content="exit status 1",
                        name="Bash", input={"command": "false"},
                        is_error=True),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "run it"})], timestamp=now())
    session.run_turn(user_msg, cancel)
    results = [b for m in session.conv.messages for b in m.content
               if b.type == "tool_result"]
    assert len(results) == 1, [m.role for m in session.conv.messages]
    block = results[0]
    assert block.data["tool_use_id"] == "t1", block.data
    assert block.data["is_error"] is True, "observed failure recorded as success"
    assert "exit status 1" in block.data["content"][0].data["text"], block.data
  finally:
    shutil.rmtree(tmp)


def exercise_observed_content_blocks_normalizes_shapes():
  """_observed_content_blocks (which builds each claude_code answer block) keeps
  text verbatim, reduces base64-bearing blocks (image / document) to a compact
  ``[<type>]`` placeholder so no blob inlines into the transcript, JSON-encodes
  small structured items, and maps empty / None to an empty text block."""
  from qttbx.widgets.chat.agent.session import _observed_content_blocks

  def _texts(blocks):
    return [b.data["text"] for b in blocks]

  assert _texts(_observed_content_blocks("hi")) == ["hi"]        # string verbatim
  blocks = _observed_content_blocks([
    {"type": "text", "text": "ok"},
    {"type": "image", "source": {"type": "base64", "data": "AAAA"}},
    {"type": "document", "source": {"type": "base64", "data": "BBBB"}},
    {"type": "tool_reference", "id": "r1"}])
  t = _texts(blocks)
  assert t[0] == "ok", t
  assert t[1] == "[image]", t
  assert t[2] == "[document]", t                                 # base64 doc, not blob
  assert "AAAA" not in "".join(t) and "BBBB" not in "".join(t), t
  assert "tool_reference" in t[3], t                             # small item -> JSON
  # empty list / empty string / None all normalize to one empty text block.
  assert _texts(_observed_content_blocks([])) == [""]
  assert _texts(_observed_content_blocks("")) == [""]
  assert _texts(_observed_content_blocks(None)) == [""]


# ---- phenix_ask_user_question (API-backend ask-user parity) --------------
# Mirrors the approval-queue parking pattern: the builtin handler emits an
# AskUserQuestionRequested and parks the worker on question_queue until the
# GUI delivers the answers via submit_question_answer (or a _Cancelled
# sentinel on cancel).


def exercise_await_question_answer_returns_preseeded_answers():
  """_await_question_answer emits the request and returns the answers the
  GUI delivered (here pre-seeded so the call doesn't actually block)."""
  session, tmp = _new_test_session([[TurnDone(stop_reason="end_turn")]])
  try:
    events = []
    session.on_event = lambda ev: events.append(ev)
    session.question_queue.put({"Q": "A"})
    answers = session._await_question_answer(
      "q1", [{"question": "Q", "options": [{"label": "A"}, {"label": "B"}]}])
    assert answers == {"Q": "A"}, answers
    asks = [e for e in events if isinstance(e, AskUserQuestionRequested)]
    assert len(asks) == 1, events
    assert asks[0].request_id == "q1", asks[0].request_id
    assert asks[0].questions[0]["question"] == "Q"
    # The in-flight id is cleared once the answer arrives.
    assert session._pending_question_id is None
  finally:
    shutil.rmtree(tmp)


def exercise_await_question_answer_cancel_raises_turn_cancelled():
  """A _Cancelled sentinel in the question queue raises TurnCancelled and
  clears the in-flight id (the finally block runs)."""
  session, tmp = _new_test_session([[TurnDone(stop_reason="end_turn")]])
  try:
    session.question_queue.put(_Cancelled())
    raised = False
    try:
      session._await_question_answer("q1", [])
    except TurnCancelled:
      raised = True
    assert raised, "TurnCancelled was not raised on a _Cancelled sentinel"
    assert session._pending_question_id is None
  finally:
    shutil.rmtree(tmp)


def exercise_submit_question_answer_routes_only_matching_id():
  """submit_question_answer accepts an answer only for the in-flight id and
  puts it on the queue; a non-matching id is rejected without queuing."""
  session, tmp = _new_test_session([[TurnDone(stop_reason="end_turn")]])
  try:
    session._pending_question_id = "q1"
    ans = {"Q": ["A", "B"]}
    assert session.submit_question_answer("q1", ans) is True
    assert session.question_queue.get_nowait() == ans
    # Wrong id -> rejected, nothing queued.
    assert session.submit_question_answer("other", ans) is False
    assert session.question_queue.empty()
  finally:
    shutil.rmtree(tmp)


def exercise_ask_user_question_tool_loop_round_trips_answers():
  """Integration: the model calls phenix_ask_user_question on turn 1; the
  registered builtin parks the worker; the GUI delivers answers; the
  tool_result message carries the answers JSON; turn 2 ends the turn."""
  from qttbx.widgets.chat.agent.tools import register_ask_user_question
  session, tmp = _new_test_session([
    # Turn 1: the model asks the user a question.
    [ToolUseRequested(
      id="t1", name="phenix_ask_user_question",
      input={"questions": [
        {"question": "Pick one", "multiSelect": False,
         "options": [{"label": "X"}, {"label": "Y"}]}]}),
     TurnDone(stop_reason="tool_use")],
    # Turn 2: final text after seeing the answers.
    [TextDelta(text="Thanks."),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    register_ask_user_question(session.tools)
    # Deliver the answers up front so the parked get() returns immediately.
    session.question_queue.put({"Pick one": "X"})
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "hi"})], timestamp=now())
    final = session.run_turn(user_msg, cancel)
    assert final.stop_reason == "end_turn"
    tr_msg = session.conv.messages[2]
    assert tr_msg.content[0].type == "tool_result"
    assert tr_msg.content[0].data["tool_use_id"] == "t1"
    assert tr_msg.content[0].data["is_error"] is False
    # The answers JSON is the tool result the model reads.
    text = tr_msg.content[0].data["content"][0].data["text"]
    assert "Pick one" in text and "X" in text, text
    import json as _json
    assert _json.loads(text) == {"Pick one": "X"}
  finally:
    shutil.rmtree(tmp)


def exercise_run_turn_drains_stale_cancelled_from_question_queue():
  """A _Cancelled left in question_queue by a previous turn's cancel must
  be drained at run_turn start, so the next turn's first question isn't
  aborted by the stale sentinel."""
  from qttbx.widgets.chat.agent.tools import register_ask_user_question
  session, tmp = _new_test_session([
    [ToolUseRequested(
      id="t1", name="phenix_ask_user_question",
      input={"questions": [
        {"question": "Q", "options": [{"label": "A"}, {"label": "B"}]}]}),
     TurnDone(stop_reason="tool_use")],
    [TextDelta(text="ok"), TurnDone(stop_reason="end_turn")],
  ])
  try:
    register_ask_user_question(session.tools)
    # Stale sentinel from a previous turn's cancel.
    session.question_queue.put(_Cancelled())
    # Real answer for this turn's question.
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "hi"})], timestamp=now())

    # Pre-seed the real answer only AFTER run_turn has drained the stale
    # sentinel. Because the answer is delivered before the question is
    # asked, we instead deliver it via on_event when the ask is surfaced.
    def _on_event(ev):
      if isinstance(ev, AskUserQuestionRequested):
        session.question_queue.put({"Q": "A"})
    session.on_event = _on_event

    final = session.run_turn(user_msg, cancel)
    assert final.stop_reason == "end_turn", final.stop_reason
    tr_msg = session.conv.messages[2]
    # If the stale sentinel had NOT been drained, the handler would have
    # raised TurnCancelled and is_error would be True.
    assert tr_msg.content[0].data["is_error"] is False
    text = tr_msg.content[0].data["content"][0].data["text"]
    assert "Q" in text and "A" in text, text
  finally:
    shutil.rmtree(tmp)


# ---- MCP result conversion -----------------------------------------------
# AgentSession._to_canonical_content_blocks is the seam between
# McpToolResult (the MCP client's typed payload) and ContentBlock (the
# chat's canonical wire shape). Tests live here because the conversion
# behaviour is owned by the session, not the MCP client.


def _new_session_for_mcp(tmp):
  from qttbx.widgets.chat.agent.conversation import Conversation
  from qttbx.widgets.chat.agent.session import AgentSession
  from qttbx.widgets.chat.agent.storage import ConversationStorage
  from qttbx.widgets.chat.agent.tools import ToolPolicy, ToolRegistry

  class _FakeAgent:
    name = "fake"
    model = "fake"
    def stream_turn(self, conversation, tools, cancel):
      yield from ()

  conv = Conversation.new(profile_name="t", model="m")
  storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
  return AgentSession(
    agent=_FakeAgent(), conversation=conv, storage=storage,
    tools=ToolRegistry(log=null_out()),
    policy=ToolPolicy(default="allow"),
    profile=None, log=null_out()), conv


def exercise_mcp_text_item_round_trips():
  from qttbx.widgets.chat.agent.mcp_client import (
    McpToolItem, McpToolResult)
  tmp = tempfile.mkdtemp()
  try:
    session, _ = _new_session_for_mcp(tmp)
    result = McpToolResult(content=[
      McpToolItem(type="text", text="hello")])
    blocks = session._to_canonical_content_blocks(result)
    assert len(blocks) == 1
    assert blocks[0].type == "text"
    assert blocks[0].data["text"] == "hello"
  finally:
    shutil.rmtree(tmp)


def exercise_mcp_image_item_uses_sha256_directly():
  from qttbx.widgets.chat.agent.mcp_client import (
    McpToolItem, McpToolResult)
  tmp = tempfile.mkdtemp()
  try:
    session, _ = _new_session_for_mcp(tmp)
    result = McpToolResult(content=[
      McpToolItem(type="image", sha256="abc",
                  mime="image/png", caption="cap")])
    blocks = session._to_canonical_content_blocks(result)
    assert len(blocks) == 1
    assert blocks[0].type == "image"
    assert blocks[0].data["attachment_sha256"] == "abc"
    assert blocks[0].data["mime"] == "image/png"
    assert blocks[0].data["caption"] == "cap"
  finally:
    shutil.rmtree(tmp)


def exercise_mcp_resource_item_renders_as_text():
  from qttbx.widgets.chat.agent.mcp_client import (
    McpToolItem, McpToolResult)
  tmp = tempfile.mkdtemp()
  try:
    session, _ = _new_session_for_mcp(tmp)
    result = McpToolResult(content=[
      McpToolItem(type="resource", uri="file:///x.pdb",
                  text_excerpt="ATOM 1 N")])
    blocks = session._to_canonical_content_blocks(result)
    assert blocks[0].type == "text"
    assert "file:///x.pdb" in blocks[0].data["text"]
    assert "ATOM 1 N" in blocks[0].data["text"]
  finally:
    shutil.rmtree(tmp)


def exercise_non_mcp_results_still_handled():
  """The non-MCP paths (bytes, str, dict) must keep working alongside
  the McpToolResult branch."""
  tmp = tempfile.mkdtemp()
  try:
    session, _ = _new_session_for_mcp(tmp)
    assert session._to_canonical_content_blocks(b"raw")[0].type == "text"
    assert session._to_canonical_content_blocks("hi")[0].type == "text"
    assert "key" in session._to_canonical_content_blocks(
      {"key": "value"})[0].data["text"]
  finally:
    shutil.rmtree(tmp)


# ---- MCP dispatch is_error propagation -----------------------------------
# The dispatch loop must carry an MCP tool's is_error flag onto the
# tool_result block: a failed MCP tool (McpToolResult.is_error=True) must
# reach the model + UI as an error, not be reported as a success.


def exercise_mcp_tool_error_result_propagates_is_error():
  """A dispatched MCP tool whose handler returns an McpToolResult flagged
  is_error=True must yield a tool_result block with is_error=True. The
  failure must reach the model + UI as an error rather than being reported
  as a success (the session must honour the handler result's is_error
  instead of hardcoding False)."""
  from qttbx.widgets.chat.agent.mcp_client import McpToolItem, McpToolResult
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="phenix_boom", input={}),
     TurnDone(stop_reason="tool_use")],
    [TextDelta(text="saw the error"),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    session.tools.register_mcp_tool(
      ToolSpec(name="phenix_boom", description="boom",
               input_schema={"type": "object"}),
      server_name="phenix",
      handler=lambda name, input, cancel: McpToolResult(
        content=[McpToolItem(type="text", text="it failed")],
        is_error=True),
      risk="write")
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "go"})], timestamp=now())
    session.run_turn(user_msg, cancel)
    tr_msg = session.conv.messages[2]
    assert tr_msg.content[0].type == "tool_result"
    assert tr_msg.content[0].data["tool_use_id"] == "t1"
    assert tr_msg.content[0].data["is_error"] is True, \
      "failed MCP tool was reported as success"
    text = tr_msg.content[0].data["content"][0].data["text"]
    assert "it failed" in text, text
  finally:
    shutil.rmtree(tmp)


def exercise_mcp_tool_success_result_keeps_is_error_false():
  """The is_error-honouring change must not flip a successful MCP tool to
  error: an McpToolResult with is_error=False yields a tool_result block
  with is_error=False."""
  from qttbx.widgets.chat.agent.mcp_client import McpToolItem, McpToolResult
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="phenix_ok", input={}),
     TurnDone(stop_reason="tool_use")],
    [TextDelta(text="done"), TurnDone(stop_reason="end_turn")],
  ])
  try:
    session.tools.register_mcp_tool(
      ToolSpec(name="phenix_ok", description="ok",
               input_schema={"type": "object"}),
      server_name="phenix",
      handler=lambda name, input, cancel: McpToolResult(
        content=[McpToolItem(type="text", text="all good")],
        is_error=False),
      risk="write")
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "go"})], timestamp=now())
    session.run_turn(user_msg, cancel)
    tr_msg = session.conv.messages[2]
    assert tr_msg.content[0].data["is_error"] is False
  finally:
    shutil.rmtree(tmp)


# ---- throttled mid-turn autosave -----------------------------------------
# AgentSession persists a non-mutating snapshot of the conversation (committed
# messages + the in-progress assistant message) at most every
# autosave_interval_s during a turn, reindex=False, top-level sessions only.
# The clock is injectable so the throttle can be crossed deterministically.


class _ManualClock(object):
  """Advanceable clock; call returns the current time in 'seconds'.

  ``start`` seeds the initial time -- a non-zero start is what distinguishes a
  throttle measured from turn start (run_turn resets ``_last_autosave`` to the
  clock) from one left at the ``__init__`` default of 0.0."""
  def __init__(self, start=0.0):
    self.t = start
  def __call__(self):
    return self.t


class _TimedFakeAgent(FakeAgent):
  """FakeAgent that advances a shared _ManualClock by a per-event delta as it
  yields, so the autosave throttle can be crossed deterministically."""
  def __init__(self, turn_scripts, clock, deltas):
    super(_TimedFakeAgent, self).__init__(turn_scripts)
    self._advance = clock
    self._deltas = list(deltas)
    self._i = 0
  def stream_turn(self, conversation, tools, cancel):
    events = self.turn_scripts[self._cursor]
    self._cursor += 1
    for ev in events:
      if self._i < len(self._deltas):
        self._advance.t += self._deltas[self._i]
      self._i += 1
      yield ev


# One recorded storage.save call: message count, trailing-assistant text length,
# the reindex flag it was called with, and whether the persisted snapshot ends
# in an assistant tool_use with no following tool_result (the un-resumable state
# a mid-turn crash must never leave on disk).
_Save = namedtuple("_Save", ["nmsgs", "text_len", "reindex", "orphan_tail"])


def _ondisk_ends_in_unanswered_tool_use(storage, conv_id):
  """True if the PERSISTED messages.json ends in an assistant tool_use with no
  following tool_result -- the un-resumable state no writer may freeze on disk.

  Reads what was actually written (after storage.save's centralized trim), so a
  test proves the real on-disk invariant regardless of which writer produced the
  snapshot or where the trim lives. Defined independently of the production
  helper so the test checks behavior, not the implementation."""
  p = storage.conv_dir(conv_id) / "messages.json"
  if not p.exists():
    return False
  with open(p, encoding="utf-8") as fh:
    msgs = json.load(fh).get("messages", [])
  return bool(msgs) and msgs[-1].get("role") == "assistant" \
      and any(b.get("type") == "tool_use" for b in msgs[-1].get("content", []))


def _install_save_recorder(session):
  """Record a _Save per storage.save call: message-count and trailing text
  length of the ARG, the reindex flag, and whether the PERSISTED messages.json
  ends in an unanswered tool_use (the crash-safety invariant, checked on disk)."""
  saves = []
  storage = session.storage
  real_save = storage.save
  def _rec(conv, reindex=True):
    text = ""
    if conv.messages:
      for b in conv.messages[-1].content:
        if b.type == "text":
          text += b.data.get("text", "")
    result = real_save(conv, reindex=reindex)
    saves.append(_Save(
      nmsgs=len(conv.messages),
      text_len=len(text),
      reindex=reindex,
      orphan_tail=_ondisk_ends_in_unanswered_tool_use(storage, conv.meta.id)))
    return result
  session.storage.save = _rec
  return saves


def _run_text_turn(session, script_text="hi"):
  cancel = CancelToken()
  user_msg = Message(role="user",
                     content=[ContentBlock(type="text",
                                           data={"text": script_text})],
                     timestamp=now())
  return session.run_turn(user_msg, cancel)


def exercise_autosave_throttled_during_stream():
  """A long content stream autosaves at throttle boundaries, and each mid-turn
  save captures the growing assistant text (proving incremental capture)."""
  clock = _ManualClock()
  script = [[TextDelta(text="x") for _ in range(11)]
            + [TurnDone(stop_reason="end_turn")]]
  agent = _TimedFakeAgent(script, clock, deltas=[2.0] * 12)
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0)
  try:
    saves = _install_save_recorder(session)
    _run_text_turn(session)
    # Content at t=2,4,...,22; last-save resets at t=6,12,18 -> 3 saves.
    assert len(saves) == 3, saves
    text_lens = [s.text_len for s in saves]
    assert text_lens == sorted(text_lens) and text_lens[0] < text_lens[-1], \
      text_lens                                    # incremental capture
    # Each snapshot carries the user message + the in-progress assistant msg.
    assert all(s.nmsgs == 2 for s in saves), saves
    # Mid-turn autosaves must skip the O(N) reindex (session passes reindex=False).
    assert all(s.reindex is False for s in saves), saves
  finally:
    shutil.rmtree(tmp)


def exercise_autosave_bursts_collapse_to_one_save():
  """Many events within a single interval collapse to one intermediate save."""
  clock = _ManualClock()
  script = [[TextDelta(text="x") for _ in range(5)]
            + [TurnDone(stop_reason="end_turn")]]
  # Only the 3rd event crosses 5.0; the rest stay inside the interval.
  agent = _TimedFakeAgent(script, clock, deltas=[1.0, 1.0, 4.0, 1.0, 1.0, 1.0])
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0)
  try:
    saves = _install_save_recorder(session)
    _run_text_turn(session)
    assert len(saves) == 1, saves
  finally:
    shutil.rmtree(tmp)


def exercise_autosave_not_fired_on_turndone():
  """The terminal TurnDone must never autosave, even when the throttle has
  elapsed -- doing so would race the GUI's turn-end save (GIL-releasing I/O
  between the queued emit and conv.append). A content event at the same time
  DOES save, proving it's the gating (not the throttle) that suppresses it."""
  # Gated case: only the TurnDone crosses the interval.
  clock = _ManualClock()
  script = [[TextDelta(text="x"), TurnDone(stop_reason="end_turn")]]
  agent = _TimedFakeAgent(script, clock, deltas=[0.0, 6.0])
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0)
  try:
    saves = _install_save_recorder(session)
    _run_text_turn(session)
    assert saves == [], saves                      # no save on TurnDone
  finally:
    shutil.rmtree(tmp)
  # Control: a content event at the same elapsed time DOES save.
  clock = _ManualClock()
  script = [[TextDelta(text="x"), TurnDone(stop_reason="end_turn")]]
  agent = _TimedFakeAgent(script, clock, deltas=[6.0, 0.0])
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0)
  try:
    saves = _install_save_recorder(session)
    _run_text_turn(session)
    assert len(saves) == 1, saves                  # content event saved
  finally:
    shutil.rmtree(tmp)


def exercise_autosave_skipped_for_subagent():
  """A depth>0 (subagent) session never autosaves; subagents persist as
  SubagentRecords, not top-level conversations."""
  clock = _ManualClock()
  script = [[TextDelta(text="x"), TextDelta(text="y"),
             TurnDone(stop_reason="end_turn")]]
  agent = _TimedFakeAgent(script, clock, deltas=[6.0, 6.0, 6.0])
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0, depth=1)
  try:
    saves = _install_save_recorder(session)
    _run_text_turn(session)
    assert saves == [], saves
  finally:
    shutil.rmtree(tmp)


def exercise_autosave_noop_without_storage():
  """With storage=None the autosave is a no-op: it neither raises NOR attempts a
  save. Capturing the log guards the `storage is None` early-return -- without
  it, _maybe_autosave would call None.save, whose AttributeError the broad except
  swallows and logs as 'mid-turn autosave failed', so a removed guard would
  otherwise ship undetected. The script emits no ImageEmitted (that would fault
  _accumulate.store_attachment on a None storage, independent of autosave)."""
  clock = _ManualClock()
  script = [[TextDelta(text="x"), TurnDone(stop_reason="end_turn")]]
  agent = _TimedFakeAgent(script, clock, deltas=[6.0, 0.0])
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0)
  try:
    session.storage = None
    log = io.StringIO()
    session.log = log
    assistant = _run_text_turn(session)            # must not raise
    assert assistant.stop_reason == "end_turn"
    assert "autosave failed" not in log.getvalue(), log.getvalue()
  finally:
    shutil.rmtree(tmp)


def _register_echo(session):
  session.tools.register_builtin(
    ToolSpec(name="echo", description="echo", input_schema={"type": "object"}),
    handler=lambda name, input, cancel, session, tool_use_id: input["text"],
    risk="write")


def exercise_autosave_over_tool_use_turn_is_crash_safe():
  """A mid-turn autosave that fires as the model streams a tool_use must NOT
  persist a snapshot ending in that unanswered tool_use: an assistant tool_use
  with no following tool_result is a non-recoverable provider 400 on reload (see
  exercise_cancel_*_does_not_orphan_tool_use). The throttle is arranged to fire
  exactly on the ToolUseRequested event; the persisted snapshot must be trimmed
  to the well-formed committed prefix instead."""
  clock = _ManualClock()
  script = [
    [TextDelta(text="calling"),
     ToolUseRequested(id="t1", name="echo", input={"text": "hi"}),
     TurnDone(stop_reason="tool_use")],
    [TextDelta(text="done"),
     TurnDone(stop_reason="end_turn")],
  ]
  # t: TextDelta=2 (<5, no save), ToolUseRequested=8 (>=5 -> the FIRST save
  # fires here, on the tool_use), then turn 2 TextDelta=14 (>=5 -> a save).
  agent = _TimedFakeAgent(script, clock, deltas=[2.0, 6.0, 0.0, 6.0, 0.0])
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0)
  try:
    _register_echo(session)
    saves = _install_save_recorder(session)
    _run_text_turn(session)
    assert saves, "expected at least one mid-turn autosave"
    # No persisted snapshot may end in an unanswered tool_use.
    assert not any(s.orphan_tail for s in saves), \
      [(s.nmsgs, s.orphan_tail) for s in saves]
  finally:
    shutil.rmtree(tmp)


def exercise_autosave_throttle_measured_from_turn_start():
  """run_turn resets _last_autosave to the current clock at turn start, so the
  throttle is measured from turn start -- NOT the __init__ default of 0.0. With
  a clock starting well past 0, a first content event inside the interval must
  NOT autosave (it would if the reset were missing and _last_autosave stayed
  0.0, since clock() - 0.0 >> the interval)."""
  clock = _ManualClock(start=1000.0)
  script = [[TextDelta(text="x"), TextDelta(text="y"),
             TurnDone(stop_reason="end_turn")]]
  # First event 2s into the turn (<5s -> no save); second 8s in (>=5s -> save).
  agent = _TimedFakeAgent(script, clock, deltas=[2.0, 6.0, 0.0])
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0)
  try:
    saves = _install_save_recorder(session)
    _run_text_turn(session)
    assert len(saves) == 1, [s.nmsgs for s in saves]
  finally:
    shutil.rmtree(tmp)


def exercise_checkpoint_fires_immediately_after_tool_result():
  """A completed tool iteration checkpoints to disk IMMEDIATELY (no 5s
  throttle): with a clock that never advances, the old throttled call could
  never fire, so the recorded nmsgs==3 save (user, assistant tool_use,
  tool_result) proves the unthrottled _checkpoint_iteration replaced it.
  reindex=False keeps the index rescan off the hot loop, and the persisted
  snapshot must not end in an unanswered tool_use (the crash-safety invariant
  that absorbs the old window-level orphan-tool_use save test)."""
  clock = _ManualClock()
  script = [
    [TextDelta(text="calling"),
     ToolUseRequested(id="t1", name="echo", input={"text": "hi"}),
     TurnDone(stop_reason="tool_use")],
    [TextDelta(text="done"),
     TurnDone(stop_reason="end_turn")],
  ]
  agent = _TimedFakeAgent(script, clock, deltas=[0.0] * 10)
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0)
  try:
    _register_echo(session)
    saves = _install_save_recorder(session)
    _run_text_turn(session)
    hits = [s for s in saves if s.nmsgs == 3]
    assert hits, [(s.nmsgs, s.reindex) for s in saves]
    assert all(s.reindex is False for s in hits), hits
    assert all(not s.orphan_tail for s in saves), saves
  finally:
    shutil.rmtree(tmp)


def exercise_checkpoint_floor_coalesces_rapid_iterations():
  """Two tool iterations 0.1s apart: the first checkpoints, the second lands
  inside the 250ms floor and is skipped -- a rapid tool burst coalesces
  instead of rewriting messages.json per iteration. (The skipped iteration
  is still committed in memory and reaches disk with the next save.)"""
  clock = _ManualClock()
  script = [
    [ToolUseRequested(id="t1", name="echo", input={"text": "a"}),
     TurnDone(stop_reason="tool_use")],
    [ToolUseRequested(id="t2", name="echo", input={"text": "b"}),
     TurnDone(stop_reason="tool_use")],
    [TextDelta(text="done"), TurnDone(stop_reason="end_turn")],
  ]
  # One 0.1s tick per event: every gap stays far under both the 5s autosave
  # interval and 250ms floor except the first checkpoint (floor starts unset).
  agent = _TimedFakeAgent(script, clock, deltas=[0.1] * 12)
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0)
  try:
    _register_echo(session)
    saves = _install_save_recorder(session)
    _run_text_turn(session)
    assert len(saves) == 1, [(s.nmsgs, s.reindex) for s in saves]
    assert saves[0].nmsgs == 3, saves
  finally:
    shutil.rmtree(tmp)


def exercise_checkpoint_floor_resets_per_turn():
  """The 250ms coalescing floor must not leak across turns: a NEW turn's
  first completed iteration checkpoints even when it lands within the floor
  of the PREVIOUS turn's last checkpoint. run_turn resets the floor sentinel
  at turn start (mirroring the autosave-clock reset), so the per-iteration
  freshness contract holds per turn, not per wall-clock window."""
  clock = _ManualClock()
  script = [
    [ToolUseRequested(id="t1", name="echo", input={"text": "a"}),
     TurnDone(stop_reason="tool_use")],
    [TextDelta(text="done"), TurnDone(stop_reason="end_turn")],
    [ToolUseRequested(id="t2", name="echo", input={"text": "b"}),
     TurnDone(stop_reason="tool_use")],
    [TextDelta(text="done2"), TurnDone(stop_reason="end_turn")],
  ]
  # 10ms per event: turn 2's iteration completes ~40ms after turn 1's
  # checkpoint -- far inside the floor, so only the per-turn reset lets it
  # write.
  agent = _TimedFakeAgent(script, clock, deltas=[0.01] * 16)
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0)
  try:
    _register_echo(session)
    saves = _install_save_recorder(session)
    _run_text_turn(session, "one")
    first_turn_saves = len(saves)
    assert first_turn_saves >= 1, [(s.nmsgs, s.reindex) for s in saves]
    _run_text_turn(session, "two")
    assert len(saves) > first_turn_saves, \
      ("turn 2's first iteration was floor-skipped",
       [(s.nmsgs, s.reindex) for s in saves])
    assert all(not s.orphan_tail for s in saves), saves
  finally:
    shutil.rmtree(tmp)


def exercise_claude_code_iterations_commit_as_alternating_pairs():
  """A claude_code-style turn (tool results observed in-stream, single
  terminal TurnDone) commits each completed iteration into the conversation
  as an (assistant tool_use, user tool_results) pair -- the SAME transcript
  shape the API backends produce -- and checkpoints it to disk mid-turn."""
  session, tmp = _new_test_session([
    [TextDelta(text="Checking. "),
     ToolUseRequested(id="t1", name="phenix_get_status", input={"job": "1"}),
     ToolResultObserved(tool_use_id="t1", content="Job 1: done",
                        name="phenix_get_status", input={"job": "1"}),
     TextDelta(text="It finished."),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    saves = _install_save_recorder(session)
    final = _run_text_turn(session, "status?")
    roles = [m.role for m in session.conv.messages]
    assert roles == ["user", "assistant", "user", "assistant"], roles
    committed = session.conv.messages[1]
    assert committed.stop_reason == "tool_use", committed.stop_reason
    assert any(b.type == "tool_use" for b in committed.content)
    results = session.conv.messages[2]
    assert [b.data["tool_use_id"] for b in results.content] == ["t1"]
    assert final.stop_reason == "end_turn", final.stop_reason
    assert final.content[0].data["text"] == "It finished."
    # The pair reached disk DURING the turn (run_turn itself never saves at
    # turn end, so any recorded save proves the mid-turn checkpoint).
    assert any(s.nmsgs == 3 and s.reindex is False and not s.orphan_tail
               for s in saves), [(s.nmsgs, s.reindex) for s in saves]
  finally:
    shutil.rmtree(tmp)


def exercise_parallel_tool_uses_commit_once_all_answered_in_order():
  """Two tool_uses in one assistant burst commit as ONE pair only after BOTH
  results are observed, results ordered by tool_use order (t1, t2) even when
  observed out of order (t2 first)."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="a", input={}),
     ToolUseRequested(id="t2", name="b", input={}),
     ToolResultObserved(tool_use_id="t2", content="B", name="b", input={}),
     ToolResultObserved(tool_use_id="t1", content="A", name="a", input={}),
     TextDelta(text="both done"),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    _run_text_turn(session, "go")
    roles = [m.role for m in session.conv.messages]
    assert roles == ["user", "assistant", "user", "assistant"], roles
    committed = session.conv.messages[1]
    assert [b.data["id"] for b in committed.content
            if b.type == "tool_use"] == ["t1", "t2"]
    results = session.conv.messages[2]
    assert [b.data["tool_use_id"] for b in results.content] == ["t1", "t2"]
    texts = [b.data["content"][0].data["text"] for b in results.content]
    assert texts == ["A", "B"], texts
  finally:
    shutil.rmtree(tmp)


def exercise_two_sequential_iterations_commit_two_pairs():
  """Sequential iterations inside one claude_code turn commit as two
  alternating pairs, then the trailing text lands as the final assistant."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="a", input={}),
     ToolResultObserved(tool_use_id="t1", content="A", name="a", input={}),
     ToolUseRequested(id="t2", name="b", input={}),
     ToolResultObserved(tool_use_id="t2", content="B", name="b", input={}),
     TextDelta(text="done"),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    _run_text_turn(session, "go")
    roles = [m.role for m in session.conv.messages]
    assert roles == ["user", "assistant", "user", "assistant", "user",
                     "assistant"], roles
    assert session.conv.messages[1].stop_reason == "tool_use"
    assert session.conv.messages[3].stop_reason == "tool_use"
    assert session.conv.messages[5].stop_reason == "end_turn"
  finally:
    shutil.rmtree(tmp)


def exercise_observed_result_for_unknown_tool_use_never_commits():
  """An observed result whose id matches no pending tool_use neither triggers
  a commit nor leaks into the committed results (same drop semantics the
  turn-end by_id assembly always had)."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="a", input={}),
     ToolResultObserved(tool_use_id="zz", content="stray", name="?", input={}),
     ToolResultObserved(tool_use_id="t1", content="A", name="a", input={}),
     TextDelta(text="done"),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    _run_text_turn(session, "go")
    results = session.conv.messages[2]
    assert [b.data["tool_use_id"] for b in results.content] == ["t1"]
  finally:
    shutil.rmtree(tmp)


def exercise_empty_residual_folds_usage_and_stop_onto_last_commit():
  """A turn ending right on a tool result (no trailing text) leaves an
  EMPTY residual assistant carrying the terminal stop_reason and usage
  (claude_code emits usage once, from the terminal ResultMessage).
  Appending it would be trimmed on save, silently dropping usage -- so it
  folds onto the last committed assistant instead, and run_turn returns
  that message."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="a", input={}),
     ToolResultObserved(tool_use_id="t1", content="A", name="a", input={}),
     TokenUsageEvent(input=7, output=3),
     TurnDone(stop_reason="end_turn")],
  ])
  try:
    final = _run_text_turn(session, "go")
    roles = [m.role for m in session.conv.messages]
    assert roles == ["user", "assistant", "user"], roles
    committed = session.conv.messages[1]
    assert final is committed, (final, committed)
    assert committed.stop_reason == "end_turn", committed.stop_reason
    assert committed.usage is not None and committed.usage.input == 7
  finally:
    shutil.rmtree(tmp)


def exercise_cancelled_empty_residual_folds_stop_without_usage_clobber():
  """The cancel flavor of the fold: no usage arrives (the agent's cancel
  path discards the ResultMessage), so only stop_reason folds -- usage on
  the committed message must not be clobbered to None."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="a", input={}),
     ToolResultObserved(tool_use_id="t1", content="A", name="a", input={}),
     TurnDone(stop_reason="cancelled")],
  ])
  try:
    final = _run_text_turn(session, "go")
    roles = [m.role for m in session.conv.messages]
    assert roles == ["user", "assistant", "user"], roles
    assert final.stop_reason == "cancelled", final.stop_reason
    assert final.usage is None
  finally:
    shutil.rmtree(tmp)


def exercise_session_propagates_allow_remember():
  """A builtin registered allow_remember=False surfaces a request with
  allow_remember False; a normal tool surfaces True. _resolve_and_approve
  must stamp the minted ToolApprovalRequest from the registry so the card
  can suppress 'Always allow this tool' for a destructive tool."""
  from qttbx.widgets.chat.agent.tools import ToolApprovalRequest
  session, tmp = _new_test_session([[TurnDone(stop_reason="end_turn")]])
  try:
    session.tools.register_builtin(
      ToolSpec(name="keepme", description="", input_schema={}),
      handler=lambda *a, **k: "", risk="destructive", allow_remember=False)
    session.tools.register_builtin(
      ToolSpec(name="normaltool", description="", input_schema={}),
      handler=lambda *a, **k: "", risk="write")
    session.policy = ToolPolicy(default="ask")
    # Answer each surfaced request from a separate thread (on_event runs under
    # the coordinator lock; inline submit would re-enter and deadlock). We only
    # care about the ToolApprovalRequests surfaced, not the answers.
    captured, on_event, stop = _answer_surfaced_approvals(session)
    session.on_event = on_event
    session._resolve_and_approve(
      ToolUseRequested(id="t1", name="keepme", input={}), batch_id=None)
    session._resolve_and_approve(
      ToolUseRequested(id="t2", name="normaltool", input={}), batch_id=None)
    stop()
    reqs = [e for e in captured if isinstance(e, ToolApprovalRequest)]
    assert len(reqs) == 2, reqs
    assert reqs[0].allow_remember is False, reqs[0]
    assert reqs[1].allow_remember is True, reqs[1]
  finally:
    shutil.rmtree(tmp)


def exercise_session_destructive_always_cards():
  """A1: a destructive / allow_remember=False builtin ALWAYS surfaces a card,
  even under a permissive tool_policy (default 'allow'); a normal tool under the
  same policy auto-approves without a card."""
  from qttbx.widgets.chat.agent.tools import ToolApprovalRequest
  session, tmp = _new_test_session([[TurnDone(stop_reason="end_turn")]])
  try:
    session.tools.register_builtin(
      ToolSpec(name="killtool", description="", input_schema={}),
      handler=lambda *a, **k: "", risk="destructive", allow_remember=False)
    session.tools.register_builtin(
      ToolSpec(name="readtool", description="", input_schema={}),
      handler=lambda *a, **k: "", risk="read")
    session.policy = ToolPolicy(default="allow")     # permissive
    # Destructive: forced to a card despite default='allow'. Answer surfaced
    # requests from a separate thread (on_event runs under the coordinator lock;
    # inline submit would re-enter and deadlock).
    captured, on_event, stop = _answer_surfaced_approvals(session)
    session.on_event = on_event
    d = session._resolve_and_approve(
      ToolUseRequested(id="t1", name="killtool", input={}), batch_id=None)
    reqs = [e for e in captured if isinstance(e, ToolApprovalRequest)]
    assert len(reqs) == 1 and d == "deny", (reqs, d)
    # Normal read tool under the same allow policy: auto-approved, no card.
    r = session._resolve_and_approve(
      ToolUseRequested(id="t2", name="readtool", input={}), batch_id=None)
    assert r == "approve", r
    assert not [e for e in captured[1:] if isinstance(e, ToolApprovalRequest)]
    # A destructive tool that STILL offers 'Always allow' (allow_remember=True)
    # is NOT force-carded -- else its shown checkbox would do nothing (#1). The
    # floor is scoped to the opt-out flag, not risk.
    session.tools.register_builtin(
      ToolSpec(name="rememberkill", description="", input_schema={}),
      handler=lambda *a, **k: "", risk="destructive", allow_remember=True)
    before = len([e for e in captured if isinstance(e, ToolApprovalRequest)])
    r2 = session._resolve_and_approve(
      ToolUseRequested(id="t3", name="rememberkill", input={}), batch_id=None)
    assert r2 == "approve", r2
    assert len([e for e in captured if isinstance(e, ToolApprovalRequest)]) \
        == before, "destructive+rememberable must not be force-carded"
    stop()
  finally:
    shutil.rmtree(tmp)


def exercise_cancel_landing_before_open_does_not_wedge():
  """A Stop landing after the dispatch cancel-check but before _await_approval
  opens its future must not park the worker forever. The runner sets the cancel
  token and then approvals.cancel_turn(), which CLOSES the coordinator; the
  subsequent open() returns a pre-cancelled future without registering or
  emitting, so _await_approval's fut.result() wakes cancelled and raises
  TurnCancelled instead of blocking. (The closed flag superseded the earlier
  post-open token re-check: cancel and open are now serialized under one lock.)"""
  import threading
  from qttbx.widgets.chat.agent.tools import ToolApprovalRequest
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="echo", input={}),
     TurnDone(stop_reason="tool_use")],
  ])
  try:
    session.policy = ToolPolicy(default="ask")
    session.tools.register_builtin(
      ToolSpec(name="echo", description="echo",
               input_schema={"type": "object"}),
      handler=lambda **kw: "x", risk="write")
    cancel = CancelToken()
    session.cancel = cancel                       # run_turn sets this per turn
    session.approvals.begin_turn()
    # Simulate a Stop landing before open(): the runner sets the token, then
    # cancel_turn() runs on the (still-empty) registry.
    cancel.set()
    session.approvals.cancel_turn()
    req = ToolApprovalRequest(
      request_id="r1", tool_name="echo", tool_source="builtin",
      input={}, risk="write", batch_id=None, allow_remember=True)
    box = []

    def _run():
      try:
        session._await_approval(req)
      except TurnCancelled:
        box.append("cancelled")

    t = threading.Thread(target=_run, daemon=True)
    t.start()
    t.join(timeout=2.0)
    assert not t.is_alive(), \
      "worker wedged: _await_approval blocked on a future a Stop already passed"
    assert box == ["cancelled"], box
    assert not session.approvals.owns("r1"), \
      "open() on the closed coordinator must not register the future"
  finally:
    shutil.rmtree(tmp)


def exercise_run_turn_cancels_pending_approval_on_exit():
  """A turn that ends with an approval still open (subprocess-death style) must
  leave the coordinator registry empty -- run_turn cancel_turn()s in a finally,
  so owns() is False when idle and no pending task is destroyed at exit."""
  from qttbx.widgets.chat.agent.tools import ToolApprovalRequest
  session, tmp = _new_test_session([])
  try:
    box = {}

    class _OpensApprovalAgent(Agent):
      name = "opens"
      model = "opens-1"

      def stream_turn(self, conversation, tools, cancel):
        # Open an approval future the way claude_code's _on_can_use_tool would,
        # then end the turn WITHOUT it resolving (as on a subprocess death).
        req = ToolApprovalRequest(
          request_id="r1", tool_name="t", tool_source="builtin",
          input={}, risk="write", batch_id=None, allow_remember=True)
        box["fut"] = session.approvals.open(req, lambda ev: None)
        yield TurnDone(stop_reason="end_turn")

      def resolve_credentials(self, cli_override=None):
        return "k"

      def credentials_dialog_class(self):
        return object

    session.agent = _OpensApprovalAgent()
    user = Message(role="user", timestamp=now(),
                   content=[ContentBlock(type="text", data={"text": "hi"})])
    session.run_turn(user, CancelToken())
    assert not session.approvals.owns("r1"), "pending approval leaked past turn end"
    assert box["fut"].cancelled(), "pending future not cancelled at turn end"
  finally:
    shutil.rmtree(tmp)


def exercise_api_remember_recorded_before_resolve_survives_cancel():
  """The API-backend remember='tool' persistence runs in record_approval_remember
  (GUI thread, before the coordinator resolves) via the request->tool map, so a
  click racing a cancel still records 'always allow' -- parity with claude_code
  (whose recording lives in submit_approval, likewise pre-resolve)."""
  from qttbx.widgets.chat.agent.tools import ToolApprovalResponse
  session, tmp = _new_test_session([])
  try:
    session.policy = ToolPolicy(default="ask")
    session.tools.register_builtin(
      ToolSpec(name="echo", description="echo",
               input_schema={"type": "object"}),
      handler=lambda **kw: "x", risk="write")
    # Simulate _resolve_and_approve having surfaced request r1 for tool 'echo'.
    session._pending_approval_calls = {"r1": "echo"}
    assert session.policy.resolve("echo") == "ask"          # not yet remembered
    session.record_approval_remember(ToolApprovalResponse(
      request_id="r1", decision="approve", remember="tool"))
    assert session.policy.resolve("echo") == "allow"        # now session-allowed
  finally:
    shutil.rmtree(tmp)


# ---- synthesized terminal TurnDone on cancel short-circuits --------------
# _run_turn's two cancel short-circuits (pre-dispatch, post-dispatch)
# historically returned with no event, leaving the stale intermediate
# finish=tool_use as the last TurnDone consumers saw. Each now synthesizes a
# terminal TurnDone(stop_reason=cancelled, finish=cancelled) AFTER the residual
# appends, restoring one terminal TurnDone per completed turn.


def exercise_cancel_before_dispatch_synthesizes_terminal_turndone():
  """Cancel landing between a tool_use-stop response and dispatch exits
  run_turn's short-circuit, which historically emitted NO event -- leaving
  GUI/headless consumers with a stale intermediate finish=TOOL_USE as the
  last word. The session now synthesizes the terminal
  TurnDone(cancelled) AFTER the residual appends, restoring the
  one-terminal-TurnDone-per-turn invariant."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="echo", input={"text": "x"}),
     TurnDone(stop_reason="tool_use")],
  ])
  try:
    _register_echo(session)
    cancel = CancelToken()
    seen = []
    at_emit = []
    def _on_event(ev):
      seen.append(ev)
      if isinstance(ev, TurnDone) and ev.stop_reason == "tool_use":
        cancel.set()                 # lands before dispatch
      elif isinstance(ev, TurnDone) and ev.stop_reason == "cancelled":
        # Snapshot, AT the synthesized terminal emit, whether this turn's
        # residual answer is already on the conversation. run_turn is driven on
        # this thread, so the closure runs synchronously inside the emit; a
        # post-hoc conv check cannot see the emit-vs-append order (conv.append
        # fires no event), so only this in-closure snapshot discriminates it.
        at_emit.append(any(
          b.type == "tool_result" and b.data.get("tool_use_id") == "t1"
          for m in session.conv.messages for b in m.content))
    session.on_event = _on_event
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "go"})], timestamp=now())
    session.run_turn(user_msg, cancel)
    cancelled = [e for e in seen if isinstance(e, TurnDone)
                 and e.stop_reason == "cancelled"]
    assert len(cancelled) == 1, [
      (e.stop_reason, e.finish) for e in seen if isinstance(e, TurnDone)]
    assert cancelled[0].finish == "cancelled", cancelled[0]
    assert isinstance(seen[-1], TurnDone) and seen[-1] is cancelled[0], \
      seen[-1]
    # The residual appends preceded the emit: at the sole cancelled emit the
    # tool_use was already answered on the conversation (at_emit == [True]).
    # This closure snapshot -- not the post-hoc conv shape below -- is what
    # proves the ordering, since reordering the emit before the appends leaves
    # the final conv shape identical and only flips at_emit to [False].
    assert at_emit == [True], at_emit
    results = [b for m in session.conv.messages for b in m.content
               if b.type == "tool_result"]
    assert results and results[0].data["tool_use_id"] == "t1"
  finally:
    shutil.rmtree(tmp)


def exercise_cancel_during_dispatch_synthesizes_terminal_turndone():
  """A tool handler that sets the cancel token (the deny_and_stop shape)
  exits through the post-dispatch short-circuit -- the second historically
  event-less return. One synthesized terminal TurnDone(cancelled), emitted
  last."""
  session, tmp = _new_test_session([
    [ToolUseRequested(id="t1", name="boom", input={}),
     TurnDone(stop_reason="tool_use")],
  ])
  try:
    session.tools.register_builtin(
      ToolSpec(name="boom", description="sets cancel",
               input_schema={"type": "object"}),
      handler=lambda name, input, cancel, session, tool_use_id:
        (cancel.set(), "stopped")[1],
      risk="write")
    seen = []
    at_emit = []
    def _on_event(ev):
      seen.append(ev)
      if isinstance(ev, TurnDone) and ev.stop_reason == "cancelled":
        # As in the pre-dispatch test: snapshot inside the synchronous emit
        # whether the post-dispatch residual (the boom tool_result) is already
        # appended. Proves the emit lands after conv.append(tool_result_msg).
        at_emit.append(any(
          b.type == "tool_result" and b.data.get("tool_use_id") == "t1"
          for m in session.conv.messages for b in m.content))
    session.on_event = _on_event
    _run_text_turn(session, "go")
    cancelled = [e for e in seen if isinstance(e, TurnDone)
                 and e.stop_reason == "cancelled"]
    assert len(cancelled) == 1 and cancelled[0].finish == "cancelled", [
      (e.stop_reason, e.finish) for e in seen if isinstance(e, TurnDone)]
    assert seen[-1] is cancelled[0], seen[-1]
    # Residual appended before the emit (post-dispatch short-circuit).
    assert at_emit == [True], at_emit
  finally:
    shutil.rmtree(tmp)


def exercise_streaming_cancel_synthesizes_no_second_terminal():
  """A cancel the AGENT detected mid-stream already carries its own terminal
  TurnDone(cancelled); the session must not add a second one (that path
  returns through the terminal branch, not a short-circuit)."""
  session, tmp = _new_test_session([
    [TextDelta(text="partial"),
     TurnDone(stop_reason="cancelled", finish="cancelled")],
  ])
  try:
    seen = []
    session.on_event = seen.append
    _run_text_turn(session, "go")
    cancelled = [e for e in seen if isinstance(e, TurnDone)
                 and e.stop_reason == "cancelled"]
    assert len(cancelled) == 1, [
      (e.stop_reason, e.finish) for e in seen if isinstance(e, TurnDone)]
  finally:
    shutil.rmtree(tmp)


def exercise():
  exercise_simple_text_turn()
  exercise_assistant_messages_stamped_with_model_and_backend()
  exercise_reconciles_meta_model_and_backend_on_continue()
  exercise_tool_use_loop_completes()
  exercise_cancel_after_tool_use_does_not_orphan_tool_use()
  exercise_cancel_during_streaming_does_not_orphan_tool_use()
  exercise_cancel_during_multi_tool_streaming_answers_all()
  exercise_cancel_before_dispatch_synthesizes_terminal_turndone()
  exercise_cancel_during_dispatch_synthesizes_terminal_turndone()
  exercise_streaming_cancel_synthesizes_no_second_terminal()
  exercise_builtin_handler_receives_tool_use_id()
  exercise_tool_denied_returns_is_error()
  exercise_cancel_while_awaiting_approval()
  exercise_cancel_landing_before_open_does_not_wedge()
  exercise_run_turn_cancels_pending_approval_on_exit()
  exercise_api_remember_recorded_before_resolve_survives_cancel()
  exercise_max_turns_cap()
  exercise_add_subagent_usage_aggregates()
  exercise_token_usage_event_carries_all_fields_to_message()
  exercise_deny_and_stop_ends_turn()
  exercise_server_tool_events_accumulate_without_dispatch()
  exercise_observed_results_answer_completed_claude_code_tools()
  exercise_partial_observed_answers_mix_real_and_cancelled()
  exercise_observed_error_result_preserves_is_error()
  exercise_observed_content_blocks_normalizes_shapes()
  exercise_await_question_answer_returns_preseeded_answers()
  exercise_await_question_answer_cancel_raises_turn_cancelled()
  exercise_submit_question_answer_routes_only_matching_id()
  exercise_ask_user_question_tool_loop_round_trips_answers()
  exercise_run_turn_drains_stale_cancelled_from_question_queue()
  exercise_mcp_text_item_round_trips()
  exercise_mcp_image_item_uses_sha256_directly()
  exercise_mcp_resource_item_renders_as_text()
  exercise_non_mcp_results_still_handled()
  exercise_mcp_tool_error_result_propagates_is_error()
  exercise_mcp_tool_success_result_keeps_is_error_false()
  exercise_autosave_throttled_during_stream()
  exercise_autosave_bursts_collapse_to_one_save()
  exercise_autosave_not_fired_on_turndone()
  exercise_autosave_skipped_for_subagent()
  exercise_autosave_noop_without_storage()
  exercise_autosave_over_tool_use_turn_is_crash_safe()
  exercise_autosave_throttle_measured_from_turn_start()
  exercise_checkpoint_fires_immediately_after_tool_result()
  exercise_checkpoint_floor_coalesces_rapid_iterations()
  exercise_checkpoint_floor_resets_per_turn()
  exercise_claude_code_iterations_commit_as_alternating_pairs()
  exercise_parallel_tool_uses_commit_once_all_answered_in_order()
  exercise_two_sequential_iterations_commit_two_pairs()
  exercise_observed_result_for_unknown_tool_use_never_commits()
  exercise_empty_residual_folds_usage_and_stop_onto_last_commit()
  exercise_cancelled_empty_residual_folds_stop_without_usage_clobber()
  exercise_session_propagates_allow_remember()
  exercise_session_destructive_always_cards()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
