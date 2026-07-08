import io
import json
import queue
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
  ToolUseRequested, TurnDone, TokenUsage as TokenUsageEvent)
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
  """The deadlock-fix path: worker parked on approval queue, GUI pushes
  _Cancelled sentinel, worker raises TurnCancelled, turn ends."""
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

    # Wire a real Queue for approvals.
    approval_queue = queue.Queue()
    session.approval_queue = approval_queue

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
    approval_queue.put(_Cancelled())

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
    approval_queue = queue.Queue()
    session.approval_queue = approval_queue
    # Pre-load the queue with a deny_and_stop response so the worker
    # doesn't actually have to block.
    approval_queue.put(ToolApprovalResponse(
      request_id="ignored",  # session uses internally-generated ids
      decision="deny_and_stop"))
    cancel = CancelToken()
    user_msg = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "hi"})], timestamp=now())
    assistant = session.run_turn(user_msg, cancel)
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


def exercise_autosave_checkpoints_after_tool_result():
  """The after-tool-result checkpoint (run_turn's _maybe_autosave() with no
  partial_msg) persists a committed-only snapshot -- tail=[] -- ending in the
  tool_result, reindex=False. This is the crash checkpoint for a long tool
  dispatch; the throttle is arranged to cross only after the tool_result is
  appended, exercising the partial_msg=None branch a text-only stream never
  reaches."""
  clock = _ManualClock()
  script = [
    [TextDelta(text="calling"),
     ToolUseRequested(id="t1", name="echo", input={"text": "hi"}),
     TurnDone(stop_reason="tool_use")],
    [TextDelta(text="done"),
     TurnDone(stop_reason="end_turn")],
  ]
  # Stay under 5s through turn-1 streaming (t=1,2); the TurnDone delta (t=8)
  # pushes past 5s so the post-tool-result checkpoint fires on the committed
  # [user, assistant(tool_use), tool_result]; turn-2 events add no time.
  agent = _TimedFakeAgent(script, clock, deltas=[1.0, 1.0, 6.0, 0.0, 0.0])
  session, tmp = _new_test_session(script, agent=agent, clock=clock,
                                   autosave_interval_s=5.0)
  try:
    _register_echo(session)
    saves = _install_save_recorder(session)
    _run_text_turn(session)
    assert any(s.nmsgs == 3 and not s.orphan_tail and s.reindex is False
               for s in saves), \
      [(s.nmsgs, s.orphan_tail, s.reindex) for s in saves]
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


def exercise():
  exercise_simple_text_turn()
  exercise_assistant_messages_stamped_with_model_and_backend()
  exercise_reconciles_meta_model_and_backend_on_continue()
  exercise_tool_use_loop_completes()
  exercise_cancel_after_tool_use_does_not_orphan_tool_use()
  exercise_cancel_during_streaming_does_not_orphan_tool_use()
  exercise_cancel_during_multi_tool_streaming_answers_all()
  exercise_builtin_handler_receives_tool_use_id()
  exercise_tool_denied_returns_is_error()
  exercise_cancel_while_awaiting_approval()
  exercise_max_turns_cap()
  exercise_add_subagent_usage_aggregates()
  exercise_token_usage_event_carries_all_fields_to_message()
  exercise_deny_and_stop_ends_turn()
  exercise_server_tool_events_accumulate_without_dispatch()
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
  exercise_autosave_checkpoints_after_tool_result()
  exercise_autosave_throttle_measured_from_turn_start()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
