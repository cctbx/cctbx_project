import queue
import shutil
import tempfile
import threading
from pathlib import Path

from libtbx.utils import format_cpu_times, null_out
from qttbx.widgets.chat.agent.base import (
  Agent, AgentCapabilities, ToolSpec)
from qttbx.widgets.chat.agent.conversation import (
  ContentBlock, Conversation, Message, TokenUsage, now)
from qttbx.widgets.chat.agent.errors import CancelToken, TurnCancelled
from qttbx.widgets.chat.agent.events import (
  ServerToolResult, ServerToolUsed, TextDelta, ToolUseRequested,
  TurnDone, TokenUsage as TokenUsageEvent)
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
  capabilities = AgentCapabilities.STREAMING | AgentCapabilities.TOOL_USE

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


def _new_test_session(turn_scripts, project_dir=None, max_depth=1):
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
    agent=FakeAgent(turn_scripts),
    conversation=conv,
    storage=storage,
    tools=registry,
    policy=policy,
    profile=profile,
    log=null_out())
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
    capabilities = 0
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


def exercise():
  exercise_simple_text_turn()
  exercise_tool_use_loop_completes()
  exercise_tool_denied_returns_is_error()
  exercise_cancel_while_awaiting_approval()
  exercise_max_turns_cap()
  exercise_add_subagent_usage_aggregates()
  exercise_deny_and_stop_ends_turn()
  exercise_server_tool_events_accumulate_without_dispatch()
  exercise_mcp_text_item_round_trips()
  exercise_mcp_image_item_uses_sha256_directly()
  exercise_mcp_resource_item_renders_as_text()
  exercise_non_mcp_results_still_handled()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
