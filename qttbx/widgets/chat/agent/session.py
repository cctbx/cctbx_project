"""``AgentSession`` — multi-iteration tool loop orchestrator.

Wraps any ``Agent``. Pure-Python; no Qt. ``QtAgentRunner`` (in
``runner.py``) wraps this to emit Qt signals on the GUI thread.
"""

import json
import queue
import sys
import uuid

from libtbx.utils import Sorry

from qttbx.widgets.chat.agent.conversation import (
  ContentBlock, Conversation, Message, TokenUsage, now)
from qttbx.widgets.chat.agent.errors import TurnCancelled
from qttbx.widgets.chat.agent.events import (
  AgentError, ImageEmitted, ServerToolResult, ServerToolUsed,
  Thinking, TextDelta, ToolResultsBatched, ToolUseRequested, TurnDone,
  TokenUsage as TokenUsageEvent)
from qttbx.widgets.chat.agent.tools import (
  ToolApprovalRequest, ToolApprovalResponse, _Cancelled)


def _new_id(prefix=""):
  return prefix + uuid.uuid4().hex[:12]


class AgentSession:
  """Pure-Python orchestrator for the chat tool loop.

  Reusable in any context (UI thread, worker thread, headless subagent).
  Wrappers (Qt or otherwise) add transport / signaling but don't change
  the loop semantics."""

  def __init__(self, agent, conversation, storage, tools, policy,
               profile, depth=0, on_event=None, log=None,
               approval_queue=None):
    self.agent = agent
    self.conv = conversation
    self.storage = storage
    self.tools = tools
    self.policy = policy
    self.profile = profile
    self.depth = depth
    self.on_event = on_event or (lambda ev: None)
    self.log = log if log is not None else sys.stdout

    # The approval queue is the worker's blocking point during 'ask' policy.
    # External callers (the GUI runner) inject a queue they also push to;
    # tests inject a queue.Queue.
    self.approval_queue = approval_queue or queue.Queue()

    # Set when run_turn starts; consulted by tool handlers.
    self.cancel = None
    self.sub_id = _new_id("sa_") if depth > 0 else None
    self.started_at = now()
    self.current_tool_use_id = None
    # Push-based rollup of nested-session token usage. Keyed by sub_id.
    self._subagent_usage_by_id = {}

  # ---- main loop -----------------------------------------------------------

  def run_turn(self, user_message, cancel, max_turns=None):
    """Append the user message and iterate the tool-use loop.

    Loop terminates on ``stop_reason != 'tool_use'``, on cancellation,
    or when the optional ``max_turns`` cap is hit. When the cap fires,
    a visible ``[Subagent stopped at turn cap (N)]`` marker is appended
    so the model (or a viewer) sees why the loop ended.

    Cancel handling: if the cancel token is set after a tool dispatch
    (``deny_and_stop``, or a ``_Cancelled`` sentinel during approval),
    we do NOT ask the model for another response — we finalize the
    in-flight assistant message with ``stop_reason='cancelled'`` and
    return.

    Parameters
    ----------
    user_message : Message
        The user's input message for this turn.
    cancel : CancelToken
        Polled between events and tool dispatches.
    max_turns : int, optional
        Safety cap for nested-agent loops. ``None`` means unbounded.

    Returns
    -------
    Message
        The final assistant message of the turn.
    """
    self.cancel = cancel
    self.conv.append(user_message)
    iterations = 0
    while True:
      iterations += 1
      assistant_msg, tool_calls = self._collect_one_response(cancel)
      self.conv.append(assistant_msg)

      if assistant_msg.stop_reason != "tool_use":
        return assistant_msg
      if cancel.is_set():
        assistant_msg.stop_reason = "cancelled"
        return assistant_msg

      if max_turns is not None and iterations >= max_turns:
        self.conv.append(Message(
          role="user", timestamp=now(),
          content=[ContentBlock(type="text", data={
            "text": "[Subagent stopped at turn cap (%d)]" % max_turns})]))
        return assistant_msg

      tool_result_msg = self._dispatch_and_build_results(tool_calls, cancel)
      self.conv.append(tool_result_msg)

      # Dispatch may have set the cancel token (deny_and_stop or sentinel
      # cancellation during approval). Don't ask the model for another
      # response — that would either re-enter a cancelled loop or, in
      # FakeAgent tests, exhaust the scripted turns.
      if cancel.is_set():
        assistant_msg.stop_reason = "cancelled"
        return assistant_msg

  def _collect_one_response(self, cancel):
    msg = Message(role="assistant", content=[], timestamp=now())
    tool_calls = []
    for event in self.agent.stream_turn(self.conv, self.tools.specs(), cancel):
      self.on_event(event)
      _accumulate(msg, event, tool_calls, self.storage, self.conv.meta.id)
    return msg, tool_calls

  # ---- tool dispatch -------------------------------------------------------

  def _dispatch_and_build_results(self, tool_calls, cancel):
    batch_id = _new_id("b_") if len(tool_calls) > 1 else None
    result_blocks = []

    for call in tool_calls:
      if cancel.is_set():
        result_blocks.append(_tool_error_block(
          call.id, "Cancelled before dispatch"))
        continue

      try:
        decision = self._resolve_and_approve(call, batch_id, cancel)
      except TurnCancelled:
        result_blocks.append(_tool_error_block(
          call.id, "Cancelled by user"))
        continue

      if decision in ("deny", "deny_and_stop"):
        result_blocks.append(_tool_error_block(
          call.id, "Denied by user policy"))
        if decision == "deny_and_stop":
          cancel.set()
        continue

      try:
        handler_result = self._invoke_tool(call, cancel)
      except TurnCancelled:
        result_blocks.append(_tool_error_block(
          call.id, "Cancelled during execution"))
        continue
      except Sorry as e:
        result_blocks.append(_tool_error_block(call.id, str(e)))
        continue
      except Exception as e:
        print("Tool '%s' raised: %s" % (call.name, e), file=self.log)
        result_blocks.append(_tool_error_block(
          call.id, "Tool error: %s" % e))
        continue

      blocks = self._to_canonical_content_blocks(handler_result)
      result_blocks.append(ContentBlock(type="tool_result", data={
        "tool_use_id": call.id,
        "content": blocks,
        "is_error": False,
      }))

    self.on_event(ToolResultsBatched(blocks=result_blocks))
    return Message(role="user", content=result_blocks, timestamp=now())

  def _resolve_and_approve(self, call, batch_id, cancel):
    policy = self.policy.resolve(call.name)
    if policy == "deny":
      return "deny"
    if policy == "allow":
      return "approve"
    # policy == "ask"
    req = ToolApprovalRequest(
      request_id=_new_id("r_"),
      tool_name=call.name,
      tool_source=self.tools.source_of(call.name) or "builtin",
      input=call.input,
      risk=self.tools.risk_of(call.name),
      summary=_summarize_call(call),
      batch_id=batch_id,
    )
    response = self._await_approval(req)
    if response.decision == "approve":
      if response.remember == "tool":
        self.policy.allow_tool_for_session(call.name)
      elif response.remember == "server":
        server = self.tools.server_of(call.name)
        if server:
          self.policy.allow_server_for_session(server)
      return "approve"
    return response.decision

  def _await_approval(self, req):
    """Park the worker on the approval queue. Wakes on either a real
    ``ToolApprovalResponse`` or a ``_Cancelled`` sentinel pushed by the
    GUI's cancel handler."""
    self.on_event(req)                              # surface to UI
    response = self.approval_queue.get()
    if isinstance(response, _Cancelled):
      raise TurnCancelled()
    return response

  def _invoke_tool(self, call, cancel):
    source = self.tools.source_of(call.name) or ""
    if source == "builtin":
      self.current_tool_use_id = call.id
      try:
        return self.tools.invoke_builtin(
          call.name, call.input,
          cancel=cancel, session=self, tool_use_id=call.id)
      finally:
        self.current_tool_use_id = None
    if source == "skill":
      return self.tools.invoke_skill(call.name, call.input)
    if source.startswith("mcp:"):
      return self.tools.invoke_mcp(call.name, call.input, cancel=cancel)
    raise Sorry("Unknown tool source for '%s'" % call.name)

  def _to_canonical_content_blocks(self, result):
    # Function-scope import: avoid pulling mcp_client into the import graph
    # for sessions that never see an MCP result (MCP is optional).
    from qttbx.widgets.chat.agent.mcp_client import (
      McpToolResult, _mcp_item_to_block)
    if isinstance(result, McpToolResult):
      return [_mcp_item_to_block(item, self.storage, self.conv.meta.id)
              for item in result.content]
    if isinstance(result, bytes):
      return [ContentBlock(type="text", data={
        "text": result.decode("utf-8", errors="replace")})]
    if isinstance(result, str):
      return [ContentBlock(type="text", data={"text": result})]
    return [ContentBlock(type="text", data={
      "text": json.dumps(result, indent=2, default=str)})]

  # ---- nested-session rollup -----------------------------------------------

  def add_subagent_usage(self, sub_id, usage):
    """Record token usage from a nested session keyed by ``sub_id``."""
    self._subagent_usage_by_id[sub_id] = usage


# ---- module helpers --------------------------------------------------------

def _accumulate(msg, event, tool_calls, storage, conv_id):
  if isinstance(event, TextDelta):
    _append_text(msg, event.text)
  elif isinstance(event, Thinking):
    _append_thinking(msg, event.text, event.signature)
  elif isinstance(event, ToolUseRequested):
    msg.content.append(ContentBlock(type="tool_use", data={
      "id": event.id, "name": event.name, "input": event.input}))
    tool_calls.append(event)
  elif isinstance(event, ServerToolUsed):
    # API-executed; no dispatch needed. Persist for the bubble + replay.
    msg.content.append(ContentBlock(type="server_tool_use", data={
      "id": event.id, "name": event.name, "input": event.input}))
  elif isinstance(event, ServerToolResult):
    msg.content.append(ContentBlock(type="server_tool_result", data={
      "tool_use_id": event.tool_use_id, "content": event.content}))
  elif isinstance(event, ImageEmitted):
    att = storage.store_attachment(conv_id, event.data, event.mime)
    msg.content.append(ContentBlock(type="image", data={
      "attachment_sha256": att.sha256, "mime": event.mime,
      "caption": event.caption}))
  elif isinstance(event, TokenUsageEvent):
    msg.usage = TokenUsage(
      input=event.input,
      output=event.output,
      cache_read=event.cache_read,
      cache_creation=event.cache_creation)
  elif isinstance(event, TurnDone):
    msg.stop_reason = event.stop_reason
  elif isinstance(event, AgentError):
    # Error events are surfaced via on_event by the caller; we also stop
    # accumulating into this assistant message by not appending content.
    pass


def _append_text(msg, text):
  if msg.content and msg.content[-1].type == "text":
    msg.content[-1].data["text"] += text
  else:
    msg.content.append(ContentBlock(type="text", data={"text": text}))


def _append_thinking(msg, text, signature):
  if msg.content and msg.content[-1].type == "thinking":
    msg.content[-1].data["text"] += text
    if signature:
      msg.content[-1].data["signature"] = signature
  else:
    msg.content.append(ContentBlock(type="thinking", data={
      "text": text, "signature": signature or ""}))


def _summarize_call(call):
  """Short human-readable summary for the approval card. Best-effort; the
  UI shows the full input when expanded."""
  return "%s(%s)" % (call.name, ", ".join(
    "%s=%s" % (k, _short(v)) for k, v in call.input.items()))


def _short(v, n=40):
  s = str(v)
  return s if len(s) <= n else s[:n] + "..."


def _tool_error_block(tool_use_id, message):
  return ContentBlock(type="tool_result", data={
    "tool_use_id": tool_use_id,
    "content": [ContentBlock(type="text", data={"text": message})],
    "is_error": True,
  })
