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
  the loop semantics.

  Parameters
  ----------
  agent : Agent
      Backend that streams turns and exposes tool specs.
  conversation : Conversation
      Conversation this session appends messages to.
  storage : object
      Attachment store used to persist emitted images.
  tools : object
      Tool registry providing specs and dispatch (builtin/skill/MCP).
  policy : object
      Approval policy resolving each tool name to allow/deny/ask.
  profile : object
      Active profile/configuration for the session.
  depth : int, optional
      Nesting depth; ``0`` is the top-level session, ``> 0`` a subagent.
  on_event : callable, optional
      Callback invoked with each emitted event. Defaults to a no-op.
  log : file-like, optional
      Stream for log output. Defaults to ``sys.stdout``.
  approval_queue : queue.Queue, optional
      Queue the worker parks on during the ``ask`` policy. External
      callers (the GUI runner) inject a queue they also push to; tests
      inject a ``queue.Queue``. A fresh queue is created when omitted.
  """

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

    # The question queue is the worker's blocking point while a
    # phenix_ask_user_question builtin is in flight (the API-backend
    # equivalent of ClaudeCodeAgent's SDK ask-user tool). It mirrors
    # approval_queue: the GUI runner flushes a _Cancelled sentinel into it
    # on cancel, and delivers the user's answers via submit_question_answer.
    self.question_queue = queue.Queue()
    self._pending_question_id = None

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
    # Drop any stale _Cancelled sentinel left in this (shared, long-lived)
    # approval queue by a Stop clicked during a PREVIOUS turn's streaming:
    # the runner pushes the sentinel while nobody is parked on
    # approval_queue.get(), so it is never consumed, and left in place it
    # would be picked up by this turn's first tool approval and cancel a
    # tool the user never denied. Real ToolApprovalResponses are preserved
    # and re-queued in order — a synchronous caller may legitimately
    # pre-seed an approval before run_turn.
    kept = []
    while True:
      try:
        item = self.approval_queue.get_nowait()
      except queue.Empty:
        break
      if not isinstance(item, _Cancelled):
        kept.append(item)
    for item in kept:
      self.approval_queue.put(item)
    # Same stale-_Cancelled drain for the question queue: a Stop clicked
    # during a previous turn's streaming flushes a sentinel here too, and
    # left in place it would abort this turn's first phenix_ask_user_question
    # the user never cancelled. Real pre-seeded answers are preserved.
    kept_q = []
    while True:
      try:
        item = self.question_queue.get_nowait()
      except queue.Empty:
        break
      if not isinstance(item, _Cancelled):
        kept_q.append(item)
    for item in kept_q:
      self.question_queue.put(item)
    self.conv.append(user_message)
    # Reconcile the conversation meta to the model/backend actually running
    # this turn. A conversation continued under a different model/backend
    # (e.g. relaunched with a new --model or --backend) then reflects the
    # current one; the per-message stamp preserves the full per-turn history.
    # Top-level only -- a subagent (depth>0) must not stamp the parent meta
    # with its own (possibly different) model/backend.
    if self.depth == 0:
      active_model = getattr(self.agent, "model", None)
      if active_model:
        self.conv.meta.model = active_model
      active_backend = getattr(self.profile, "backend", None)
      if active_backend:
        self.conv.meta.backend = active_backend
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
    # Stamp each assistant message with the model (from the agent) and the
    # backend (from the profile) that produced it. getattr guards keep this
    # robust for agents without a public .model (handled by exposing one)
    # and for sessions built without a profile (synchronous tests).
    msg = Message(role="assistant", content=[], timestamp=now(),
                  model=getattr(self.agent, "model", None),
                  backend=getattr(self.profile, "backend", None))
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
    """Park the worker on the approval queue until a decision arrives.

    Wakes on either a real ``ToolApprovalResponse`` or a ``_Cancelled``
    sentinel pushed by the GUI's cancel handler.

    Parameters
    ----------
    req : ToolApprovalRequest
        Request surfaced to the UI before blocking.

    Returns
    -------
    ToolApprovalResponse
        The user's approval decision.

    Raises
    ------
    TurnCancelled
        If a ``_Cancelled`` sentinel is received instead of a response.
    """
    self.on_event(req)                              # surface to UI
    response = self.approval_queue.get()
    if isinstance(response, _Cancelled):
      raise TurnCancelled()
    return response

  def _await_question_answer(self, request_id, questions):
    """Emit an ``AskUserQuestionRequested`` and park the worker for answers.

    The API-backend counterpart of ``ClaudeCodeAgent``'s SDK ask-user
    tool. Mirrors ``_await_approval``: surface the request to the UI, then
    block on ``question_queue`` until the GUI delivers the user's answers
    (via ``submit_question_answer``) or pushes a ``_Cancelled`` sentinel.

    Parameters
    ----------
    request_id : str
        Identifier correlating this request with the answers delivered
        back through ``submit_question_answer``.
    questions : list of dict
        The model's question specs (``question`` / ``options`` / ...).

    Returns
    -------
    object
        The answers object the GUI delivered.

    Raises
    ------
    TurnCancelled
        If a ``_Cancelled`` sentinel is received instead of answers.
    """
    from qttbx.widgets.chat.agent.events import AskUserQuestionRequested
    self._pending_question_id = request_id
    try:
      self.on_event(AskUserQuestionRequested(
        request_id=request_id, questions=list(questions or [])))
      answers = self.question_queue.get()
      if isinstance(answers, _Cancelled):
        raise TurnCancelled()
      return answers
    finally:
      self._pending_question_id = None

  def submit_question_answer(self, request_id, answers):
    """Deliver the user's answers to a parked ``phenix_ask_user_question``.

    Called from the GUI thread (via the runner) when the user submits a
    QuestionCard. Wakes the worker parked in ``_await_question_answer``.

    Parameters
    ----------
    request_id : str
        Identifier of the in-flight question.
    answers : object
        The user's answers to deliver to the parked handler.

    Returns
    -------
    bool
        ``True`` iff a question with this id is currently in flight (so the
        answers were queued); ``False`` otherwise.
    """
    if request_id != self._pending_question_id:
      return False
    self.question_queue.put(answers)
    return True

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
    """Record token usage from a nested session.

    Parameters
    ----------
    sub_id : str
        Identifier of the nested session the usage belongs to.
    usage : TokenUsage
        Token usage rolled up from that nested session.
    """
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
  """Build a short human-readable summary for the approval card.

  Best-effort; the UI shows the full input when expanded.

  Parameters
  ----------
  call : ToolUseRequested
      The tool-use request to summarize.

  Returns
  -------
  str
      A ``name(arg=value, ...)`` rendering with truncated values.
  """
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
