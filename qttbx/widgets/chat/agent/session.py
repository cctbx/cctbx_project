"""``AgentSession`` — multi-iteration tool loop orchestrator.

Wraps any ``Agent``. Pure-Python; no Qt. ``QtAgentRunner`` (in
``runner.py``) wraps this to emit Qt signals on the GUI thread.
"""

import json
import queue
import sys
import time
import uuid

from libtbx.utils import Sorry

from qttbx.widgets.chat.agent.conversation import (
  Conversation, ContentBlock, Message, now)
from qttbx.widgets.chat.agent.errors import TurnCancelled
from qttbx.widgets.chat.agent.events import (
  AgentError, ImageEmitted, ServerToolResult, ServerToolUsed,
  Thinking, TextDelta, ToolResultObserved, ToolResultsBatched,
  ToolUseRequested, TurnDone, TokenUsage as TokenUsageEvent)
from qttbx.widgets.chat.agent.tools import (
  ToolApprovalRequest, _Cancelled)


# Events that add content to the in-progress assistant message. The mid-turn
# autosave fires only on these -- never on the terminal TurnDone, whose
# GIL-releasing save I/O would race the GUI turn-end save (see
# _collect_one_response).
_CONTENT_EVENT_TYPES = (
  TextDelta, Thinking, ToolUseRequested, ServerToolUsed, ServerToolResult,
  ImageEmitted)


def _new_id(prefix=""):
  return prefix + uuid.uuid4().hex[:12]


def _drain_stale_cancelled(q):
  """Drop stale ``_Cancelled`` sentinels from ``q``, keeping real items.

  A Stop clicked during a PREVIOUS turn's streaming pushes a sentinel
  while nobody is parked on the (shared, long-lived) queue, so it is
  never consumed; left in place it would cancel this turn's first park
  (a tool approval or a ``phenix_ask_user_question``) that the user
  never denied. Real items are preserved and re-queued in order -- a
  synchronous caller may legitimately pre-seed one before ``run_turn``.
  """
  kept = []
  while True:
    try:
      item = q.get_nowait()
    except queue.Empty:
      break
    if not isinstance(item, _Cancelled):
      kept.append(item)
  for item in kept:
    q.put(item)


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
      Queue the worker parks on during the ``ask`` policy. A fresh queue is
      created here and exposed as ``self.approval_queue``; the GUI runner and
      tests push the user's decision onto it. No caller injects one today,
      though the parameter permits it.
  """

  def __init__(self, agent, conversation, storage, tools, policy,
               profile, depth=0, on_event=None, log=None,
               approval_queue=None, autosave_interval_s=5.0, clock=None):
    self.agent = agent
    self.conv = conversation
    self.storage = storage
    self.tools = tools
    self.policy = policy
    self.profile = profile
    self.depth = depth
    self.autosave_interval_s = autosave_interval_s
    self._clock = clock or time.monotonic
    self._last_autosave = 0.0
    self.on_event = on_event or (lambda ev: None)
    self.log = log if log is not None else sys.stdout

    # The approval queue is the worker's blocking point during 'ask' policy.
    # A fresh queue is created here (no caller injects one today); the GUI
    # runner and tests push the user's decision onto self.approval_queue.
    self.approval_queue = approval_queue or queue.Queue()

    # The question queue is the worker's blocking point while a
    # phenix_ask_user_question builtin is in flight (the API-backend
    # equivalent of ClaudeCodeAgent's SDK ask-user tool). It mirrors
    # approval_queue: the GUI runner flushes a _Cancelled sentinel into it
    # on cancel, and delivers the user's answers via submit_question_answer.
    self.question_queue = queue.Queue()
    self._pending_question_id = None
    self._pending_approval_id = None

    # Set when run_turn starts; consulted by tool handlers.
    self.cancel = None
    self.sub_id = _new_id("sa_") if depth > 0 else None
    self.started_at = now()
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
    self._last_autosave = self._clock()            # reset mid-turn autosave
    # Drop any stale _Cancelled sentinel a previous turn's Stop left in
    # the shared queues (see _drain_stale_cancelled).
    _drain_stale_cancelled(self.approval_queue)
    _drain_stale_cancelled(self.question_queue)
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
      assistant_msg, tool_calls, answered = self._collect_one_response(cancel)
      self.conv.append(assistant_msg)

      if assistant_msg.stop_reason != "tool_use":
        # A tool_use can ride a non-"tool_use" stop. Two backends reach here:
        #   - claude_code runs its OWN tool loop in the SDK subprocess, so its
        #     tool_use blocks arrive as ToolUseRequested and each real result as
        #     a ToolResultObserved that _accumulate has ALREADY turned into an
        #     answering tool_result in ``answered`` -- so those tools are never
        #     pending here; they are answered with their real result (not a bogus
        #     cancel that used to surface on resume). A tool_use with NO answer
        #     genuinely did not run (the turn was cancelled between its tool_use
        #     and its result) and is answered "Cancelled by user.".
        #   - Anthropic streams the tool_use block and only THEN surfaces a
        #     cancel (trailing TurnDone(stop_reason='cancelled')): nothing was
        #     observed, so every tool_use is answered "Cancelled by user.".
        # Either way the tool_uses must be answered, or the saved transcript
        # orphans them -> non-recoverable replay 400. (Backends that batch tool
        # calls only after their stream loop never reach this with tool_calls.)
        self._answer_pending_tool_uses(
          tool_calls, "Cancelled by user.", answered=answered)
        return assistant_msg
      if cancel.is_set():
        # Cancel set AFTER a tool_use-stop response but before dispatch (the
        # narrower window where stop_reason stays 'tool_use'): answer the
        # pending tool_uses so the saved transcript isn't orphaned.
        self._answer_pending_tool_uses(tool_calls, "Cancelled by user.")
        assistant_msg.stop_reason = "cancelled"
        return assistant_msg

      if max_turns is not None and iterations >= max_turns:
        # Stopped at the turn cap before dispatching this turn's tool_uses:
        # answer them so the tool_use isn't orphaned, then the human-visible
        # cap marker.
        self._answer_pending_tool_uses(
          tool_calls, "[Subagent stopped at turn cap (%d)]" % max_turns)
        self.conv.append(Message(
          role="user", timestamp=now(),
          content=[ContentBlock(type="text", data={
            "text": "[Subagent stopped at turn cap (%d)]" % max_turns})]))
        return assistant_msg

      tool_result_msg = self._dispatch_and_build_results(tool_calls, cancel)
      self.conv.append(tool_result_msg)
      self._maybe_autosave()                       # checkpoint after tool result

      # Dispatch may have set the cancel token (deny_and_stop or sentinel
      # cancellation during approval). Don't ask the model for another
      # response — that would either re-enter a cancelled loop or, in
      # FakeAgent tests, exhaust the scripted turns.
      if cancel.is_set():
        assistant_msg.stop_reason = "cancelled"
        return assistant_msg

  def _maybe_autosave(self, partial_msg=None):
    """Throttled mid-turn persistence (top-level sessions only).

    Writes a non-mutating snapshot -- the committed messages plus the
    in-progress assistant message -- at most once per ``autosave_interval_s``.
    ``storage.save`` trims a snapshot that would end in an unanswered ``tool_use``
    to the committed prefix, so a crash never freezes an un-resumable transcript
    on disk. ``reindex=False`` keeps the O(N) index rescan (and the worker's
    ``index.json`` writes) out of the hot loop. Never call this on the terminal
    ``TurnDone`` event; see ``_collect_one_response``.

    Parameters
    ----------
    partial_msg : Message, optional
        The assistant message currently being streamed (not yet appended to the
        conversation). Included in the snapshot when it has content.
    """
    if self.depth != 0 or self.storage is None:
      return
    now_t = self._clock()
    if now_t - self._last_autosave < self.autosave_interval_s:
      return
    self._last_autosave = now_t                    # throttle even if save raises
    tail = [partial_msg] if (partial_msg is not None and partial_msg.content) \
        else []
    # storage.save trims a trailing unanswered tool_use, so a snapshot caught
    # mid-tool_use (before dispatch appends the tool_result) is persisted as a
    # resumable prefix rather than an un-resumable provider-400 transcript.
    snapshot = Conversation(meta=self.conv.meta,
                            messages=self.conv.messages + tail)
    try:
      self.storage.save(snapshot, reindex=False)
    except Exception as exc:
      print("chat: mid-turn autosave failed: %s" % exc, file=self.log)

  def _collect_one_response(self, cancel):
    # Stamp each assistant message with the model (from the agent) and the
    # backend (from the profile) that produced it. getattr guards keep this
    # robust for agents without a public .model (handled by exposing one)
    # and for sessions built without a profile (synchronous tests).
    msg = Message(role="assistant", content=[], timestamp=now(),
                  model=getattr(self.agent, "model", None),
                  backend=getattr(self.profile, "backend", None))
    tool_calls = []
    answered = []
    for event in self.agent.stream_turn(self.conv, self.tools.specs(), cancel):
      self.on_event(event)
      _accumulate(msg, event, tool_calls, answered, self.storage,
                  self.conv.meta.id)
      # Throttled mid-turn checkpoint -- ONLY on content events. Never on the
      # terminal TurnDone: _collect_one_response returns straight into
      # conv.append(assistant_msg), and autosave's blocking I/O (which releases
      # the GIL) between the queued TurnDone emit and that append would let the
      # GUI's _on_turn_done save a conv.messages still missing this message.
      if isinstance(event, _CONTENT_EVENT_TYPES):
        self._maybe_autosave(msg)
    return msg, tool_calls, answered

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
        decision = self._resolve_and_approve(call, batch_id)
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
      # Honor a handler result that carries its own is_error (an MCP
      # McpToolResult flags failures with is_error=True instead of raising,
      # so a failed MCP tool must reach the model + UI as an error, not a
      # success). Builtin/skill results are str/bytes/dict with no is_error
      # attribute, so they default to False -- their failures already arrive
      # as Sorry/Exception and were turned into _tool_error_block above.
      is_error = bool(getattr(handler_result, "is_error", False))
      result_blocks.append(ContentBlock(type="tool_result", data={
        "tool_use_id": call.id,
        "content": blocks,
        "is_error": is_error,
      }))

    self.on_event(ToolResultsBatched(blocks=result_blocks))
    return Message(role="user", content=result_blocks, timestamp=now())

  def _answer_pending_tool_uses(self, tool_calls, reason, answered=None):
    """Append one user message answering every tool_use of the just-ended turn.

    A turn that ends right after an assistant tool_use block -- Stop pressed
    before dispatch, the subagent turn-cap, or a claude_code turn (which runs
    tools in its SDK subprocess) -- must not leave an UNMATCHED tool_use in the
    saved transcript: on the next turn that transcript is replayed verbatim, and
    an assistant tool_use with no following tool_result is a non-recoverable
    provider 400 (anthropic / openai / portkey / gemini) that permanently breaks
    the conversation. (claude_code resumes via its own SDK session and is immune
    to the 400, but an unanswered tool_use would still be trimmed off the saved
    transcript by ``persistable_prefix`` and lost from the rebuilt view.)

    ``answered`` is the ordered list of tool_result blocks ``_accumulate``
    already built from this turn's observed results (claude_code tools that
    actually ran, each carrying its real content + is_error). A tool_use with no
    answer genuinely did not run -- a pre-dispatch cancel, a mid-turn cancel
    before its result arrived, or the turn cap -- and gets an is_error result
    carrying ``reason``. Blocks are emitted in tool_use order. No-op when there
    is nothing to answer.
    """
    by_id = {}
    for b in answered or []:
      by_id[b.data.get("tool_use_id")] = b
    blocks = [by_id.get(c.id) or _tool_error_block(c.id, reason)
              for c in tool_calls]
    if not blocks:
      return
    self.conv.append(Message(role="user", timestamp=now(), content=blocks))

  def _resolve_and_approve(self, call, batch_id):
    policy = self.policy.resolve(call.name)
    # A tool that opted out of standing approval (allow_remember=False -- the
    # coot force-close) ALWAYS confirms via a card, regardless of a permissive
    # tool_policy. Scoped to the opt-out flag, NOT risk, so a destructive tool
    # that still offers "Always allow" keeps it. Only the blanket --auto-approve
    # (handled in the chat window) may skip the card.
    if policy == "allow" and not self.tools.allow_remember_of(call.name):
      policy = "ask"
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
      batch_id=batch_id,
      allow_remember=self.tools.allow_remember_of(call.name),
    )
    response = self._await_approval(req)
    if response.decision == "approve":
      # Re-check allow_remember_of, not just the card's remember flag: keep the
      # "never per-tool auto-approve" guarantee in the data layer, not only the
      # UI that suppresses the checkbox.
      if (response.remember == "tool"
          and self.tools.allow_remember_of(call.name)):
        self.policy.allow_tool_for_session(call.name)
      return "approve"
    return response.decision

  def _await_approval(self, req):
    """Park the worker on the approval queue until a decision arrives.

    Wakes on either a real ``ToolApprovalResponse`` or a ``_Cancelled``
    sentinel pushed by the GUI's cancel handler. A response whose
    ``request_id`` does not match ``req`` is a stale click from an abandoned
    earlier turn's card and is discarded, so it can't be misapplied to this
    request (mirrors ``submit_question_answer``'s id check).

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
    self._pending_approval_id = req.request_id
    try:
      self.on_event(req)                              # surface to UI
      while True:
        response = self.approval_queue.get()
        if isinstance(response, _Cancelled):
          raise TurnCancelled()
        if getattr(response, "request_id", None) == req.request_id:
          return response
        # A response for a DIFFERENT request is a stale card click from an
        # abandoned earlier turn (the approval path is is_busy-gated, so it can
        # be queued while this later turn runs). Drop it -- returning it would
        # misapply an unrelated decision to THIS tool. Batch siblings are queued
        # and consumed in the same dispatch order, so the only non-matching item
        # reached here is a genuinely stale one, never a sibling still needed.
    finally:
      self._pending_approval_id = None

  def owns_pending_approval(self, request_id):
    """Whether the worker is currently parked in ``_await_approval`` waiting for
    this exact ``request_id``.

    Lets the runner gate a ``deny_and_stop`` cancel on a matching request so a
    stale card's stop can't cancel an unrelated later turn.
    """
    return request_id is not None and request_id == self._pending_approval_id

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
      return self.tools.invoke_builtin(
        call.name, call.input,
        cancel=cancel, session=self, tool_use_id=call.id)
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
      return [_mcp_item_to_block(item) for item in result.content]
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

def _accumulate(msg, event, tool_calls, answered, storage, conv_id):
  if isinstance(event, TextDelta):
    _append_text(msg, event.text)
  elif isinstance(event, Thinking):
    _append_thinking(msg, event.text, event.signature)
  elif isinstance(event, ToolUseRequested):
    msg.content.append(ContentBlock(type="tool_use", data={
      "id": event.id, "name": event.name, "input": event.input}))
    tool_calls.append(event)
  elif isinstance(event, ToolResultObserved):
    # A backend that runs its OWN tool loop (claude_code, in the SDK subprocess)
    # reports each executed result here instead of the session dispatching it.
    # The matching tool_use rode in as a ToolUseRequested (pending in tool_calls)
    # but the turn ends 'end_turn' (not 'tool_use'), so build its answering
    # tool_result NOW -- carrying the observed is_error -- into ``answered``.
    # run_turn appends these as the tool_use's real answer, so the tool is never
    # left pending to be fabricated as a "Cancelled by user." error on resume.
    answered.append(ContentBlock(type="tool_result", data={
      "tool_use_id": event.tool_use_id,
      "content": _observed_content_blocks(event.content),
      "is_error": bool(getattr(event, "is_error", False))}))
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
    msg.usage = event.to_stored()
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


def _tool_error_block(tool_use_id, message):
  return ContentBlock(type="tool_result", data={
    "tool_use_id": tool_use_id,
    "content": [ContentBlock(type="text", data={"text": message})],
    "is_error": True,
  })


def _is_binary_block(item):
  """True for a content block reduced to a compact ``[<type>]`` placeholder
  rather than JSON-inlined: any image (regardless of source), or any block with
  a base64 ``source`` (a document / PDF, ...). Keeps a multi-MB base64 blob out
  of the persisted transcript; an image is placeholdered even when url-sourced,
  since its bytes are not needed in a text transcript."""
  if item.get("type") == "image":
    return True
  source = item.get("source")
  return isinstance(source, dict) and source.get("type") == "base64"


def _observed_content_blocks(content):
  """Canonical tool_result content blocks from a claude_code observed result.

  The SDK surfaces a tool result's content as a plain string, a list of
  Anthropic content dicts (text / image / document / tool_reference / ...), or
  None. Text is preserved verbatim; a base64-bearing block (image, document,
  ...) is reduced to a compact ``[<type>]`` placeholder -- so the persisted
  transcript stays text-sized and never inlines a base64 blob (the bytes are
  re-derivable from the claude_code transcript on resume; images are also
  surfaced live via ImageEmitted). Any other structured item is JSON-encoded,
  and an empty or None result becomes a single empty text block.
  """
  if content is None:
    return [ContentBlock(type="text", data={"text": ""})]
  if isinstance(content, str):
    return [ContentBlock(type="text", data={"text": content})]
  if isinstance(content, list):
    blocks = []
    for item in content:
      if isinstance(item, dict) and item.get("type") == "text":
        blocks.append(ContentBlock(type="text",
                                   data={"text": item.get("text", "")}))
      elif isinstance(item, dict) and _is_binary_block(item):
        blocks.append(ContentBlock(
          type="text", data={"text": "[%s]" % (item.get("type") or "binary")}))
      else:
        blocks.append(ContentBlock(type="text",
                                   data={"text": json.dumps(item, default=str)}))
    return blocks or [ContentBlock(type="text", data={"text": ""})]
  return [ContentBlock(type="text",
                       data={"text": json.dumps(content, default=str)})]
