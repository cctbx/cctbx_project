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
  ToolApprovalRequest, record_tool_remember, _Cancelled)


# Events that add content to the in-progress assistant message. The mid-turn
# autosave fires only on these -- never on the terminal TurnDone, whose
# GIL-releasing save I/O would race the GUI turn-end save (see
# _collect_one_response).
_CONTENT_EVENT_TYPES = (
  TextDelta, Thinking, ToolUseRequested, ServerToolUsed, ServerToolResult,
  ImageEmitted)


# Two rapid iteration checkpoints inside this window coalesce into one write:
# storage.save is a full-file rewrite, so a burst of sub-250ms tool calls
# must not rewrite messages.json per call. Freshness target is per-iteration,
# not per-write; a hard kill inside the window loses only those sub-250ms
# iterations until the next save heals.
_CHECKPOINT_FLOOR_S = 0.25


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
  approvals : ApprovalCoordinator, optional
      Shared approval coordinator; a fresh one is created when not injected.
  """

  def __init__(self, agent, conversation, storage, tools, policy,
               profile, depth=0, on_event=None, log=None,
               approvals=None, autosave_interval_s=5.0, clock=None):
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
    self._last_checkpoint = None       # None: first checkpoint always writes
    self.on_event = on_event or (lambda ev: None)
    self.log = log if log is not None else sys.stdout

    # Cancel-safe approval bookkeeping shared with the claude_code backend
    # (chat_window builds one coordinator and injects it into both). A fresh
    # one is created when no caller injects it (tests, headless subagents).
    from qttbx.widgets.chat.agent.approval import ApprovalCoordinator
    self.approvals = approvals or ApprovalCoordinator()
    # request_id -> tool name for approvals surfaced by _resolve_and_approve, so
    # record_approval_remember can persist a remember='tool' choice on the GUI
    # thread. Cleared at each turn's start.
    self._pending_approval_calls = {}

    # The question queue is the worker's blocking point while a
    # phenix_ask_user_question builtin is in flight (the API-backend
    # equivalent of ClaudeCodeAgent's SDK ask-user tool). The GUI runner
    # flushes a _Cancelled sentinel into it on cancel, and delivers the
    # user's answers via submit_question_answer.
    self.question_queue = queue.Queue()
    self._pending_question_id = None

    # Set when run_turn starts; consulted by tool handlers.
    self.cancel = None
    self.sub_id = _new_id("sa_") if depth > 0 else None
    self.started_at = now()
    # Push-based rollup of nested-session token usage. Keyed by sub_id.
    self._subagent_usage_by_id = {}

  # ---- main loop -----------------------------------------------------------

  def run_turn(self, user_message, cancel, max_turns=None):
    """Run one turn, then guarantee the approval registry is left clean.

    A turn that ends with an approval still open -- e.g. a subprocess-death
    claude_code turn whose ``_on_can_use_tool`` future never resolved -- must
    not leave a pending future in the coordinator: it would keep
    ``owns_pending_approval`` True while idle (so a stale Deny+Stop could set the
    idle cancel token) and destroy a still-pending SDK task at exit.
    ``cancel_turn`` is a no-op on a clean turn."""
    try:
      return self._run_turn(user_message, cancel, max_turns)
    finally:
      self.approvals.cancel_turn()

  def _run_turn(self, user_message, cancel, max_turns=None):
    """Append the user message and iterate the tool-use loop.

    Loop terminates on ``stop_reason != 'tool_use'``, on cancellation,
    or when the optional ``max_turns`` cap is hit. When the cap fires,
    a visible ``[Subagent stopped at turn cap (N)]`` marker is appended
    so the model (or a viewer) sees why the loop ended.

    Cancel handling: if the cancel token is set after a tool dispatch
    (``deny_and_stop``, or the coordinator cancelling the pending approval
    future during a Stop), we do NOT ask the model for another response
    — we finalize the
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
    self._last_checkpoint = None                   # floor is per-turn: a new
                                                   # turn's first iteration
                                                   # always writes
    # Turn boundary: abandon anything a previous (cancelled) turn left pending.
    # Replaces the old _drain_stale_cancelled(approval_queue). Clear the
    # remember map at ENTRY (not turn end) so a late click after the prior turn
    # still recorded its choice.
    self.approvals.begin_turn()
    self._pending_approval_calls.clear()
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
      (assistant_msg, tool_calls, answered,
       last_committed) = self._collect_one_response(cancel)
      # An empty residual after a per-iteration commit (the turn ended right
      # on a tool result) must NOT be appended: persistable_prefix would trim
      # it on save, silently dropping the terminal stop_reason/usage riding
      # it. Fold both onto the last committed assistant (same object already
      # in conv.messages) instead -- usage only when present, never
      # clobbering with None. A tool_use-bearing or text residual, and any
      # turn with no commit (API backends; error-before-content), keeps the
      # existing append. Mutating an already-appended message after the
      # agent's terminal TurnDone was queued shares the race discipline of
      # _collect_one_response's autosave note: no blocking I/O intervenes
      # before these fields land, exactly as the bulk append below always
      # relied on -- and the at-risk window here is two fields, not a whole
      # message.
      folded = (last_committed is not None
                and not assistant_msg.content
                and assistant_msg.stop_reason != "tool_use")
      if folded:
        if assistant_msg.stop_reason:
          last_committed.stop_reason = assistant_msg.stop_reason
        if assistant_msg.usage is not None:
          last_committed.usage = assistant_msg.usage
      else:
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
        return last_committed if folded else assistant_msg
      if cancel.is_set():
        # Cancel set AFTER a tool_use-stop response but before dispatch (the
        # narrower window where stop_reason stays 'tool_use'): answer the
        # pending tool_uses so the saved transcript isn't orphaned.
        self._answer_pending_tool_uses(tool_calls, "Cancelled by user.")
        assistant_msg.stop_reason = "cancelled"
        # This return path emits no agent TurnDone (the last one said
        # finish='tool_use'), so consumers keying turn end on the terminal
        # event -- the GUI composer unlock, headless disposition -- would
        # wait forever. Synthesize it. Literal "cancelled": the canonical
        # CANCELLED value documented on TurnDone (qttbx must not import the
        # phenix finish module). MUST stay AFTER the appends above -- the
        # GUI's turn-end save races the worker (see _collect_one_response's
        # autosave note), and emitting first would let it persist a
        # conversation missing this turn's answers.
        self.on_event(TurnDone(stop_reason="cancelled", finish="cancelled"))
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
      self._checkpoint_iteration()                 # per-iteration disk cadence

      # Dispatch may have set the cancel token (deny_and_stop or sentinel
      # cancellation during approval). Don't ask the model for another
      # response — that would either re-enter a cancelled loop or, in
      # FakeAgent tests, exhaust the scripted turns.
      if cancel.is_set():
        assistant_msg.stop_reason = "cancelled"
        # Same synthesized terminal as the pre-dispatch short-circuit above,
        # same ordering constraint: after this iteration's appends.
        self.on_event(TurnDone(stop_reason="cancelled", finish="cancelled"))
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

  def _checkpoint_iteration(self):
    """Immediate (floored, unthrottled) persistence of a completed tool
    iteration -- the per-iteration disk cadence every backend shares. Same
    failure contract and depth gate as ``_maybe_autosave``; ``reindex=False``
    for the same reason. Also advances the autosave clock so the throttled
    path doesn't double-write the snapshot it just wrote."""
    if self.depth != 0 or self.storage is None:
      return
    now_t = self._clock()
    if (self._last_checkpoint is not None
        and now_t - self._last_checkpoint < _CHECKPOINT_FLOOR_S):
      return
    self._last_checkpoint = now_t
    self._last_autosave = now_t
    snapshot = Conversation(meta=self.conv.meta,
                            messages=list(self.conv.messages))
    try:
      self.storage.save(snapshot, reindex=False)
    except Exception as exc:
      print("chat: iteration checkpoint failed: %s" % exc, file=self.log)

  def _new_assistant_msg(self):
    """A fresh in-progress assistant message stamped with the model (from the
    agent) and backend (from the profile) producing this turn."""
    return Message(role="assistant", content=[], timestamp=now(),
                   model=getattr(self.agent, "model", None),
                   backend=getattr(self.profile, "backend", None))

  def _collect_one_response(self, cancel):
    msg = self._new_assistant_msg()
    tool_calls = []
    answered = []
    last_committed = None
    for event in self.agent.stream_turn(self.conv, self.tools.specs(), cancel):
      self.on_event(event)
      _accumulate(msg, event, tool_calls, answered, self.storage,
                  self.conv.meta.id)
      # A backend that runs its OWN tool loop (claude_code) streams the whole
      # multi-iteration turn through one stream_turn call. Commit each
      # completed iteration -- every pending tool_use answered by an observed
      # result -- as the same (assistant, tool_results) message pair the API
      # loop appends, so the persisted shape is backend-uniform and the
      # checkpoint below gives claude_code the same per-iteration disk
      # cadence. API backends never emit ToolResultObserved, so this branch
      # cannot fire for them.
      if isinstance(event, ToolResultObserved) and tool_calls:
        answered_ids = set(b.data.get("tool_use_id") for b in answered)
        if all(c.id in answered_ids for c in tool_calls):
          msg.stop_reason = "tool_use"
          self.conv.append(msg)
          self.conv.append(Message(
            role="user", timestamp=now(),
            content=_ordered_tool_result_blocks(
              tool_calls, answered, "Result not observed.")))
          last_committed = msg
          msg = self._new_assistant_msg()
          tool_calls = []
          answered = []
          self._checkpoint_iteration()
      # Throttled mid-turn checkpoint -- ONLY on content events. Never on the
      # terminal TurnDone: _collect_one_response returns straight into
      # conv.append(assistant_msg), and autosave's blocking I/O (which releases
      # the GIL) between the queued TurnDone emit and that append would let the
      # GUI's _on_turn_done save a conv.messages still missing this message.
      if isinstance(event, _CONTENT_EVENT_TYPES):
        self._maybe_autosave(msg)
    return msg, tool_calls, answered, last_committed

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
    blocks = _ordered_tool_result_blocks(tool_calls, answered, reason)
    if not blocks:
      return
    self.conv.append(Message(role="user", timestamp=now(), content=blocks))

  def _resolve_and_approve(self, call, batch_id):
    policy = self.policy.resolve(call.name)
    # A tool that opted out of standing approval (allow_remember=False -- the
    # coot force-close) ALWAYS confirms via a card, regardless of a permissive
    # tool_policy. Scoped to the opt-out flag, NOT risk, so a destructive tool
    # that still offers "Always allow" keeps it. Even the blanket --auto-approve
    # honors this floor (the chat window shows the card for an allow_remember=
    # False tool rather than auto-approving it).
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
    # Map request -> tool so record_approval_remember (called by the runner on
    # the GUI thread, BEFORE the coordinator resolves this future) can persist a
    # remember='tool' choice even when a click races a cancel that wakes this
    # _await_approval cancelled. Cleared at the next turn's start.
    self._pending_approval_calls[req.request_id] = call.name
    response = self._await_approval(req)
    if response.decision == "approve":
      return "approve"
    return response.decision

  def record_approval_remember(self, response):
    """Persist a remember='tool' choice for an API-backend approval by wiring the
    session's pending-calls map, the registry opt-out recheck, and the policy's
    session-allow set into the shared ``record_tool_remember`` contract (the same
    one ``ClaudeCodeAgent.submit_approval`` uses). No-op for claude_code, whose
    SDK tools never populate ``_pending_approval_calls`` (so ``tool_name_of``
    returns None). Called by the runner on the GUI thread BEFORE the coordinator
    resolves the future, so a click racing a cancel still records it."""
    record_tool_remember(
      response,
      tool_name_of=self._pending_approval_calls.get,
      allow_remember=self.tools.allow_remember_of,
      remember=self.policy.allow_tool_for_session)

  def _await_approval(self, req):
    """Block the worker on the shared coordinator until a decision arrives.

    Surfaces ``req`` via ``on_event``, then parks on the coordinator's future.
    A stale click from an abandoned earlier turn is dropped by the coordinator
    (its request already left the registry at the turn boundary); a cancel
    (runner Stop / next turn) cancels the future.

    Raises
    ------
    TurnCancelled
        If the turn is cancelled before a decision arrives.
    """
    import concurrent.futures
    # A Stop that lands between the dispatch loop's cancel check and here is
    # handled inside open(): a closed coordinator (cancel_turn ran) returns a
    # pre-cancelled future without surfacing a card, so fut.result() unwinds
    # into TurnCancelled at once instead of parking on a future a Stop passed.
    fut = self.approvals.open(req, self.on_event)
    try:
      return fut.result()
    except concurrent.futures.CancelledError:
      raise TurnCancelled()

  def owns_pending_approval(self, request_id):
    """Whether the worker is currently parked awaiting this exact
    ``request_id`` (lets the runner gate a deny_and_stop cancel)."""
    return self.approvals.owns(request_id)

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


def _ordered_tool_result_blocks(tool_calls, answered, fallback_reason):
  """The answering tool_result blocks for ``tool_calls``, in tool_use order:
  each call's observed block from ``answered`` when present, else an is_error
  block carrying ``fallback_reason``. Single assembly shared by the
  per-iteration commit and ``_answer_pending_tool_uses`` so the two can't
  drift."""
  by_id = {}
  for b in answered or []:
    by_id[b.data.get("tool_use_id")] = b
  return [by_id.get(c.id) or _tool_error_block(c.id, fallback_reason)
          for c in tool_calls]


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
