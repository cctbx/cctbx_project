"""Wrap ``AgentSession`` to emit Qt signals on the GUI thread.

The runner owns:

- the ``AgentSession`` (no Qt dependency itself)
- a worker ``QThread`` that runs ``session.run_turn``
- the ``CancelToken`` for the current turn
- the question queue (shared with the session) so the GUI Stop button can
  flush a ``_Cancelled`` sentinel into it; tool approvals instead cancel via
  ``session.approvals.cancel_turn()``
"""

import queue
import traceback

from qttbx.qt import QtCore

from qttbx.widgets.chat.agent.errors import AgentError, CancelToken
from qttbx.widgets.chat.agent.events import (
  AskUserQuestionRequested, ImageEmitted, ServerToolResult,
  ServerToolUsed, TextDelta, Thinking,
  ToolResultObserved, ToolResultsBatched,
  ToolUseRequested, TurnDone)
from qttbx.widgets.chat.agent.tools import (
  ToolApprovalRequest, _Cancelled)


class _Worker(QtCore.QObject):
  """Run a turn on a worker QThread and re-emit its events as Qt signals.

  Lives on a worker ``QThread``. Runs ``session.run_turn`` and re-emits its
  ``AgentEvent``\\ s as a single ``agent_event`` signal, then ``finished``.

  Notes
  -----
  The signal is named ``agent_event`` (not ``event``) to avoid shadowing
  ``QObject.event(QEvent)`` — PySide2 emits a spurious ``TypeError: native
  Qt signal is not callable`` from C++ when a Signal shares the name of an
  inherited slot, even when ``emit()`` succeeds.

  Parameters
  ----------
  session : AgentSession
      Session whose ``run_turn`` this worker drives.
  user_message : Message
      User message passed to ``run_turn`` for the turn.
  cancel : CancelToken
      Cancel token for the turn.
  """

  agent_event = QtCore.Signal(object)
  finished = QtCore.Signal()

  def __init__(self, session, user_message, cancel):
    super().__init__()
    self.session = session
    self.user_message = user_message
    self.cancel = cancel

  @QtCore.Slot()
  def run(self):
    """Drive the turn on the worker thread, emitting ``finished`` at end."""
    try:
      # The runner owns the session for the turn; override on_event to
      # forward through the worker's Qt signal (queued to the GUI thread).
      self.session.on_event = self.agent_event.emit
      try:
        self.session.run_turn(self.user_message, self.cancel)
      except Exception as exc:
        # Last-ditch: surface unexpected exceptions as an AgentError event
        # so the GUI sees them instead of crashing the worker silently.
        self.agent_event.emit(AgentError(
          message="Worker crash: %s" % exc,
          recoverable=False, kind="internal"))
        # The traceback is printed to the session log via session.log.
        try:
          print(traceback.format_exc(), file=self.session.log)
        except Exception:
          pass
    finally:
      self.finished.emit()


class QtAgentRunner(QtCore.QObject):
  """GUI-side wrapper for ``AgentSession``.

  Re-emits session events as Qt signals on the GUI thread. One runner per
  ``ChatWindow``, reused across turns. Only one turn may run at a time;
  ``start_turn`` while busy is a no-op (returns ``False``).

  Parameters
  ----------
  session : AgentSession
      Session this runner drives and forwards events from.
  parent : QtCore.QObject, optional
      Qt parent object.

  Notes
  -----
  Signals (all GUI-thread):

  - ``text_delta(str)``
  - ``thinking_delta(str)``
  - ``tool_use_requested(object)`` — ``ToolUseRequested`` or
    ``ToolApprovalRequest``
  - ``tool_results_batched(object)`` — ``ToolResultsBatched``
  - ``tool_result_observed(object)`` — ``ToolResultObserved``
  - ``server_tool_used(object)`` — ``ServerToolUsed``
  - ``server_tool_result(object)`` — ``ServerToolResult``
  - ``image_emitted(object)``
  - ``ask_user_question_requested(object)`` — ``AskUserQuestionRequested``
  - ``turn_done(object)`` — the full ``TurnDone`` event; carries
    ``stop_reason`` and the canonical ``finish`` (the window routes
    iteration-boundary vs terminal on ``finish``)
  - ``error(str, bool, str)`` — message, recoverable, kind
  - ``idle()`` — the turn's worker/thread have been torn down and
    ``is_busy()`` is False again. Emitted once per turn, after the turn's
    terminal event (``turn_done`` / ``error``); a slot may synchronously
    start the next turn. ``shutdown()`` does NOT emit it (teardown callers
    manage their own state).

  Thread affinity: ``start_turn``, ``cancel``, ``submit_approval``,
  ``submit_question_answer``, and ``wait_for_idle`` must be called from the
  GUI thread; they touch ``_thread`` / ``_worker`` / ``_cancel`` which are
  GUI-thread-owned.
  """

  text_delta = QtCore.Signal(str)
  thinking_delta = QtCore.Signal(str)
  tool_use_requested = QtCore.Signal(object)
  tool_results_batched = QtCore.Signal(object)
  tool_result_observed = QtCore.Signal(object)
  server_tool_used = QtCore.Signal(object)
  server_tool_result = QtCore.Signal(object)
  image_emitted = QtCore.Signal(object)
  ask_user_question_requested = QtCore.Signal(object)
  turn_done = QtCore.Signal(object)
  error = QtCore.Signal(str, bool, str)
  idle = QtCore.Signal()

  def __init__(self, session, parent=None):
    super().__init__(parent)
    self.session = session
    self._thread = None
    self._worker = None
    self._cancel = CancelToken()
    # Only the question queue still uses the sentinel-flush wake-up; approvals
    # now wake via ApprovalCoordinator.cancel_turn().
    self._question_cancel_queues = [self.session.question_queue]

  # ---- public API ----------------------------------------------------------

  def start_turn(self, user_message):
    """Kick off a new turn on the worker thread.

    Must be called from the GUI thread.

    Parameters
    ----------
    user_message : Message
        User message to run the turn with.

    Returns
    -------
    bool
        ``True`` if the turn was started, ``False`` if one is already
        running.
    """
    if self.is_busy():
      return False
    self._cancel = CancelToken()
    self._worker = _Worker(self.session, user_message, self._cancel)
    self._thread = QtCore.QThread()
    self._worker.moveToThread(self._thread)
    self._thread.started.connect(self._worker.run)
    self._worker.agent_event.connect(self._dispatch_event)
    self._worker.finished.connect(self._on_worker_finished)
    self._thread.start()
    return True

  def cancel(self):
    """Cancel the in-flight turn.

    Sets the cancel token, cancels any pending tool approval via
    ``session.approvals.cancel_turn()``, and flushes a ``_Cancelled``
    sentinel into the question queue so a worker blocked on ``queue.get()``
    wakes up. No-op when no turn is running. Must be called from the GUI
    thread.
    """
    if not self.is_busy():
      # No worker to cancel; pushing a sentinel now would mis-cancel the
      # next turn (the queue is shared with the session).
      return
    self._cancel.set()
    self.session.approvals.cancel_turn()
    for q in self._question_cancel_queues:      # question queue only
      # An unbounded queue.Queue; put_nowait never blocks.
      q.put_nowait(_Cancelled())

  def submit_approval(self, response):
    """Route a ``ToolApprovalResponse`` from the GUI to whoever asked.

    First persists any remember='tool' choice on the GUI thread --
    ``session.record_approval_remember`` (API backends) and
    ``agent.submit_approval`` (claude_code) -- BEFORE resolving, so a click
    that races a cancel still records "always allow". Then resolves the
    decision through the session's ApprovalCoordinator, but only while a turn
    is in flight: a late click after the turn ended is dropped (the is_busy
    gate) so it can't leak into the next turn. ``agent.submit_approval`` now
    returns False for every backend, so the coordinator is always the resolver.

    Must be called from the GUI thread.

    Parameters
    ----------
    response : ToolApprovalResponse
        The user's approval decision to route.
    """
    handled = False
    if self.is_busy():
      # A live turn: persist any remember='tool' on the GUI thread BEFORE the
      # coordinator resolves (API backends via the session, claude_code via
      # agent.submit_approval), so a click racing a cancel still records "always
      # allow", then resolve. An idle click on a stale card records NOTHING and
      # is dropped -- it must not grant a standing session auto-approve for a
      # tool that never ran (recording above this gate was finding-1's bug).
      self.session.record_approval_remember(response)
      self.session.agent.submit_approval(response)
      handled = self.session.approvals.submit(response)
    # A deny_and_stop cancels the in-flight turn -- but only when this response
    # belongs to it: the agent owned it (claude_code's can_use_tool callback),
    # or it answers the approval the session worker is parked on. A stale card's
    # stop from an abandoned earlier turn must NOT cancel an unrelated later
    # turn (the same hazard cancel()'s is_busy() guard addresses).
    if response.decision == "deny_and_stop" and (
        handled or self.session.owns_pending_approval(response.request_id)):
      self._cancel.set()

  def submit_question_answer(self, request_id, answers):
    """Route the user's answers to whoever asked the question.

    Tries ``agent.submit_question_answer`` first: the Claude Code backend
    runs its own loop and owns its own ask-user futures, so it resolves
    them itself. The API backends inherit the ``Agent`` base default
    (returns ``False``) and so fall through to the session's
    ``submit_question_answer`` -- the path a ``phenix_ask_user_question``
    builtin parks on. Every backend exposes the ``Agent`` contract (the API
    backends subclass it; claude_code duck-types it), so the method is always
    present (a direct call, mirroring ``submit_approval`` above).

    Must be called from the GUI thread.

    Parameters
    ----------
    request_id : str
        Identifier of the in-flight ``AskUserQuestionRequested``.
    answers : object
        The user's answers to forward.

    Returns
    -------
    bool
        ``True`` if either the agent or the session owned and handled the
        request, ``False`` otherwise.
    """
    if self.session.agent.submit_question_answer(request_id, answers):
      return True
    return self.session.submit_question_answer(request_id, answers)

  def is_busy(self):
    return self._thread is not None and self._thread.isRunning()

  def wait_for_idle(self, timeout_ms=5000):
    """Block until the current turn finishes.

    Called by the GUI on the conversation-switch and close teardown paths
    (ChatWindow) and by tests. Must be called from the GUI thread.

    Parameters
    ----------
    timeout_ms : int, optional
        Maximum time to wait, in milliseconds.
    """
    if self._thread is None:
      return
    self._thread.wait(timeout_ms)

  def shutdown(self, timeout_ms=2000):
    """Synchronously stop and join the worker thread from the GUI thread.

    Unlike ``wait_for_idle`` -- which only ``wait()``s and relies on the queued
    ``_on_worker_finished`` slot to ``quit()`` the thread's ``exec()`` loop, a
    slot that cannot run while the GUI thread is blocked in that ``wait()`` --
    this quits the thread's event loop itself, then joins it. Call it on
    teardown (``closeEvent`` / conversation switch, after ``cancel()``) so the
    ``QThread`` is actually finished before the runner/window is destroyed;
    otherwise ``wait()`` blocks the full timeout and the still-running thread
    aborts with "QThread: Destroyed while thread is still running". Idempotent;
    must be called from the GUI thread.
    """
    thread = self._thread
    if thread is None:
      return
    thread.quit()               # ask the thread's exec() loop to exit
    thread.wait(timeout_ms)     # then join it (bounded)
    thread.deleteLater()
    worker = self._worker
    if worker is not None:
      worker.deleteLater()
    self._thread = None
    self._worker = None

  # ---- event routing -------------------------------------------------------

  def _dispatch_event(self, ev):
    if isinstance(ev, TextDelta):
      self.text_delta.emit(ev.text)
    elif isinstance(ev, Thinking):
      self.thinking_delta.emit(ev.text)
    elif isinstance(ev, ToolApprovalRequest):
      self.tool_use_requested.emit(ev)
    elif isinstance(ev, ToolUseRequested):
      self.tool_use_requested.emit(ev)
    elif isinstance(ev, ToolResultsBatched):
      self.tool_results_batched.emit(ev)
    elif isinstance(ev, ToolResultObserved):
      self.tool_result_observed.emit(ev)
    elif isinstance(ev, ServerToolUsed):
      self.server_tool_used.emit(ev)
    elif isinstance(ev, ServerToolResult):
      self.server_tool_result.emit(ev)
    elif isinstance(ev, ImageEmitted):
      self.image_emitted.emit(ev)
    elif isinstance(ev, AskUserQuestionRequested):
      self.ask_user_question_requested.emit(ev)
    elif isinstance(ev, TurnDone):
      self.turn_done.emit(ev)
    elif isinstance(ev, AgentError):
      self.error.emit(ev.message, ev.recoverable, ev.kind or "")

  def _on_worker_finished(self):
    """Tear down the finished worker/thread and drain the question queue.

    Runs on the GUI thread when the turn's worker emits ``finished``.
    Draining here removes any question answer left by a click during
    teardown -- submitted while ``is_busy()`` was still true (so the
    submit_approval guard let it through) but after the worker had stopped
    consuming. The next turn cannot start until this slot clears
    ``is_busy()``, so the stray item is gone before it could be applied to
    the next turn's first approval.
    """
    if self._thread is not None:
      self._thread.quit()
      self._thread.wait()
      self._thread.deleteLater()
    if self._worker is not None:
      self._worker.deleteLater()
    self._thread = None
    self._worker = None
    for q in self._question_cancel_queues:
      while True:
        try:
          q.get_nowait()
        except queue.Empty:
          break
    # Announce end-of-turn teardown LAST -- is_busy() is False and the stale
    # question answers are drained, so a slot may synchronously start the next
    # turn (the chat window fires its deferred auth retry from here instead of
    # polling is_busy() on a timer).
    self.idle.emit()
