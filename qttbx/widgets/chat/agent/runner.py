"""QtAgentRunner — wrap AgentSession to emit Qt signals on the GUI thread.

Section 4.5 of the design spec. The runner owns:
  - the AgentSession (no Qt dependency itself)
  - a worker QThread that runs session.run_turn
  - the CancelToken for the current turn
  - the approval queue (shared with the session) so the GUI Stop button can
    flush a _Cancelled sentinel into it (Section 10.4)
"""

import traceback

from qttbx.qt import QtCore

from qttbx.widgets.chat.agent.errors import AgentError, CancelToken
from qttbx.widgets.chat.agent.events import (
  ImageEmitted, TextDelta, Thinking, TokenUsage as TokenUsageEvent,
  ToolResultsBatched, ToolUseRequested, TurnDone)
from qttbx.widgets.chat.agent.tools import (
  ToolApprovalRequest, _Cancelled)


class _Worker(QtCore.QObject):
  """Lives on a worker QThread. Runs session.run_turn and re-emits its
  AgentEvents as a single 'agent_event' signal, then 'finished'.

  Note: the signal is named ``agent_event`` (not ``event``) to avoid
  shadowing ``QObject.event(QEvent)`` — PySide2 emits a spurious
  ``TypeError: native Qt signal is not callable`` from C++ when a Signal
  shares the name of an inherited slot, even when emit() succeeds.
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
  """GUI-side wrapper for AgentSession.

  Signals (all GUI-thread):
    text_delta(str)
    thinking_delta(str)
    tool_use_requested(object)   # ToolUseRequested or ToolApprovalRequest
    tool_results_batched(object) # ToolResultsBatched
    image_emitted(object)
    usage(object)                # TokenUsage event
    turn_done(str)               # stop_reason
    error(str, bool, str)        # message, recoverable, kind

  Lifecycle: one runner per ChatWindow, reused across turns. Only one turn
  may run at a time; start_turn while busy is a no-op (returns False).

  Thread affinity: start_turn, cancel, submit_approval, and wait_for_idle
  must be called from the GUI thread; they touch _thread / _worker / _cancel
  which are GUI-thread-owned."""

  text_delta = QtCore.Signal(str)
  thinking_delta = QtCore.Signal(str)
  tool_use_requested = QtCore.Signal(object)
  tool_results_batched = QtCore.Signal(object)
  image_emitted = QtCore.Signal(object)
  usage = QtCore.Signal(object)
  turn_done = QtCore.Signal(str)
  error = QtCore.Signal(str, bool, str)

  def __init__(self, session, parent=None):
    super().__init__(parent)
    self.session = session
    self._thread = None
    self._worker = None
    self._cancel = CancelToken()
    # Approval queues parked on by the session worker — the GUI flushes
    # _Cancelled sentinels into these on cancel (Section 10.4).
    self._pending_approval_queues = [self.session.approval_queue]

  # ---- public API ----------------------------------------------------------

  def start_turn(self, user_message):
    """Kick off a new turn on the worker thread. Must be called from the
    GUI thread."""
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
    """Set the cancel token AND flush _Cancelled sentinels into every
    parked approval queue so a worker blocked on queue.get() wakes up
    (Section 10.4). Must be called from the GUI thread."""
    if not self.is_busy():
      # No worker to cancel; pushing a sentinel now would mis-cancel the
      # next turn (the queue is shared with the session).
      return
    self._cancel.set()
    for q in self._pending_approval_queues:
      # Approval queues are unbounded queue.Queue; put_nowait never blocks.
      q.put_nowait(_Cancelled())

  def submit_approval(self, response):
    """Push a ToolApprovalResponse into the session's approval queue.
    Must be called from the GUI thread."""
    self.session.approval_queue.put(response)
    if response.decision == "deny_and_stop":
      self._cancel.set()

  def is_busy(self):
    return self._thread is not None and self._thread.isRunning()

  def wait_for_idle(self, timeout_ms=5000):
    """Block until the current turn finishes. Used by tests; the GUI never
    calls this. Must be called from the GUI thread."""
    if self._thread is None:
      return
    self._thread.wait(timeout_ms)

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
    elif isinstance(ev, ImageEmitted):
      self.image_emitted.emit(ev)
    elif isinstance(ev, TokenUsageEvent):
      self.usage.emit(ev)
    elif isinstance(ev, TurnDone):
      self.turn_done.emit(ev.stop_reason)
    elif isinstance(ev, AgentError):
      self.error.emit(ev.message, ev.recoverable, ev.kind or "")

  def _on_worker_finished(self):
    if self._thread is not None:
      self._thread.quit()
      self._thread.wait()
      self._thread.deleteLater()
    if self._worker is not None:
      self._worker.deleteLater()
    self._thread = None
    self._worker = None
