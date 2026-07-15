"""Shared, cancel-safe tool-approval bookkeeping.

One ``ApprovalCoordinator`` is the single owner of the pending-approval
registry: it drops every pending request at each turn boundary
(``begin_turn`` / ``cancel_turn``), so a stale response from an abandoned turn
is no longer in the registry and ``submit`` rejects it, and it wakes waiters on
cancel. Both
``AgentSession`` (the API backends, which block a worker thread on
``future.result()``) and ``ClaudeCodeAgent`` (the SDK backend, whose async
``_on_can_use_tool`` does ``await asyncio.wrap_future(future)``) delegate to it,
so neither carries its own copy of the invariant that a cancelled turn must not
wedge the next turn's approvals. Callers supply only how to *emit* a request to
the UI; the coordinator owns everything cancel-sensitive.

Pure sync primitive: no Qt, no I/O. ``concurrent.futures.Future`` serves both a
blocking ``.result()`` (worker thread) and ``asyncio.wrap_future`` (event loop).
"""

import concurrent.futures
import threading


class ApprovalCoordinator:
  """Cancel-safe registry of in-flight tool approvals, keyed by request id.

  ``begin_turn`` / ``cancel_turn`` abandon (cancel + drop) every pending
  request, so a request from an abandoned turn is no longer in the registry
  and a late response for it is dropped by ``submit``. ``cancel_turn`` also
  CLOSES the coordinator (level-triggered; the next ``begin_turn`` reopens it):
  an ``open`` that races a Stop is judged against the closed flag under the SAME
  lock that ``cancel_turn`` sets it, so it registers nothing and surfaces no
  card. Cancellation state and register+emit no longer live under different
  locks -- that retires the cancel-vs-open race class (and no epoch counter is
  needed).
  """

  def __init__(self):
    self._pending = {}                 # request_id -> Future
    self._closed = False               # True between cancel_turn and begin_turn
    self._lock = threading.Lock()

  def begin_turn(self):
    """Start a turn: REOPEN the coordinator and abandon anything a previous
    (cancelled) turn left pending, waking its waiter with CancelledError.
    Usually a no-op on the registry (a completed/cancelled turn already emptied
    it); always clears the closed flag."""
    self._abandon_pending(closed=False)

  def open(self, req, emit):
    """Register ``req``, surface it via ``emit``, and return the future the
    caller waits on in its own idiom.

    Parameters
    ----------
    req : ToolApprovalRequest
        The request to surface; ``req.request_id`` keys the registry.
    emit : callable
        ``emit(req)`` pushes the request toward the UI (the session passes its
        ``on_event``; claude_code passes its per-turn event queue's ``put``).

    Returns
    -------
    concurrent.futures.Future
        Resolves to a ``ToolApprovalResponse`` on ``submit``; is cancelled on
        ``cancel_turn`` / the next ``begin_turn``. If the coordinator is closed
        (a Stop already cancelled this turn), returns a pre-cancelled future and
        does NOT emit -- no card is surfaced.
    """
    with self._lock:
      fut = concurrent.futures.Future()
      if self._closed:
        # A Stop closed the coordinator this turn (cancel_turn ran). Register
        # nothing and do NOT emit: surfacing a card here is the finding-1/2
        # phantom -- it would render after the turn's cards were finalized and
        # never be disabled. The pre-cancelled future makes the caller (worker
        # fut.result() / SDK await) unwind into TurnCancelled at once.
        fut.cancel()
        return fut
      self._pending[req.request_id] = fut
      # Emit UNDER the lock: register + emit is one critical section, so a
      # cancel_turn() (Stop) can't interleave between them and strand a card it
      # already finalized past (the last cancel-vs-open crack). Safe to hold the
      # lock here -- both emit targets are non-blocking (on_event queues a Qt
      # signal; the SDK path does q.put), never UI work that could block.
      emit(req)
    return fut

  def submit(self, response):
    """Resolve the future for ``response.request_id``.

    Returns ``True`` if it matched a pending request; ``False`` for an unknown
    id or a stale click from an abandoned earlier turn (already dropped from
    the registry at the turn boundary), which is discarded, never misapplied.
    """
    with self._lock:
      fut = self._pending.pop(response.request_id, None)
    if fut is None:
      return False
    if not fut.done():
      try:
        fut.set_result(response)
      except concurrent.futures.InvalidStateError:
        # Concurrently cancelled/resolved between the done() check and here
        # (e.g. SDK-loop teardown racing a GUI-thread submit). Nothing to do.
        pass
    return True

  def cancel_turn(self):
    """Cancel every approval the current turn is waiting on AND close the
    coordinator (the runner's Stop handler). Waiters wake with CancelledError;
    a subsequent ``open`` this turn returns a pre-cancelled future without
    surfacing a card, so a Stop racing a tool-approval ``open`` leaves no
    phantom card. Reopened by the next ``begin_turn``."""
    self._abandon_pending(closed=True)

  def owns(self, request_id):
    """Whether ``request_id`` is pending (lets the runner gate a
    deny_and_stop cancel on a matching request)."""
    with self._lock:
      return request_id in self._pending

  def _abandon_pending(self, closed):
    with self._lock:
      self._closed = closed
      pending = list(self._pending.values())
      self._pending.clear()
    for fut in pending:
      if not fut.done():
        fut.cancel()
