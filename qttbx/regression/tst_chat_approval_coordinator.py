import concurrent.futures
import threading

from libtbx.utils import format_cpu_times
from qttbx.widgets.chat.agent.approval import ApprovalCoordinator
from qttbx.widgets.chat.agent.tools import (
  ToolApprovalRequest, ToolApprovalResponse)


def _req(request_id):
  return ToolApprovalRequest(
    request_id=request_id, tool_name="t", tool_source="builtin",
    input={}, risk="write", batch_id=None, allow_remember=True)


def _resp(request_id, decision="approve", remember="none"):
  return ToolApprovalResponse(
    request_id=request_id, decision=decision, remember=remember)


def exercise_open_then_submit_resolves_future():
  """A submitted response resolves the future the opener is waiting on."""
  coord = ApprovalCoordinator()
  coord.begin_turn()
  emitted = []
  fut = coord.open(_req("r1"), emitted.append)
  assert [e.request_id for e in emitted] == ["r1"]
  assert not fut.done()
  assert coord.submit(_resp("r1")) is True
  assert fut.result(timeout=1).decision == "approve"


def exercise_submit_unknown_id_returns_false():
  """A response for an id we never opened is dropped, not applied."""
  coord = ApprovalCoordinator()
  coord.begin_turn()
  assert coord.submit(_resp("nope")) is False


def exercise_begin_turn_cancels_prior_pending():
  """A new turn abandons an approval left pending by the previous turn:
  the old waiter wakes with CancelledError, and the new turn's opens work.
  This is the log's cancel-then-resend scenario at the unit level."""
  coord = ApprovalCoordinator()
  coord.begin_turn()
  stale = coord.open(_req("old"), lambda r: None)
  coord.begin_turn()                                  # resend / next turn
  try:
    stale.result(timeout=1)
    raise AssertionError("stale future should have been cancelled")
  except concurrent.futures.CancelledError:
    pass
  fresh = coord.open(_req("new"), lambda r: None)
  assert coord.submit(_resp("new")) is True
  assert fresh.result(timeout=1).decision == "approve"


def exercise_submit_after_new_turn_is_dropped():
  """A late click on an abandoned earlier turn's card is dropped (its request
  left the registry at the turn boundary), so it can never be misapplied to
  the current turn."""
  coord = ApprovalCoordinator()
  coord.begin_turn()
  coord.open(_req("old"), lambda r: None)
  coord.begin_turn()
  assert coord.submit(_resp("old")) is False


def exercise_cancel_turn_cancels_all_pending():
  """Stop wakes every approval the turn is waiting on (parallel tool calls)."""
  coord = ApprovalCoordinator()
  coord.begin_turn()
  a = coord.open(_req("a"), lambda r: None)
  b = coord.open(_req("b"), lambda r: None)
  coord.cancel_turn()
  for fut in (a, b):
    try:
      fut.result(timeout=1)
      raise AssertionError("cancel_turn should have cancelled the future")
    except concurrent.futures.CancelledError:
      pass


def exercise_concurrent_opens_resolve_independently():
  """Two approvals pending at once (the SDK can request parallel tool calls)
  each resolve to their own response -- the case the session's scalar
  _pending_approval_id could not represent."""
  coord = ApprovalCoordinator()
  coord.begin_turn()
  fa = coord.open(_req("a"), lambda r: None)
  fb = coord.open(_req("b"), lambda r: None)
  assert coord.submit(_resp("b", decision="deny")) is True
  assert coord.submit(_resp("a", decision="approve")) is True
  assert fa.result(timeout=1).decision == "approve"
  assert fb.result(timeout=1).decision == "deny"


def exercise_owns_tracks_pending():
  coord = ApprovalCoordinator()
  coord.begin_turn()
  coord.open(_req("a"), lambda r: None)
  assert coord.owns("a") is True
  assert coord.owns("b") is False
  coord.begin_turn()
  assert coord.owns("a") is False


def exercise_blocking_waiter_wakes_on_submit_from_another_thread():
  """result() on the opener's thread wakes when submit() lands on another --
  the AgentSession worker-thread blocking model."""
  coord = ApprovalCoordinator()
  coord.begin_turn()
  fut = coord.open(_req("a"), lambda r: None)
  box = []

  def _wait():
    box.append(fut.result(timeout=2).decision)

  t = threading.Thread(target=_wait)
  t.start()
  assert coord.submit(_resp("a")) is True
  t.join(timeout=2)
  assert box == ["approve"]


def exercise_open_on_closed_coordinator_precancels_without_emitting():
  """After cancel_turn (Stop) the coordinator is CLOSED: open() must NOT surface
  a card and must return a pre-cancelled future -- so a Stop racing a
  tool-approval open() leaves no phantom card and the caller unwinds at once.
  The next begin_turn reopens it."""
  coord = ApprovalCoordinator()
  coord.begin_turn()
  coord.cancel_turn()                          # Stop closed the coordinator
  emitted = []
  fut = coord.open(_req("r1"), emitted.append)
  assert fut.cancelled(), "closed open() must return a pre-cancelled future"
  assert emitted == [], "closed open() must not emit a card"
  assert not coord.owns("r1"), "closed open() must not register"
  coord.begin_turn()                           # a fresh turn reopens it
  fut2 = coord.open(_req("r2"), emitted.append)
  assert not fut2.cancelled()
  assert len(emitted) == 1, emitted            # now emits normally


def exercise_open_emits_under_the_lock():
  """open() must surface the card while STILL holding the coordinator lock, so a
  cancel_turn() (Stop) on another thread cannot interleave between registering
  the future and emitting. Emitting after the lock is released leaves a window
  where cancel_turn() cancels the future and clears the registry, yet the card
  still surfaces -- a phantom that renders after the turn's cards were finalized
  and is never disabled. Probe it directly: a non-blocking acquire from inside
  emit fails iff open() still holds the lock. (The closed-flag test covers cancel
  BEFORE open; this covers cancel racing DURING open.)"""
  coord = ApprovalCoordinator()
  coord.begin_turn()
  seen = {}

  def probe_emit(req):
    got = coord._lock.acquire(blocking=False)
    seen["lock_free_during_emit"] = got
    if got:
      coord._lock.release()

  coord.open(_req("r1"), probe_emit)
  assert seen["lock_free_during_emit"] is False, \
    "open() emitted outside the lock -- a cancel_turn() can interleave between " \
    "register and emit and strand a phantom card"


def exercise():
  exercise_open_then_submit_resolves_future()
  exercise_open_on_closed_coordinator_precancels_without_emitting()
  exercise_open_emits_under_the_lock()
  exercise_submit_unknown_id_returns_false()
  exercise_begin_turn_cancels_prior_pending()
  exercise_submit_after_new_turn_is_dropped()
  exercise_cancel_turn_cancels_all_pending()
  exercise_concurrent_opens_resolve_independently()
  exercise_owns_tracks_pending()
  exercise_blocking_waiter_wakes_on_submit_from_another_thread()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
