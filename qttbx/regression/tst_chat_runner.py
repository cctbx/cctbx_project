"""QtAgentRunner unit tests. Uses a FakeAgent that yields scripted events,
exercises the QObject-side wiring (signals, start_turn, cancel) without
hitting any real network or Qt event-loop machinery beyond a brief
processEvents tick."""

import os
import shutil
import sys
import tempfile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from libtbx.utils import format_cpu_times, null_out

try:
  from qttbx.qt import QtCore, QtWidgets
except ImportError:
  print("PySide2/PySide6 not available; skipping")
  print("OK")
  sys.exit(0)

from qttbx.widgets.chat.agent.base import Agent, ToolSpec
from qttbx.widgets.font_init import init_default_app_font
from qttbx.widgets.chat.agent.conversation import (
  ContentBlock, Conversation, Message, now)
from qttbx.widgets.chat.agent.events import (
  TextDelta, TokenUsage as TokenUsageEvent, ToolUseRequested, TurnDone)
from qttbx.widgets.chat.agent.runner import QtAgentRunner
from qttbx.widgets.chat.agent.session import AgentSession
from qttbx.widgets.chat.agent.storage import ConversationStorage
from qttbx.widgets.chat.agent.tools import (
  ToolApprovalRequest, ToolPolicy, ToolRegistry)


class _ScriptedAgent(Agent):
  name = "scripted"
  model = "scripted-1"

  def __init__(self, events):
    self._events = events

  def stream_turn(self, conversation, tools, cancel):
    for ev in self._events:
      if cancel.is_set():
        return
      yield ev

  def resolve_credentials(self, cli_override=None):
    return "scripted-key"

  def credentials_dialog_class(self):
    return object


def _make_session(events, tmp):
  storage = ConversationStorage(project_dir=tmp, log=null_out())
  conv = Conversation.new(profile_name="t", model="scripted-1")
  registry = ToolRegistry(log=null_out())
  policy = ToolPolicy(default="ask")
  agent = _ScriptedAgent(events)
  return AgentSession(
    agent=agent, conversation=conv, storage=storage, tools=registry,
    policy=policy, profile=None, log=null_out()), conv


def _pump(app, ms=200):
  """Spin the Qt loop briefly so worker thread signals deliver."""
  deadline = QtCore.QElapsedTimer()
  deadline.start()
  while deadline.elapsed() < ms:
    app.processEvents(QtCore.QEventLoop.AllEvents, 25)


def _req(request_id):
  return ToolApprovalRequest(
    request_id=request_id, tool_name="t", tool_source="builtin",
    input={}, risk="write", batch_id=None, allow_remember=True)


def exercise_emits_text_delta_then_turn_done():
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    # A TokenUsage event in the stream must be tolerated even though the
    # runner no longer re-emits it as a Qt signal: token usage flows via
    # msg.usage (set on the message and persisted by the session), so the
    # runner just ignores the event without breaking dispatch.
    events = [TextDelta(text="hello "), TextDelta(text="world"),
              TokenUsageEvent(input=10, output=2),
              TurnDone(stop_reason="end_turn")]
    session, _ = _make_session(events, tmp)
    runner = QtAgentRunner(session)
    text_deltas, turn_dones = [], []
    runner.text_delta.connect(lambda s: text_deltas.append(s))
    runner.turn_done.connect(lambda r: turn_dones.append(r))

    user = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "hi"})], timestamp=now())
    runner.start_turn(user)
    runner.wait_for_idle(timeout_ms=2000)
    _pump(app)

    assert text_deltas == ["hello ", "world"], text_deltas
    assert turn_dones == ["end_turn"], turn_dones
  finally:
    shutil.rmtree(tmp)


def exercise_cancel_stops_the_turn():
  """A streamed delta drives a cancel that propagates onto the runner's token.

  The text delta is delivered to the GUI slot, which calls ``runner.cancel()``
  -- that must set the cancel token the session worker polls. The weak form
  only checked the runner went idle, which happens anyway once this short
  scripted turn ends; pinning the delta delivery + the set token catches a
  broken event dispatch or a ``cancel()`` that never reaches the token."""
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    events = [TextDelta(text="hello"), TurnDone(stop_reason="end_turn")]
    session, _ = _make_session(events, tmp)
    runner = QtAgentRunner(session)
    received = []
    runner.text_delta.connect(lambda s: (received.append(s), runner.cancel()))
    runner.start_turn(Message(role="user", content=[
      ContentBlock(type="text", data={"text": "x"})], timestamp=now()))
    runner.wait_for_idle(timeout_ms=2000)
    _pump(app)
    assert not runner.is_busy()
    # The streamed delta reached the GUI slot (so the cancel could fire) ...
    assert received == ["hello"], received
    # ... and runner.cancel() set the token the worker polls. Going idle alone
    # does not prove the cancel path ran.
    assert runner._cancel.is_set()
  finally:
    shutil.rmtree(tmp)


def exercise_error_event_emits_error_signal():
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    from qttbx.widgets.chat.agent.errors import AgentError
    events = [AgentError(message="boom", recoverable=False, kind="auth")]
    session, _ = _make_session(events, tmp)
    runner = QtAgentRunner(session)
    errs = []
    runner.error.connect(lambda msg, recoverable, kind:
                         errs.append((msg, recoverable, kind)))
    runner.start_turn(Message(role="user", content=[
      ContentBlock(type="text", data={"text": "x"})], timestamp=now()))
    runner.wait_for_idle(timeout_ms=2000)
    _pump(app)
    assert len(errs) == 1, errs
    assert errs[0][0] == "boom"
    assert errs[0][1] is False
    assert errs[0][2] == "auth"
  finally:
    shutil.rmtree(tmp)


def exercise_idle_fires_after_worker_teardown():
  """The runner announces end-of-turn teardown with ``idle``: emitted once per
  turn, AFTER the turn's terminal event (error / turn_done) has been dispatched
  and AFTER _on_worker_finished has torn down the worker/thread (is_busy()
  False) and drained the question queue -- so a GUI slot may synchronously
  start the next turn from the signal. shutdown() must NOT emit it (teardown
  paths manage their own state)."""
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    from qttbx.widgets.chat.agent.errors import AgentError
    events = [AgentError(message="boom", recoverable=False, kind="auth")]
    session, _ = _make_session(events, tmp)
    runner = QtAgentRunner(session)
    seq = []
    runner.error.connect(lambda *a: seq.append("error"))
    runner.idle.connect(lambda: seq.append(("idle", runner.is_busy())))
    runner.start_turn(Message(role="user", content=[
      ContentBlock(type="text", data={"text": "x"})], timestamp=now()))
    runner.wait_for_idle(timeout_ms=2000)
    _pump(app)
    assert seq == ["error", ("idle", False)], seq
    runner.shutdown()                        # already torn down: no re-emit
    _pump(app, 50)
    assert seq == ["error", ("idle", False)], seq
  finally:
    shutil.rmtree(tmp)


def exercise_cancel_when_idle_does_not_poison_next_turn():
  """Clicking Stop after a turn finishes must not cancel an approval a later
  turn has open on the coordinator -- cancel() no-ops while idle."""
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    events = [TextDelta(text="x"), TurnDone(stop_reason="end_turn")]
    session, _ = _make_session(events, tmp)
    runner = QtAgentRunner(session)
    runner.start_turn(Message(role="user", content=[
      ContentBlock(type="text", data={"text": "x"})], timestamp=now()))
    runner.wait_for_idle(timeout_ms=2000)
    _pump(app)
    assert not runner.is_busy()
    session.approvals.begin_turn()
    fut = session.approvals.open(_req("live"), lambda r: None)
    runner.cancel()                          # idle — must NOT cancel it
    assert not fut.done(), \
      "cancel() while idle must not cancel a pending approval"
  finally:
    shutil.rmtree(tmp)


def exercise_idle_click_records_no_remember_and_is_dropped():
  """A remember='tool' click on a stale card AFTER the turn ended must record
  NOTHING and reach neither the agent recorder nor the coordinator:
  runner.submit_approval gates record_approval_remember / agent.submit_approval
  / coordinator.submit behind is_busy(), so an idle click can't grant a standing
  session auto-approve for a tool that never ran."""
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    from qttbx.widgets.chat.agent.tools import ToolApprovalResponse
    session, _ = _make_session([], tmp)          # no turn running -> idle
    runner = QtAgentRunner(session)
    assert not runner.is_busy()
    recorded, agent_calls, submitted = [], [], []
    session.record_approval_remember = lambda r: recorded.append(r)
    session.agent.submit_approval = lambda r: agent_calls.append(r)
    session.approvals.submit = lambda r: submitted.append(r) or False
    runner.submit_approval(ToolApprovalResponse(
      request_id="stray", decision="approve", remember="tool"))
    assert recorded == [], "idle click must not record a remember"
    assert agent_calls == [], "idle click must not reach the agent recorder"
    assert submitted == [], "idle click must not reach the coordinator"
  finally:
    shutil.rmtree(tmp)


def exercise_submit_approval_resolves_parked_worker():
  """The GUI approval path end to end: a turn parks in _await_approval (the
  coordinator surfaces a ToolApprovalRequest via the runner's
  tool_use_requested signal); submit_approval resolves it and the turn ends."""
  from qttbx.widgets.chat.agent.tools import (
    ToolApprovalRequest, ToolApprovalResponse)
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    class _AskThenDone(Agent):
      name = "askdone"
      model = "askdone-1"
      def __init__(self):
        self._n = 0
      def stream_turn(self, conversation, tools, cancel):
        self._n += 1
        if self._n == 1:
          yield ToolUseRequested(id="t1", name="echo", input={})
          yield TurnDone(stop_reason="tool_use")
        else:
          yield TurnDone(stop_reason="end_turn")
      def resolve_credentials(self, cli_override=None):
        return "k"
      def credentials_dialog_class(self):
        return object

    storage = ConversationStorage(project_dir=tmp, log=null_out())
    conv = Conversation.new(profile_name="t", model="askdone-1")
    session = AgentSession(
      agent=_AskThenDone(), conversation=conv, storage=storage,
      tools=ToolRegistry(log=null_out()),
      policy=ToolPolicy(default="ask"), profile=None, log=null_out())
    session.tools.register_builtin(
      ToolSpec(name="echo", description="echo",
               input_schema={"type": "object"}),
      handler=lambda **kw: "x", risk="write")
    runner = QtAgentRunner(session)
    reqs = []
    runner.tool_use_requested.connect(
      lambda ev: reqs.append(ev) if isinstance(ev, ToolApprovalRequest)
      else None)
    runner.start_turn(Message(
      role="user", timestamp=now(),
      content=[ContentBlock(type="text", data={"text": "hi"})]))
    deadline = QtCore.QElapsedTimer()
    deadline.start()
    while not reqs and deadline.elapsed() < 3000:
      _pump(app, 20)
    assert reqs, "no ToolApprovalRequest surfaced"
    runner.submit_approval(ToolApprovalResponse(
      request_id=reqs[0].request_id, decision="approve"))
    runner.wait_for_idle(timeout_ms=3000)
    _pump(app, 50)
    assert not runner.is_busy()
  finally:
    shutil.rmtree(tmp)


def exercise_submit_question_answer_falls_through_to_session():
  """When the agent inherits the base submit_question_answer (the API
  backends return the default False), runner.submit_question_answer routes
  the answers to the session's submit_question_answer. The agent path wins
  when it owns the id; the session path is the fallback the API agents rely
  on."""
  import queue as _queue

  class _AgentWithoutQuestions(Agent):
    """An API-style agent that does NOT override submit_question_answer; it
    inherits the Agent base default (returns False), so the runner falls
    through to the session."""
    name = "no-questions"
    model = "no-questions-1"

    def stream_turn(self, conversation, tools, cancel):
      return iter(())

    def resolve_credentials(self, cli_override=None):
      return "k"

    def credentials_dialog_class(self):
      return object

  class _AgentOwning:
    """A claude_code-style agent that owns the id."""
    def submit_question_answer(self, request_id, answers):
      return request_id == "agent_owned"

  class _FakeSession:
    def __init__(self, agent):
      self.agent = agent
      self.approval_queue = _queue.Queue()
      self.question_queue = _queue.Queue()
      self.routed = []
    def submit_question_answer(self, request_id, answers):
      self.routed.append((request_id, answers))
      return True

  # Agent lacks the method -> falls through to the session.
  sess = _FakeSession(_AgentWithoutQuestions())
  runner = QtAgentRunner(sess)
  assert runner.submit_question_answer("q1", {"Q": "A"}) is True
  assert sess.routed == [("q1", {"Q": "A"})], sess.routed

  # Agent owns the id -> handled by the agent, session not consulted.
  sess2 = _FakeSession(_AgentOwning())
  runner2 = QtAgentRunner(sess2)
  assert runner2.submit_question_answer("agent_owned", {"Q": "A"}) is True
  assert sess2.routed == [], sess2.routed

  # Agent has the method but doesn't own the id -> falls through to session.
  sess3 = _FakeSession(_AgentOwning())
  runner3 = QtAgentRunner(sess3)
  assert runner3.submit_question_answer("session_owned", {"Q": "B"}) is True
  assert sess3.routed == [("session_owned", {"Q": "B"})], sess3.routed


def exercise_cancel_flushes_question_queue():
  """cancel() must push a _Cancelled into the session's question_queue too,
  so a worker parked on a pending phenix_ask_user_question wakes up."""
  import threading
  from qttbx.widgets.chat.agent.tools import _Cancelled
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)

  class _BlockingAgent(Agent):
    name = "blocking-q"
    model = "blocking-q-1"

    def __init__(self):
      self.release = threading.Event()

    def stream_turn(self, conversation, tools, cancel):
      self.release.wait(timeout=5)
      yield TurnDone(stop_reason="end_turn")

    def resolve_credentials(self, cli_override=None):
      return "k"

    def credentials_dialog_class(self):
      return object

  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=tmp, log=null_out())
    conv = Conversation.new(profile_name="t", model="blocking-q-1")
    agent = _BlockingAgent()
    session = AgentSession(
      agent=agent, conversation=conv, storage=storage,
      tools=ToolRegistry(log=null_out()),
      policy=ToolPolicy(default="ask"), profile=None, log=null_out())
    runner = QtAgentRunner(session)
    runner.start_turn(Message(
      role="user", timestamp=now(),
      content=[ContentBlock(type="text", data={"text": "hi"})]))
    deadline = QtCore.QElapsedTimer()
    deadline.start()
    while not runner.is_busy() and deadline.elapsed() < 2000:
      _pump(app, 20)
    assert runner.is_busy()
    runner.cancel()
    # The cancel must have flushed a sentinel into the question queue.
    item = session.question_queue.get_nowait()
    assert isinstance(item, _Cancelled)
    agent.release.set()
    runner.wait_for_idle(timeout_ms=5000)
    _pump(app, 50)
  finally:
    shutil.rmtree(tmp)


def exercise_stray_approval_after_turn_end_is_dropped():
  """A response routed through the runner after the worker finished must be
  dropped by the is_busy gate -- it must never reach the coordinator (spied
  below), where it could resolve the next turn's first approval."""
  from qttbx.widgets.chat.agent.tools import ToolApprovalResponse
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    events = [TextDelta(text="x"), TurnDone(stop_reason="end_turn")]
    session, _ = _make_session(events, tmp)
    runner = QtAgentRunner(session)
    runner.start_turn(Message(role="user", timestamp=now(),
      content=[ContentBlock(type="text", data={"text": "x"})]))
    runner.wait_for_idle(timeout_ms=2000)
    _pump(app)
    assert not runner.is_busy()
    # Route the stray click through the runner so the is_busy drop-gate is
    # exercised: an idle-time submit must NOT reach the coordinator (where it
    # could resolve the next turn's first approval). Spying on submit catches a
    # regression that deletes `and self.is_busy()` from runner.submit_approval.
    calls = []
    real_submit = session.approvals.submit

    def _spy(r):
      calls.append(r)
      return real_submit(r)

    session.approvals.submit = _spy
    runner.submit_approval(ToolApprovalResponse(
      request_id="stray", decision="approve"))
    assert calls == [], "stray click after turn end must not reach the coordinator"
  finally:
    shutil.rmtree(tmp)


def exercise_stale_deny_and_stop_does_not_cancel_unrelated_turn():
  """runner.submit_approval sets the turn's cancel token on a deny_and_stop ONLY
  when the response belongs to the in-flight turn (agent-owned, or the approval
  the session worker is parked on). A stale card's stop from an abandoned turn
  must not cancel an unrelated later turn (root fix for the unguarded cancel)."""
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    from qttbx.widgets.chat.agent.tools import ToolApprovalResponse
    session, _ = _make_session([], tmp)
    runner = QtAgentRunner(session)
    # Simulate the worker parked on THIS turn's approval "r_current".
    # was: session._pending_approval_id = "r_current"
    session.approvals.begin_turn()
    session.approvals.open(_req("r_current"), lambda r: None)
    assert not runner._cancel.is_set()
    # A stale deny_and_stop for a DIFFERENT (abandoned) request must NOT cancel.
    runner.submit_approval(ToolApprovalResponse(
      request_id="r_stale", decision="deny_and_stop"))
    assert not runner._cancel.is_set(), "stale stop must not cancel the turn"
    # The in-flight request's own deny_and_stop DOES cancel.
    runner.submit_approval(ToolApprovalResponse(
      request_id="r_current", decision="deny_and_stop"))
    assert runner._cancel.is_set(), "the in-flight request's stop must cancel"
  finally:
    shutil.rmtree(tmp)


def exercise_shutdown_joins_worker_without_event_loop():
  """closeEvent / conversation-switch tears down a still-running turn from the
  GUI thread. It must fully stop AND join the worker QThread WITHOUT relying on
  the queued _on_worker_finished slot (which can't run while the GUI thread is
  blocked in the join) -- otherwise the QThread outlives the window and aborts
  with 'QThread: Destroyed while thread is still running'. shutdown() quits the
  thread's exec loop itself then joins it, so cancel()+shutdown() with NO
  event-loop pump (mimicking the blocked GUI thread) leaves the runner idle and
  the thread reference dropped."""
  import threading           # noqa: F401  (parity with the other blocking tests)
  import time
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    class _CancellableBlockingAgent(Agent):
      name = "cancellable"
      model = "cancellable-1"

      def stream_turn(self, conversation, tools, cancel):
        # Stay in flight (worker parked in next()) until the user cancels --
        # the closeEvent "Stop and exit" path.
        while not cancel.is_set():
          time.sleep(0.01)
        return
        yield  # pragma: no cover -- makes this a generator

      def resolve_credentials(self, cli_override=None):
        return "k"

      def credentials_dialog_class(self):
        return object

    storage = ConversationStorage(project_dir=tmp, log=null_out())
    conv = Conversation.new(profile_name="t", model="cancellable-1")
    session = AgentSession(
      agent=_CancellableBlockingAgent(), conversation=conv, storage=storage,
      tools=ToolRegistry(log=null_out()),
      policy=ToolPolicy(default="ask"), profile=None, log=null_out())
    runner = QtAgentRunner(session)
    runner.start_turn(Message(role="user", timestamp=now(),
      content=[ContentBlock(type="text", data={"text": "hi"})]))
    deadline = QtCore.QElapsedTimer()
    deadline.start()
    while not runner.is_busy() and deadline.elapsed() < 2000:
      _pump(app, 20)
    assert runner.is_busy(), "runner should be busy during the blocked turn"
    # closeEvent teardown: cancel + shutdown, with NO event-loop pump after --
    # in the real closeEvent the GUI thread is blocked here, so the queued
    # _on_worker_finished slot cannot run.
    runner.cancel()
    t0 = QtCore.QElapsedTimer()
    t0.start()
    runner.shutdown(timeout_ms=3000)
    assert not runner.is_busy(), "shutdown must stop + join the worker thread"
    assert runner._thread is None, runner._thread
    # The thread really finished -- shutdown did not just spin the full timeout.
    assert t0.elapsed() < 2500, \
      "shutdown blocked the full timeout (thread never joined): %dms" \
      % t0.elapsed()
  finally:
    shutil.rmtree(tmp)


def exercise():
  exercise_emits_text_delta_then_turn_done()
  exercise_cancel_stops_the_turn()
  exercise_error_event_emits_error_signal()
  exercise_idle_fires_after_worker_teardown()
  exercise_cancel_when_idle_does_not_poison_next_turn()
  exercise_idle_click_records_no_remember_and_is_dropped()
  exercise_submit_approval_resolves_parked_worker()
  exercise_submit_question_answer_falls_through_to_session()
  exercise_cancel_flushes_question_queue()
  exercise_stray_approval_after_turn_end_is_dropped()
  exercise_stale_deny_and_stop_does_not_cancel_unrelated_turn()
  exercise_shutdown_joins_worker_without_event_loop()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
