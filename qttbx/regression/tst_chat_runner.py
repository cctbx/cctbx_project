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

from qttbx.widgets.chat.agent.base import Agent, AgentCapabilities
from qttbx.widgets.font_init import init_default_app_font
from qttbx.widgets.chat.agent.conversation import (
  ContentBlock, Conversation, Message, now)
from qttbx.widgets.chat.agent.events import (
  TextDelta, TokenUsage as TokenUsageEvent, TurnDone)
from qttbx.widgets.chat.agent.runner import QtAgentRunner
from qttbx.widgets.chat.agent.session import AgentSession
from qttbx.widgets.chat.agent.storage import ConversationStorage
from qttbx.widgets.chat.agent.tools import ToolPolicy, ToolRegistry


class _ScriptedAgent(Agent):
  name = "scripted"
  model = "scripted-1"
  capabilities = AgentCapabilities.STREAMING

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


def exercise_emits_text_delta_then_turn_done():
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    events = [TextDelta(text="hello "), TextDelta(text="world"),
              TokenUsageEvent(input=10, output=2),
              TurnDone(stop_reason="end_turn")]
    session, _ = _make_session(events, tmp)
    runner = QtAgentRunner(session)
    text_deltas, turn_dones, usages = [], [], []
    runner.text_delta.connect(lambda s: text_deltas.append(s))
    runner.turn_done.connect(lambda r: turn_dones.append(r))
    runner.usage.connect(lambda u: usages.append(u))

    user = Message(role="user", content=[
      ContentBlock(type="text", data={"text": "hi"})], timestamp=now())
    runner.start_turn(user)
    runner.wait_for_idle(timeout_ms=2000)
    _pump(app)

    assert text_deltas == ["hello ", "world"], text_deltas
    assert turn_dones == ["end_turn"], turn_dones
    assert len(usages) == 1
    assert usages[0].input == 10
  finally:
    shutil.rmtree(tmp)


def exercise_cancel_stops_the_turn():
  """Cancel mid-stream — the scripted agent honors cancel between events
  so the turn ends with stop_reason='cancelled'."""
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


def exercise_cancel_when_idle_does_not_poison_next_turn():
  """Clicking Stop after a turn finishes must not push a sentinel that the
  next turn's first approval.get() would consume."""
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
    runner.cancel()                          # idle — must NOT enqueue sentinel
    assert session.approval_queue.empty(), \
      "cancel() while idle should not push to approval queue"
  finally:
    shutil.rmtree(tmp)


def exercise_submit_approval_routes_to_agent_when_agent_owns_request_id():
  """When the agent owns a request_id (e.g. Claude Code's can_use_tool
  callback emitted it), runner.submit_approval forwards to
  agent.submit_approval and does NOT push to the session queue (the
  session might be parked on its own unrelated approval). A non-owned
  response submitted while no turn is in flight is dropped, so a late
  click after the turn ended can't leak into the next turn's first
  approval."""
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    from qttbx.widgets.chat.agent.tools import ToolApprovalResponse

    class _OwningAgent:
      """Pretends to own request_id 'agent_owned'."""
      received = []

      def submit_approval(self, response):
        _OwningAgent.received.append(response)
        return response.request_id == "agent_owned"

    events = []   # no streamed events
    session, _ = _make_session(events, tmp)
    session.agent = _OwningAgent()
    runner = QtAgentRunner(session)
    # request_id the agent owns -> session queue stays empty.
    runner.submit_approval(ToolApprovalResponse(
      request_id="agent_owned", decision="approve"))
    assert session.approval_queue.empty(), \
      "agent-owned response should NOT land in session queue"
    assert _OwningAgent.received[-1].request_id == "agent_owned"
    # request_id the agent doesn't own AND no turn in flight -> dropped
    # (queueing it would leak into the next turn's first approval).
    runner.submit_approval(ToolApprovalResponse(
      request_id="session_owned", decision="approve"))
    assert session.approval_queue.empty(), \
      "non-owned response submitted while idle must be dropped"
  finally:
    shutil.rmtree(tmp)


def exercise_submit_approval_queues_for_session_while_turn_in_flight():
  """While a turn IS running, a non-owned response is queued for the
  session worker -- the normal GUI approval path (the worker is parked on
  approval_queue.get() waiting for it)."""
  import threading
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    from qttbx.widgets.chat.agent.tools import ToolApprovalResponse

    class _BlockingAgent(Agent):
      name = "blocking"
      model = "blocking-1"
      capabilities = AgentCapabilities.STREAMING

      def __init__(self):
        self.release = threading.Event()

      def stream_turn(self, conversation, tools, cancel):
        # Hold the turn open (runner stays busy) until released.
        self.release.wait(timeout=5)
        yield TurnDone(stop_reason="end_turn")

      def resolve_credentials(self, cli_override=None):
        return "blocking-key"

      def credentials_dialog_class(self):
        return object

    storage = ConversationStorage(project_dir=tmp, log=null_out())
    conv = Conversation.new(profile_name="t", model="blocking-1")
    agent = _BlockingAgent()
    session = AgentSession(
      agent=agent, conversation=conv, storage=storage,
      tools=ToolRegistry(log=null_out()),
      policy=ToolPolicy(default="ask"), profile=None, log=null_out())
    runner = QtAgentRunner(session)
    runner.start_turn(Message(
      role="user", timestamp=now(),
      content=[ContentBlock(type="text", data={"text": "hi"})]))
    # Wait until the worker thread is actually running.
    deadline = QtCore.QElapsedTimer()
    deadline.start()
    while not runner.is_busy() and deadline.elapsed() < 2000:
      _pump(app, 20)
    assert runner.is_busy(), "runner should be busy during the blocked turn"
    runner.submit_approval(ToolApprovalResponse(
      request_id="session_owned", decision="approve"))
    assert not session.approval_queue.empty(), \
      "non-owned response while a turn is in flight should be queued"
    assert session.approval_queue.get_nowait().request_id == "session_owned"
    agent.release.set()
    runner.wait_for_idle(timeout_ms=5000)
    _pump(app, 50)
  finally:
    shutil.rmtree(tmp)


def exercise_submit_question_answer_falls_through_to_session():
  """When the agent has no submit_question_answer (the API backends),
  runner.submit_question_answer routes the answers to the session's
  submit_question_answer. The agent path wins when it owns the id; the
  session path is the fallback the API agents rely on."""
  import queue as _queue

  class _AgentWithoutQuestions:
    """An API-style agent: no submit_question_answer override at all (it
    is not even an Agent subclass, so getattr returns None)."""

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
    capabilities = AgentCapabilities.STREAMING

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


def exercise_worker_finish_drains_stray_approval_response():
  """A response left in the approval queue when the turn ends (e.g. a late
  click during teardown, submitted while is_busy() was still true) is
  drained, so it can't leak into the next turn's first approval."""
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
  init_default_app_font(app)
  tmp = tempfile.mkdtemp()
  try:
    from qttbx.widgets.chat.agent.tools import ToolApprovalResponse
    session, _ = _make_session([TurnDone(stop_reason="end_turn")], tmp)
    runner = QtAgentRunner(session)
    # Stray response with no consumer: run_turn's start drain keeps real
    # responses, so only the turn-end drain in _on_worker_finished removes
    # it.
    session.approval_queue.put(ToolApprovalResponse(
      request_id="stray", decision="approve"))
    runner.start_turn(Message(
      role="user", timestamp=now(),
      content=[ContentBlock(type="text", data={"text": "hi"})]))
    deadline = QtCore.QElapsedTimer()
    deadline.start()
    while runner.is_busy() and deadline.elapsed() < 5000:
      _pump(app, 20)
    assert not runner.is_busy(), "turn should have finished"
    assert session.approval_queue.empty(), \
      "a stray approval response must be drained when the turn finishes"
  finally:
    shutil.rmtree(tmp)


def exercise():
  exercise_emits_text_delta_then_turn_done()
  exercise_cancel_stops_the_turn()
  exercise_error_event_emits_error_signal()
  exercise_cancel_when_idle_does_not_poison_next_turn()
  exercise_submit_approval_routes_to_agent_when_agent_owns_request_id()
  exercise_submit_approval_queues_for_session_while_turn_in_flight()
  exercise_submit_question_answer_falls_through_to_session()
  exercise_cancel_flushes_question_queue()
  exercise_worker_finish_drains_stray_approval_response()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
