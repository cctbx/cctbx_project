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
  """When the agent owns a request_id (e.g. Claude Code's
  can_use_tool callback emitted it), runner.submit_approval must
  forward to agent.submit_approval and NOT push to the session
  queue. The session might be parked on its own unrelated approval,
  and grabbing the SDK-side response would mis-route."""
  app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
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
    # request_id the agent doesn't own -> falls through to session.
    runner.submit_approval(ToolApprovalResponse(
      request_id="session_owned", decision="approve"))
    assert not session.approval_queue.empty()
    pulled = session.approval_queue.get_nowait()
    assert pulled.request_id == "session_owned", pulled
  finally:
    shutil.rmtree(tmp)


def exercise():
  exercise_emits_text_delta_then_turn_done()
  exercise_cancel_stops_the_turn()
  exercise_error_event_emits_error_signal()
  exercise_cancel_when_idle_does_not_poison_next_turn()
  exercise_submit_approval_routes_to_agent_when_agent_owns_request_id()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
