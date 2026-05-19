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


def exercise():
  exercise_emits_text_delta_then_turn_done()
  exercise_cancel_stops_the_turn()
  exercise_error_event_emits_error_signal()
  exercise_cancel_when_idle_does_not_poison_next_turn()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
