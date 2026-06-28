"""McpServerConnection tests against an in-process fastmcp server
fixture: tool listing + dispatch round-trip, error paths, empty-command
clarity, and risk-annotation derivation. Subprocess coverage of the
same client lives in tst_chat_mcp_e2e.py (phenix repo)."""

import shutil
import sys
import tempfile
from pathlib import Path

from libtbx.utils import format_cpu_times, null_out, Sorry

try:
  import fastmcp                                          # noqa: F401
except ImportError:
  print("fastmcp not available; skipping")
  print("OK")
  sys.exit(0)


from qttbx.widgets.chat.agent.errors import CancelToken
from qttbx.widgets.chat.agent.mcp_client import (
  McpToolResult, McpServerConnection)
from qttbx.widgets.chat.agent.storage import ConversationStorage


def _build_app():
  from fastmcp import FastMCP
  app = FastMCP("test")

  @app.tool()
  def echo(message: str) -> str:
    """Echo back."""
    return "echo:" + message

  @app.tool()
  def boom() -> str:
    raise RuntimeError("intentional")

  return app


def _build_risk_app():
  from fastmcp import FastMCP
  app = FastMCP("risk-test")

  @app.tool(annotations={"readOnlyHint": True})
  def look() -> str:
    """Read-only."""
    return "ok"

  @app.tool(annotations={"destructiveHint": True})
  def nuke() -> str:
    """Destructive."""
    return "ok"

  @app.tool()
  def write() -> str:
    """Default risk."""
    return "ok"

  return app


class _FakeServerConfig:
  def __init__(self, name="t"):
    self.name = name
    self.command = None
    self.args = []
    self.env = {}
    self.tool_policy = {}


def _conn(tmp, in_process_app):
  """Build a connection whose transport talks to the in-process app."""
  storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
  conn = McpServerConnection(
    config=_FakeServerConfig(name="t"),
    project_dir=Path(tmp), storage=storage, conv_id="conv-1",
    log=null_out())
  conn._inject_in_process_app(in_process_app)
  return conn


def _bare_conn(tmp):
  """A connection with no in-process app, for exercising the event-loop
  helpers (_start_loop / _run_async / _stop_loop) in isolation -- without
  spawning a server or calling start()."""
  storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
  return McpServerConnection(
    config=_FakeServerConfig(name="t"),
    project_dir=Path(tmp), storage=storage, conv_id="conv-1",
    log=null_out())


def exercise_start_lists_tools():
  tmp = tempfile.mkdtemp()
  try:
    conn = _conn(tmp, _build_app())
    conn.start()
    try:
      assert conn.state == McpServerConnection.STATE_READY
      names = {t.name for t in conn.tools}
      assert "echo" in names
      assert "boom" in names
    finally:
      conn.stop()
  finally:
    shutil.rmtree(tmp)


def exercise_call_tool_text_round_trips():
  tmp = tempfile.mkdtemp()
  try:
    conn = _conn(tmp, _build_app())
    conn.start()
    try:
      cancel = CancelToken()
      result = conn.call_tool("echo", {"message": "hi"}, cancel=cancel)
      assert isinstance(result, McpToolResult)
      assert result.is_error is False
      assert any("echo:hi" in (i.text or "") for i in result.content)
    finally:
      conn.stop()
  finally:
    shutil.rmtree(tmp)


def exercise_call_tool_failure_returns_error_result():
  tmp = tempfile.mkdtemp()
  try:
    conn = _conn(tmp, _build_app())
    conn.start()
    try:
      cancel = CancelToken()
      result = conn.call_tool("boom", {}, cancel=cancel)
      assert result.is_error is True
      assert any("intentional" in (i.text or "") for i in result.content)
    finally:
      conn.stop()
  finally:
    shutil.rmtree(tmp)


def exercise_cancelled_tool_names_cancellation_in_error():
  """[F6] A cancelled tool call must surface an error_result whose message NAMES
  the cancellation -- not an empty 'Tool call failed: '. In the cancel branch of
  _async_call, `await task` raises asyncio.CancelledError, which on py3.8+ is a
  BaseException and so slips past a bare `except Exception`; the intended
  `raise Sorry("Cancelled")` was dead code, and the cancellation surfaced to
  call_tool as an empty-message concurrent.futures.CancelledError -> error_result
  'Tool call failed: '. Catching asyncio.CancelledError first makes the
  Sorry("Cancelled") reachable, so the message names the cause.

  Revert-proof: with the bare `except Exception` (no asyncio.CancelledError
  branch), the result message is 'Tool call failed: ' with no cause and the
  'cancel' substring assert fails. The fix does not touch
  KeyboardInterrupt/SystemExit (neither is an asyncio.CancelledError)."""
  tmp = tempfile.mkdtemp()
  try:
    conn = _conn(tmp, _build_app())
    conn.start()
    try:
      cancel = CancelToken()
      cancel.set()                       # already cancelled before dispatch
      result = conn.call_tool("echo", {"message": "hi"}, cancel=cancel)
      assert isinstance(result, McpToolResult)
      assert result.is_error is True, result
      msg = " ".join((i.text or "") for i in result.content).lower()
      assert "cancel" in msg, repr(msg)  # names the cancellation, non-empty
    finally:
      conn.stop()
  finally:
    shutil.rmtree(tmp)


def exercise_call_tool_when_not_ready_returns_error():
  tmp = tempfile.mkdtemp()
  try:
    conn = _conn(tmp, _build_app())
    # Skip start(); state stays STOPPED.
    result = conn.call_tool("echo", {"message": "x"}, cancel=CancelToken())
    assert result.is_error is True
    assert "not ready" in str(result.content[0].text)
  finally:
    shutil.rmtree(tmp)


def exercise_call_tool_convert_failure_returns_error_result():
  """call_tool's never-raise contract must cover result conversion too. If
  _convert_result raises (e.g. an attachment-store OSError when persisting
  an image on a full disk), call_tool must return an error result rather
  than letting the exception escape and abort the turn."""
  tmp = tempfile.mkdtemp()
  try:
    conn = _conn(tmp, _build_app())
    conn.start()
    try:
      def _boom(raw):
        raise OSError("disk full")
      # Instance attribute shadows the bound method; called as _boom(raw).
      conn._convert_result = _boom
      cancel = CancelToken()
      result = conn.call_tool("echo", {"message": "hi"}, cancel=cancel)
      assert isinstance(result, McpToolResult)
      assert result.is_error is True
      assert any("disk full" in (i.text or "") for i in result.content), \
        result.content
    finally:
      conn.stop()
  finally:
    shutil.rmtree(tmp)


def exercise_start_failure_after_connect_closes_client():
  """Partial start: __aenter__ connects (spawns the subprocess) but
  list_tools() raises. start() must close the client (await __aexit__) so
  the subprocess is told to exit -- and it must do so BEFORE _stop_loop()
  tears down the loop the __aexit__ coroutine needs. Otherwise the
  subprocess leaks and stop() can't recover (the loop is already None)."""
  import fastmcp
  exited = []

  class _PartialStartClient:
    """Connects on __aenter__ but fails in list_tools, recording whether
    __aexit__ (the subprocess-exit signal) was awaited."""
    def __init__(self, transport):
      pass
    async def __aenter__(self):
      return self
    async def __aexit__(self, *exc):
      exited.append(True)
      return False
    async def list_tools(self):
      raise RuntimeError("list_tools failed after connect")

  tmp = tempfile.mkdtemp()
  saved = fastmcp.Client
  fastmcp.Client = _PartialStartClient
  try:
    conn = _conn(tmp, _build_app())          # in-process branch: skip cmd checks
    raised = False
    try:
      conn.start()
    except Sorry:
      raised = True
    assert raised, "start() must surface a Sorry when list_tools fails"
    assert exited == [True], (
      "client __aexit__ must be awaited on the partial-start failure path "
      "(subprocess told to exit) before the loop is torn down; got %r"
      % exited)
    assert conn.state == McpServerConnection.STATE_FAILED
    # Client handles cleared and the loop torn down after the close.
    assert conn._client_cm is None
    assert conn._loop is None
  finally:
    fastmcp.Client = saved
    shutil.rmtree(tmp)


def exercise_empty_command_raises_clear_sorry():
  """Regression: if the profile's command field is empty (e.g. a
  variable expansion resolved to nothing, or the field was omitted),
  start() must surface a clear error rather than silently shifting
  args[0] into the command slot."""

  class _EmptyCmdCfg:
    name = "phenix"
    command = ""                                 # the bug condition
    args = ["-m", "phenix.mcp"]
    env = {}
    tool_policy = {}

  tmp = tempfile.mkdtemp()
  try:
    storage = ConversationStorage(project_dir=Path(tmp), log=null_out())
    conn = McpServerConnection(
      config=_EmptyCmdCfg(), project_dir=Path(tmp),
      storage=storage, conv_id="c1", log=null_out())
    try:
      conn.start()
    except Sorry as exc:
      # Error must clearly name the cause and point at the profile.
      msg = str(exc).lower()
      assert "empty command" in msg, msg
      assert "profile" in msg, msg
    else:
      raise AssertionError("expected Sorry on empty command")
    assert conn.state == McpServerConnection.STATE_FAILED
    # The except-Sorry branch must also run _stop_loop()'s teardown -- the
    # loop thread started by start() is joined and its handle nulled -- so a
    # rejected start leaks no background thread (today the except-Exception
    # branch asserts this; pin it here for the Sorry branch too).
    assert conn._loop is None
  finally:
    shutil.rmtree(tmp)


def exercise_risk_derived_from_annotations():
  """The fastmcp tool annotations (readOnlyHint / destructiveHint)
  must propagate into McpServerConnection.risk_for() so the approval
  card can colour-code by risk."""
  tmp = tempfile.mkdtemp()
  try:
    conn = _conn(tmp, _build_risk_app())
    conn.start()
    try:
      assert conn.risk_for("look") == "read"
      assert conn.risk_for("nuke") == "destructive"
      assert conn.risk_for("write") == "write"
    finally:
      conn.stop()
  finally:
    shutil.rmtree(tmp)


def exercise_subprocess_env_injects_phenix_when_opted_in():
  from qttbx.widgets.chat.agent.mcp_client import McpServerConnection
  from qttbx.widgets.chat.agent.profile import McpServerConfig
  cfg = McpServerConfig(name="phenix", command="x")  # inject_phenix_env defaults True
  conn = McpServerConnection(config=cfg, project_dir="/tmp/proj",
                             storage=None, conv_id="c1", log=null_out())
  env = conn._subprocess_env()
  assert env["PHENIX_PROJECT_DIR"] == "/tmp/proj"
  assert "PHENIX_CHAT_HOME" in env and env["PHENIX_CHAT_HOME"]


def exercise_subprocess_env_omits_phenix_for_foreign_server():
  """A foreign server (inject_phenix_env=False, e.g. Coot) must not receive the
  Phenix scoping vars -- even when they are EXPORTED in the parent environment
  they are stripped from the subprocess env, so Phenix scoping never leaks. The
  export also makes the test hermetic: it can't pass vacuously on a machine that
  happens not to export them, nor fail spuriously on one that does."""
  import os
  from qttbx.widgets.chat.agent.mcp_client import McpServerConnection
  from qttbx.widgets.chat.agent.profile import McpServerConfig
  saved = {k: os.environ.get(k)
           for k in ("PHENIX_PROJECT_DIR", "PHENIX_CHAT_HOME")}
  os.environ["PHENIX_PROJECT_DIR"] = "/exported/proj"
  os.environ["PHENIX_CHAT_HOME"] = "/exported/home"
  try:
    cfg = McpServerConfig(name="coot", command="x", inject_phenix_env=False)
    conn = McpServerConnection(config=cfg, project_dir="/tmp/proj",
                               storage=None, conv_id="c1", log=null_out())
    env = conn._subprocess_env()
    assert "PHENIX_PROJECT_DIR" not in env, env.get("PHENIX_PROJECT_DIR")
    assert "PHENIX_CHAT_HOME" not in env, env.get("PHENIX_CHAT_HOME")
  finally:
    for k, v in saved.items():
      if v is None:
        os.environ.pop(k, None)
      else:
        os.environ[k] = v


def exercise_run_async_times_out_instead_of_hanging():
  """_run_async must be bounded: a coroutine that never returns (a wedged
  MCP server) raises Sorry within the timeout rather than blocking the
  caller forever. start()/stop() run on the GUI thread, so an unbounded
  wait would freeze the whole chat."""
  import asyncio
  import time
  tmp = tempfile.mkdtemp()
  try:
    conn = _bare_conn(tmp)
    conn._start_loop()
    try:
      async def _never():
        await asyncio.sleep(30)
      t0 = time.time()
      raised = False
      try:
        conn._run_async(_never(), timeout=0.2)
      except Sorry:
        raised = True
      elapsed = time.time() - t0
      assert raised, "expected Sorry when the coroutine exceeds the timeout"
      assert elapsed < 5.0, \
        "a timed-out call must return promptly, not hang: %s" % elapsed
      # Pump the loop once so the cancelled _never task is reaped before we
      # close the loop -- otherwise it is gc'd while pending (a noisy but
      # harmless asyncio warning).
      conn._run_async(asyncio.sleep(0), timeout=2.0)
    finally:
      conn._stop_loop()
  finally:
    shutil.rmtree(tmp)


def exercise_stop_loop_survives_unstoppable_loop_thread():
  """_stop_loop must not raise when the loop thread won't stop in time. If a
  coroutine wedges the loop thread (so the join times out and the loop is
  still running), calling loop.close() would raise 'Cannot close a running
  event loop'. _stop_loop must skip the close and null out its handles."""
  import asyncio
  import threading
  import time
  tmp = tempfile.mkdtemp()
  try:
    conn = _bare_conn(tmp)
    conn._start_loop()
    started = threading.Event()

    async def _wedge():
      started.set()
      time.sleep(3.0)          # blocks the loop thread (no await point)

    asyncio.run_coroutine_threadsafe(_wedge(), conn._loop)
    assert started.wait(timeout=2.0), "wedge coroutine did not start"
    # The join inside _stop_loop times out (the thread is stuck in
    # time.sleep), so the loop is still running -- closing it would raise.
    # _stop_loop must not raise and must clear its handles.
    conn._stop_loop()
    assert conn._loop is None
    assert conn._loop_thread is None
  finally:
    shutil.rmtree(tmp)


def exercise():
  exercise_start_lists_tools()
  exercise_call_tool_text_round_trips()
  exercise_call_tool_failure_returns_error_result()
  exercise_cancelled_tool_names_cancellation_in_error()
  exercise_call_tool_when_not_ready_returns_error()
  exercise_call_tool_convert_failure_returns_error_result()
  exercise_start_failure_after_connect_closes_client()
  exercise_empty_command_raises_clear_sorry()
  exercise_risk_derived_from_annotations()
  exercise_subprocess_env_injects_phenix_when_opted_in()
  exercise_subprocess_env_omits_phenix_for_foreign_server()
  exercise_run_async_times_out_instead_of_hanging()
  exercise_stop_loop_survives_unstoppable_loop_thread()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
