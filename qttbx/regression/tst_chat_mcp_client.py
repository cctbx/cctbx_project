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


def exercise():
  exercise_start_lists_tools()
  exercise_call_tool_text_round_trips()
  exercise_call_tool_failure_returns_error_result()
  exercise_call_tool_when_not_ready_returns_error()
  exercise_empty_command_raises_clear_sorry()
  exercise_risk_derived_from_annotations()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
