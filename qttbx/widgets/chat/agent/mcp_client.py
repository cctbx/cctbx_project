"""MCP client types, content converters, and ``McpServerConnection``.

The dataclasses live here (not in ``session.py``) so a session that
doesn't use MCP doesn't pull ``fastmcp`` into its import graph; only
the function-scope ``fastmcp`` imports inside ``McpServerConnection``
touch the SDK.
"""

import asyncio
import os
import shutil
import sys
import threading
from dataclasses import dataclass, field

from libtbx.utils import Sorry

from qttbx.widgets.chat.agent.base import ToolSpec
from qttbx.widgets.chat.agent.conversation import ContentBlock


@dataclass
class McpToolItem:
  type: str                        # "text" | "image" | "resource"
  text: str = None                 # for text and resource excerpts
  uri: str = None                  # for resource
  sha256: str = None               # for image: populated by the client
                                   # after storing the bytes via storage
  mime: str = None                 # for image
  caption: str = None              # for image (optional)
  text_excerpt: str = None         # for resource (optional short snippet)


@dataclass
class McpToolResult:
  content: list = field(default_factory=list)         # list[McpToolItem]
  is_error: bool = False


def error_result(message):
  """Construct an McpToolResult that the dispatch loop will surface as
  tool_result is_error=true content."""
  return McpToolResult(
    content=[McpToolItem(type="text", text=message)],
    is_error=True)


def _derive_risk(raw):
  ann = getattr(raw, "annotations", None) or {}
  # Pydantic models vs plain dicts both supported.
  if not isinstance(ann, dict):
    ann = getattr(ann, "model_dump", lambda: {})() or {}
  if ann.get("destructiveHint"):
    return "destructive"
  if ann.get("readOnlyHint"):
    return "read"
  return "write"


def _mcp_item_to_block(item, storage, conv_id):
  """Translate one MCP content item into a canonical ContentBlock. Image
  bytes were already stored by McpServerConnection._convert_result, which
  populated item.sha256 — we just reference it here. The storage and
  conv_id arguments are kept for symmetry with future paths that may need
  to write at conversion time."""
  if item.type == "text":
    return ContentBlock(type="text", data={"text": item.text or ""})
  if item.type == "image":
    return ContentBlock(type="image", data={
      "attachment_sha256": item.sha256,
      "mime": item.mime,
      "caption": item.caption,
    })
  if item.type == "resource":
    excerpt = item.text_excerpt or ""
    return ContentBlock(type="text", data={
      "text": "Resource %s\n%s" % (item.uri, excerpt)})
  return ContentBlock(type="text", data={
    "text": "Unsupported MCP content type: %s" % item.type})


class McpServerConnection:
  """One MCP server subprocess and its tools.

  Lifecycle: ``STOPPED`` -> ``STARTING`` -> ``READY`` -> (``STOPPED`` |
  ``STOPPED_UNEXPECTEDLY`` | ``FAILED``). Public surface: ``start()``,
  ``stop()``, ``call_tool(name, input, cancel)``. Tools are exposed as
  ``ToolSpec`` instances after ``start()`` so the ``ToolRegistry`` can
  register them with ``register_mcp_tool``.

  Async-sync bridge: ``fastmcp``'s ``Client`` is async but
  ``AgentSession`` (which calls ``call_tool``) is sync. A dedicated
  daemon thread runs an asyncio event loop; calls are dispatched via
  ``asyncio.run_coroutine_threadsafe``.
  """

  STATE_STOPPED = "stopped"
  STATE_STARTING = "starting"
  STATE_READY = "ready"
  STATE_FAILED = "failed"
  STATE_STOPPED_UNEXPECTEDLY = "stopped_unexpectedly"

  def __init__(self, config, project_dir, storage, conv_id, log=None):
    self.config = config
    self.project_dir = project_dir
    self.storage = storage
    self.conv_id = conv_id
    self.log = log if log is not None else sys.stdout
    self.state = self.STATE_STOPPED
    self.tools = []
    self._risk_by_name = {}
    self._client = None
    self._loop = None
    self._loop_thread = None
    self._in_process_app = None
    self._client_cm = None                            # async context mgr

  # ---- public API ----------------------------------------------------------

  def start(self):
    """Spawn the subprocess (or use the injected in-process app for tests),
    list tools, and mark READY. Raises Sorry on hard failures."""
    self.state = self.STATE_STARTING
    self._start_loop()
    try:
      self._run_async(self._async_start())
      self.state = self.STATE_READY
      print("MCP %s: ready (%d tools)" %
            (self.config.name, len(self.tools)), file=self.log)
    except Sorry:
      self.state = self.STATE_FAILED
      raise
    except Exception as exc:
      self.state = self.STATE_FAILED
      raise Sorry("MCP %s: start failed: %s" %
                  (self.config.name, exc))

  def stop(self):
    if self._client_cm is not None:
      try:
        self._run_async(self._client_cm.__aexit__(None, None, None))
      except Exception:
        # Best-effort shutdown; subprocess may already be gone.
        pass
      self._client_cm = None
      self._client = None
    self._stop_loop()
    self.state = self.STATE_STOPPED

  def call_tool(self, name, input, cancel, timeout=60.0):
    """Returns McpToolResult. On failure (server down, tool raised,
    cancel) the result has is_error=True so the dispatch loop in
    AgentSession surfaces it to the model rather than aborting the turn."""
    if self.state != self.STATE_READY:
      return error_result("MCP server '%s' not ready" % self.config.name)
    try:
      raw = self._run_async(
        self._async_call(name, input, timeout=timeout, cancel=cancel))
    except Exception as exc:
      return error_result("Tool call failed: %s" % exc)
    return self._convert_result(raw)

  # ---- test hooks ----------------------------------------------------------

  def _inject_in_process_app(self, app):
    """Test seam: swap in a FastMCP app to be served via FastMCPTransport
    instead of spawning a subprocess. Used by unit tests only."""
    self._in_process_app = app

  # ---- async surface -------------------------------------------------------

  async def _async_start(self):
    from fastmcp import Client
    if self._in_process_app is not None:
      from fastmcp.client.transports import FastMCPTransport
      transport = FastMCPTransport(self._in_process_app)
    else:
      from fastmcp.client.transports import StdioTransport
      cmd = self._expanded_command()
      if not cmd or not cmd[0]:
        raise Sorry(
          "MCP server '%s': empty command. The profile's `command` field "
          "resolved to an empty string — check the profile JSON for an "
          "unset variable expansion or a missing command entry."
          % self.config.name)
      if not shutil.which(cmd[0]):
        raise Sorry("MCP server '%s': command not found: %s" %
                    (self.config.name, cmd[0]))
      env = dict(os.environ)
      env["PHENIX_PROJECT_DIR"] = str(self.project_dir)
      if getattr(self.config, "env", None):
        env.update(self.config.env)
      transport = StdioTransport(command=cmd[0], args=cmd[1:], env=env)
    self._client_cm = Client(transport)
    self._client = await self._client_cm.__aenter__()
    raw_tools = await self._client.list_tools()
    self.tools = [self._tool_to_spec(t) for t in raw_tools]

  async def _async_call(self, name, input, timeout, cancel):
    """Run the call as a Task so we can cancel it when the chat cancels.
    raise_on_error=False makes fastmcp return CallToolResult with
    is_error=True instead of raising ToolError on tool-raised exceptions."""
    coro = self._client.call_tool(
      name, input, timeout=timeout, raise_on_error=False)
    task = asyncio.ensure_future(coro)
    while not task.done():
      if cancel.is_set():
        task.cancel()
        try:
          await task
        except Exception:
          pass
        raise Sorry("Cancelled")
      try:
        return await asyncio.wait_for(asyncio.shield(task), timeout=0.1)
      except asyncio.TimeoutError:
        continue
    return task.result()

  # ---- conversion ----------------------------------------------------------

  def _expanded_command(self):
    # Preserve the leading element even when empty so callers can detect
    # the "command unset" case (a previous filter dropped it and silently
    # promoted args[0] to the command slot). _async_start validates.
    cmd = [self.config.command or ""] + list(getattr(self.config, "args", []))
    return cmd

  def _tool_to_spec(self, raw):
    name = raw.name
    self._risk_by_name[name] = _derive_risk(raw)
    return ToolSpec(
      name=name,
      description=getattr(raw, "description", "") or "",
      input_schema=getattr(raw, "inputSchema", None)
        or getattr(raw, "input_schema", {}) or {})

  def risk_for(self, name):
    return self._risk_by_name.get(name, "write")

  def _convert_result(self, raw):
    items = []
    is_error = bool(getattr(raw, "is_error", False))
    content = getattr(raw, "content", None) or []
    for entry in content:
      t = getattr(entry, "type", None)
      if t == "text":
        items.append(McpToolItem(type="text",
                                 text=getattr(entry, "text", "")))
      elif t == "image":
        import base64
        data_b64 = getattr(entry, "data", "")
        mime = getattr(entry, "mimeType", None) or \
               getattr(entry, "mime_type", "image/png")
        try:
          data = base64.b64decode(data_b64)
        except Exception:
          data = b""
        att = self.storage.store_attachment(
          self.conv_id, data, mime)
        items.append(McpToolItem(type="image",
                                 sha256=att.sha256, mime=mime))
      elif t == "resource":
        items.append(McpToolItem(type="resource",
                                 uri=getattr(entry, "uri", ""),
                                 text_excerpt=getattr(entry, "text", "")))
    return McpToolResult(content=items, is_error=is_error)

  # ---- event loop helpers --------------------------------------------------

  def _start_loop(self):
    if self._loop is not None:
      return
    self._loop = asyncio.new_event_loop()
    self._loop_thread = threading.Thread(
      target=self._loop.run_forever,
      name="mcp-loop-%s" % self.config.name, daemon=True)
    self._loop_thread.start()

  def _stop_loop(self):
    if self._loop is None:
      return
    self._loop.call_soon_threadsafe(self._loop.stop)
    self._loop_thread.join(timeout=2.0)
    self._loop.close()
    self._loop = None
    self._loop_thread = None

  def _run_async(self, coro):
    future = asyncio.run_coroutine_threadsafe(coro, self._loop)
    return future.result()
