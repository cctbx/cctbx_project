"""MCP client types, content converters, and ``McpServerConnection``.

The dataclasses live here (not in ``session.py``) so a session that
doesn't use MCP doesn't pull ``fastmcp`` into its import graph; only
the function-scope ``fastmcp`` imports inside ``McpServerConnection``
touch the SDK.
"""

import asyncio
import shutil
import sys
import threading
from dataclasses import dataclass, field

from libtbx.utils import Sorry

from qttbx.widgets.chat.agent.base import ToolSpec
from qttbx.widgets.chat.agent.conversation import ContentBlock


def phenix_server_env(project_dir):
  """The PHENIX_* env a phenix-aware MCP server subprocess should receive.

  Centralizes the project/chat-home env so the launcher's claude_code SDK
  loop and McpServerConnection inject identical values.

  Parameters
  ----------
  project_dir : str or pathlib.Path
      The active project directory.

  Returns
  -------
  dict
      ``{"PHENIX_PROJECT_DIR": ..., "PHENIX_CHAT_HOME": ...}``.
  """
  from qttbx.widgets.chat.agent.paths import chat_root_for
  return {
    "PHENIX_PROJECT_DIR": str(project_dir),
    "PHENIX_CHAT_HOME": str(chat_root_for(project_dir)),
  }


@dataclass
class McpToolItem:
  """One content item from an MCP tool result.

  Which fields are populated depends on ``type``: ``text`` uses ``text``;
  ``image`` uses ``sha256``, ``mime``, and optional ``caption``;
  ``resource`` uses ``uri`` and optional ``text_excerpt``.
  """
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
  """Result of an MCP tool call: content items plus an error flag."""
  content: list = field(default_factory=list)         # list[McpToolItem]
  is_error: bool = False


def error_result(message):
  """Construct an ``McpToolResult`` flagged as an error.

  The dispatch loop surfaces it as ``tool_result`` content with
  ``is_error=True``.

  Parameters
  ----------
  message : str
      Error text to wrap as the single text item of the result.

  Returns
  -------
  McpToolResult
      A result with ``is_error=True`` carrying ``message``.
  """
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
  """Translate one MCP content item into a canonical ``ContentBlock``.

  Image bytes were already stored by
  ``McpServerConnection._convert_result``, which populated ``item.sha256``
  -- we just reference it here.

  Parameters
  ----------
  item : McpToolItem
      The MCP content item to translate.
  storage : object
      Attachment store. Kept for symmetry with future paths that may need
      to write at conversion time.
  conv_id : str
      Conversation id. Kept for symmetry with future paths that may need
      to write at conversion time.

  Returns
  -------
  ContentBlock
      The canonical block; a text block describing the type for
      unsupported content types.
  """
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

  Parameters
  ----------
  config : McpServerConfig
      Server entry (name, command, args, env) from the profile.
  project_dir : str or pathlib.Path
      Project directory. Exported to the subprocess as
      ``PHENIX_PROJECT_DIR`` only for phenix-aware servers
      (``config.inject_phenix_env`` True); foreign servers do not see it.
  storage : object
      Attachment store used to persist image content from tool results.
  conv_id : str
      Conversation id under which attachments are stored.
  log : file-like, optional
      Destination for status messages. Defaults to ``sys.stdout``.
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
    """Spawn the subprocess, list tools, and mark the server ready.

    Uses the injected in-process app instead of a subprocess when one was
    set via ``_inject_in_process_app`` (tests).

    Raises
    ------
    libtbx.utils.Sorry
        On hard start failures (e.g. missing command, connection error).
    """
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
    """Close the client, stop the event loop, and mark the server stopped.

    Best-effort: a subprocess that has already exited does not raise.
    """
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
    """Invoke an MCP tool and return its converted result.

    On failure (server down, tool raised, cancel) the result has
    ``is_error=True`` so the dispatch loop in ``AgentSession`` surfaces it
    to the model rather than aborting the turn.

    Parameters
    ----------
    name : str
        Tool name as exposed by the server.
    input : dict
        Tool arguments matching the tool's input schema.
    cancel : CancelToken
        Polled while the call is in flight; cancels the task when set.
    timeout : float, optional
        Per-call timeout in seconds. Defaults to ``60.0``.

    Returns
    -------
    McpToolResult
        The converted result, possibly flagged ``is_error=True``.
    """
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
    """Swap in a FastMCP app served via ``FastMCPTransport`` (test seam).

    Avoids spawning a subprocess. Used by unit tests only.

    Parameters
    ----------
    app : object
        A FastMCP app instance to serve in-process.
    """
    self._in_process_app = app

  # ---- async surface -------------------------------------------------------

  def _subprocess_env(self):
    """Build the subprocess environment for this server.

    Starts from os.environ, merges the server's own sanitized env, then --
    only when the server opts in via inject_phenix_env (default True) --
    overlays the authoritative PHENIX_PROJECT_DIR / PHENIX_CHAT_HOME.

    Returns
    -------
    dict
        The environment to hand to the spawned MCP-server subprocess.
    """
    import os
    env = dict(os.environ)
    if getattr(self.config, "env", None):
      from qttbx.widgets.chat.agent.profile import sanitize_server_env
      env.update(sanitize_server_env(self.config.env))
    if getattr(self.config, "inject_phenix_env", True):
      env.update(phenix_server_env(self.project_dir))
    else:
      # A foreign server (inject_phenix_env=False, e.g. Coot) must not see the
      # Phenix scoping vars. The flag skips the overlay above, but they can
      # still leak in from os.environ if the user exported them -- strip those
      # so scoping never reaches a non-Phenix subprocess.
      for k in phenix_server_env(self.project_dir):
        env.pop(k, None)
    return env

  async def _async_start(self):
    """Open the client connection and populate ``self.tools``.

    Validates the resolved command and raises before connecting.

    Raises
    ------
    libtbx.utils.Sorry
        When the command is empty or not found on ``PATH``.
    """
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
      env = self._subprocess_env()
      transport = StdioTransport(command=cmd[0], args=cmd[1:], env=env)
    self._client_cm = Client(transport)
    self._client = await self._client_cm.__aenter__()
    raw_tools = await self._client.list_tools()
    self.tools = [self._tool_to_spec(t) for t in raw_tools]

  async def _async_call(self, name, input, timeout, cancel):
    """Run the tool call as a cancellable asyncio task.

    Runs the call as a ``Task`` so it can be cancelled when the chat
    cancels. ``raise_on_error=False`` makes fastmcp return a
    ``CallToolResult`` with ``is_error=True`` instead of raising
    ``ToolError`` on tool-raised exceptions.

    Parameters
    ----------
    name : str
        Tool name.
    input : dict
        Tool arguments.
    timeout : float
        Per-call timeout in seconds.
    cancel : CancelToken
        Polled while the task runs; cancels it when set.

    Returns
    -------
    object
        The raw ``CallToolResult`` from fastmcp.

    Raises
    ------
    libtbx.utils.Sorry
        When ``cancel`` is set before the call completes.
    """
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
    """Build the command vector, preserving a possibly-empty command slot.

    The leading element is kept even when empty so callers can detect the
    "command unset" case; ``_async_start`` validates it.

    Returns
    -------
    list of str
        ``[command] + args`` with ``command`` coerced to ``""`` if unset.
    """
    # An empty leading element guards against a previous filter dropping it
    # and silently promoting args[0] to the command slot.
    cmd = [self.config.command or ""] + list(getattr(self.config, "args", []))
    return cmd

  def _tool_to_spec(self, raw):
    """Build a ``ToolSpec`` from a raw MCP tool and record its risk."""
    name = raw.name
    self._risk_by_name[name] = _derive_risk(raw)
    return ToolSpec(
      name=name,
      description=getattr(raw, "description", "") or "",
      input_schema=getattr(raw, "inputSchema", None)
        or getattr(raw, "input_schema", {}) or {})

  def risk_for(self, name):
    """Return the risk class for ``name`` (``"write"`` if unknown)."""
    return self._risk_by_name.get(name, "write")

  def _convert_result(self, raw):
    """Convert a raw fastmcp result into an ``McpToolResult``.

    Text and resource items are copied through; image items are
    base64-decoded and persisted via ``self.storage``, with the resulting
    sha256 referenced on the item.

    Returns
    -------
    McpToolResult
        The converted result, preserving the raw ``is_error`` flag.
    """
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
