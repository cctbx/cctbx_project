"""Tool policy, registry, and approval message types.

The registry unifies built-in tools, skill-wrapped tools, and MCP tools
under one dispatch surface. ``ToolPolicy`` resolves allow / ask / deny
decisions with session-scoped overrides.
"""

import sys
from dataclasses import dataclass

from libtbx.utils import Sorry

from qttbx.widgets.chat.agent.errors import AgentEvent


# ---- approval message types ------------------------------------------------

@dataclass
class ToolApprovalRequest(AgentEvent):
  """Event asking the user to approve a pending tool call.

  Parameters
  ----------
  request_id : str
      Unique id correlating this request with its response.
  tool_name : str
      Name of the tool awaiting approval.
  tool_source : str
      Origin of the tool: ``'builtin'``, ``'skill'``, or
      ``'mcp:<server>'``.
  input : dict
      The arguments the tool would be invoked with.
  risk : str
      Risk level: ``'read'``, ``'write'``, or ``'destructive'``.
  summary : str, optional
      Short human-readable summary of the call for the approval card.
  batch_id : str, optional
      Identifier grouping requests issued together in one batch.
  """
  request_id: str
  tool_name: str
  tool_source: str                 # 'builtin' | 'skill' | 'mcp:<server>'
  input: dict
  risk: str                        # 'read' | 'write' | 'destructive'
  summary: str = None
  batch_id: str = None


@dataclass
class ToolApprovalResponse:
  """User's decision on a pending tool call.

  Parameters
  ----------
  request_id : str
      The id of the ``ToolApprovalRequest`` this answers.
  decision : str
      The user's choice: ``'approve'``, ``'deny'``, or
      ``'deny_and_stop'``.
  remember : str, optional
      Scope to remember the choice for this session: ``'none'``
      (default), ``'tool'``, or ``'server'``.
  """
  request_id: str
  decision: str                    # 'approve' | 'deny' | 'deny_and_stop'
  remember: str = "none"           # 'none' | 'tool' | 'server'


class _Cancelled:
  """Sentinel pushed into the approval queue when the turn is cancelled.

  ``queue.Queue.get()`` with no timeout doesn't observe a ``CancelToken``
  — we must push something to wake it.
  """


# ---- policy ----------------------------------------------------------------

class ToolPolicy:
  """allow / ask / deny policy with session-scoped remembered choices.

  Resolution order, highest first: a per-tool entry (from profile or
  session memory); then a per-server entry (from profile
  ``mcp_servers[].tool_policy['*']``, which requires the
  ``tool_to_source`` mapping to know which server owns a tool); then the
  default (``'ask'`` unless overridden by ``profile.tool_policy_default``).

  Parameters
  ----------
  default : str, optional
      Fallback decision when no per-tool or per-server entry matches.
      Defaults to ``'ask'``.
  per_tool : dict, optional
      Mapping of tool name to ``allow`` / ``ask`` / ``deny``.
  per_server : dict, optional
      Mapping of server name to ``allow`` / ``ask`` / ``deny``.
  tool_to_source : dict, optional
      Mapping of tool name to its source (``'mcp:<server>'`` or
      ``'builtin'`` / ``'skill'``), used to resolve per-server entries.
  """

  def __init__(self, default="ask",
               per_tool=None, per_server=None,
               tool_to_source=None):
    self.default = default
    self.per_tool = dict(per_tool or {})           # name -> allow|ask|deny
    self.per_server = dict(per_server or {})       # server -> allow|ask|deny
    self.tool_to_source = dict(tool_to_source or {})  # tool -> 'mcp:server' or 'builtin'/'skill'

  def resolve(self, tool_name):
    """Return the ``allow`` / ``ask`` / ``deny`` decision for a tool."""
    if tool_name in self.per_tool:
      return self.per_tool[tool_name]
    source = self.tool_to_source.get(tool_name, "")
    if source.startswith("mcp:"):
      server = source.split(":", 1)[1]
      if server in self.per_server:
        return self.per_server[server]
    return self.default

  def allow_tool_for_session(self, tool_name):
    """Remember ``allow`` for this tool for the rest of the session."""
    self.per_tool[tool_name] = "allow"

  def allow_server_for_session(self, server_name):
    """Remember ``allow`` for this server for the rest of the session."""
    self.per_server[server_name] = "allow"


# ---- registry --------------------------------------------------------------

@dataclass
class _ToolEntry:
  """Internal registry record pairing a tool spec with its handler."""
  spec: object                     # ToolSpec
  source: str                      # 'builtin' | 'skill' | 'mcp:<server>'
  handler: object = None           # callable; signature varies by source
  risk: str = "write"


class ToolRegistry:
  """Union of built-in, skill, and MCP tools.

  Collision resolution: built-in beats skill beats MCP. Within MCP,
  the first-registered tool keeps the bare name and later ones get
  ``<server>:<tool>`` namespacing.
  """

  def __init__(self, log=None):
    self.log = log if log is not None else sys.stdout
    self._entries = {}             # name -> _ToolEntry

  # ---- registration --------------------------------------------------------

  def register_builtin(self, spec, handler, risk="write"):
    """Register a built-in tool under ``spec.name``."""
    self._add(spec.name, _ToolEntry(
      spec=spec, source="builtin", handler=handler, risk=risk))

  def register_skill_tool(self, spec, handler):
    """Register a skill-wrapped tool (always ``read`` risk)."""
    self._add(spec.name, _ToolEntry(
      spec=spec, source="skill", handler=handler, risk="read"))

  def register_mcp_tool(self, spec, server_name, handler, risk="write"):
    """Register an MCP-sourced tool.

    On name collision with an existing entry, the new MCP tool is
    registered as ``<server>:<name>`` instead of the bare ``name``.
    """
    bare = spec.name
    if bare in self._entries and self._entries[bare].source != "mcp:" + server_name:
      namespaced = "%s:%s" % (server_name, bare)
      print("tool '%s' collides; registering MCP server '%s' tool as '%s'"
            % (bare, server_name, namespaced), file=self.log)
      # Replace name in a copy of the spec.
      from qttbx.widgets.chat.agent.base import ToolSpec
      spec = ToolSpec(name=namespaced,
                      description=spec.description,
                      input_schema=spec.input_schema)
    self._add(spec.name, _ToolEntry(
      spec=spec, source="mcp:" + server_name, handler=handler, risk=risk))

  def _add(self, name, entry):
    """Insert an entry under ``name``, skipping if already registered."""
    if name in self._entries:
      existing = self._entries[name].source
      print("tool '%s' already registered from %s; skipping new %s"
            % (name, existing, entry.source), file=self.log)
      return
    self._entries[name] = entry

  # ---- queries -------------------------------------------------------------

  def specs(self):
    """Return the list of registered tool specs."""
    return [e.spec for e in self._entries.values()]

  def source_of(self, name):
    """Return a tool's source string, or ``None`` if unregistered."""
    return self._entries[name].source if name in self._entries else None

  def server_of(self, name):
    """Return the MCP server owning a tool, or ``None`` if not MCP."""
    src = self.source_of(name) or ""
    if src.startswith("mcp:"):
      return src.split(":", 1)[1]
    return None

  def risk_of(self, name):
    """Return a tool's risk level, defaulting to ``write`` if unknown."""
    return self._entries[name].risk if name in self._entries else "write"

  # ---- invocation ----------------------------------------------------------

  def invoke_builtin(self, name, input, cancel, session, tool_use_id):
    """Invoke a built-in tool's handler.

    Raises
    ------
    libtbx.utils.Sorry
        If ``name`` is not a registered built-in tool.
    """
    e = self._entries.get(name)
    if e is None or e.source != "builtin":
      raise Sorry("Built-in tool not found: %s" % name)
    return e.handler(name=name, input=input, cancel=cancel,
                     session=session, tool_use_id=tool_use_id)

  def invoke_skill(self, name, input):
    """Invoke a skill tool's handler.

    Raises
    ------
    libtbx.utils.Sorry
        If ``name`` is not a registered skill tool.
    """
    e = self._entries.get(name)
    if e is None or e.source != "skill":
      raise Sorry("Skill tool not found: %s" % name)
    return e.handler(name=name, input=input)

  def invoke_mcp(self, name, input, cancel):
    """Invoke an MCP tool's handler.

    Raises
    ------
    libtbx.utils.Sorry
        If ``name`` is not a registered MCP tool.
    """
    e = self._entries.get(name)
    if e is None or not e.source.startswith("mcp:"):
      raise Sorry("MCP tool not found: %s" % name)
    return e.handler(name=name, input=input, cancel=cancel)

  # ---- session memory snapshot ---------------------------------------------

  def tool_to_source_map(self):
    """Return a name-to-source snapshot for ``ToolPolicy`` construction."""
    return {name: e.source for name, e in self._entries.items()}
