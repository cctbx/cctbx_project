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
  request_id: str
  tool_name: str
  tool_source: str                 # 'builtin' | 'skill' | 'mcp:<server>'
  input: dict
  risk: str                        # 'read' | 'write' | 'destructive'
  summary: str = None
  batch_id: str = None


@dataclass
class ToolApprovalResponse:
  request_id: str
  decision: str                    # 'approve' | 'deny' | 'deny_and_stop'
  remember: str = "none"           # 'none' | 'tool' | 'server'


class _Cancelled:
  """Sentinel pushed into the approval queue when the turn is cancelled.
  ``queue.Queue.get()`` with no timeout doesn't observe a
  ``CancelToken`` — we must push something to wake it."""


# ---- policy ----------------------------------------------------------------

class ToolPolicy:
  """allow / ask / deny policy with session-scoped remembered choices.

  Resolution order:
    1. Per-tool entry (from profile or session memory) — highest.
    2. Per-server entry (from profile mcp_servers[].tool_policy['*']) —
       requires tool_to_source mapping to know which server owns a tool.
    3. Default ('ask' unless overridden by profile.tool_policy_default).
  """

  def __init__(self, default="ask",
               per_tool=None, per_server=None,
               tool_to_source=None):
    self.default = default
    self.per_tool = dict(per_tool or {})           # name -> allow|ask|deny
    self.per_server = dict(per_server or {})       # server -> allow|ask|deny
    self.tool_to_source = dict(tool_to_source or {})  # tool -> 'mcp:server' or 'builtin'/'skill'

  def resolve(self, tool_name):
    if tool_name in self.per_tool:
      return self.per_tool[tool_name]
    source = self.tool_to_source.get(tool_name, "")
    if source.startswith("mcp:"):
      server = source.split(":", 1)[1]
      if server in self.per_server:
        return self.per_server[server]
    return self.default

  def allow_tool_for_session(self, tool_name):
    self.per_tool[tool_name] = "allow"

  def allow_server_for_session(self, server_name):
    self.per_server[server_name] = "allow"


# ---- registry --------------------------------------------------------------

@dataclass
class _ToolEntry:
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
    self._add(spec.name, _ToolEntry(
      spec=spec, source="builtin", handler=handler, risk=risk))

  def register_skill_tool(self, spec, handler):
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
    if name in self._entries:
      existing = self._entries[name].source
      print("tool '%s' already registered from %s; skipping new %s"
            % (name, existing, entry.source), file=self.log)
      return
    self._entries[name] = entry

  # ---- queries -------------------------------------------------------------

  def specs(self):
    return [e.spec for e in self._entries.values()]

  def source_of(self, name):
    return self._entries[name].source if name in self._entries else None

  def server_of(self, name):
    src = self.source_of(name) or ""
    if src.startswith("mcp:"):
      return src.split(":", 1)[1]
    return None

  def risk_of(self, name):
    return self._entries[name].risk if name in self._entries else "write"

  # ---- invocation ----------------------------------------------------------

  def invoke_builtin(self, name, input, cancel, session, tool_use_id):
    e = self._entries.get(name)
    if e is None or e.source != "builtin":
      raise Sorry("Built-in tool not found: %s" % name)
    return e.handler(name=name, input=input, cancel=cancel,
                     session=session, tool_use_id=tool_use_id)

  def invoke_skill(self, name, input):
    e = self._entries.get(name)
    if e is None or e.source != "skill":
      raise Sorry("Skill tool not found: %s" % name)
    return e.handler(name=name, input=input)

  def invoke_mcp(self, name, input, cancel):
    e = self._entries.get(name)
    if e is None or not e.source.startswith("mcp:"):
      raise Sorry("MCP tool not found: %s" % name)
    return e.handler(name=name, input=input, cancel=cancel)

  # ---- session memory snapshot ---------------------------------------------

  def tool_to_source_map(self):
    """Snapshot for ToolPolicy construction."""
    return {name: e.source for name, e in self._entries.items()}
