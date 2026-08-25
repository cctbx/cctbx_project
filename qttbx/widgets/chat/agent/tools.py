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
  batch_id : str, optional
      Identifier grouping requests issued together in one batch.
  allow_remember : bool, optional
      Whether the approval card may offer 'Always allow this tool'
      (standing per-tool auto-approval). Defaults to ``True``; a
      destructive tool sets ``False`` to force a decision every time.
  """
  request_id: str
  tool_name: str
  tool_source: str                 # 'builtin' | 'skill' | 'mcp:<server>'
  input: dict
  risk: str                        # 'read' | 'write' | 'destructive'
  batch_id: str = None
  allow_remember: bool = True      # False -> card omits "Always allow this tool"


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
      (default) or ``'tool'``.
  """
  request_id: str
  decision: str                    # 'approve' | 'deny' | 'deny_and_stop'
  remember: str = "none"           # 'none' | 'tool'


def record_tool_remember(response, tool_name_of, allow_remember, remember):
  """Apply a remember='tool' approval to a backend's session-allow sink.

  The single implementation of the remember='tool' contract, shared by both
  backends -- the API path (``AgentSession``) and the SDK path
  (``ClaudeCodeAgent``), which differ only in where the request_id->tool_name
  map, the opt-out recheck, and the allow set live. Records nothing unless the
  user APPROVED with remember='tool' AND the tool still permits being remembered
  (the data-layer opt-out floor, enforced here rather than trusting the card's
  checkbox). Called on the GUI thread BEFORE the coordinator resolves the
  future, so a click that races a cancel still records "always allow".

  Parameters
  ----------
  response : ToolApprovalResponse
      The user's decision.
  tool_name_of : callable
      ``request_id -> tool_name or None`` (the per-request map recorded at
      ``open``). Returns None for a backend that never surfaced this request, so
      that backend's call is a no-op -- which is how the runner can hand every
      response to both backends and have exactly one record it.
  allow_remember : callable
      ``tool_name -> bool``; a False vetoes recording (e.g. a force-close tool
      registered ``allow_remember=False``).
  remember : callable
      ``tool_name -> None``; the sink that marks the tool session-allowed.
  """
  if (getattr(response, "decision", None) != "approve"
      or getattr(response, "remember", "none") != "tool"):
    return
  tool_name = tool_name_of(response.request_id)
  if tool_name is None or not allow_remember(tool_name):
    return
  remember(tool_name)


class _Cancelled:
  """Sentinel pushed into the approval queue when the turn is cancelled.

  ``queue.Queue.get()`` with no timeout doesn't observe a ``CancelToken``
  — we must push something to wake it.
  """


# ---- policy ----------------------------------------------------------------

class ToolPolicy:
  """allow / ask / deny policy with session-scoped remembered choices.

  Resolution order, highest first: an exact per-tool entry (a session
  ``allow`` remembered under the registered name); then the profile's
  per-server-scoped per-tool entry for the tool's owning server (so a
  collision-renamed ``<server>:<tool>`` still resolves to its own server's
  decision); then the per-server (``'*'``) entry; then the default (``'ask'``
  unless overridden by ``profile.tool_policy_default``). The per-server and
  per-server-scoped lookups require the ``tool_to_source`` mapping to know
  which server owns a tool. There is deliberately no bare-name fallback for a
  renamed tool, so one server's allow can never bleed onto another's.

  Parameters
  ----------
  default : str, optional
      Fallback decision when nothing more specific matches. Defaults to
      ``'ask'``.
  per_tool : dict, optional
      Mapping of registered tool name to ``allow`` / ``ask`` / ``deny`` --
      session-remembered allows and directly-constructed (non-server-scoped)
      entries.
  per_server : dict, optional
      Mapping of server name to ``allow`` / ``ask`` / ``deny`` (the
      ``tool_policy['*']`` entry).
  tool_to_source : dict, optional
      Mapping of registered tool name to its source (``'mcp:<server>'`` or
      ``'builtin'`` / ``'skill'``), used to resolve per-server entries.
  per_server_tool : dict, optional
      Mapping of ``(server_name, bare_tool_name)`` to ``allow`` / ``ask`` /
      ``deny`` -- the profile's per-tool policy, kept server-scoped so a
      cross-server name collision can't make one server's decision clobber
      another's. Built by :meth:`from_server_configs`.
  """

  def __init__(self, default="ask",
               per_tool=None, per_server=None,
               tool_to_source=None, per_server_tool=None):
    self.default = default
    self.per_tool = dict(per_tool or {})           # name -> allow|ask|deny
    self.per_server = dict(per_server or {})       # server -> allow|ask|deny
    # (server, bare tool) -> decision. Profile per-tool policy is per-server,
    # so it must stay keyed by its owning server: an MCP tool whose name
    # collides across servers is registered as '<server>:<tool>', and this
    # keying lets resolve() apply each server's own decision to the right tool
    # regardless of the rename (instead of flattening both into one bare key).
    self.per_server_tool = dict(per_server_tool or {})
    self.tool_to_source = dict(tool_to_source or {})  # tool -> 'mcp:server' or 'builtin'/'skill'

  def _server_and_bare(self, tool_name):
    """Return ``(server, bare_tool)`` for *tool_name*, else ``(None, tool_name)``.

    A collision-renamed MCP tool is registered as ``<server>:<tool>``; recover
    the owning server (from ``tool_to_source``) and strip its prefix so the
    per-server-scoped policy keys match whether or not the tool was renamed.
    """
    source = self.tool_to_source.get(tool_name, "")
    if not source.startswith("mcp:"):
      return None, tool_name
    server = source.split(":", 1)[1]
    prefix = server + ":"
    bare = tool_name[len(prefix):] if tool_name.startswith(prefix) else tool_name
    return server, bare

  def resolve(self, tool_name):
    """Return the ``allow`` / ``ask`` / ``deny`` decision for a tool.

    Order, highest first: an exact ``per_tool`` entry (a session-remembered
    allow, keyed by the registered name); the profile's per-server-scoped
    per-tool entry for the tool's owning server; the server-wide (``'*'``)
    entry; then the default.

    There is deliberately NO bare-name ``per_tool`` fallback for a
    collision-renamed ``<server>:<tool>``: it would let one server's session
    allow (or a seeded bare builtin/skill name) override a DIFFERENT server's
    decision -- a cross-server deny->allow bypass. A session allow only ever
    matches the exact registered name it was remembered under (step 1).
    """
    if tool_name in self.per_tool:
      return self.per_tool[tool_name]
    server, bare = self._server_and_bare(tool_name)
    if server is not None and (server, bare) in self.per_server_tool:
      return self.per_server_tool[(server, bare)]
    if server is not None and server in self.per_server:
      return self.per_server[server]
    return self.default

  def allow_tool_for_session(self, tool_name):
    """Remember ``allow`` for this tool for the rest of the session."""
    self.per_tool[tool_name] = "allow"

  @classmethod
  def from_server_configs(cls, configs, default="ask", tool_to_source=None):
    """Build a policy from per-server profile ``tool_policy`` dicts.

    Each config carries ``name`` and a ``tool_policy`` mapping of tool
    name -> decision, where the ``"*"`` key applies to the whole server.
    Shared by the launcher (claude_code path) and ChatWindow (API-backend
    path) so both derive identical policy from one profile.

    Parameters
    ----------
    configs : iterable
        Server config objects (``McpServerConfig`` or equivalent).
    default : str, optional
        Fallback decision; defaults to ``'ask'``.
    tool_to_source : dict, optional
        Forwarded to the constructor (resolves per-server entries).
    """
    per_server, per_server_tool = {}, {}
    for cfg in (configs or []):
      name = getattr(cfg, "name", "")
      for tool_name, decision in (getattr(cfg, "tool_policy", {}) or {}).items():
        if tool_name == "*":
          per_server[name] = decision
        else:
          # Keep per-tool policy scoped to its server so a cross-server name
          # collision (and the registry's '<server>:<tool>' rename) can't make
          # one server's decision clobber another's or be silently bypassed.
          per_server_tool[(name, tool_name)] = decision
    return cls(default=default, per_server=per_server,
               per_server_tool=per_server_tool, tool_to_source=tool_to_source)


# ---- registry --------------------------------------------------------------

@dataclass
class _ToolEntry:
  """Internal registry record pairing a tool spec with its handler."""
  spec: object                     # ToolSpec
  source: str                      # 'builtin' | 'skill' | 'mcp:<server>'
  handler: object = None           # callable; signature varies by source
  risk: str = "write"
  allow_remember: bool = True      # False -> approval card hides "Always allow"


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

  def register_builtin(self, spec, handler, risk="write", allow_remember=True):
    """Register a built-in tool under ``spec.name``.

    Built-ins take precedence: a same-named non-builtin tool (skill or MCP)
    already registered is replaced, so a same-named MCP tool can't shadow a
    trusted built-in and inherit its pre-authorization. This enforces the
    registry's built-in > skill > MCP collision order regardless of
    registration order (production registers MCP before these built-ins).

    ``allow_remember=False`` marks a tool whose approval card must not offer
    'Always allow this tool' -- a destructive tool that should be re-confirmed
    on every call rather than granted standing per-tool auto-approval.
    """
    self._add(spec.name, _ToolEntry(
      spec=spec, source="builtin", handler=handler, risk=risk,
      allow_remember=allow_remember),
      overwrite=True)

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

  def _add(self, name, entry, overwrite=False):
    """Insert an entry under ``name``.

    First-wins by default: a duplicate name is skipped. ``overwrite=True``
    (built-ins) replaces an existing entry of a DIFFERENT source so a
    trusted built-in always wins a name collision; a same-source duplicate
    still stays first-wins.
    """
    existing = self._entries.get(name)
    if existing is not None:
      if overwrite and existing.source != entry.source:
        print("tool '%s' from %s overridden by %s"
              % (name, existing.source, entry.source), file=self.log)
        self._entries[name] = entry
        return
      print("tool '%s' already registered from %s; skipping new %s"
            % (name, existing.source, entry.source), file=self.log)
      return
    self._entries[name] = entry

  # ---- queries -------------------------------------------------------------

  def specs(self):
    """Return the list of registered tool specs."""
    return [e.spec for e in self._entries.values()]

  def source_of(self, name):
    """Return a tool's source string, or ``None`` if unregistered."""
    return self._entries[name].source if name in self._entries else None

  def risk_of(self, name):
    """Return a tool's risk level, defaulting to ``write`` if unknown."""
    return self._entries[name].risk if name in self._entries else "write"

  def allow_remember_of(self, name):
    """Whether the approval card may offer 'Always allow this tool' for
    ``name``. Defaults to True for unknown names."""
    return self._entries[name].allow_remember if name in self._entries else True

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


# ---- builtin: phenix_ask_user_question -------------------------------------

# Input schema for the phenix_ask_user_question tool. Module-level so the
# Claude Code backend (which registers its own SDK copy of the tool) imports
# the SAME schema: the GUI QuestionCard renders these fields, so the two
# backends' ask-user contracts must not drift.
ASK_USER_QUESTION_SCHEMA = {
  "type": "object",
  "properties": {
    "questions": {
      "type": "array", "minItems": 1,
      "items": {
        "type": "object",
        "properties": {
          "question": {"type": "string"},
          "header": {"type": "string"},
          "multiSelect": {"type": "boolean"},
          "options": {
            "type": "array", "minItems": 2,
            "items": {
              "type": "object",
              "properties": {"label": {"type": "string"},
                             "description": {"type": "string"}},
              "required": ["label"]}}},
        "required": ["question", "options"]}}},
  "required": ["questions"]}


def register_ask_user_question(registry):
  """Register the ``phenix_ask_user_question`` builtin on ``registry``.

  The API-backend equivalent of ``ClaudeCodeAgent``'s SDK ask-user tool
  (claude_code runs its own loop and owns its own copy; this is for the
  anthropic / openai / portkey / google backends that drive
  ``AgentSession``'s tool loop). The handler emits an
  ``AskUserQuestionRequested`` and parks the ``AgentSession`` worker on its
  ``question_queue`` until the user answers, then returns the user's
  selections as the tool result the model reads on its next turn.

  Parameters
  ----------
  registry : ToolRegistry
      Registry the builtin is added to.
  """
  import json
  import uuid
  from qttbx.widgets.chat.agent.base import ToolSpec

  description = (
    "Ask the user one or more multiple-choice questions. Use this whenever "
    "you need user input to disambiguate, pick between options, or confirm a "
    "destructive choice. Returns the user's selections as JSON keyed by "
    "question text: single-select -> a string (chosen label or free-form "
    "'Other' text); multi-select -> a list of strings.")
  def _handler(name, input, cancel, session, tool_use_id):
    request_id = "q_" + uuid.uuid4().hex[:12]
    questions = (input or {}).get("questions", [])
    answers = session._await_question_answer(request_id, questions)
    # Return a plain JSON string: AgentSession._to_canonical_content_blocks
    # renders a str as a single text block, so the model reads exactly this
    # JSON. (A dict would be re-wrapped via json.dumps, which is also fine
    # but adds no value here.)
    return json.dumps(answers, indent=2, default=str)

  spec = ToolSpec(name="phenix_ask_user_question", description=description,
                  input_schema=ASK_USER_QUESTION_SCHEMA)
  registry.register_builtin(spec=spec, handler=_handler, risk="read")
