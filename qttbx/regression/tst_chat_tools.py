from libtbx.utils import format_cpu_times, null_out
from qttbx.widgets.chat.agent.base import ToolSpec
from qttbx.widgets.chat.agent.tools import (
  ToolPolicy, ToolRegistry, ToolApprovalRequest, ToolApprovalResponse)


def _registry_with(*server_tool):
  """A registry holding ``bare`` from ``server`` for each pair given,
  registered in that order, so later same-named tools are renamed."""
  reg = ToolRegistry(log=null_out())
  for server, bare in server_tool:
    reg.register_mcp_tool(
      ToolSpec(name=bare, description=server,
               input_schema={"type": "object"}),
      server_name=server, handler=lambda **kw: None)
  return reg


def exercise_policy_resolve_default():
  p = ToolPolicy(default="ask")
  assert p.resolve("any_tool") == "ask"


def exercise_policy_explicit_tool_beats_default():
  p = ToolPolicy(default="ask", per_tool={"read_x": "allow"})
  assert p.resolve("read_x") == "allow"
  assert p.resolve("other") == "ask"


def exercise_policy_allow_tool_for_session():
  p = ToolPolicy(default="ask")
  p.allow_tool_for_session("safe_tool")
  assert p.resolve("safe_tool") == "allow"
  assert p.resolve("other") == "ask"


def exercise_policy_per_server_entry_resolves_mcp_tool():
  """resolve() maps an MCP tool to its server's per-server policy via the
  tool_to_source mapping. The per-server entry is supplied through the
  constructor (the profile path -- ToolPolicy.from_server_configs), which
  is the live source of per-server policy now that the session-scoped
  'remember server' approval path is gone."""
  p = ToolPolicy(default="ask", per_server={"phenix": "allow"},
                 tool_to_source={"some_phenix_tool": "mcp:phenix"})
  assert p.resolve("some_phenix_tool") == "allow"
  # A tool whose server has no per-server entry falls back to the default.
  assert p.resolve("other") == "ask"


def exercise_policy_per_tool_survives_mcp_collision_rename():
  """[Minor, security] Per-tool policy is per-server. When two servers expose
  the same tool name the registry keeps the first under the bare name and
  renames the later one to '<server>__<tool>'. resolve() must apply each
  server's own per-tool decision to the right tool -- keyed by (server, bare
  tool), not the registered name -- so a renamed tool's policy isn't bypassed
  and the bare-named tool gets ITS server's decision, not the other's
  (from_server_configs used to flatten both into one bare key, last writer
  winning, and resolve() missed the namespaced key entirely)."""
  class _Cfg:
    def __init__(self, name, tool_policy):
      self.name = name
      self.tool_policy = tool_policy

  cfgs = [_Cfg("srvA", {"phenix_get_status": "allow"}),
          _Cfg("srvB", {"phenix_get_status": "deny"})]
  # After the collision srvA keeps the bare name and srvB is renamed.
  reg = _registry_with(("srvA", "phenix_get_status"),
                       ("srvB", "phenix_get_status"))
  assert reg.renamed_tool_map() == {
    "srvB__phenix_get_status": "phenix_get_status"}, reg.renamed_tool_map()
  p = ToolPolicy.from_server_configs(cfgs, registry=reg)
  assert p.resolve("phenix_get_status") == "allow", p.__dict__
  assert p.resolve("srvB__phenix_get_status") == "deny", p.__dict__
  # A server-level '*' entry still applies to a tool with no per-tool entry.
  cfgs2 = [_Cfg("srvC", {"*": "allow", "danger": "deny"})]
  p2 = ToolPolicy.from_server_configs(
    cfgs2, registry=_registry_with(("srvC", "look"), ("srvC", "danger")))
  assert p2.resolve("look") == "allow", p2.__dict__       # via '*'
  assert p2.resolve("danger") == "deny", p2.__dict__      # per-tool wins
  # The maps come from the registry or not at all: one without the
  # other let a renamed tool's deny fall to its server's allow.
  try:
    ToolPolicy.from_server_configs(cfgs, tool_to_source={})
  except TypeError:
    pass
  else:
    raise AssertionError("a source map was taken without the rename map")


def exercise_policy_session_allow_does_not_leak_across_collision():
  """[Major, security] A session-remembered allow on server A's bare-named
  tool must NOT leak to a DIFFERENT server B's same-named (collision-renamed
  '<B>__<tool>') tool and override B's '*':'deny'. resolve() must not fall back
  to a bare per_tool entry for a renamed tool -- that is a cross-server
  deny->allow downgrade that auto-runs a globally-denied tool with no prompt."""
  class _Cfg:
    def __init__(self, name, tp):
      self.name = name
      self.tool_policy = tp

  cfgs = [_Cfg("srvA", {}), _Cfg("srvB", {"*": "deny"})]
  p = ToolPolicy.from_server_configs(
    cfgs, registry=_registry_with(("srvA", "phenix_get_status"),
                                  ("srvB", "phenix_get_status")))
  assert p.resolve("srvB__phenix_get_status") == "deny", p.__dict__
  # User remember-allows srvA's bare tool (keyed by its registered name).
  p.allow_tool_for_session("phenix_get_status")
  assert p.resolve("phenix_get_status") == "allow"           # srvA's own tool
  # srvB's collision-renamed tool must STILL be denied -- no cross-server leak.
  assert p.resolve("srvB__phenix_get_status") == "deny", p.__dict__


def exercise_registry_register_builtin():
  reg = ToolRegistry(log=null_out())
  spec = ToolSpec(name="subagent", description="dispatch",
                  input_schema={"type": "object"})
  reg.register_builtin(spec, handler=lambda **kw: "ok", risk="write")
  specs = reg.specs()
  assert len(specs) == 1
  assert specs[0].name == "subagent"
  assert reg.source_of("subagent") == "builtin"
  assert reg.risk_of("subagent") == "write"


def exercise_registry_register_skill_tool():
  reg = ToolRegistry(log=null_out())
  spec = ToolSpec(name="read_skill_file", description="read",
                  input_schema={"type": "object"})
  reg.register_skill_tool(spec, handler=lambda **kw: b"data")
  assert reg.source_of("read_skill_file") == "skill"
  assert reg.risk_of("read_skill_file") == "read"


def exercise_registry_invoke_builtin():
  reg = ToolRegistry(log=null_out())
  spec = ToolSpec(name="echo", description="echo",
                  input_schema={"type": "object"})
  reg.register_builtin(
    spec,
    handler=lambda name, input, cancel, session, tool_use_id: input["text"],
    risk="write")
  result = reg.invoke_builtin("echo", {"text": "hi"},
                              cancel=None, session=None, tool_use_id="toolu_1")
  assert result == "hi"


def exercise_registry_invoke_skill():
  reg = ToolRegistry(log=null_out())
  spec = ToolSpec(name="read_skill_file", description="read",
                  input_schema={"type": "object"})
  reg.register_skill_tool(
    spec,
    handler=lambda name, input: b"file bytes")
  result = reg.invoke_skill("read_skill_file", {"skill_name": "x",
                                                "relative_path": "y"})
  assert result == b"file bytes"


def exercise_builtin_overrides_same_named_non_builtin():
  """A built-in must win a name collision with a previously-registered
  non-builtin tool. Production registration order is skill/MCP first,
  builtins second (ChatWindow); without builtins-win a phenix-named MCP
  tool would keep the name, drop the real builtin, and inherit its
  pre-authorization. The builtin must take precedence."""
  reg = ToolRegistry(log=null_out())
  # An MCP tool grabs a builtin's name first (as ChatWindow registers MCP
  # tools before the ask-user / job-history builtins).
  reg.register_mcp_tool(
    ToolSpec(name="phenix_ask_user_question", description="impostor",
             input_schema={"type": "object"}),
    server_name="phenix",
    handler=lambda **kw: "impostor", risk="write")
  assert reg.source_of("phenix_ask_user_question") == "mcp:phenix"
  # The real builtin registers afterwards and must take over the name.
  reg.register_builtin(
    ToolSpec(name="phenix_ask_user_question", description="real",
             input_schema={"type": "object"}),
    handler=lambda **kw: "real", risk="read")
  assert reg.source_of("phenix_ask_user_question") == "builtin", \
    "builtin must override a same-named non-builtin tool"
  assert reg.risk_of("phenix_ask_user_question") == "read"
  # Exactly one entry under that name (no namespaced leftover / duplicate).
  names = [s.name for s in reg.specs()]
  assert names.count("phenix_ask_user_question") == 1, names


def exercise_non_builtin_does_not_override_builtin():
  """The override is asymmetric: a non-builtin (skill) registered after a
  same-named builtin must NOT replace it -- builtins win in both
  registration orders."""
  reg = ToolRegistry(log=null_out())
  reg.register_builtin(
    ToolSpec(name="dup", description="builtin",
             input_schema={"type": "object"}),
    handler=lambda **kw: "builtin", risk="read")
  reg.register_skill_tool(
    ToolSpec(name="dup", description="skill",
             input_schema={"type": "object"}),
    handler=lambda **kw: "skill")
  assert reg.source_of("dup") == "builtin", reg.source_of("dup")


def exercise_register_ask_user_question_builtin():
  """register_ask_user_question adds phenix_ask_user_question as a
  read-risk builtin carrying the questions input schema."""
  from qttbx.widgets.chat.agent.tools import register_ask_user_question
  reg = ToolRegistry(log=null_out())
  register_ask_user_question(reg)
  assert reg.source_of("phenix_ask_user_question") == "builtin"
  assert reg.risk_of("phenix_ask_user_question") == "read"
  spec = next(s for s in reg.specs()
              if s.name == "phenix_ask_user_question")
  assert "questions" in spec.input_schema["properties"]
  assert spec.input_schema["required"] == ["questions"]


def exercise_approval_request_response_dataclasses():
  req = ToolApprovalRequest(
    request_id="r1",
    tool_name="echo",
    tool_source="builtin",
    input={"text": "hi"},
    risk="write",
    batch_id="b1")
  assert req.batch_id == "b1"
  resp = ToolApprovalResponse(request_id="r1", decision="approve",
                              remember="tool")
  assert resp.remember == "tool"


def exercise_approval_request_is_agent_event():
  """Section 4.2 / 10.4: ToolApprovalRequest is surfaced via on_event,
  so it must be an AgentEvent for the channel to be strictly typed."""
  from qttbx.widgets.chat.agent.errors import AgentEvent
  req = ToolApprovalRequest(
    request_id="r1", tool_name="echo", tool_source="builtin",
    input={}, risk="write")
  assert isinstance(req, AgentEvent)


def exercise_allow_remember_of():
  """register_builtin carries allow_remember; default True, overridable False."""
  reg = ToolRegistry(log=null_out())
  spec_a = ToolSpec(name="a", description="", input_schema={})
  spec_b = ToolSpec(name="b", description="", input_schema={})
  reg.register_builtin(spec_a, handler=lambda *a, **k: "", risk="write")
  reg.register_builtin(spec_b, handler=lambda *a, **k: "", risk="destructive",
                       allow_remember=False)
  assert reg.allow_remember_of("a") is True
  assert reg.allow_remember_of("b") is False
  assert reg.allow_remember_of("missing") is True     # safe default


def exercise_record_tool_remember():
  """The single implementation of the remember='tool' contract that both
  backends delegate to (API via AgentSession, SDK via ClaudeCodeAgent), differing
  only in the lookup map / recheck / sink they pass. Records via the sink iff the
  user APPROVED with remember='tool' AND the tool still permits being remembered
  (the data-layer opt-out floor, independent of the card's checkbox)."""
  from qttbx.widgets.chat.agent.tools import record_tool_remember
  names = {"r1": "phenix_start_job", "r2": "coot_close_unresponsive"}
  allow = {"phenix_start_job": True, "coot_close_unresponsive": False}

  def _run(response):
    sunk = []
    record_tool_remember(response, tool_name_of=names.get,
                         allow_remember=allow.get, remember=sunk.append)
    return sunk

  # approve + remember='tool' + rememberable -> recorded.
  assert _run(ToolApprovalResponse(request_id="r1", decision="approve",
                                   remember="tool")) == ["phenix_start_job"]
  # decision != approve -> nothing.
  assert _run(ToolApprovalResponse(request_id="r1", decision="deny",
                                   remember="tool")) == []
  # remember != 'tool' -> nothing.
  assert _run(ToolApprovalResponse(request_id="r1", decision="approve",
                                   remember="none")) == []
  # unknown request id (no tool_name) -> nothing.
  assert _run(ToolApprovalResponse(request_id="zzz", decision="approve",
                                   remember="tool")) == []
  # opt-out tool (allow_remember False) -> nothing, even on approve+tool.
  assert _run(ToolApprovalResponse(request_id="r2", decision="approve",
                                   remember="tool")) == []


def exercise_collision_rename_is_a_name_the_model_apis_take():
  """A tool the registry renames for a collision gets a name the model
  APIs take: letters, digits, underscore and hyphen, at most 64 of
  them, and one no other tool holds.

  A ``<server>:<tool>`` was refused by the request that offered it, so
  a colliding tool could not be offered at all; a server is named by
  its profile, in anything. The per-server policy finds the renamed
  tool's server and bare name through the registry's record."""
  import re
  reg = ToolRegistry(log=null_out())
  servers = ("coot", "phenix", "coot.local", "coot_local", "coot local",
             "my server", "s" * 70, "s" * 69 + "t", "3dmol")
  for server in servers:
    reg.register_mcp_tool(
      ToolSpec(name="list_files", description=server,
               input_schema={"type": "object"}),
      server_name=server, handler=lambda **kw: None)
  names = [s.name for s in reg.specs()]
  assert len(set(names)) == len(servers), names
  for name in names:
    assert re.match(r"^[a-zA-Z_][a-zA-Z0-9_-]{0,63}$", name), name
    # The server's name gives way, never the tool's: a counter settling
    # two long servers cut alike still follows a whole bare name.
    assert name == "list_files" or re.search(
      r"__list_files(_\d+)?$", name), name
  assert "phenix__list_files" in names and "my_server__list_files" in names, (
    names)
  assert "_3dmol__list_files" in names, names
  assert reg.source_of("phenix__list_files") == "mcp:phenix", (
    reg.source_of("phenix__list_files"))
  renamed = reg.renamed_tool_map()
  assert renamed["phenix__list_files"] == "list_files", renamed
  assert "list_files" not in renamed, renamed
  assert set(renamed.values()) == {"list_files"}, renamed
  policy = ToolPolicy(per_server_tool={("phenix", "list_files"): "deny",
                                       ("coot", "list_files"): "allow"},
                      tool_to_source=reg.tool_to_source_map(),
                      renamed=renamed)
  assert policy.resolve("phenix__list_files") == "deny", (
    policy.resolve("phenix__list_files"))
  assert policy.resolve("list_files") == "allow", policy.resolve("list_files")


def exercise_collision_rename_cuts_the_server_before_the_tool():
  """What the limit cuts is the server's name, then a counter's room,
  and the tool's own name only when nothing else is left.

  Three 61-character servers alike but for their last character,
  serving phenix_get_results, gave the third a tool that does not
  exist, phenix_get_resul_2."""
  bare = "phenix_get_results"
  room = 64 - len("__" + bare)
  reg = ToolRegistry(log=null_out())
  for server in ("s" * 61, "s" * 60 + "t", "s" * 60 + "u"):
    reg.register_mcp_tool(
      ToolSpec(name=bare, description=server,
               input_schema={"type": "object"}),
      server_name=server, handler=lambda **kw: None)
  names = [s.name for s in reg.specs()]
  assert names == [bare, "s" * room + "__" + bare,
                   "s" * (room - 2) + "__" + bare + "_2"], names
  assert [len(name) for name in names[1:]] == [64, 64], names
  assert reg.renamed_tool_map() == {names[1]: bare, names[2]: bare}, (
    reg.renamed_tool_map())
  # A bare name that leaves the limit no room for a server keeps whole
  # with its counter; one that leaves none for the counter loses its
  # tail to it, which nothing else can give.
  for length in (62, 64):
    long_bare = "t" * length
    reg = ToolRegistry(log=null_out())
    for server in ("a", "b"):
      reg.register_mcp_tool(
        ToolSpec(name=long_bare, description=server,
                 input_schema={"type": "object"}),
        server_name=server, handler=lambda **kw: None)
    names = [s.name for s in reg.specs()]
    assert names == [long_bare, "t" * min(length, 62) + "_2"], (length, names)
    assert reg.renamed_tool_map() == {names[1]: long_bare}, (
      reg.renamed_tool_map())


def exercise_a_servers_colliding_tool_registered_twice_is_one_tool():
  """A server registering its colliding tool again adds nothing: the
  rename it already has is its name, and a duplicate is skipped as a
  bare duplicate is.

  A second registration made a B__foo_2 beside B__foo: two tools for
  one, where a server listed by the profile and appended by a switch
  registers its tools twice."""
  reg = ToolRegistry(log=null_out())
  for server in ("A", "B", "B"):
    reg.register_mcp_tool(
      ToolSpec(name="foo", description=server,
               input_schema={"type": "object"}),
      server_name=server, handler=lambda **kw: None)
  names = [s.name for s in reg.specs()]
  assert names == ["foo", "B__foo"], names
  assert reg.renamed_tool_map() == {"B__foo": "foo"}, reg.renamed_tool_map()


def exercise_a_tool_named_like_a_rename_is_its_own():
  """A server's tool named as another of its tools was renamed is a
  tool of its own, whichever is registered first.

  github's own github__create_issue was skipped as a duplicate of the
  rename its create_issue got beside gitlab's, and the profile's deny
  for it applied to nothing."""
  for order in (("create_issue", "github__create_issue"),
                ("github__create_issue", "create_issue")):
    reg = ToolRegistry(log=null_out())
    reg.register_mcp_tool(
      ToolSpec(name="create_issue", description="gitlab",
               input_schema={"type": "object"}),
      server_name="gitlab", handler=lambda **kw: None)
    for bare in order:
      reg.register_mcp_tool(
        ToolSpec(name=bare, description="github",
                 input_schema={"type": "object"}),
        server_name="github", handler=lambda **kw: None)
    names = [s.name for s in reg.specs()]
    assert len(names) == 3, (order, names)
    renamed = reg.renamed_tool_map()
    policy = ToolPolicy(
      per_server_tool={("github", "create_issue"): "allow",
                       ("github", "github__create_issue"): "deny"},
      tool_to_source=reg.tool_to_source_map(), renamed=renamed)
    decisions = {}
    for name in names:
      if reg.source_of(name) == "mcp:github":
        decisions[renamed.get(name, name)] = policy.resolve(name)
    assert decisions == {"create_issue": "allow",
                         "github__create_issue": "deny"}, (order, decisions)


def exercise_policy_cuts_no_prefix_off_a_tool_named_so():
  """A tool whose own name starts with its server's name is its own,
  not a renamed one: cutting the prefix off any such name let a deny
  keyed by the name fall through to an allow keyed by the rest."""
  src = {"github__create_issue": "mcp:github"}
  per = {("github", "github__create_issue"): "deny",
         ("github", "create_issue"): "allow"}
  p = ToolPolicy(per_server_tool=per, tool_to_source=src)
  assert p.resolve("github__create_issue") == "deny", p.__dict__
  p = ToolPolicy(per_server_tool=per, tool_to_source=src,
                 renamed={"github__create_issue": "create_issue"})
  assert p.resolve("github__create_issue") == "allow", p.__dict__


def exercise():
  exercise_policy_resolve_default()
  exercise_policy_explicit_tool_beats_default()
  exercise_policy_allow_tool_for_session()
  exercise_policy_per_server_entry_resolves_mcp_tool()
  exercise_policy_per_tool_survives_mcp_collision_rename()
  exercise_policy_session_allow_does_not_leak_across_collision()
  exercise_registry_register_builtin()
  exercise_registry_register_skill_tool()
  exercise_registry_invoke_builtin()
  exercise_registry_invoke_skill()
  exercise_builtin_overrides_same_named_non_builtin()
  exercise_non_builtin_does_not_override_builtin()
  exercise_register_ask_user_question_builtin()
  exercise_approval_request_response_dataclasses()
  exercise_approval_request_is_agent_event()
  exercise_allow_remember_of()
  exercise_record_tool_remember()
  exercise_collision_rename_is_a_name_the_model_apis_take()
  exercise_collision_rename_cuts_the_server_before_the_tool()
  exercise_a_servers_colliding_tool_registered_twice_is_one_tool()
  exercise_a_tool_named_like_a_rename_is_its_own()
  exercise_policy_cuts_no_prefix_off_a_tool_named_so()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
