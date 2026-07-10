from libtbx.utils import format_cpu_times, null_out
from qttbx.widgets.chat.agent.base import ToolSpec
from qttbx.widgets.chat.agent.tools import (
  ToolPolicy, ToolRegistry, ToolApprovalRequest, ToolApprovalResponse)


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
  renames the later one to '<server>:<tool>'. resolve() must apply each
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
  # Registry state after the collision: srvA keeps the bare name, srvB renamed.
  tool_to_source = {
    "phenix_get_status": "mcp:srvA",
    "srvB:phenix_get_status": "mcp:srvB",
  }
  p = ToolPolicy.from_server_configs(cfgs, tool_to_source=tool_to_source)
  assert p.resolve("phenix_get_status") == "allow", p.__dict__
  assert p.resolve("srvB:phenix_get_status") == "deny", p.__dict__
  # A server-level '*' entry still applies to a tool with no per-tool entry.
  cfgs2 = [_Cfg("srvC", {"*": "allow", "danger": "deny"})]
  src2 = {"look": "mcp:srvC", "danger": "mcp:srvC"}
  p2 = ToolPolicy.from_server_configs(cfgs2, tool_to_source=src2)
  assert p2.resolve("look") == "allow", p2.__dict__       # via '*'
  assert p2.resolve("danger") == "deny", p2.__dict__      # per-tool wins


def exercise_policy_session_allow_does_not_leak_across_collision():
  """[Major, security] A session-remembered allow on server A's bare-named
  tool must NOT leak to a DIFFERENT server B's same-named (collision-renamed
  '<B>:<tool>') tool and override B's '*':'deny'. resolve() must not fall back
  to a bare per_tool entry for a renamed tool -- that is a cross-server
  deny->allow downgrade that auto-runs a globally-denied tool with no prompt."""
  class _Cfg:
    def __init__(self, name, tp):
      self.name = name
      self.tool_policy = tp

  cfgs = [_Cfg("srvA", {}), _Cfg("srvB", {"*": "deny"})]
  src = {"phenix_get_status": "mcp:srvA",
         "srvB:phenix_get_status": "mcp:srvB"}
  p = ToolPolicy.from_server_configs(cfgs, tool_to_source=src)
  assert p.resolve("srvB:phenix_get_status") == "deny", p.__dict__
  # User remember-allows srvA's bare tool (keyed by its registered name).
  p.allow_tool_for_session("phenix_get_status")
  assert p.resolve("phenix_get_status") == "allow"           # srvA's own tool
  # srvB's collision-renamed tool must STILL be denied -- no cross-server leak.
  assert p.resolve("srvB:phenix_get_status") == "deny", p.__dict__


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


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
