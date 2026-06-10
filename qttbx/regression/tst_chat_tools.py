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


def exercise_policy_allow_server_for_session():
  p = ToolPolicy(default="ask",
                 tool_to_source={"some_phenix_tool": "mcp:phenix"})
  p.allow_server_for_session("phenix")
  assert p.resolve("some_phenix_tool") == "allow"


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
  assert reg.server_of("subagent") is None


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
    summary="echo hi",
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
  exercise_policy_allow_server_for_session()
  exercise_registry_register_builtin()
  exercise_registry_register_skill_tool()
  exercise_registry_invoke_builtin()
  exercise_registry_invoke_skill()
  exercise_register_ask_user_question_builtin()
  exercise_approval_request_response_dataclasses()
  exercise_approval_request_is_agent_event()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
