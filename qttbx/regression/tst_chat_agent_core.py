"""Agent foundation types: the Agent ABC + capability flags + ToolSpec
(base.py), the event dataclasses streamed from stream_turn (events.py),
the cancel / error sentinels (errors.py), and the conversation data
model (conversation.py). These are the types every backend touches; the
tests are tiny, so they live together."""

from datetime import datetime

from libtbx.test_utils import Exception_expected
from libtbx.utils import format_cpu_times

from qttbx.widgets.chat.agent.base import (
  Agent, AgentCapabilities, ToolSpec)
from qttbx.widgets.chat.agent.conversation import (
  Attachment, ContentBlock, Conversation, ConversationMeta, Message,
  SubagentRecord, TokenUsage, now)
from qttbx.widgets.chat.agent.errors import (
  AgentError, AgentEvent, CancelToken, TurnCancelled)
from qttbx.widgets.chat.agent.events import (
  ImageEmitted, TextDelta, Thinking, ToolResultsBatched,
  ToolUseRequested, TurnDone)
from qttbx.widgets.chat.agent.events import TokenUsage as TokenUsageEvent


# ---- Agent ABC + ToolSpec + AgentCapabilities ----------------------------


def exercise_agent_capabilities_flags():
  caps = (AgentCapabilities.STREAMING
          | AgentCapabilities.TOOL_USE
          | AgentCapabilities.VISION_INPUT)
  assert AgentCapabilities.STREAMING in caps
  assert AgentCapabilities.TOOL_USE in caps
  assert AgentCapabilities.IMAGE_OUTPUT not in caps


def exercise_tool_spec_basic():
  spec = ToolSpec(
    name="echo",
    description="Echo back the input text.",
    input_schema={"type": "object",
                  "properties": {"text": {"type": "string"}},
                  "required": ["text"]})
  assert spec.name == "echo"
  assert "text" in spec.input_schema["properties"]


def exercise_agent_is_abstract():
  try:
    Agent()
  except TypeError:
    pass
  else:
    raise Exception_expected


def exercise_subclass_must_implement_methods():
  # Concrete subclass that only declares the attributes -- instantiation
  # should still fail because abstract methods aren't implemented.
  class IncompleteAgent(Agent):
    name = "incomplete"
    model = "x"
    capabilities = AgentCapabilities.STREAMING
  try:
    IncompleteAgent()
  except TypeError:
    pass
  else:
    raise Exception_expected


def exercise_subclass_with_methods_instantiable():
  class GoodAgent(Agent):
    name = "good"
    model = "x"
    capabilities = AgentCapabilities.STREAMING

    def stream_turn(self, conversation, tools, cancel):
      return iter([])

    def resolve_credentials(self, cli_override=None):
      return "key"

    def credentials_dialog_class(self):
      return object

  g = GoodAgent()
  assert g.name == "good"


# ---- event dataclasses ---------------------------------------------------


def exercise_text_delta():
  e = TextDelta(text="hello")
  assert e.text == "hello"


def exercise_thinking_with_signature():
  e = Thinking(text="reasoning", signature="sig-abc")
  assert e.text == "reasoning"
  assert e.signature == "sig-abc"


def exercise_thinking_default_signature():
  e = Thinking(text="reasoning")
  assert e.signature == ""


def exercise_tool_use_requested():
  e = ToolUseRequested(id="toolu_1", name="echo",
                       input={"text": "hi"})
  assert e.id == "toolu_1"
  assert e.name == "echo"
  assert e.input == {"text": "hi"}


def exercise_image_emitted():
  e = ImageEmitted(id="img_1", mime="image/png", data=b"\x89PNG",
                   caption="plot")
  assert e.mime == "image/png"
  assert e.data == b"\x89PNG"
  assert e.caption == "plot"


def exercise_image_emitted_default_caption():
  e = ImageEmitted(id="img_1", mime="image/png", data=b"")
  assert e.caption is None


def exercise_token_usage_event_defaults():
  e = TokenUsageEvent()
  assert e.input == 0
  assert e.output == 0
  assert e.cache_read == 0
  assert e.cache_creation == 0


def exercise_token_usage_event_with_values():
  e = TokenUsageEvent(
    input=100, output=50, cache_read=80, cache_creation=20)
  assert e.input == 100
  assert e.cache_creation == 20


def exercise_turn_done():
  e = TurnDone(stop_reason="end_turn")
  assert e.stop_reason == "end_turn"


def exercise_agent_error_is_agent_event():
  # AgentError lives in errors.py but inherits from AgentEvent for the
  # session loop to dispatch it uniformly.
  e = AgentError(message="x")
  assert isinstance(e, AgentEvent)


def exercise_tool_results_batched():
  e = ToolResultsBatched(blocks=[])
  assert e.blocks == []


# ---- CancelToken / AgentError / TurnCancelled ----------------------------


def exercise_cancel_token_starts_clear():
  tok = CancelToken()
  assert not tok.is_set()


def exercise_cancel_token_set_and_check():
  tok = CancelToken()
  tok.set()
  assert tok.is_set()


def exercise_agent_error_defaults_recoverable():
  e = AgentError(message="hello")
  assert e.message == "hello"
  assert e.recoverable is True


def exercise_agent_error_non_recoverable():
  e = AgentError(message="auth", recoverable=False)
  assert e.recoverable is False


def exercise_turn_cancelled_is_exception():
  assert issubclass(TurnCancelled, Exception)
  try:
    raise TurnCancelled()
  except TurnCancelled:
    pass


# ---- ContentBlock / Message / Conversation data model --------------------


def exercise_content_block_text():
  b = ContentBlock(type="text", data={"text": "hello"})
  assert b.type == "text"
  assert b.data["text"] == "hello"


def exercise_message_basic():
  ts = now()
  m = Message(role="user", content=[
    ContentBlock(type="text", data={"text": "hi"})], timestamp=ts)
  assert m.role == "user"
  assert m.timestamp == ts
  assert m.stop_reason is None
  assert m.usage is None


def exercise_message_assistant_with_usage():
  m = Message(role="assistant", content=[], timestamp=now(),
              stop_reason="end_turn",
              usage=TokenUsage(input=100, output=50))
  assert m.stop_reason == "end_turn"
  assert m.usage.input == 100


def exercise_conversation_meta_defaults():
  m = ConversationMeta(id="01HX",
                       title="Test",
                       profile_name="phenix_expert",
                       model="claude-opus-4-7",
                       created_at=now(),
                       updated_at=now())
  assert m.archived is False
  assert m.pinned is False
  assert m.summary == ""
  assert m.schema_version == "1.0"


def exercise_conversation_new():
  c = Conversation.new(profile_name="phenix_expert",
                       model="claude-opus-4-7",
                       title="New chat")
  assert c.meta.profile_name == "phenix_expert"
  assert c.meta.model == "claude-opus-4-7"
  assert c.meta.title == "New chat"
  assert isinstance(c.meta.id, str) and len(c.meta.id) > 0
  assert c.messages == []
  assert c.attachments == {}
  assert c.subagents == []


def exercise_conversation_append():
  c = Conversation.new(profile_name="p", model="m")
  m = Message(role="user", content=[
    ContentBlock(type="text", data={"text": "hi"})], timestamp=now())
  c.append(m)
  assert len(c.messages) == 1
  assert c.messages[0] is m


def exercise_attachment():
  a = Attachment(sha256="abc123", mime="image/png", path="sha256-abc123.png")
  assert a.sha256 == "abc123"
  assert a.path == "sha256-abc123.png"


def exercise_subagent_record():
  r = SubagentRecord(
    sub_id="sa_01",
    parent_conversation_id="01HX",
    parent_tool_use_id="toolu_1",
    task="monitor job",
    profile_name="phenix_expert_subagent",
    model="claude-opus-4-7",
    started_at=now(),
    finished_at=now(),
    final_text="done",
    token_usage=TokenUsage(input=10),
    messages=[])
  assert r.sub_id == "sa_01"
  assert r.token_usage.input == 10


def exercise_now_returns_datetime():
  assert isinstance(now(), datetime)


def exercise():
  # base
  exercise_agent_capabilities_flags()
  exercise_tool_spec_basic()
  exercise_agent_is_abstract()
  exercise_subclass_must_implement_methods()
  exercise_subclass_with_methods_instantiable()
  # events
  exercise_text_delta()
  exercise_thinking_with_signature()
  exercise_thinking_default_signature()
  exercise_tool_use_requested()
  exercise_image_emitted()
  exercise_image_emitted_default_caption()
  exercise_token_usage_event_defaults()
  exercise_token_usage_event_with_values()
  exercise_turn_done()
  exercise_agent_error_is_agent_event()
  exercise_tool_results_batched()
  # errors
  exercise_cancel_token_starts_clear()
  exercise_cancel_token_set_and_check()
  exercise_agent_error_defaults_recoverable()
  exercise_agent_error_non_recoverable()
  exercise_turn_cancelled_is_exception()
  # conversation
  exercise_content_block_text()
  exercise_message_basic()
  exercise_message_assistant_with_usage()
  exercise_conversation_meta_defaults()
  exercise_conversation_new()
  exercise_conversation_append()
  exercise_attachment()
  exercise_subagent_record()
  exercise_now_returns_datetime()


if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
