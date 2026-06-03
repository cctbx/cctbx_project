"""Agent event types yielded by ``Agent.stream_turn``.

Consumers (``AgentSession``, ``QtAgentRunner``) dispatch on isinstance.
"""

from dataclasses import dataclass, field

# Re-export AgentEvent + AgentError from errors.py so consumers can import
# everything they need from events.py uniformly. The base class lives in
# errors.py because errors and the marker base are tightly coupled.
from qttbx.widgets.chat.agent.errors import AgentEvent, AgentError  # noqa: F401


@dataclass
class TextDelta(AgentEvent):
  text: str = ""


@dataclass
class Thinking(AgentEvent):
  """Claude extended-thinking block.

  ``signature`` is preserved on the block and sent back on subsequent
  turns (an Anthropic requirement for tool-use continuation).
  """
  text: str = ""
  signature: str = ""


@dataclass
class ToolUseRequested(AgentEvent):
  id: str = ""
  name: str = ""
  input: dict = field(default_factory=dict)


@dataclass
class ServerToolUsed(AgentEvent):
  """Claude invoked a server-side tool.

  Examples include ``web_search``, ``web_fetch``, and
  ``code_execution``. The Anthropic API runs these on the model's behalf
  — no client dispatch is needed; the corresponding ``ServerToolResult``
  arrives later in the same assistant turn. Surfaced to the UI so the
  user can see what the model ran.
  """
  id: str = ""
  name: str = ""
  input: dict = field(default_factory=dict)


@dataclass
class ServerToolResult(AgentEvent):
  """Result payload of a server-side tool call.

  ``content`` is the raw API dict (opaque at this layer); ``tool_use_id``
  pairs the result with a prior ``ServerToolUsed`` event.
  """
  tool_use_id: str = ""
  content: dict = field(default_factory=dict)


@dataclass
class ImageEmitted(AgentEvent):
  """Emitted when an image lands in the assistant stream.

  Sources include MCP tool results carrying image content and (for the
  Claude Code backend) tool-result ``UserMessage`` payloads from
  claude's built-in tools (e.g. ``Read`` on a PNG). The session stores
  the raw bytes via ``storage.store_attachment`` and appends an image
  ``ContentBlock`` referencing the resulting ``sha256``.
  """
  id: str = ""
  mime: str = ""
  data: bytes = b""
  caption: str = None


@dataclass
class TokenUsage(AgentEvent):
  input: int = 0
  output: int = 0
  cache_read: int = 0
  cache_creation: int = 0


@dataclass
class TurnDone(AgentEvent):
  stop_reason: str = ""


@dataclass
class ToolResultsBatched(AgentEvent):
  """Emitted after dispatching the tool_use blocks from one assistant turn.

  Carries the list of ``tool_result`` ``ContentBlock`` instances the
  session appended as the next user message. The UI uses this to render
  the batch as a single visual unit.
  """
  blocks: list = field(default_factory=list)


@dataclass
class AskUserQuestionRequested(AgentEvent):
  """The agent needs the user to answer one or more multiple-choice questions.

  Emitted by the backend's in-process MCP tool handler when the model
  calls the question-asking tool; the GUI renders a QuestionCard for the
  user. The answers come back via the runner's ``submit_question_answer``,
  with which the agent fulfills the matching pending future.

  ``questions`` is a list of dicts shaped like the model's tool input::

      [{
        "question": str,
        "header": str | None,
        "options": [{"label": str, "description": str | None}, ...],
        "multiSelect": bool,
      }, ...]
  """
  request_id: str = ""
  questions: list = field(default_factory=list)
