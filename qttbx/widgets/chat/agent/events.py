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
  """Claude extended-thinking block. signature is preserved on the block
  and sent back on subsequent turns (Anthropic requirement for tool-use
  continuation)."""
  text: str = ""
  signature: str = ""


@dataclass
class ToolUseRequested(AgentEvent):
  id: str = ""
  name: str = ""
  input: dict = field(default_factory=dict)


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
  """Emitted by the session after dispatching the tool_use blocks from one
  assistant turn. Carries the list of tool_result ContentBlocks the session
  appended as the next user message. The UI uses this to render the batch
  as a single visual unit."""
  blocks: list = field(default_factory=list)
