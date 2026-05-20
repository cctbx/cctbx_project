"""Agent ABC and supporting types.

Provider-agnostic. ``phenix.gui.chat.anthropic_agent.AnthropicAgent``
and ``phenix.gui.chat.claude_code_agent.ClaudeCodeAgent`` are the two
concrete backends; both are selected by
``phenix.gui.chat.agent_factory.build_agent``.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import IntFlag, auto


class AgentCapabilities(IntFlag):
  STREAMING = auto()
  VISION_INPUT = auto()
  IMAGE_OUTPUT = auto()
  TOOL_USE = auto()
  EXTENDED_THINKING = auto()


@dataclass
class ToolSpec:
  """Provider-agnostic tool specification. Agents convert to/from their
  provider's tool shape at the API boundary."""
  name: str
  description: str
  input_schema: dict           # JSON Schema


class Agent(ABC):
  """Provider-agnostic chat agent.

  Implementations live in ``phenix.gui.chat``. ``AgentSession``
  (``session.py``) wraps any ``Agent`` to orchestrate the multi-turn
  tool dispatch loop.
  """

  # Subclasses set these as class attributes (or instance attributes in
  # __init__).
  name: str
  model: str
  capabilities: AgentCapabilities

  @abstractmethod
  def stream_turn(self, conversation, tools, cancel):
    """Yield ``AgentEvent`` instances as the agent processes one API call.

    Tool execution is performed by ``AgentSession``, not by the agent.
    The agent yields ``ToolUseRequested`` events; the session collects
    them after ``TurnDone`` and runs the tool dispatch loop separately.

    Parameters
    ----------
    conversation : Conversation
        Full conversation so far. The latest message is the user input
        for this turn.
    tools : list of ToolSpec
        Tools the registry has surfaced for this turn.
    cancel : CancelToken
        Polled between events; honor it to abort cleanly.
    """

  @abstractmethod
  def resolve_credentials(self, cli_override=None):
    """Resolve credentials for this provider.

    Parameters
    ----------
    cli_override : object, optional
        A credential value the launcher accepted on the command line.
        Implementations may use it verbatim or normalize first.

    Returns
    -------
    object or None
        A value (string or dict) suitable for the provider's SDK
        constructor, or ``None`` when user input is still needed.
    """

  @abstractmethod
  def credentials_dialog_class(self):
    """Return the Qt dialog class to prompt for this provider's
    credentials. Each provider gets its own dialog matching its auth
    model."""
