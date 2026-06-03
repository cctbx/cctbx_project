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
  """Provider-agnostic tool specification.

  Agents convert to and from their provider's tool shape at the API
  boundary.
  """
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
    """Return the Qt dialog class to prompt for this provider's credentials.

    Each provider gets its own dialog matching its auth model.
    """

  def submit_approval(self, response):
    """Forward a user approval decision to the agent if it owns the request.

    The default implementation returns ``False`` -- the agent didn't
    originate the request, so the runner should fall through to the
    session's approval queue.

    Backends that gate tool execution via a provider-side callback
    (e.g. the Claude Code SDK's ``can_use_tool``) override this to
    fulfill a pending future when ``response.request_id`` matches one
    they emitted. Returning ``True`` tells the runner the request was
    handled here and should NOT also land in the session queue (which
    is parked on a different tool's approval).

    Parameters
    ----------
    response : ToolApprovalResponse
        The user's decision delivered from the GUI.

    Returns
    -------
    bool
        ``True`` if the agent owned ``response.request_id``, else
        ``False``.
    """
    return False

  def submit_question_answer(self, request_id, answers):
    """Forward the user's answers to a question asked by the agent.

    The question was asked via ``AskUserQuestionRequested``. The default
    returns ``False`` -- the agent never asks structured questions.

    Backends that expose an "ask the user" MCP tool override this to
    fulfill the pending future when ``request_id`` matches one they
    emitted.

    Parameters
    ----------
    request_id : str
        The same id the agent put on the ``AskUserQuestionRequested``.
    answers : dict
        ``{question_text: selected_label}`` (single-select) or
        ``{question_text: [labels...]}`` (multi-select), plus an
        optional ``"_notes"`` sub-dict for free-form additions.

    Returns
    -------
    bool
        ``True`` if the agent owned ``request_id``, else ``False``.
    """
    return False
