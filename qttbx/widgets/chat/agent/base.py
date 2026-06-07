"""Agent ABC and supporting types.

Provider-agnostic. The concrete backends all live in ``phenix.gui.chat``
(``AnthropicAgent``, ``ClaudeCodeAgent``, ``OpenAIAgent``, ``PortkeyAgent``,
``GeminiAgent``); each is selected by
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
        ``{question_text: label}`` (single-select; ``label`` is the
        chosen option or the user's free-form "Other" text) or
        ``{question_text: [labels...]}`` (multi-select; any "Other"
        text is appended to the list).

    Returns
    -------
    bool
        ``True`` if the agent owned ``request_id``, else ``False``.
    """
    return False

  def close(self):
    """Release any resources the agent holds (HTTP client connection
    pools, a subprocess, an asyncio loop).

    Default: best-effort close of ``self.client`` -- the HTTP/SDK client
    every API backend stores there -- via ``close_client``; a no-op for
    agents without one. Backends holding more than a client (subprocess,
    asyncio loop) override. Called once at window teardown
    (``ChatWindow.closeEvent``); it must be safe to call even if the agent
    never opened anything, and safe to call more than once.
    """
    close_client(getattr(self, "client", None))

  def _write_debug_request_record(self, label, count, tool_names):
    """One-line redacted request record to ``self._debug_log``, if set.

    Logs only the backend name, model, item count, and tool names --
    never api keys, never message content (which can carry image bytes).
    Shared by the API backends' ``_write_debug_request`` wrappers so the
    redaction policy lives in one place; must never break a turn.
    """
    log = getattr(self, "_debug_log", None)
    if log is None:
      return
    try:
      log.write("[%s] request model=%s %s=%d tools=[%s]\n" % (
        self.name, self.model, label, count, ",".join(tool_names)))
      log.flush()
    except Exception:
      pass


def close_client(client):
  """Best-effort close of an HTTP/SDK client, releasing its connection pool.

  Safe to pass ``None`` or a client without a ``close()`` method, and
  swallows any error ``close()`` raises -- teardown and key-rotation paths
  must never fail on cleanup. Shared by the HTTP-client agents' ``close()``
  and ``set_api_key`` methods.

  Parameters
  ----------
  client : object or None
      The provider client to close.
  """
  closer = getattr(client, "close", None)
  if callable(closer):
    try:
      closer()
    except Exception:
      pass
