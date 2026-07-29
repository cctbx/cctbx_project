"""Agent event types yielded by ``Agent.stream_turn``.

Consumers (``AgentSession``, ``QtAgentRunner``) dispatch on isinstance.
"""

from dataclasses import dataclass, field

# Re-export AgentEvent + AgentError from errors.py so consumers can import
# everything they need from events.py uniformly. The base class lives in
# errors.py because errors and the marker base are tightly coupled.
from qttbx.widgets.chat.agent.errors import AgentEvent, AgentError  # noqa: F401


# ---- canonical TurnDone.finish vocabulary ----------------------------------
#
# The values ``TurnDone.finish`` may take, defined HERE -- beside the field
# they describe -- rather than only in ``phenix.gui.chat.finish``, which is
# where the per-backend MAPPERS onto them live. Two consumers need the
# vocabulary and only one of them may import phenix: ``AgentSession`` stamps
# CANCELLED on the terminal event it synthesizes, and ``run_child`` branches on
# the whole set to tell a cut-off child from a clean one. phenix's finish
# module re-exports these names, so the two repos share ONE definition and a
# value cannot drift between the side that writes it and the side that reads
# it. ``ERROR`` is a finish value, not the ``AgentError`` event.
CLEAN = "clean"                # complete answer
TRUNCATED = "truncated"        # output-token limit hit mid-answer
TOOL_USE = "tool_use"          # ended wanting tools; FINAL means loop capped
CANCELLED = "cancelled"        # Ctrl-C / Stop / --timeout
CAP = "cap"                    # the backend's own max-turns cap
ERROR = "error"                # backend error result


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
  """Token accounting for one assistant turn.

  Field meanings -- in particular that ``context_tokens`` is a peak, not a
  total, and must never be summed -- are documented on the stored
  ``conversation.TokenUsage``, which is the canonical definition.
  """
  input: int = 0
  output: int = 0
  cache_read: int = 0
  cache_creation: int = 0
  context_tokens: int = 0

  def to_stored(self):
    """Convert this usage event to the canonical stored ``TokenUsage``.

    ``conversation.TokenUsage`` is the single source of truth for the
    persisted usage field list. Copying field-by-field off it (rather than
    naming each field here) means a usage field added there flows through
    automatically and cannot be silently dropped when an event becomes a
    stored message.

    Returns
    -------
    qttbx.widgets.chat.agent.conversation.TokenUsage
        A stored ``TokenUsage`` carrying every field the two types share.
    """
    from dataclasses import fields
    from qttbx.widgets.chat.agent.conversation import TokenUsage
    return TokenUsage(**{f.name: getattr(self, f.name)
                         for f in fields(TokenUsage)
                         if hasattr(self, f.name)})


@dataclass
class TurnDone(AgentEvent):
  """One assistant round finished.

  ``stop_reason`` is the raw, load-bearing backend string (``AgentSession``
  loops while it stays ``"tool_use"``; the GUI/storage special-case
  ``"cancelled"``) and MUST keep its native per-backend vocabulary.

  ``finish`` is the backend-agnostic disposition, normalized AT THE AGENT
  BOUNDARY so a headless consumer never has to interpret per-backend
  stop-reason strings itself. Each backend maps its native finish state onto
  one of the canonical values defined above (and re-exported, with the
  per-backend mappers, by ``phenix.gui.chat.finish``): ``CLEAN`` /
  ``TRUNCATED`` / ``TOOL_USE`` / ``CANCELLED`` / ``CAP`` / ``ERROR``. A
  mapper that FAILS CLOSED may also put a raw provider reason here
  (``content_filter``, ``pause_turn``, ``SAFETY``), which is non-clean by
  construction. Empty string means "unset" -- a consumer that requires a
  definite disposition treats it as no information, not as success.
  """
  stop_reason: str = ""
  finish: str = ""


@dataclass
class ToolResultsBatched(AgentEvent):
  """Emitted after dispatching the tool_use blocks from one assistant turn.

  Carries the list of ``tool_result`` ``ContentBlock`` instances the
  session appended as the next user message. The UI uses this to render
  the batch as a single visual unit.
  """
  blocks: list = field(default_factory=list)


@dataclass
class ToolResultObserved(AgentEvent):
  """A tool result observed from an agent that runs its own tool loop.

  The Claude Code backend dispatches MCP tools inside the SDK subprocess, so
  the ``ToolResultsBatched`` path (which the session emits only for tools it
  dispatches itself) never fires for it. The backend instead surfaces one of
  these per in-subprocess tool result, carrying ``tool_use_id`` + ``content`` +
  ``is_error``. ``name``/``input`` are the originating tool call, resolved by
  the agent (it sees the tool_use before the result, in order) so the client
  need not correlate against the conversation -- the current turn's tool_use is
  not yet persisted when this fires. ``AgentSession`` builds the answering
  tool_result from this event AS IT ARRIVES, so the tool is never left pending
  and the observed result (with its error flag) is persisted to the transcript;
  the client also captures it into job history.
  """
  tool_use_id: str = ""
  content: object = None
  name: str = ""
  input: dict = field(default_factory=dict)
  is_error: bool = False


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
