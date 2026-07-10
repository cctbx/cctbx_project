"""Agent ABC and supporting types.

Provider-agnostic. The concrete backends all live in ``phenix.gui.chat``
(``AnthropicAgent``, ``ClaudeCodeAgent``, ``OpenAIAgent``, ``PortkeyAgent``,
``GeminiAgent``); each is selected by
``phenix.gui.chat.agent_factory.build_agent``.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass


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

  def resolve_api_key(self, api_key, missing_message):
    """Return an explicit ``api_key`` else the resolved credential.

    The shared constructor credential guard for the API backends (Anthropic /
    OpenAI-compatible / Gemini), each of which previously inlined the identical
    ``api_key if api_key is not None else self.resolve_credentials()`` lookup
    and ``Sorry`` raise. Raises ``Sorry(missing_message)`` when neither yields a
    usable key.

    Parameters
    ----------
    api_key : object or None
        An explicit key passed to the constructor (wins when not ``None``).
    missing_message : str
        The provider-specific message for the ``Sorry`` raised when no key is
        available.
    """
    key = api_key if api_key is not None else self.resolve_credentials()
    if not key:
      from libtbx.utils import Sorry
      raise Sorry(missing_message)
    return key

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


def flatten_tool_result_text(blocks):
  """Flatten a canonical ``tool_result`` block list into a single text string.

  ``text`` blocks contribute their text verbatim; ``image`` blocks contribute
  a placeholder (the OpenAI and Gemini tool_result wire formats are text-only,
  so an image can't ride inside one). Shared by the OpenAI-compatible and
  Gemini backends, whose tool_result flattening is byte-identical. Returns
  ``""`` for an empty or imageless list.

  Parameters
  ----------
  blocks : iterable
      Canonical content blocks (each with ``.type`` and ``.data``).

  Returns
  -------
  str
      The flattened text, newline-joined.
  """
  chunks = []
  for b in blocks:
    if b.type == "text":
      chunks.append(b.data.get("text", ""))
    elif b.type == "image":
      chunks.append("[image omitted from tool result]")
  return "\n".join(chunks)


def anthropic_image_block(mime, data):
  """Build an Anthropic-API base64 image content block from raw bytes.

  Returns the ``{"type": "image", "source": {"type": "base64", ...}}`` dict
  the Anthropic messages API expects. Shared by the Anthropic backend (its
  canonical->wire image conversion) and the Claude Code backend (its outbound
  image envelope), which both speak the Anthropic wire image shape.

  Parameters
  ----------
  mime : str
      The image MIME type (e.g. ``"image/png"``) -> ``source.media_type``.
  data : bytes
      The raw image bytes; base64-encoded into ``source.data``.

  Returns
  -------
  dict
      The Anthropic image content block.
  """
  import base64
  return {
    "type": "image",
    "source": {
      "type": "base64",
      "media_type": mime,
      "data": base64.b64encode(data).decode("ascii"),
    },
  }


# Shared exception-text classifier signals. The Gemini and Claude Code backends
# both sniff a code-less provider/SDK error's message for auth vs. rate/quota
# signals; their phrase lists had DRIFTED, so they live here ONCE (the union of
# what each had) to stay in lockstep. Match whole signal phrases, NOT bare
# "auth"/"rate" substrings -- a plain `"rate" in text` false-matches the "rate"
# inside "generate"/"moderate" and `"auth" in text` matches "author".
_AUTH_SIGNALS = ("authentication", "authenticated", "unauthorized",
                 "authorization", "permission", "api key", "api_key",
                 "oauth", "login")
_RATE_SIGNALS = ("rate limit", "rate-limit", "rate_limit", "ratelimit",
                 "429", "quota", "resource_exhausted", "resource exhausted")


def classify_error_text(text):
  """Classify a provider/SDK error message as ``"auth"``, ``"rate_limited"``,
  or ``None`` (no recognized signal).

  ``text`` is matched case-insensitively against whole signal phrases. Rate is
  checked first: a subscription/quota rate-limit is recoverable (the user just
  waits), so a message carrying both signals favors the retryable read. Shared
  by the Gemini and Claude Code backends so their phrase lists can't drift apart
  again; each caller maps the returned category to its own ``AgentError`` kind +
  recoverable policy (Gemini: ``"auth"`` / ``"rate_limited"``; Claude Code:
  ``"auth_failed"`` / ``"rate_limited"`` / ``"sdk_error"``).

  Parameters
  ----------
  text : str
      The raw exception text (need not be pre-lowercased).

  Returns
  -------
  str or None
      ``"rate_limited"`` / ``"auth"`` when a signal matches, else ``None``.
  """
  t = (text or "").lower()
  if any(s in t for s in _RATE_SIGNALS):
    return "rate_limited"
  if any(s in t for s in _AUTH_SIGNALS):
    return "auth"
  return None


def load_image_attachment(storage, conv_id, block, error_cls):
  """Load + validate a canonical ``image`` block's bytes from a
  ``ConversationStorage``.

  Shared by the OpenAI-compatible and Gemini backends, whose image
  load+validate step was identical apart from the exception type and what they
  do with the bytes afterward. Returns ``(mime, data)``; raises ``error_cls``
  when there is no storage, or the attachment is missing/empty, or its MIME is
  not an image. The caller turns the bytes into its provider's image shape (a
  base64 data URI for OpenAI; raw bytes for Gemini).

  Parameters
  ----------
  storage : ConversationStorage or None
      The store to load the attachment from; ``None`` means image blocks can't
      be sent.
  conv_id : str or None
      Conversation id the attachment belongs to.
  block : object
      The canonical ``image`` block (``block.data['attachment_sha256']`` and
      ``block.data.get('mime')``).
  error_cls : type
      Exception class raised on any failure; each backend passes its own
      (OpenAI's ``_Unsendable`` marker, Gemini's ``Sorry``).

  Returns
  -------
  tuple
      ``(mime, data)`` -- the MIME string and the raw image bytes.
  """
  if storage is None:
    raise error_cls("cannot send image blocks without a ConversationStorage")
  data = storage.load_attachment(conv_id, block.data["attachment_sha256"])
  mime = block.data.get("mime", "")
  if not (mime.startswith("image/") and data):
    raise error_cls("invalid or empty image attachment (mime=%r)" % mime)
  return mime, data


def drop_null_tool_args(args):
  """Drop null-valued keys from a tool-call argument dict.

  Both the OpenAI-compatible and Gemini models emit an explicit ``null`` for an
  optional parameter they chose to skip, and MCP servers reject ``null`` against
  a bare-type annotation -- ``null`` means 'omitted' at this boundary. A
  non-dict ``args`` is returned unchanged (only dicts are stripped), preserving
  the OpenAI path where a non-object tool argument passes through verbatim.
  Shared so both backends strip nulls identically.

  Parameters
  ----------
  args : object
      The parsed tool-call arguments (typically a dict).

  Returns
  -------
  object
      ``args`` with null-valued keys removed when it is a dict, else ``args``
      unchanged.
  """
  if not isinstance(args, dict):
    return args
  return {k: v for k, v in args.items() if v is not None}


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
