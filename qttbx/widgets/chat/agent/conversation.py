"""Conversation data model.

Provider-agnostic content blocks. Each backend (``AnthropicAgent``,
``ClaudeCodeAgent``) converts to and from this shape at its API
boundary.
"""

import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone


def now():
  """Return the current UTC datetime.

  This is the canonical timestamp for messages and meta.
  """
  return datetime.now(timezone.utc)


_MIN_DT = datetime.min.replace(tzinfo=timezone.utc)


def recency_key(meta):
  """tz-aware ``updated_at`` for sorting conversations most-recent-first.

  Robust to a hand-edited/legacy meta whose ``updated_at`` is ``None`` (returns
  ``_MIN_DT``, so it sorts last) or naive (coerced to UTC) -- comparing naive
  and aware datetimes raises ``TypeError``, the same threat the null case
  guards against. Shared by the launcher's startup restore and the chat
  window's sidebar sort.
  """
  dt = meta.updated_at
  if dt is None:
    return _MIN_DT
  return dt if dt.tzinfo is not None else dt.replace(tzinfo=timezone.utc)


def _new_id():
  return uuid.uuid4().hex


# A truthy ``data`` entry under this key marks a content block as EPHEMERAL:
# it is sent to the model as part of the live conversation, but is never
# written to disk (``storage.save`` drops it via ``is_ephemeral_block``). Used
# for transient injections that must reach the model exactly for the turn they
# ride on and neither persist nor replay -- e.g. the context-pressure handoff
# note. Filtering at the serialization choke point (rather than at one caller's
# turn-end) is what makes "never on disk" hold across every save path: mid-turn
# autosave, turn-end, errored-turn, and close saves alike.
EPHEMERAL_BLOCK_KEY = "ephemeral"


def is_ephemeral_block(b):
  """True when content block ``b`` is tagged ephemeral (see
  ``EPHEMERAL_BLOCK_KEY``) and so must not be persisted."""
  return bool(getattr(b, "data", None)) and bool(b.data.get(EPHEMERAL_BLOCK_KEY))


@dataclass
class ContentBlock:
  """One element of a ``Message.content`` list.

  ``type`` is one of: ``'text'``, ``'image'``, ``'thinking'``,
  ``'tool_use'``, ``'tool_result'``, ``'server_tool_use'``,
  ``'server_tool_result'``. The ``data`` shape is type-specific::

      text                {'text': str}
      image               {'attachment_sha256': str, 'mime': str,
                           'caption': str | None}
      thinking            {'text': str, 'signature': str | None}
      tool_use            {'id': str, 'name': str, 'input': dict}
      tool_result         {'tool_use_id': str,
                           'content': list[ContentBlock],
                           'is_error': bool}
      server_tool_use     {'id': str, 'name': str, 'input': dict}
                          -- API-executed (e.g. ``web_search``); no
                          client dispatch.
      server_tool_result  {'tool_use_id': str, 'content': dict}
                          -- opaque provider payload paired with the
                          matching ``server_tool_use``.
  """
  type: str
  data: dict


@dataclass
class TokenUsage:
  """Per-message token accounting.

  The first four are COST fields: they accumulate over everything the message
  cost to produce, and summing them across messages gives a conversation's
  total. Backends that run their own tool loop report them summed over every
  API call in the turn.

  ``context_tokens`` is NOT one of them and must never be summed. It is the
  peak PROMPT size (input + cache_read + cache_creation) of a single API call
  in that turn -- i.e. how full the context window was, which the cost fields
  cannot answer because they aggregate many calls (real sessions record cost
  sums of 3M-114M against a 1M window). To read a conversation's context
  usage, take the MAXIMUM of this field across messages, not the total.
  0 means "not measured" -- an older saved conversation, or a backend that
  does not report per-call usage -- never "context empty".
  """
  input: int = 0
  output: int = 0
  cache_read: int = 0
  cache_creation: int = 0
  context_tokens: int = 0


@dataclass
class Message:
  role: str                          # 'user' or 'assistant'
  content: list                      # list[ContentBlock]
  timestamp: datetime
  stop_reason: str = None            # assistant only
  usage: TokenUsage = None           # assistant only
  model: str = None                  # assistant only: model that produced it
  backend: str = None                # assistant only: backend that produced it
  # Assistant only: whether the session's verbose mode was on when this turn
  # was produced (phenix.agent --verbose). Provenance, like model/backend --
  # never read back to configure a session. None = unstamped: a pre-stamp
  # message, or a host app without the concept.
  verbose: bool = None


@dataclass
class Attachment:
  """Content-addressed binary stored under ``<conv_dir>/attachments/``."""
  sha256: str
  mime: str
  path: str                          # relative to attachments/, e.g. 'sha256-abc.png'


@dataclass
class SubagentRecord:
  """Persisted sub-conversation under ``<conv_dir>/subagents/<sub_id>.json``.

  Loaded on demand when the user inspects the parent's subagent tool-use
  cell.

  ``incomplete_reason`` is the machine-readable answer to "did this child
  actually finish, having measured something?": ``''`` when it ran to its own
  end AND something it could measure with answered, ``'turn_cap'`` when its
  turn cap cut it off, ``'cancelled'`` when a Stop did, ``'tools_denied'`` /
  ``'no_measurement'`` when it measured nothing (refused every tool it tried,
  or held measurement tools and got a result back from none of them), and
  otherwise the canonical finish that ended it (``'error'``, ``'truncated'``,
  ``'no_terminal_event'``, ...). A child that stopped early -- or that ran to
  its own end having checked nothing -- reads exactly like one that finished
  and found nothing, so any consumer deciding what a child's output MEANS must
  branch on this field rather than on the prose in ``final_text``.
  """
  sub_id: str
  parent_conversation_id: str
  parent_tool_use_id: str
  task: str
  profile_name: str
  model: str
  started_at: datetime
  finished_at: datetime
  final_text: str
  token_usage: TokenUsage
  messages: list                     # list[Message]
  # Defaulted, hence last: every field above is default-less.
  # '' | 'turn_cap' | 'cancelled' | 'tools_denied' | 'no_measurement' | ...
  incomplete_reason: str = ""
  # {absolute path -> sha256} of the artifacts the child was spawned to
  # examine, hashed BY THIS LAYER at spawn time -- never supplied as prose by
  # the agent that benefits from the answer. It is what lets a later reader ask
  # "did this review measure the file as it now stands?" by re-hashing, instead
  # of believing a brief. A path that could not be read maps to "" (recorded as
  # unverifiable, which is not the same as unchanged).
  subject_digests: dict = field(default_factory=dict)


@dataclass
class ConversationMeta:
  id: str
  title: str
  profile_name: str
  model: str
  created_at: datetime
  updated_at: datetime
  backend: str = ""                  # backend of the most recent turn
  archived: bool = False
  pinned: bool = False
  summary: str = ""
  agent_session_id: str = None       # session-based backends (claude_code) only
  # A read-only view of something that is NOT a conversation directory --
  # today, a finished subagent's transcript surfaced in the sidebar. It is
  # synthesized from the parent's subagents/ record and owns no directory, so
  # it cannot be saved, renamed or continued, and the window shows it with no
  # input box.
  # The parent is not stored: a view's id IS "<parent>~<sub_id>", so
  # split_subagent_view_id re-derives it, and the one consumer does exactly
  # that. A second copy could only drift from the id it duplicates.
  readonly: bool = False
  schema_version: str = "1.0"


@dataclass
class Conversation:
  """Runtime conversation.

  Persisted by ``ConversationStorage`` as a directory containing
  ``meta.json``, ``messages.json``, ``attachments/``, and ``subagents/``.
  """
  meta: ConversationMeta
  messages: list = field(default_factory=list)        # list[Message]
  attachments: dict = field(default_factory=dict)     # sha256 -> Attachment
  subagents: list = field(default_factory=list)       # list[SubagentRecord]

  @classmethod
  def new(cls, profile_name, model, title="", backend=""):
    """Create an empty conversation with freshly generated meta."""
    ts = now()
    meta = ConversationMeta(
      id=_new_id(),
      title=title or "New conversation",
      profile_name=profile_name,
      model=model,
      created_at=ts,
      updated_at=ts,
      backend=backend,
    )
    return cls(meta=meta)

  def append(self, message):
    """Append a message and refresh the conversation's update timestamp."""
    self.messages.append(message)
    self.meta.updated_at = now()
