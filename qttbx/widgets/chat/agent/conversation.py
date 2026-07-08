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
  input: int = 0
  output: int = 0
  cache_read: int = 0
  cache_creation: int = 0


@dataclass
class Message:
  role: str                          # 'user' or 'assistant'
  content: list                      # list[ContentBlock]
  timestamp: datetime
  stop_reason: str = None            # assistant only
  usage: TokenUsage = None           # assistant only
  model: str = None                  # assistant only: model that produced it
  backend: str = None                # assistant only: backend that produced it


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
