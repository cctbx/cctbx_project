"""Conversation data model.

Provider-agnostic content blocks. Each backend (``AnthropicAgent``,
``ClaudeCodeAgent``) converts to and from this shape at its API
boundary.
"""

import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone


def now():
  """UTC datetime; the canonical timestamp for messages and meta."""
  return datetime.now(timezone.utc)


def _new_id():
  return uuid.uuid4().hex


@dataclass
class ContentBlock:
  """One element of a ``Message.content`` list.

  ``type`` is one of: ``'text'``, ``'image'``, ``'thinking'``,
  ``'tool_use'``, ``'tool_result'``. The ``data`` shape is
  type-specific::

      text         {'text': str}
      image        {'attachment_sha256': str, 'mime': str,
                    'caption': str | None}
      thinking     {'text': str, 'signature': str | None}
      tool_use     {'id': str, 'name': str, 'input': dict}
      tool_result  {'tool_use_id': str, 'content': list[ContentBlock],
                    'is_error': bool}
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


@dataclass
class Attachment:
  """Content-addressed binary stored under <conv_dir>/attachments/."""
  sha256: str
  mime: str
  path: str                          # relative to attachments/, e.g. 'sha256-abc.png'


@dataclass
class SubagentRecord:
  """Persisted sub-conversation under <conv_dir>/subagents/<sub_id>.json
  (Section 11.5). Loaded on demand when the user inspects the parent's
  subagent tool-use cell."""
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
  archived: bool = False
  pinned: bool = False
  summary: str = ""
  schema_version: str = "1.0"


@dataclass
class Conversation:
  """Runtime conversation. Persisted by ConversationStorage as a directory
  containing meta.json, messages.json, attachments/, subagents/."""
  meta: ConversationMeta
  messages: list = field(default_factory=list)        # list[Message]
  attachments: dict = field(default_factory=dict)     # sha256 -> Attachment
  subagents: list = field(default_factory=list)       # list[SubagentRecord]

  @classmethod
  def new(cls, profile_name, model, title=""):
    ts = now()
    meta = ConversationMeta(
      id=_new_id(),
      title=title or "New conversation",
      profile_name=profile_name,
      model=model,
      created_at=ts,
      updated_at=ts,
    )
    return cls(meta=meta)

  def append(self, message):
    self.messages.append(message)
    self.meta.updated_at = now()
