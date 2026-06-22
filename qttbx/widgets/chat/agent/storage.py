"""Per-project conversation storage.

Layout::

    <project_dir>/.phenix_chat/
      index.json
      conversations/<conv_id>/
        meta.json
        messages.json
        attachments/sha256-<hex>.<ext>
        subagents/<sub_id>.json

Atomic writes via tmp + rename. Lazy directory creation. Content-addressed
attachments (sha256-keyed; automatic dedup within a conversation).
"""

import hashlib
import json
import mimetypes
import os
import sys
from dataclasses import asdict
from datetime import datetime
from pathlib import Path

from qttbx.widgets.chat.agent.conversation import (
  Attachment, ContentBlock, Conversation, ConversationMeta, Message,
  SubagentRecord, TokenUsage)
from qttbx.widgets.chat.agent.paths import chat_root_for


_SCHEMA_VERSION = "1.0"


class ConversationStorage:
  """Read and write conversations under a project's ``.phenix_chat/`` dir.

  Construction does not touch the filesystem. Directories appear lazily on
  first write.

  Parameters
  ----------
  project_dir : str or pathlib.Path
      Project directory whose chat root holds the conversations.
  log : file-like, optional
      Stream for diagnostic messages. Defaults to ``sys.stdout``.
  """

  def __init__(self, project_dir, log=None):
    self.project_dir = Path(project_dir)
    self.root = chat_root_for(self.project_dir)
    self.log = log if log is not None else sys.stdout

  # ---- conversations -------------------------------------------------------

  def list_conversations(self):
    """Return conversation metadata, rebuilding the index if needed.

    Reads the cached ``index.json``, or rebuilds it from the on-disk
    conversations when the index is missing or corrupt.

    Returns
    -------
    list of ConversationMeta
        Metadata for the stored conversations.
    """
    index_path = self.root / "index.json"
    if index_path.exists():
      try:
        with open(index_path) as fh:
          data = json.load(fh)
        return [_meta_from_dict(d) for d in data.get("conversations", [])]
      except Exception as e:
        print("storage: index.json unreadable (%s); rebuilding"
              % e, file=self.log)
    return self._rebuild_index()

  def load(self, conv_id):
    """Load a conversation's meta and messages.

    Attachments and subagents are loaded on demand, not here.

    Parameters
    ----------
    conv_id : str
        Identifier of the conversation to load.

    Returns
    -------
    Conversation
        The conversation with ``meta`` and ``messages`` populated and
        empty ``attachments`` / ``subagents``.

    Raises
    ------
    libtbx.utils.Sorry
        If the conversation directory does not exist, or a document's
        ``schema_version`` is unsupported.
    """
    conv_dir = self._conv_dir(conv_id)
    if not conv_dir.exists():
      from libtbx.utils import Sorry
      raise Sorry("Conversation not found: %s" % conv_id)
    meta_doc = _read_json(conv_dir / "meta.json")
    _check_schema_version(meta_doc, str(conv_dir / "meta.json"))
    meta = _meta_from_dict(meta_doc)
    messages_doc = _read_json(conv_dir / "messages.json")
    _check_schema_version(messages_doc, str(conv_dir / "messages.json"))
    messages = [_message_from_dict(m) for m in messages_doc.get("messages", [])]
    # Attachments and subagents loaded on demand.
    return Conversation(meta=meta, messages=messages,
                        attachments={}, subagents=[])

  def save(self, conv):
    """Write ``meta.json`` and ``messages.json`` atomically, then reindex.

    Parameters
    ----------
    conv : Conversation
        The conversation to persist.
    """
    self._ensure_root()
    conv_dir = self._conv_dir(conv.meta.id)
    conv_dir.mkdir(parents=True, exist_ok=True)
    _atomic_write_json(conv_dir / "meta.json", _meta_to_dict(conv.meta))
    _atomic_write_json(conv_dir / "messages.json", {
      "schema_version": _SCHEMA_VERSION,
      "messages": [_message_to_dict(m) for m in conv.messages],
    })
    self._refresh_index()

  # ---- attachments ---------------------------------------------------------

  def store_attachment(self, conv_id, data, mime):
    """Store attachment bytes content-addressed by their sha256.

    Idempotent: writing the same bytes twice produces one file.

    Parameters
    ----------
    conv_id : str
        Conversation the attachment belongs to.
    data : bytes
        The attachment bytes.
    mime : str
        MIME type, used to choose the file extension.

    Returns
    -------
    Attachment
        Reference to the stored bytes, keyed by their sha256.
    """
    self._ensure_root()
    sha = hashlib.sha256(data).hexdigest()
    ext = mimetypes.guess_extension(mime) or ".bin"
    fname = "sha256-%s%s" % (sha, ext)
    att_dir = self._conv_dir(conv_id) / "attachments"
    att_dir.mkdir(parents=True, exist_ok=True)
    fpath = att_dir / fname
    if not fpath.exists():
      tmp = fpath.with_suffix(fpath.suffix + ".tmp")
      with open(tmp, "wb") as fh:
        fh.write(data)
      os.replace(tmp, fpath)
    return Attachment(sha256=sha, mime=mime, path=fname)

  def load_attachment(self, conv_id, sha256):
    """Load attachment bytes by their sha256.

    Parameters
    ----------
    conv_id : str
        Conversation the attachment belongs to.
    sha256 : str
        Hex sha256 of the attachment bytes.

    Returns
    -------
    bytes
        The stored attachment bytes.

    Raises
    ------
    libtbx.utils.Sorry
        If ``sha256`` is unsafe as a path segment, or no matching
        attachment exists.
    """
    import glob as _glob
    sha256 = _safe_segment(sha256, "attachment")
    att_dir = self._conv_dir(conv_id) / "attachments"
    # _glob.escape neutralizes any glob metacharacters in the (validated)
    # sha so it is matched literally, not as a pattern.
    matches = sorted(att_dir.glob("sha256-%s.*" % _glob.escape(sha256)))
    if not matches:
      from libtbx.utils import Sorry
      raise Sorry("Attachment not found: sha256=%s" % sha256)
    with open(matches[0], "rb") as fh:
      return fh.read()

  # ---- subagents -----------------------------------------------------------

  def store_subagent(self, conv_id, record):
    """Write a subagent record atomically under the conversation."""
    self._ensure_root()
    sub_dir = self._conv_dir(conv_id) / "subagents"
    sub_dir.mkdir(parents=True, exist_ok=True)
    path = sub_dir / ("%s.json" % _safe_segment(record.sub_id, "subagent"))
    _atomic_write_json(path, _subagent_to_dict(record))

  def load_subagent(self, conv_id, sub_id):
    """Load a subagent record by its id.

    Parameters
    ----------
    conv_id : str
        Conversation the subagent belongs to.
    sub_id : str
        Identifier of the subagent record.

    Returns
    -------
    SubagentRecord
        The loaded subagent record.

    Raises
    ------
    libtbx.utils.Sorry
        If ``sub_id`` is unsafe as a path segment, or no record exists.
    """
    path = (self._conv_dir(conv_id) / "subagents"
            / ("%s.json" % _safe_segment(sub_id, "subagent")))
    if not path.exists():
      from libtbx.utils import Sorry
      raise Sorry("Subagent record not found: %s" % sub_id)
    return _subagent_from_dict(_read_json(path))

  def conv_dir(self, conv_id):
    """Public path to a conversation's directory (creates nothing)."""
    return self._conv_dir(conv_id)

  # ---- internal ------------------------------------------------------------

  def _ensure_root(self):
    self.root.mkdir(parents=True, exist_ok=True)
    (self.root / "conversations").mkdir(exist_ok=True)

  def _conv_dir(self, conv_id):
    return self.root / "conversations" / _safe_segment(conv_id, "conversation")

  def _refresh_index(self):
    """Re-derive the index from on-disk conversations and write it.

    Done as a full rewrite rather than an incremental update; the index
    is small enough that this is simpler.
    """
    self._ensure_root()
    metas = self._scan_conversation_metas()
    _atomic_write_json(self.root / "index.json", {
      "schema_version": _SCHEMA_VERSION,
      "conversations": [_meta_to_dict(m) for m in metas],
    })

  def _rebuild_index(self):
    metas = self._scan_conversation_metas()
    if metas:
      self._ensure_root()
      _atomic_write_json(self.root / "index.json", {
        "schema_version": _SCHEMA_VERSION,
        "conversations": [_meta_to_dict(m) for m in metas],
      })
    return metas

  def _scan_conversation_metas(self):
    conv_root = self.root / "conversations"
    if not conv_root.exists():
      return []
    out = []
    for conv_dir in sorted(conv_root.iterdir()):
      meta_path = conv_dir / "meta.json"
      if not meta_path.exists():
        continue
      try:
        out.append(_meta_from_dict(_read_json(meta_path)))
      except Exception as e:
        print("storage: skipping unreadable %s (%s)"
              % (meta_path, e), file=self.log)
    return out


# ---- serialization helpers -------------------------------------------------

def _safe_segment(value, kind):
  """Validate an externally-supplied identifier before joining it into a path.

  Conversation ids, attachment sha256s and subagent ids are normally
  uuids / hex, but they reach storage from on-disk conversation files that
  may have been shared or hand-edited. Rejecting empty / ``.`` / ``..`` /
  path separators / drive markers / absolute paths keeps a crafted id from
  escaping the chat root (``pathlib`` resolves ``..`` and drops the left
  operand entirely when the right operand is absolute).

  Parameters
  ----------
  value : str
      Identifier to validate (a conversation id, attachment sha256, or
      subagent id).
  kind : str
      Human-readable category used in the error message, e.g.
      ``"conversation"``, ``"attachment"``, or ``"subagent"``.

  Returns
  -------
  str
      ``str(value)`` unchanged, once validated as a single safe path
      segment.

  Raises
  ------
  libtbx.utils.Sorry
      If ``value`` is empty, ``.`` / ``..``, absolute, or contains a path
      separator, drive marker, or embedded NUL byte.
  """
  s = str(value)
  if (not s or s in (".", "..")
      or "/" in s or "\\" in s or ":" in s or "\x00" in s
      or os.path.isabs(s) or s != os.path.basename(s)):
    from libtbx.utils import Sorry
    raise Sorry("Unsafe %s identifier: %r" % (kind, value))
  return s


def _read_json(path):
  with open(path) as fh:
    return json.load(fh)


def _atomic_write_json(path, obj):
  tmp = path.with_suffix(path.suffix + ".tmp")
  with open(tmp, "w") as fh:
    json.dump(obj, fh, indent=2, default=_json_default)
  os.replace(tmp, path)


def _json_default(o):
  if isinstance(o, datetime):
    return o.isoformat()
  raise TypeError("Not JSON serializable: %r" % o)


def _parse_dt(s):
  if s is None:
    return None
  # fromisoformat handles offset-aware ISO timestamps in 3.11+.
  return datetime.fromisoformat(s)


def _check_schema_version(doc, source):
  """Validate a document's ``schema_version`` against the supported one.

  Currently v1 only (no migrations needed); this is the seam future
  migrations plug into.

  Parameters
  ----------
  doc : dict
      The loaded JSON document. Non-dict values are accepted and skipped.
  source : str
      Path or label of the document, used in the error message.

  Raises
  ------
  libtbx.utils.Sorry
      If ``schema_version`` is from a version this client cannot migrate.
  """
  if not isinstance(doc, dict):
    return
  version = doc.get("schema_version", _SCHEMA_VERSION)
  if version != _SCHEMA_VERSION:
    from libtbx.utils import Sorry
    raise Sorry(
      "%s: schema_version '%s' is not supported by this client "
      "(expected '%s'). Update PhenixChat or use a compatible version."
      % (source, version, _SCHEMA_VERSION))


def _meta_to_dict(m):
  return {
    "id": m.id,
    "title": m.title,
    "profile_name": m.profile_name,
    "model": m.model,
    "backend": m.backend,
    "created_at": m.created_at,
    "updated_at": m.updated_at,
    "archived": m.archived,
    "pinned": m.pinned,
    "summary": m.summary,
    "agent_session_id": m.agent_session_id,
    "schema_version": m.schema_version,
  }


def _meta_from_dict(d):
  return ConversationMeta(
    id=d["id"],
    title=d.get("title", ""),
    profile_name=d.get("profile_name", ""),
    model=d.get("model", ""),
    backend=d.get("backend", ""),
    created_at=_parse_dt(d["created_at"]),
    updated_at=_parse_dt(d["updated_at"]),
    archived=d.get("archived", False),
    pinned=d.get("pinned", False),
    summary=d.get("summary", ""),
    agent_session_id=d.get("agent_session_id"),
    schema_version=d.get("schema_version", _SCHEMA_VERSION),
  )


def _content_block_to_dict(b):
  # Recursively serialize nested ContentBlocks inside tool_result.content.
  data = {}
  for k, v in b.data.items():
    if isinstance(v, list) and v and isinstance(v[0], ContentBlock):
      data[k] = [_content_block_to_dict(inner) for inner in v]
    else:
      data[k] = v
  return {"type": b.type, "data": data}


def _content_block_from_dict(d):
  data = {}
  for k, v in d.get("data", {}).items():
    if isinstance(v, list) and v and isinstance(v[0], dict) and "type" in v[0]:
      data[k] = [_content_block_from_dict(inner) for inner in v]
    else:
      data[k] = v
  return ContentBlock(type=d["type"], data=data)


def _message_to_dict(m):
  out = {
    "role": m.role,
    "timestamp": m.timestamp,
    "content": [_content_block_to_dict(b) for b in m.content],
  }
  if m.stop_reason is not None:
    out["stop_reason"] = m.stop_reason
  if m.usage is not None:
    out["usage"] = asdict(m.usage)
  if m.model is not None:
    out["model"] = m.model
  if m.backend is not None:
    out["backend"] = m.backend
  return out


def _message_from_dict(d):
  usage = None
  if "usage" in d and d["usage"] is not None:
    u = d["usage"]
    usage = TokenUsage(
      input=u.get("input", 0),
      output=u.get("output", 0),
      cache_read=u.get("cache_read", 0),
      cache_creation=u.get("cache_creation", 0),
    )
  return Message(
    role=d["role"],
    content=[_content_block_from_dict(b) for b in d.get("content", [])],
    timestamp=_parse_dt(d["timestamp"]),
    stop_reason=d.get("stop_reason"),
    usage=usage,
    model=d.get("model"),
    backend=d.get("backend"),
  )


def _subagent_to_dict(r):
  return {
    "schema_version": _SCHEMA_VERSION,
    "sub_id": r.sub_id,
    "parent_conversation_id": r.parent_conversation_id,
    "parent_tool_use_id": r.parent_tool_use_id,
    "task": r.task,
    "profile_name": r.profile_name,
    "model": r.model,
    "started_at": r.started_at,
    "finished_at": r.finished_at,
    "final_text": r.final_text,
    "token_usage": asdict(r.token_usage),
    "messages": [_message_to_dict(m) for m in r.messages],
  }


def _subagent_from_dict(d):
  tu = d.get("token_usage", {})
  return SubagentRecord(
    sub_id=d["sub_id"],
    parent_conversation_id=d.get("parent_conversation_id", ""),
    parent_tool_use_id=d.get("parent_tool_use_id", ""),
    task=d.get("task", ""),
    profile_name=d.get("profile_name", ""),
    model=d.get("model", ""),
    started_at=_parse_dt(d["started_at"]),
    finished_at=_parse_dt(d["finished_at"]),
    final_text=d.get("final_text", ""),
    token_usage=TokenUsage(
      input=tu.get("input", 0),
      output=tu.get("output", 0),
      cache_read=tu.get("cache_read", 0),
      cache_creation=tu.get("cache_creation", 0),
    ),
    messages=[_message_from_dict(m) for m in d.get("messages", [])],
  )
