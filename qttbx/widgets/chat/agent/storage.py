"""Per-project conversation storage.

Section 6 of the design spec. Layout:
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
  """Reads and writes conversations under a project's .phenix_chat/ dir.

  Construction does not touch the filesystem. Directories appear lazily on
  first write."""

  def __init__(self, project_dir, log=None):
    self.project_dir = Path(project_dir)
    self.root = chat_root_for(self.project_dir)
    self.log = log if log is not None else sys.stdout

  # ---- conversations -------------------------------------------------------

  def list_conversations(self):
    """Read the cached index, or rebuild from on-disk conversations if
    missing/corrupt."""
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
    """Write meta.json + messages.json atomically, then update index."""
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
    """Content-addressed store. Returns an Attachment referencing the
    sha256 of the bytes. Idempotent: writing the same bytes twice
    produces one file."""
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
    att_dir = self._conv_dir(conv_id) / "attachments"
    matches = sorted(att_dir.glob("sha256-%s.*" % sha256))
    if not matches:
      from libtbx.utils import Sorry
      raise Sorry("Attachment not found: sha256=%s" % sha256)
    with open(matches[0], "rb") as fh:
      return fh.read()

  # ---- subagents -----------------------------------------------------------

  def store_subagent(self, conv_id, record):
    self._ensure_root()
    sub_dir = self._conv_dir(conv_id) / "subagents"
    sub_dir.mkdir(parents=True, exist_ok=True)
    path = sub_dir / ("%s.json" % record.sub_id)
    _atomic_write_json(path, _subagent_to_dict(record))

  def load_subagent(self, conv_id, sub_id):
    path = self._conv_dir(conv_id) / "subagents" / ("%s.json" % sub_id)
    if not path.exists():
      from libtbx.utils import Sorry
      raise Sorry("Subagent record not found: %s" % sub_id)
    return _subagent_from_dict(_read_json(path))

  # ---- internal ------------------------------------------------------------

  def _ensure_root(self):
    self.root.mkdir(parents=True, exist_ok=True)
    (self.root / "conversations").mkdir(exist_ok=True)

  def _conv_dir(self, conv_id):
    return self.root / "conversations" / conv_id

  def _refresh_index(self):
    """Re-derive the index from on-disk conversations and write atomically.
    Simpler than incremental update; index is small enough."""
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
  """Raise Sorry if the document's schema_version is from a future version
  we don't know how to migrate. Currently v1 only (no migrations needed);
  this is the seam future migrations plug into (Section 6.6)."""
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
    "created_at": m.created_at,
    "updated_at": m.updated_at,
    "archived": m.archived,
    "pinned": m.pinned,
    "summary": m.summary,
    "schema_version": m.schema_version,
  }


def _meta_from_dict(d):
  return ConversationMeta(
    id=d["id"],
    title=d.get("title", ""),
    profile_name=d.get("profile_name", ""),
    model=d.get("model", ""),
    created_at=_parse_dt(d["created_at"]),
    updated_at=_parse_dt(d["updated_at"]),
    archived=d.get("archived", False),
    pinned=d.get("pinned", False),
    summary=d.get("summary", ""),
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
