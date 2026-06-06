"""Serialize a Conversation to a markdown string.

Used by the chat window's 'Save chat' button. Pure function — takes a
Conversation and (optionally) a ConversationStorage so attachment paths
resolve to real on-disk locations the rendered markdown can link to.
Without storage, image cells fall back to a sha256-tagged placeholder."""

import json


def conversation_to_markdown(conv, storage=None):
  """Return a markdown string representing ``conv``.

  The layout is::

      # <title>
      *<created_at> · model: <model> · profile: <profile>*

      ---

      ## You
      <text content>

      ## <assistant>
      <text content, possibly multiple paragraphs>

  Each role block stitches the message's content blocks in order:

  - text -> emitted verbatim (already markdown from the model)
  - image -> ``![caption](path)`` or placeholder when storage absent
  - thinking -> skipped (extended-thinking is internal, not the chat)
  - tool_use -> fenced code block tagged 'tool-use' with JSON input
  - tool_result -> fenced code block tagged 'tool-result' with the text
    content; nested image blocks render as image links

  Errors / malformed blocks fall back to a one-line ``*[unknown]*``
  placeholder so a single bad block doesn't sink the export.

  Parameters
  ----------
  conv : Conversation
      The conversation to serialize.
  storage : ConversationStorage, optional
      Attachment store used so image cells resolve to real on-disk
      paths the rendered markdown can link to. Without it, image cells
      fall back to a sha256-tagged placeholder.

  Returns
  -------
  str
      The conversation rendered as markdown, ending with a newline.
  """
  meta = conv.meta
  out = []
  out.append("# %s" % (meta.title or "Conversation"))
  out.append("")
  out.append("*%s · model: %s · profile: %s*" % (
    _fmt_ts(meta.created_at), meta.model or "(unset)",
    meta.profile_name or "(unset)"))
  out.append("")
  out.append("---")
  out.append("")
  for msg in conv.messages or []:
    role = _role_label(msg)
    out.append("## %s" % role)
    out.append("")
    for block in msg.content or []:
      rendered = _render_block(block, storage, conv.meta.id)
      if rendered:
        out.append(rendered)
        out.append("")
  # Trim a single trailing blank line; leave the final newline so the
  # file ends with one.
  while len(out) > 1 and out[-1] == "" and out[-2] == "":
    out.pop()
  return "\n".join(out) + "\n"


_ROLE_LABELS = {
  "user": "You",
  "system": "System",
}


def _role_label(msg):
  """Section label for a message: the backend's display name for an assistant
  turn (so each section names what produced it), else the role word."""
  role = getattr(msg, "role", None)
  if role == "assistant":
    backend = getattr(msg, "backend", None)
    if backend:
      from qttbx.widgets.chat.agent.profile import backend_display_name
      return backend_display_name(backend)
    return "Assistant"
  return _ROLE_LABELS.get(role, (role or "?").title())


def _fmt_ts(ts):
  if ts is None:
    return "(no timestamp)"
  try:
    return ts.strftime("%Y-%m-%d %H:%M UTC")
  except Exception:
    return str(ts)


def _render_block(block, storage, conv_id):
  t = getattr(block, "type", None)
  data = getattr(block, "data", None) or {}
  if t == "text":
    return (data.get("text", "") or "").rstrip()
  if t == "image":
    return _render_image(data, storage, conv_id)
  if t == "thinking":
    return ""
  if t == "tool_use":
    name = data.get("name", "?")
    try:
      payload = json.dumps(data.get("input", {}), indent=2, default=str)
    except Exception:
      payload = repr(data.get("input"))
    return "```tool-use %s\n%s\n```" % (name, payload)
  if t == "tool_result":
    return _render_tool_result(data, storage, conv_id)
  return "*[unknown block: %s]*" % (t or "")


def _render_image(data, storage, conv_id):
  sha = data.get("attachment_sha256", "") or ""
  caption = (data.get("caption") or "image").strip()
  mime = data.get("mime", "") or ""
  path = _resolve_attachment_path(storage, conv_id, sha)
  if path:
    return "![%s](%s)" % (caption, path)
  short = sha[:12] if sha else "(no sha)"
  return "*[image %s %s]*" % (short, mime)


def _resolve_attachment_path(storage, conv_id, sha):
  if storage is None or not sha or not conv_id:
    return None
  try:
    att_dir = storage._conv_dir(conv_id) / "attachments"
    matches = sorted(att_dir.glob("sha256-%s.*" % sha))
  except Exception:
    return None
  if not matches:
    return None
  return str(matches[0])


def _render_tool_result(data, storage, conv_id):
  content = data.get("content", []) or []
  if isinstance(content, str):
    return "```tool-result\n%s\n```" % content.rstrip()
  parts = []
  images = []
  for inner in content:
    it = getattr(inner, "type", None) or (
      inner.get("type") if isinstance(inner, dict) else None)
    idata = getattr(inner, "data", None) or (
      inner.get("data", {}) if isinstance(inner, dict) else {}) or {}
    if it == "text":
      parts.append((idata.get("text", "") or "").rstrip())
    elif it == "image":
      images.append(_render_image(idata, storage, conv_id))
    else:
      parts.append("[%s]" % (it or "block"))
  out = []
  if parts:
    out.append("```tool-result\n%s\n```" % "\n".join(
      p for p in parts if p))
  out.extend(images)
  return "\n\n".join(out)
