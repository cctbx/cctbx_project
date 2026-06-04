"""``MessageBubble`` — one ``Message`` rendered as a Qt widget.

Renders text, thinking, tool-use (via ``ToolCallDisclosure``),
tool-result, and image cells. Flattened layout: no frame border, role
bold-prefixed onto the first text cell.
"""

import json

from qttbx.qt import QtCore, QtGui, QtWidgets

from qttbx.widgets.chat.agent.conversation import ContentBlock, Message, now
from qttbx.widgets.chat.markdown_view import MarkdownView


class _ToolResultCell(QtWidgets.QFrame):
  def __init__(self, tool_use_id, content_blocks, is_error, parent=None):
    super().__init__(parent)
    self.setFrameShape(QtWidgets.QFrame.StyledPanel)
    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(6, 4, 6, 4)
    header = QtWidgets.QLabel(
      "result%s" % (" (error)" if is_error else ""), self)
    if is_error:
      header.setStyleSheet("color: #c0392b;")
    layout.addWidget(header)
    text = _flatten_result_text(content_blocks)
    if text:
      body = QtWidgets.QLabel(text, self)
      # Tool-result text is model/tool-controlled: render literally so an
      # embedded <img src="file://..."> can't load a local file.
      body.setTextFormat(QtCore.Qt.PlainText)
      body.setWordWrap(True)
      layout.addWidget(body)


class _ThinkingCell(QtWidgets.QFrame):
  def __init__(self, text, parent=None):
    super().__init__(parent)
    self.setFrameShape(QtWidgets.QFrame.NoFrame)
    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(6, 2, 6, 2)
    label = QtWidgets.QLabel("[thinking] %s" % (text or ""), self)
    # Thinking text is model-controlled: render literally (no rich text).
    label.setTextFormat(QtCore.Qt.PlainText)
    label.setWordWrap(True)
    label.setStyleSheet("color: palette(mid); font-style: italic;")
    layout.addWidget(label)
    self._label = label

  def append(self, text):
    cur = self._label.text()
    self._label.setText(cur + (text or ""))


class _ImageCell(QtWidgets.QFrame):
  """Inline image cell rendering a QPixmap as a thumbnail.

  Holds a QPixmap, renders it as a height-capped thumbnail
  (aspect-preserved), and opens ImageLightbox on click. The caption is
  rendered below the image when present; hidden when absent.

  ``clicked(conv_id, sha256)`` is the only public signal -- it
  propagates up through MessageBubble.image_clicked and
  ConversationView.image_clicked to ChatWindow._on_image_clicked, which
  uses the two strings to look up the artifact in the side panel.

  Parameters
  ----------
  pixmap : QtGui.QPixmap
      The full-resolution image; thumbnailed for display and shown
      at full size in the lightbox.
  caption : str, optional
      Text rendered below the thumbnail; the caption is hidden when
      absent.
  sha256 : str, optional
      Attachment hash; travels with the ``clicked`` signal.
  mime : str, optional
      Attachment MIME type; stored as metadata.
  conv_id : str, optional
      Conversation id; travels with the ``clicked`` signal.
  parent : QtWidgets.QWidget, optional
      Parent widget.
  """

  MAX_HEIGHT = 240

  clicked = QtCore.Signal(str, str)                # conv_id, sha256

  def __init__(self, pixmap, caption=None, sha256=None, mime=None,
               conv_id=None, parent=None):
    super().__init__(parent)
    self._pixmap = pixmap
    self._last_lightbox = None
    # Metadata travels with the clicked signal so the chat window can
    # focus / fetch the corresponding attachment.
    self.sha256 = sha256 or ""
    self.mime = mime or ""
    self.conv_id = conv_id or ""
    self.setFrameShape(QtWidgets.QFrame.NoFrame)
    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(0, 4, 0, 4)
    layout.setSpacing(2)
    if pixmap.height() > self.MAX_HEIGHT:
      scaled = pixmap.scaledToHeight(
        self.MAX_HEIGHT, QtCore.Qt.SmoothTransformation)
    else:
      scaled = pixmap
    self.thumbnail = QtWidgets.QLabel(self)
    self.thumbnail.setPixmap(scaled)
    self.thumbnail.setCursor(QtCore.Qt.PointingHandCursor)
    # Click-to-open. Override mousePressEvent rather than installing an
    # event filter -- simpler and the label has no other use for the
    # mouse press.
    self.thumbnail.mousePressEvent = self._thumb_clicked
    layout.addWidget(self.thumbnail)
    self.caption_label = QtWidgets.QLabel(caption or "", self)
    # Caption text is model/tool-controlled: render literally.
    self.caption_label.setTextFormat(QtCore.Qt.PlainText)
    self.caption_label.setStyleSheet(
      "color: palette(mid); font-style: italic;")
    self.caption_label.setWordWrap(True)
    layout.addWidget(self.caption_label)
    if not caption:
      self.caption_label.hide()

  def click(self):
    """Fire ``clicked`` programmatically without opening the lightbox.

    Used by tests so the assertion can run synchronously.
    """
    self.clicked.emit(self.conv_id or "", self.sha256 or "")

  def _thumb_clicked(self, _event):
    self.clicked.emit(self.conv_id or "", self.sha256 or "")
    self._open_lightbox()

  def _open_lightbox(self):
    from qttbx.widgets.chat.image_lightbox import ImageLightbox
    dlg = ImageLightbox(pixmap=self._pixmap, parent=self.window())
    # Record before showing so tests can introspect even if the dialog
    # is closed asynchronously.
    self._last_lightbox = dlg
    # `.open()` is non-blocking (vs `.exec()` which spins the modal
    # event loop). The dialog still acts modally and closes on click
    # via its own mousePressEvent. Using .open() keeps the call
    # testable under QT_QPA_PLATFORM=offscreen.
    dlg.open()


_ROLE_LABELS = {"user": "You", "assistant": "Claude"}


class MessageBubble(QtWidgets.QFrame):
  """One chat turn (user or assistant) rendered as a widget.

  The bubble accumulates child cells in vertical order matching
  ``message.content``.

  Flattened layout: no QFrame border / background; the role word is
  bold-prefixed onto the first text cell (e.g. ``<b>You:</b> refine
  1yjp``) rather than rendered as a separate header label. Turn
  separation comes from the 12 px top margin.

  The bubble has two construction modes: pass ``message`` (the legacy
  ConversationView path, with content blocks rendered immediately), or
  pass ``role`` (the lightweight path for tests and future direct
  callers, with no Message wired in until needed). Exactly one of the
  two must be supplied.

  Parameters
  ----------
  message : Message, optional
      The message to render. Mutually exclusive with ``role``.
  parent : QtWidgets.QWidget, optional
      Parent widget.
  storage : object, optional
      Attachment store used to decode image blocks; without it image
      cells fall back to a placeholder pixmap.
  conv_id : str, optional
      Conversation id used to resolve image attachments.
  role : str, optional
      ``"user"`` or ``"assistant"`` to build an empty bubble. Mutually
      exclusive with ``message``.

  Raises
  ------
  TypeError
      If neither ``message`` nor ``role`` is supplied.
  """

  image_clicked = QtCore.Signal(str, str)        # conv_id, sha256

  def __init__(self, message=None, parent=None, storage=None, conv_id=None,
               role=None):
    super().__init__(parent)
    # Two construction modes:
    #   * legacy: MessageBubble(message=<Message>) -- ConversationView path,
    #     content blocks rendered immediately.
    #   * lightweight: MessageBubble(role="user"|"assistant") -- tests and
    #     future direct callers; no Message wired in until needed.
    if message is None:
      if role is None:
        raise TypeError(
          "MessageBubble requires either `message` or `role`")
      message = Message(role=role, timestamp=now(), content=[])
    self.message = message
    self.storage = storage
    self.conv_id = conv_id
    self._role = message.role
    # Flattened: no frame styling, no background.
    self.setFrameStyle(QtWidgets.QFrame.NoFrame)
    self._layout = QtWidgets.QVBoxLayout(self)
    # 12 px top margin per bubble separates turns; no padding on sides
    # (the ConversationView owns the horizontal margins).
    self._layout.setContentsMargins(0, 12, 0, 0)
    self._layout.setSpacing(4)

    self._text_view = None
    self._thinking_cell = None
    self._first_text_cell = None
    self._first_text_cell_added = False
    # tool_id -> ToolCallDisclosure. Populated by both add_tool_use_cell()
    # (runner-driven, live tool calls) and by _add_block() when replaying
    # a stored conversation -- so set_tool_use_finished() can fold a later
    # tool_result block into the right disclosure widget regardless of
    # whether the bubble was built live or loaded from storage.
    self._tool_cells_by_id = {}

    for block in message.content:
      self._add_block(block)

  # ---- block rendering -----------------------------------------------------

  def _add_block(self, block):
    if block.type == "text":
      self._append_text_markdown(block.data.get("text", ""))
    elif block.type == "thinking":
      cell = _ThinkingCell(block.data.get("text", ""), self)
      self._thinking_cell = cell
      # Non-text cell: drop _text_view so subsequent text creates a fresh
      # MarkdownView (otherwise text-after-thinking would merge into the
      # pre-thinking view, breaking visual ordering).
      self._text_view = None
      self._layout.addWidget(cell)
    elif block.type == "tool_use":
      # Route legacy message-replay tool_use blocks through the same
      # ToolCallDisclosure widget the runner-driven path uses, so the UI
      # is consistent whether the bubble was loaded from storage or built
      # live by streaming.
      self.add_tool_use_cell(
        tool_id=block.data.get("id", ""),
        name=block.data.get("name", "?"),
        args=block.data.get("input", {}))
    elif block.type == "tool_result":
      # If the matching tool_use cell is in this bubble, fold the result
      # into it; otherwise (orphan tool_result -- e.g. result block that
      # arrived in a follow-up message) fall back to the legacy cell.
      tool_use_id = block.data.get("tool_use_id", "")
      is_error = bool(block.data.get("is_error", False))
      result_text = _flatten_result_text(block.data.get("content", []))
      if tool_use_id and tool_use_id in self._tool_cells_by_id:
        if is_error:
          self.set_tool_use_finished(
            tool_id=tool_use_id, error=result_text or "error")
        else:
          self.set_tool_use_finished(
            tool_id=tool_use_id, result=result_text)
      else:
        cell = _ToolResultCell(
          tool_use_id,
          block.data.get("content", []),
          is_error,
          self)
        self._text_view = None
        self._layout.addWidget(cell)
    elif block.type == "server_tool_use":
      # API-executed tool call (e.g., web_search, web_fetch,
      # code_execution). Reuse ToolCallDisclosure; the result lands
      # later in the same turn as a server_tool_result block and folds
      # into this cell via set_tool_use_finished.
      self.add_tool_use_cell(
        tool_id=block.data.get("id", ""),
        name="server: " + (block.data.get("name", "?") or "?"),
        args=block.data.get("input", {}))
    elif block.type == "server_tool_result":
      tool_use_id = block.data.get("tool_use_id", "")
      result_text = _format_server_tool_result(block.data.get("content", {}))
      if tool_use_id and tool_use_id in self._tool_cells_by_id:
        self.set_tool_use_finished(
          tool_id=tool_use_id, result=result_text)
      else:
        # Orphan (no matching use cell in this bubble): render as a
        # plain text cell so the result still appears.
        self._append_text_markdown(
          "*server tool result:* `" + tool_use_id + "`\n\n```\n"
          + result_text + "\n```")
    elif block.type == "image":
      # Load the attachment into a QPixmap via image_cache (decoded once
      # per session) and render as an inline thumbnail. Falls back to a
      # placeholder pixmap when storage/sha256 are missing.
      sha256 = block.data.get("attachment_sha256", "")
      mime = block.data.get("mime", "")
      caption = block.data.get("caption")
      if self.storage is not None and self.conv_id is not None and sha256:
        from qttbx.widgets.chat.image_cache import get_image
        img = get_image(self.storage, self.conv_id, sha256)
        pixmap = QtGui.QPixmap.fromImage(img)
      else:
        # No storage context (legacy/tests). A 1x1 placeholder keeps the
        # cell renderable.
        pixmap = QtGui.QPixmap(1, 1)
        pixmap.fill(QtCore.Qt.lightGray)
      cell = self.add_image_cell(
        pixmap=pixmap, caption=caption,
        sha256=sha256, mime=mime, conv_id=self.conv_id)
      cell.clicked.connect(self.image_clicked)
    # Unknown block types are silently dropped - they shouldn't reach the
    # bubble in normal operation (storage rejects unknown types on load).

  def _ensure_text_view(self):
    if self._text_view is None:
      self._text_view = MarkdownView(self)
      self._layout.addWidget(self._text_view)
      if self._first_text_cell is None:
        self._first_text_cell = self._text_view

  def _role_label_text(self):
    return _ROLE_LABELS.get(self._role, self._role.title())

  def _append_text_markdown(self, text):
    """Append markdown to the shared text view.

    Prepends the bold role marker on the first text cell.
    """
    self._ensure_text_view()
    if not self._first_text_cell_added:
      prefix = "**%s:** " % self._role_label_text()
      self._text_view.append_markdown(prefix + (text or ""))
      self._first_text_cell_added = True
    else:
      self._text_view.append_markdown(text or "")

  def add_text(self, text):
    """Append a text block to the bubble.

    The first call gets the bold-prefixed role marker (``You:`` /
    ``Claude:``); subsequent calls render as plain markdown. Also
    mirrors into ``message.content`` so ``combined_text()`` stays
    consistent.

    Returns
    -------
    MarkdownView
        The shared text view the block was appended to.
    """
    self._append_text_markdown(text)
    last = self.message.content[-1] if self.message.content else None
    if last is not None and last.type == "text":
      last.data["text"] = (last.data.get("text", "") or "") + (text or "")
    else:
      self.message.content.append(
        ContentBlock(type="text", data={"text": text or ""}))
    return self._text_view

  def first_text_cell_html(self):
    """Return the HTML of the first text cell.

    Used by tests and any future inline assertions.

    Returns
    -------
    str
        The first text cell's HTML, or an empty string if no text yet.
    """
    if self._first_text_cell is None:
      return ""
    if hasattr(self._first_text_cell, "toHtml"):
      return self._first_text_cell.toHtml()
    return self._first_text_cell.text()

  # ---- image cells --------------------------------------------------------

  def add_image_cell(self, pixmap, caption=None, sha256=None, mime=None,
                     conv_id=None):
    """Insert an inline image cell.

    The cell renders a height-capped thumbnail and opens ImageLightbox
    on click. The ``sha256`` / ``mime`` / ``conv_id`` kwargs travel
    with the cell's ``clicked(conv_id, sha256)`` signal so the
    ConversationView can re-emit them to the chat window's
    artifact-focus handler.

    Parameters
    ----------
    pixmap : QtGui.QPixmap
        The image to render as a thumbnail.
    caption : str, optional
        Text rendered below the thumbnail.
    sha256 : str, optional
        Attachment hash; travels with the cell's ``clicked`` signal.
    mime : str, optional
        Attachment MIME type.
    conv_id : str, optional
        Conversation id; travels with the cell's ``clicked`` signal.

    Returns
    -------
    _ImageCell
        The inserted image cell.
    """
    cell = _ImageCell(
      pixmap=pixmap, caption=caption,
      sha256=sha256, mime=mime, conv_id=conv_id, parent=self)
    # Reset _text_view so subsequent add_text creates a fresh
    # MarkdownView -- same invariant as add_tool_use_cell.
    self._text_view = None
    self._layout.addWidget(cell)
    return cell

  # ---- tool cells ---------------------------------------------------------

  def add_tool_use_cell(self, tool_id, name, args):
    """Insert a collapsed ToolCallDisclosure row keyed by ``tool_id``.

    The runner calls ``set_tool_use_finished(tool_id, ...)`` to
    transition the cell to a terminal state once the tool returns.

    Parameters
    ----------
    tool_id : str
        Identifier the cell is keyed by for later finishing.
    name : str
        Tool name shown in the disclosure header.
    args : dict
        Tool input arguments; shown in the disclosure when non-empty.

    Returns
    -------
    ToolCallDisclosure
        The inserted disclosure row.
    """
    from qttbx.widgets.chat.tool_call_disclosure import ToolCallDisclosure
    cell = ToolCallDisclosure(name=name, status="running", parent=self)
    if args:
      cell.set_args(args)
    if tool_id:
      self._tool_cells_by_id[tool_id] = cell
    # Reset _text_view so subsequent add_text creates a fresh MarkdownView
    # -- without this, text-after-tool would merge into the pre-tool view
    # (out-of-order visual rendering).
    self._text_view = None
    self._layout.addWidget(cell)
    return cell

  def set_tool_use_finished(self, tool_id, elapsed=None, result=None,
                            error=None, cancelled=False):
    """Transition a previously-added tool cell to a terminal state.

    No-op for an unknown ``tool_id`` (so the runner can call this
    safely even if it never added the cell).

    Parameters
    ----------
    tool_id : str
        Identifier of the cell to finish.
    elapsed : optional
        Elapsed-time suffix appended to the finished status when set.
    result : optional
        Result payload rendered into the cell when set.
    error : optional
        Error message; transitions the cell to a failed state when set.
    cancelled : bool, optional
        When true, transitions the cell to a cancelled state.
    """
    cell = self._tool_cells_by_id.get(tool_id)
    if cell is None:
      return
    if cancelled:
      cell.set_status("cancelled", color="muted")
    elif error is not None:
      cell.set_status("failed: %s" % error, color="error")
    else:
      suffix = ", %s" % elapsed if elapsed else ""
      cell.set_status("finished%s" % suffix, color="default")
    if result is not None:
      cell.set_result(result)

  # ---- streaming -----------------------------------------------------------

  def append_text_delta(self, text):
    """Append a streamed text delta to the bubble and message content."""
    self._append_text_markdown(text)
    # Keep message.content in sync so combined_text() reflects what was
    # streamed in (the bubble is the source of truth for on-screen content).
    last = self.message.content[-1] if self.message.content else None
    if last is not None and last.type == "text":
      last.data["text"] = (last.data.get("text", "") or "") + (text or "")
    else:
      self.message.content.append(
        ContentBlock(type="text", data={"text": text or ""}))

  def append_thinking_delta(self, text):
    """Append a streamed thinking delta to the bubble and message content."""
    if self._thinking_cell is None:
      self._thinking_cell = _ThinkingCell("", self)
      self._layout.addWidget(self._thinking_cell)
    self._thinking_cell.append(text)
    # Mirror into message.content (same shape as append_text_delta): extend
    # the LAST block only if it is thinking; otherwise append a new one.
    last = self.message.content[-1] if self.message.content else None
    if last is not None and last.type == "thinking":
      last.data["text"] = (last.data.get("text", "") or "") + (text or "")
    else:
      self.message.content.append(
        ContentBlock(type="thinking", data={"text": text or ""}))

  def append_block(self, block):
    """Append a new ContentBlock to the in-progress bubble.

    Post-stream tool_use and tool_result blocks land here.

    Parameters
    ----------
    block : ContentBlock
        The content block to render and mirror into ``message.content``.
    """
    self.message.content.append(block)
    self._add_block(block)

  # ---- introspection (for tests + sidebar previews) -----------------------

  def combined_text(self):
    """Return a flat-text representation of the bubble's content.

    The bubble is the source of truth for the on-screen content.

    Returns
    -------
    str
        The content blocks flattened to newline-joined text.
    """
    parts = []
    for block in self.message.content:
      if block.type == "text":
        parts.append(block.data.get("text", ""))
      elif block.type == "thinking":
        parts.append("[thinking] " + block.data.get("text", ""))
      elif block.type == "tool_use":
        parts.append("[tool_use %s %s]" % (
          block.data.get("name", ""),
          _short_json(block.data.get("input", {}))))
      elif block.type == "tool_result":
        parts.append("[tool_result]")
      elif block.type == "image":
        parts.append("[image %s %s]" % (
          (block.data.get("attachment_sha256", "") or "")[:8],
          block.data.get("mime", "")))
    return "\n".join(parts)


def _short_json(d, limit=80):
  try:
    s = json.dumps(d, default=str)
  except Exception:
    s = repr(d)
  if len(s) > limit:
    s = s[:limit - 1] + "..."
  return s


def _flatten_result_text(blocks):
  out = []
  for b in blocks or []:
    if isinstance(b, dict):
      t = b.get("type")
      data = b.get("data", {})
    else:
      t = getattr(b, "type", None)
      data = getattr(b, "data", {})
    if t == "text":
      out.append(data.get("text", ""))
  return "\n".join(out)


def _format_server_tool_result(content):
  """Render a server-tool result payload for the disclosure body.

  The Anthropic API returns an opaque dict per server tool. Pretty-print
  it as JSON; for ``web_search`` specifically, extract the result list
  into a more readable bullet list so the disclosure isn't a wall of
  JSON for the most common case.
  """
  if not isinstance(content, dict):
    return str(content)
  # web_search: content has "content": [{"type":"web_search_result", ...}, ...]
  inner = content.get("content")
  if isinstance(inner, list) and inner and isinstance(inner[0], dict) \
      and inner[0].get("type") == "web_search_result":
    lines = []
    for r in inner:
      title = r.get("title", "(no title)")
      url = r.get("url", "")
      lines.append("- %s\n  %s" % (title, url))
    return "\n".join(lines)
  try:
    return json.dumps(content, indent=2, default=str)
  except Exception:
    return repr(content)
