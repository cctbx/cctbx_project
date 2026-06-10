"""Multi-line chat input with attachments, drag-drop, paste, and buttons.

The button row holds Save chat / Auto-approve / attach / Send.

The Send button's label flips to 'Stop' when ``set_busy(True)`` is
called; ``click_send`` emits ``stop`` in that mode. Attachments are
added via ``attach_bytes`` (also fed by drag-drop, Ctrl/Cmd+V paste,
or the ``@`` attach button) and travel with the next ``send`` signal
as a ``list[dict]`` of ``{"bytes": bytes, "mime": str, "filename":
str}``. ``auto_approve_changed`` carries the auto-approve toggle's
checked state; ``save_chat`` is a parameterless signal the chat window
listens to so it can prompt for a destination and write the
markdown export.
"""

from qttbx.qt import QtCore, QtGui, QtWidgets


class MessageInput(QtWidgets.QWidget):
  """Multi-line chat input with attachments and a button row."""

  send = QtCore.Signal(str, list)                  # text, attachments (list)
  stop = QtCore.Signal()
  attachment_rejected = QtCore.Signal(str)
  save_chat = QtCore.Signal()                      # 'Save chat' button click
  auto_approve_changed = QtCore.Signal(bool)       # checked state

  # Idle placeholder text. The assistant name defaults to "Claude" but is
  # rewritten per session by ChatWindow via set_assistant_name() so the box
  # names the active backend (GPT / Gemini / …). ChatWindow also swaps in a
  # cycling 'Thinking...' placeholder while a turn is in flight and restores
  # the idle one on turn_done.
  _PLACEHOLDER_FMT = "Message %s...  (Ctrl/Cmd+Enter to send)"
  DEFAULT_PLACEHOLDER = _PLACEHOLDER_FMT % "Claude"

  _MAX_IMAGE_BYTES = 20 * 1024 * 1024
  _ALLOWED_MIMES = frozenset((
    "image/png", "image/jpeg", "image/webp", "image/gif"))

  def __init__(self, parent=None):
    super().__init__(parent)
    self._busy = False
    self._attachments = []          # list[dict]
    self._max_image_bytes = self._MAX_IMAGE_BYTES
    # Up/Down input history. The recall list is supplied per conversation
    # by the chat window (set_history) and appended to on every send.
    # _history_index is None while editing live text, else an index into
    # _history while navigating; _history_draft preserves the live text so
    # Down can restore it past the newest entry.
    self._history = []
    self._history_index = None
    self._history_draft = ""
    self.setAcceptDrops(True)
    self._build_ui()

  def _build_ui(self):
    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(4, 4, 4, 4)
    self._chip_bar = QtWidgets.QWidget(self)
    self._chip_layout = QtWidgets.QHBoxLayout(self._chip_bar)
    self._chip_layout.setContentsMargins(0, 0, 0, 0)
    self._chip_bar.setVisible(False)
    layout.addWidget(self._chip_bar)
    # Text edit gets the full width; the VBox owns horizontal expansion
    # so the edit grows with the panel automatically.
    self._edit = QtWidgets.QPlainTextEdit(self)
    # Per-session idle placeholder; set_assistant_name() rewrites it.
    self._idle_placeholder = self.DEFAULT_PLACEHOLDER
    self._edit.setPlaceholderText(self._idle_placeholder)
    self._edit.installEventFilter(self)
    # Cache the theme's stock placeholder colour (dim grey on most
    # palettes). set_placeholder(..., dim=False) swaps in the regular
    # text colour so 'Thinking...' verbs read at full contrast;
    # reset_placeholder restores this dim value.
    self._dim_placeholder_color = (
      self._edit.palette().color(QtGui.QPalette.PlaceholderText))
    layout.addWidget(self._edit)
    # Button row below the edit. Save chat sits on the left; the
    # auto-approve toggle sits in the centre (flanked by stretches so
    # it stays centred as the row grows); Attach + Send are on the
    # right.
    button_row = QtWidgets.QHBoxLayout()
    button_row.setContentsMargins(0, 0, 0, 0)
    self._save_chat_btn = QtWidgets.QPushButton("Save chat", self)
    self._save_chat_btn.setToolTip(
      "Save the current conversation as a markdown file")
    # QPushButton.clicked emits (checked: bool); relay through a slot
    # so the public save_chat signal stays parameterless across both
    # PySide2 and PySide6 (signal-to-signal arity bridging is allowed
    # but spells less clearly than an explicit emit).
    self._save_chat_btn.clicked.connect(self._on_save_chat_clicked)
    button_row.addWidget(self._save_chat_btn)
    button_row.addStretch(1)
    self._auto_approve_btn = QtWidgets.QPushButton("Auto-approve", self)
    self._auto_approve_btn.setCheckable(True)
    self._auto_approve_btn.setToolTip(
      "Approve every tool request automatically for this session. "
      "Skips the approval card -- use with care.")
    self._auto_approve_btn.toggled.connect(self._on_auto_approve_toggled)
    button_row.addWidget(self._auto_approve_btn)
    button_row.addStretch(1)
    self._attach_btn = QtWidgets.QToolButton(self)
    self._attach_btn.setText("@")    # ASCII-safe; UI later replaces with icon
    self._attach_btn.setToolTip("Attach an image")
    self._attach_btn.clicked.connect(self._on_attach_clicked)
    button_row.addWidget(self._attach_btn)
    self._button = QtWidgets.QPushButton("Send", self)
    self._button.clicked.connect(self.click_send)
    button_row.addWidget(self._button)
    layout.addLayout(button_row)

  # ---- text helpers --------------------------------------------------------

  def text(self):
    return self._edit.toPlainText()

  def set_text(self, value):
    self._edit.setPlainText(value)

  def clear(self):                                       # noqa: A003
    self._edit.clear()

  # ---- input history (Up / Down recall) ------------------------------------

  def set_history(self, texts):
    """Replace the Up/Down recall list and reset navigation.

    Called by the chat window when a conversation is loaded or switched,
    with that conversation's user-message texts (oldest first), so recall
    is per-conversation and survives a reload.

    Parameters
    ----------
    texts : iterable of str
        Previously-sent input texts, oldest first. Empty or
        whitespace-only entries are dropped.
    """
    self._history = [t for t in texts if t and t.strip()]
    self._history_index = None
    self._history_draft = ""

  def _history_recall_prev(self):
    """Step one entry back (older). Returns ``True`` when handled."""
    if not self._history:
      return False
    if self._history_index is None:
      self._history_draft = self.text()
      self._history_index = len(self._history)
    if self._history_index > 0:
      self._history_index -= 1
      self._set_edit_text(self._history[self._history_index])
    return True

  def _history_recall_next(self):
    """Step one entry forward (newer), restoring the saved draft past the
    newest entry. Returns ``True`` only when we were navigating."""
    if self._history_index is None:
      return False
    self._history_index += 1
    if self._history_index >= len(self._history):
      self._history_index = None
      self._set_edit_text(self._history_draft)
    else:
      self._set_edit_text(self._history[self._history_index])
    return True

  def _set_edit_text(self, text):
    """Replace the edit contents and put the cursor at the end."""
    self._edit.setPlainText(text or "")
    cursor = self._edit.textCursor()
    cursor.movePosition(QtGui.QTextCursor.End)
    self._edit.setTextCursor(cursor)

  def _cursor_at_first_line(self):
    """True when the cursor is on the first visual line -- the only place
    Up recalls history instead of moving the cursor (wrap-aware)."""
    probe = QtGui.QTextCursor(self._edit.textCursor())
    return not probe.movePosition(QtGui.QTextCursor.Up)

  def _cursor_at_last_line(self):
    """True when the cursor is on the last visual line -- the only place
    Down walks history forward instead of moving the cursor."""
    probe = QtGui.QTextCursor(self._edit.textCursor())
    return not probe.movePosition(QtGui.QTextCursor.Down)

  def set_placeholder(self, text, dim=True):
    """Update the edit's placeholder text.

    Parameters
    ----------
    text : str
        Placeholder string. Empty / None clears it.
    dim : bool, optional
        When ``True`` (default) the placeholder is rendered in the
        theme's ``PlaceholderText`` palette colour (the usual dim
        grey). When ``False`` it uses the regular ``Text`` colour so
        the placeholder reads at full contrast -- used by ChatWindow
        for the 'Thinking...' verb cycle so the verbs aren't visually
        muted.
    """
    palette = self._edit.palette()
    target_color = (
      self._dim_placeholder_color if dim
      else palette.color(QtGui.QPalette.Text))
    palette.setColor(QtGui.QPalette.PlaceholderText, target_color)
    self._edit.setPalette(palette)
    self._edit.setPlaceholderText(text or "")
    # QPlainTextEdit only repaints the cursor area on placeholder
    # change when focused, so the user sees a stale tail (just the
    # first letter of the new verb updates). Force the viewport to
    # repaint the whole placeholder region.
    self._edit.viewport().update()

  def reset_placeholder(self):
    """Restore the idle placeholder (per-session text + dim colour)."""
    self.set_placeholder(self._idle_placeholder, dim=True)

  def set_assistant_name(self, name):
    """Set the assistant display name shown in the idle placeholder.

    ChatWindow calls this with the active backend's display name (Claude /
    GPT / Gemini / …). Applied immediately when the box is idle; the
    'Thinking...' verb cycle reverts to it via :meth:`reset_placeholder`.

    Parameters
    ----------
    name : str
        The assistant display name. Empty / ``None`` falls back to a
        generic label.
    """
    self._idle_placeholder = self._PLACEHOLDER_FMT % (name or "the assistant")
    if not getattr(self, "_busy", False):
      self.set_placeholder(self._idle_placeholder, dim=True)

  # ---- busy / button state -------------------------------------------------

  def set_busy(self, busy):
    self._busy = bool(busy)
    self._button.setText("Stop" if self._busy else "Send")

  def is_busy(self):
    return self._busy

  # ---- attachments ---------------------------------------------------------

  def attach_bytes(self, data, mime, filename=""):
    """Add bytes to the pending attachment list.

    Validates mime type and size; emits ``attachment_rejected`` on
    failure.

    Parameters
    ----------
    data : bytes
        Raw image bytes to attach.
    mime : str
        MIME type; must be one of the allowed image types.
    filename : str, optional
        Display name for the attachment. Defaults to ``"image"``.

    Returns
    -------
    bool
        ``True`` if the attachment was accepted, ``False`` if rejected.
    """
    if mime not in self._ALLOWED_MIMES:
      self.attachment_rejected.emit(
        "Unsupported attachment type: %s" % mime)
      return False
    data = self._maybe_resample(data, mime)
    if len(data) > self._max_image_bytes:
      self.attachment_rejected.emit(
        "Attachment too large (%d bytes)" % len(data))
      return False
    self._attachments.append({
      "bytes": data, "mime": mime, "filename": filename or "image"})
    self._refresh_chip_bar()
    return True

  def attachment_count(self):
    return len(self._attachments)

  def remove_attachment(self, index):
    if 0 <= index < len(self._attachments):
      self._attachments.pop(index)
      self._refresh_chip_bar()

  def _maybe_resample(self, data, mime):
    """Shrink an oversized image to fit under the byte cap.

    Repeatedly halves dimensions until under the byte cap, up to 5
    iterations.

    Parameters
    ----------
    data : bytes
        Raw image bytes.
    mime : str
        MIME type; selects the re-encode format (PNG vs JPG).

    Returns
    -------
    bytes
        The smaller bytes, or the original if Qt can't read them
        (which lets the size check catch it later).
    """
    if len(data) <= self._max_image_bytes:
      return data
    img = QtGui.QImage()
    if not img.loadFromData(data):
      return data
    target_fmt = "PNG" if mime == "image/png" else "JPG"
    out = data
    for _ in range(5):
      img = img.scaled(max(1, img.width() // 2),
                       max(1, img.height() // 2),
                       QtCore.Qt.KeepAspectRatio,
                       QtCore.Qt.SmoothTransformation)
      buf = QtCore.QBuffer()
      buf.open(QtCore.QBuffer.WriteOnly)
      img.save(buf, target_fmt)
      out = bytes(buf.data())
      if len(out) <= self._max_image_bytes:
        return out
    return out

  def _refresh_chip_bar(self):
    while self._chip_layout.count():
      item = self._chip_layout.takeAt(0)
      w = item.widget()
      if w is not None:
        w.setParent(None)
        w.deleteLater()
    for i, att in enumerate(self._attachments):
      chip = QtWidgets.QFrame(self)
      chip.setFrameShape(QtWidgets.QFrame.StyledPanel)
      h = QtWidgets.QHBoxLayout(chip)
      h.setContentsMargins(4, 1, 4, 1)
      h.addWidget(QtWidgets.QLabel(att["filename"], chip))
      btn = QtWidgets.QToolButton(chip)
      btn.setText("x")
      btn.clicked.connect(lambda _checked=False, idx=i:
                          self.remove_attachment(idx))
      h.addWidget(btn)
      self._chip_layout.addWidget(chip)
    self._chip_bar.setVisible(bool(self._attachments))

  # ---- save chat -----------------------------------------------------------

  def _on_save_chat_clicked(self, _checked=False):
    self.save_chat.emit()

  # ---- auto-approve --------------------------------------------------------

  def _on_auto_approve_toggled(self, checked):
    """Reflect the toggled state and propagate it.

    Reflects the toggled state in the button's label and propagates it
    via the public ``auto_approve_changed`` signal. The button is
    checkable; Qt's native pressed-in 'checked' rendering plus the
    label flip ('Auto-approve' -> 'Auto-approve: ON') signals the
    state without a colour override (a hardcoded colour would read
    poorly on one theme or the other).
    """
    self._auto_approve_btn.setText(
      "Auto-approve: ON" if checked else "Auto-approve")
    self.auto_approve_changed.emit(bool(checked))

  def set_auto_approve(self, on):
    """Programmatically sync the auto-approve toggle state.

    For tests and direct callers. Avoids a feedback loop by checking
    the button's state before re-toggling.
    """
    on = bool(on)
    if self._auto_approve_btn.isChecked() != on:
      self._auto_approve_btn.setChecked(on)

  # ---- send / stop ---------------------------------------------------------

  def click_send(self):
    if self._busy:
      self.stop.emit()
      return
    msg = self.text().strip()
    if not msg and not self._attachments:
      return
    atts = list(self._attachments)
    self._attachments = []
    self._refresh_chip_bar()
    self.send.emit(msg, atts)
    self.clear()
    # Record the sent text for Up/Down recall; reset navigation so the
    # next Up starts from this newest entry.
    if msg:
      self._history.append(msg)
    self._history_index = None
    self._history_draft = ""

  # ---- drag-drop / paste / file picker ------------------------------------

  def dragEnterEvent(self, event):
    if event.mimeData().hasUrls() or event.mimeData().hasImage():
      event.acceptProposedAction()

  def dropEvent(self, event):
    md = event.mimeData()
    if md.hasImage():
      self._handle_qimage(md.imageData())
    if md.hasUrls():
      for url in md.urls():
        if url.isLocalFile():
          self._handle_file(url.toLocalFile())
    event.acceptProposedAction()

  def _handle_qimage(self, qimage_or_obj):
    img = qimage_or_obj if isinstance(qimage_or_obj, QtGui.QImage) \
          else QtGui.QImage(qimage_or_obj)
    if img.isNull():
      return
    buf = QtCore.QBuffer()
    buf.open(QtCore.QBuffer.WriteOnly)
    img.save(buf, "PNG")
    self.attach_bytes(bytes(buf.data()), "image/png",
                      filename="pasted.png")

  def _handle_file(self, path):
    import mimetypes
    mime, _ = mimetypes.guess_type(path)
    if mime is None:
      return
    try:
      with open(path, "rb") as fh:
        data = fh.read()
    except OSError:
      return
    import os
    self.attach_bytes(data, mime, filename=os.path.basename(path))

  def _on_attach_clicked(self):
    paths, _ = QtWidgets.QFileDialog.getOpenFileNames(
      self, "Attach image", "",
      "Images (*.png *.jpg *.jpeg *.webp *.gif)")
    for p in paths:
      self._handle_file(p)

  def eventFilter(self, obj, event):
    if obj is self._edit and event.type() == QtCore.QEvent.KeyPress:
      key = event.key()
      mods = event.modifiers()
      if key in (QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter):
        if mods & (QtCore.Qt.ControlModifier | QtCore.Qt.MetaModifier):
          self.click_send()
          return True
      # Plain Up / Down recall the input history (shell-style), but only at
      # the first / last line so multi-line editing and Shift-selection are
      # left alone.
      nav_mods = mods & (QtCore.Qt.ShiftModifier | QtCore.Qt.ControlModifier
                         | QtCore.Qt.AltModifier | QtCore.Qt.MetaModifier)
      if key == QtCore.Qt.Key_Up and not nav_mods \
         and self._cursor_at_first_line():
        if self._history_recall_prev():
          return True
      if key == QtCore.Qt.Key_Down and not nav_mods \
         and self._cursor_at_last_line():
        if self._history_recall_next():
          return True
      if key == QtCore.Qt.Key_V and \
         (mods & (QtCore.Qt.ControlModifier | QtCore.Qt.MetaModifier)):
        cb = QtWidgets.QApplication.clipboard()
        md = cb.mimeData()
        if md is not None and md.hasImage():
          self._handle_qimage(md.imageData())
          return True
    return super().eventFilter(obj, event)
