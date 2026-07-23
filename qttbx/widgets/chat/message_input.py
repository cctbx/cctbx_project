"""Multi-line chat input with attachments, drag-drop, paste, and buttons.

The button row holds Save chat / Auto-approve / search / attach / Send.

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


class _DropTextEdit(QtWidgets.QPlainTextEdit):
  """Entry box that hands dropped/pasted files and images to its owner.

  Qt delivers drag-and-drop to a QPlainTextEdit's viewport, so the parent
  ``MessageInput.dropEvent`` never fires for a drop on the visible text
  area -- the default handling would insert the raw ``file://`` URI as
  text. Both drop and paste funnel through ``insertFromMimeData``, so
  intercepting it here is the single chokepoint that covers both: a
  URL/image payload is routed to the owner (becoming an attachment), while
  ordinary text falls through to the normal editor behaviour.
  """

  def __init__(self, owner):
    super().__init__(owner)
    self._owner = owner

  def insertFromMimeData(self, source):
    # Intercept ONLY payloads that actually become attachments: an image, or
    # urls that include a LOCAL file. A remote http(s) link (hasUrls() but no
    # local file) must paste as text, not be swallowed -- so it falls through
    # to the normal editor behaviour.
    if source is not None and (
        source.hasImage()
        or (source.hasUrls()
            and any(u.isLocalFile() for u in source.urls()))):
      self._owner._handle_dropped_mime(source)
      return
    super().insertFromMimeData(source)


class MessageInput(QtWidgets.QWidget):
  """Multi-line chat input with attachments and a button row."""

  send = QtCore.Signal(str, list)                  # text, attachments (list)
  stop = QtCore.Signal()
  attachment_rejected = QtCore.Signal(str)
  save_chat = QtCore.Signal()                      # 'Save chat' button click
  auto_approve_changed = QtCore.Signal(bool)       # checked state
  search_clicked = QtCore.Signal()                 # 🔍 button click

  # Idle placeholder text. The assistant name defaults to "Claude" but is
  # rewritten per session by ChatWindow via set_assistant_name() so the box
  # names the active backend (GPT / Gemini / …). ChatWindow also swaps in a
  # cycling 'Thinking...' placeholder while a turn is in flight. turn_done
  # now carries the full TurnDone event; only a TERMINAL finish returns the
  # idle placeholder and unlocks the composer.
  _PLACEHOLDER_FMT = "Message %s...  (Ctrl/Cmd+Enter to send)"
  DEFAULT_PLACEHOLDER = _PLACEHOLDER_FMT % "Claude"

  _MAX_IMAGE_BYTES = 20 * 1024 * 1024
  _ALLOWED_MIMES = frozenset((
    "image/png", "image/jpeg", "image/webp", "image/gif"))

  def __init__(self, parent=None):
    super().__init__(parent)
    self._busy = False
    # Optional 0-arg callable -> bool consulted by click_send BEFORE the draft
    # is consumed. Returning False refuses the send and KEEPS the typed text
    # and attachments in place: the `send` signal is fire-and-forget (the
    # widget clears itself right after emitting), so a host window that can
    # reject a send (e.g. a turn still in flight) must gate it here -- a
    # rejection inside its slot would come after the draft was already wiped.
    # The gate owns any user feedback (status banner).
    self.send_gate = None
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
    self._edit = _DropTextEdit(self)
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
    # Cache the last-applied placeholder text + dim so set_placeholder can skip
    # the palette re-polish (dim unchanged) and the viewport repaint (text
    # unchanged) -- the in-flight 'Thinking...' spinner drives it ~8x/s. Seeded
    # to the dim idle text set above.
    self._placeholder_text = self._idle_placeholder
    self._placeholder_dim = True
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
    self._search_btn = QtWidgets.QToolButton(self)
    self._search_btn.setText("🔍")
    self._search_btn.setToolTip("Search conversation (Ctrl+F / ⌘F)")
    self._search_btn.clicked.connect(self._on_search_clicked)
    button_row.addWidget(self._search_btn)
    self._attach_btn = QtWidgets.QToolButton(self)
    self._attach_btn.setText("@")    # ASCII-safe; UI later replaces with icon
    self._attach_btn.setToolTip("Attach an image")
    self._attach_btn.clicked.connect(self._on_attach_clicked)
    button_row.addWidget(self._attach_btn)
    self._button = QtWidgets.QPushButton("Send", self)
    self._button.clicked.connect(self.click_send)
    button_row.addWidget(self._button)
    layout.addLayout(button_row)
    # Focusing the composite lands in the editor -- the search bar's
    # close path refocuses the input via setFocus().
    self.setFocusProxy(self._edit)

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
    text = text or ""
    changed = False
    # Only rebuild the palette when the dim state actually changed: the palette
    # depends solely on `dim`, so re-applying it (a style re-polish) on every
    # spinner tick -- which holds dim constant -- is pure waste.
    if dim != self._placeholder_dim:
      palette = self._edit.palette()
      target_color = (
        self._dim_placeholder_color if dim
        else QtWidgets.QApplication.palette().color(QtGui.QPalette.Text))
      palette.setColor(QtGui.QPalette.PlaceholderText, target_color)
      self._edit.setPalette(palette)
      self._placeholder_dim = dim
      changed = True
    # Only re-set the text (and repaint) when it actually changed.
    if text != self._placeholder_text:
      self._edit.setPlaceholderText(text)
      self._placeholder_text = text
      changed = True
    if changed:
      # QPlainTextEdit only repaints the cursor area on placeholder change
      # when focused, so the user sees a stale tail (just the first letter of
      # the new verb updates). Force the viewport to repaint the whole
      # placeholder region.
      self._edit.viewport().update()

  def changeEvent(self, event):
    """Re-resolve the placeholder colour on an app theme (palette) switch.

    set_placeholder caches the applied dim state and skips re-writing the
    explicitly-set PlaceholderText palette role while it's unchanged. That role
    does not follow a theme switch on its own, so without this a 'Thinking...'
    cue begun before a light->dark switch would keep its old (near-black) colour
    and render dark-on-dark for the rest of the turn. Invalidate the cached dim
    state and re-apply against the new palette (a rare event, so the ~8x/s
    fast-path stays cheap)."""
    super().changeEvent(event)
    et = event.type()
    # Guard construction: a PaletteChange can arrive before _build_ui seeds the
    # placeholder cache.
    if et in (QtCore.QEvent.PaletteChange, QtCore.QEvent.ThemeChange) \
       and hasattr(self, "_placeholder_text"):
      self._dim_placeholder_color = (
        QtWidgets.QApplication.palette().color(QtGui.QPalette.PlaceholderText))
      dim = self._placeholder_dim
      self._placeholder_dim = None                 # force set_placeholder to re-apply
      self.set_placeholder(
        self._placeholder_text, dim=True if dim is None else dim)

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
    if not self._busy:
      self.set_placeholder(self._idle_placeholder, dim=True)

  # ---- busy / button state -------------------------------------------------

  def set_busy(self, busy):
    self._busy = bool(busy)
    self._button.setText("Stop" if self._busy else "Send")

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
    # _maybe_resample re-encodes a non-PNG type (webp/gif/jpeg) to JPEG when
    # it shrinks, so it returns the mime that matches the bytes it hands back.
    # The attachment must carry that mime -- shipping webp/gif bytes that are
    # really JPEG makes the provider reject the request (400).
    data, mime = self._maybe_resample(data, mime)
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
    iterations. PNG re-encodes as PNG; every other allowed type
    (jpeg/webp/gif) re-encodes as JPEG, so the returned mime is
    ``"image/jpeg"`` in that case -- the caller must advertise the bytes'
    real format or the provider rejects the attachment.

    Parameters
    ----------
    data : bytes
        Raw image bytes.
    mime : str
        MIME type; selects the re-encode format (PNG vs JPG).

    Returns
    -------
    (bytes, str)
        The (possibly smaller) bytes and the mime matching them. The
        original bytes + original mime are returned unchanged when the
        input is already under the cap or Qt can't read it (which lets the
        size check catch an unreadable oversized blob later).
    """
    if len(data) <= self._max_image_bytes:
      return data, mime
    img = QtGui.QImage()
    if not img.loadFromData(data):
      return data, mime
    if mime == "image/png":
      target_fmt, out_mime = "PNG", "image/png"
    else:
      target_fmt, out_mime = "JPG", "image/jpeg"
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
        return out, out_mime
    return out, out_mime

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

  def _on_search_clicked(self, _checked=False):
    self.search_clicked.emit()

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
    if self.send_gate is not None and not self.send_gate():
      # Refused (e.g. a turn in flight): the draft stays in the composer.
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
    self._handle_dropped_mime(event.mimeData())
    event.acceptProposedAction()

  def _handle_dropped_mime(self, md):
    """Route a dropped/pasted image or local-file payload to attachments.

    Shared by the entry box (:class:`_DropTextEdit`, which intercepts
    drops/pastes on the visible text area) and this widget's own
    ``dropEvent`` (drops on the chip bar / button row). A dropped file
    becomes an attachment instead of a raw ``file://`` URI in the prompt;
    non-local URLs are ignored rather than inserted as text.
    """
    if md is None:
      return
    # elif, not a second if: a payload that carries BOTH image data and a
    # local file:// URL (some apps/clipboards populate both) must yield ONE
    # attachment. Two independent ifs attached the image twice -- once as
    # 'pasted.png' here, once by reading the dropped file below.
    if md.hasImage():
      self._handle_qimage(md.imageData())
    elif md.hasUrls():
      for url in md.urls():
        if url.isLocalFile():
          self._handle_file(url.toLocalFile())

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
      # Extension mimetypes can't map -> we can't confirm an allowed image
      # type, so reject via the same attachment_rejected path a known-
      # unsupported mime uses rather than returning silently (no signal,
      # no status).
      self.attachment_rejected.emit(
        "Could not determine attachment type: %s" % path)
      return
    # Validate the guessed mime BEFORE opening/reading the file. Reading first
    # (as this method used to) meant a dropped multi-GB file with a non-None
    # but disallowed mime (e.g. video/quicktime) was pulled entirely into
    # memory on the GUI thread -- freezing the event loop, risking MemoryError
    # -- only for attach_bytes to reject it afterwards. Reject here instead,
    # via the same attachment_rejected path attach_bytes uses for a disallowed
    # mime, so behaviour is unchanged for the user -- just moved earlier, and
    # the bytes are never read.
    if mime not in self._ALLOWED_MIMES:
      self.attachment_rejected.emit(
        "Unsupported attachment type: %s" % mime)
      return
    try:
      with open(path, "rb") as fh:
        data = fh.read()
    except OSError as exc:
      # Surface the read failure via attachment_rejected instead of
      # swallowing it silently; carry the OSError detail so it isn't
      # discarded bare.
      self.attachment_rejected.emit(
        "Could not read attachment %s: %s" % (path, exc))
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
    return super().eventFilter(obj, event)
