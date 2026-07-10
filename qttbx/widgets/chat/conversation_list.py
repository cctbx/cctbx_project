"""Sidebar listing of saved conversations."""

from qttbx.qt import QtCore, QtGui, QtWidgets


# Extra item-data role: True on a row whose conversation is open in another
# process. The delegate reads it to keep such a row dimmed even when selected --
# the default view paints a selected row's text in HighlightedText, which would
# otherwise hide the 'locked' cue.
_LOCKED_ROLE = QtCore.Qt.UserRole + 1


class _LockAwareDelegate(QtWidgets.QStyledItemDelegate):
  """Paints a locked conversation's name dimmed in EVERY state -- including when
  it is the selected row -- so selecting a locked conversation still shows it as
  locked."""

  def initStyleOption(self, option, index):
    super().initStyleOption(option, index)
    if index.data(_LOCKED_ROLE):
      pal = option.palette
      dim = pal.color(QtGui.QPalette.Disabled, QtGui.QPalette.Text)
      pal.setColor(QtGui.QPalette.Text, dim)
      pal.setColor(QtGui.QPalette.HighlightedText, dim)
      option.palette = pal


class ConversationList(QtWidgets.QWidget):
  """Sidebar widget listing saved conversations with new/rename/delete."""

  selected = QtCore.Signal(str)                    # conv_id (drives lock re-check)
  activated = QtCore.Signal(str)                    # conv_id (double-click / Enter)
  new_requested = QtCore.Signal()
  delete_requested = QtCore.Signal(str)            # conv_id
  rename_requested = QtCore.Signal(str, str)       # conv_id, new_title
  unlock_requested = QtCore.Signal(str)            # conv_id (clear a held lock)

  # After the shared Delete/Unlock button's label flips, it is disabled for this
  # many ms so a click aimed at the old label (e.g. the lock-poll relabelled it
  # under the cursor) drops harmlessly. See _debounce_action_button.
  _ACTION_DEBOUNCE_MS = 400

  def __init__(self, parent=None):
    super().__init__(parent)
    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(4, 4, 4, 4)
    self._list = QtWidgets.QListWidget(self)
    # A delegate that keeps a locked row dimmed even when it's the selected row.
    self._list.setItemDelegate(_LockAwareDelegate(self._list))
    self._list.currentRowChanged.connect(self._on_row_changed)
    # F2 starts an in-place rename editor; itemChanged fires once the
    # editor commits. The default trigger set on QListWidget is
    # NoEditTriggers, so we have to opt in explicitly for the items'
    # ItemIsEditable flag (set in set_conversations) to take effect.
    # NOT DoubleClicked: double-click activates a non-active conversation (rename
    # only the already-active one, via _on_item_double_clicked). F2 / the Rename
    # button still open the editor via EditKeyPressed + click_rename.
    self._list.setEditTriggers(
      QtWidgets.QAbstractItemView.EditKeyPressed)
    self._list.itemDoubleClicked.connect(self._on_item_double_clicked)
    self._list.installEventFilter(self)
    # itemDoubleClicked fires for ANY mouse button; filter the viewport to gate
    # activation/rename on a LEFT double-click only (see eventFilter).
    self._list.viewport().installEventFilter(self)
    self._suppress_next_return = False               # bounded post-rename guard
    self._list.itemDelegate().closeEditor.connect(
      self._arm_return_suppression)
    self._list.itemChanged.connect(self._on_item_changed)
    layout.addWidget(self._list, stretch=1)
    button_row = QtWidgets.QHBoxLayout()
    self._new_btn = QtWidgets.QPushButton("New", self)
    self._rename_btn = QtWidgets.QPushButton("Rename", self)
    self._del_btn = QtWidgets.QPushButton("Delete", self)
    self._new_btn.clicked.connect(self.click_new)
    self._rename_btn.clicked.connect(self.click_rename)
    self._del_btn.clicked.connect(self._on_del_button)
    # Reserve room for the wider of the two labels this button toggles between
    # (Delete / Unlock) so switching text never resizes the button -- and thus
    # the sidebar -- by a pixel.
    self._del_btn.setText("Unlock")
    _unlock_w = self._del_btn.sizeHint().width()
    self._del_btn.setText("Delete")
    self._del_btn.setMinimumWidth(
      max(self._del_btn.sizeHint().width(), _unlock_w))
    # Single-shot timer (child of self, so teardown can't fire it on a dead
    # button) that re-enables the Delete/Unlock button after a debounce.
    self._action_debounce = QtCore.QTimer(self)
    self._action_debounce.setSingleShot(True)
    self._action_debounce.setInterval(self._ACTION_DEBOUNCE_MS)
    self._action_debounce.timeout.connect(self._enable_action_button)
    button_row.addWidget(self._new_btn)
    button_row.addWidget(self._rename_btn)
    button_row.addWidget(self._del_btn)
    button_row.addStretch(1)
    layout.addLayout(button_row)
    self._metas = []
    # id -> ConversationMeta, rebuilt with self._metas in set_conversations, so
    # the per-row title/meta lookups (_cached_title, _set_row_locked, the rename
    # sites) are O(1) instead of each re-scanning self._metas. Holds the same
    # meta objects, so an in-place title mutation is visible through both.
    self._meta_by_id = {}
    self._active_id = None
    self._locked_ids = set()          # conversations open in another process
    # Fixed, small icon gutter so the checkmark can't shift the conversation
    # name: pin a small iconSize, then center the standard check on an EXACTLY
    # icon_px-square transparent canvas (the standard pixmap keeps its own
    # aspect, e.g. 11x12) so the checkmark and the inactive-row blank share the
    # same footprint -- the name starts at the same x either way, and the
    # checkmark is small (matched to the text). Built here (QApplication exists
    # once the widget is constructed).
    icon_px = 12
    self._list.setIconSize(QtCore.QSize(icon_px, icon_px))
    glyph = self.style().standardIcon(
      QtWidgets.QStyle.SP_DialogApplyButton).pixmap(
        QtCore.QSize(icon_px, icon_px))
    # glyph.width()/height() are DEVICE pixels (e.g. 24 for a 12px logical icon
    # on a 2x display); centre using the glyph's LOGICAL size and give the canvas
    # the same devicePixelRatio -- otherwise the offset goes negative on HiDPI
    # and the check is clipped into the corner (the dpr-1 offscreen test can't
    # catch it).
    dpr = glyph.devicePixelRatio() or 1.0
    gw, gh = glyph.width() / dpr, glyph.height() / dpr           # logical size
    check = QtGui.QPixmap(int(round(icon_px * dpr)), int(round(icon_px * dpr)))
    check.setDevicePixelRatio(dpr)
    check.fill(QtCore.Qt.transparent)
    painter = QtGui.QPainter(check)
    painter.drawPixmap(int(round((icon_px - gw) / 2)),
                       int(round((icon_px - gh) / 2)), glyph)
    painter.end()
    self._check_icon = QtGui.QIcon(check)
    blank = QtGui.QPixmap(int(round(icon_px * dpr)), int(round(icon_px * dpr)))
    blank.setDevicePixelRatio(dpr)
    blank.fill(QtCore.Qt.transparent)
    self._blank_icon = QtGui.QIcon(blank)
    # Locked rows are dimmed (this palette's disabled text colour) but stay
    # selectable, so a click can re-check the lock and offer Unlock.
    self._dim_brush = QtGui.QBrush(
      self.palette().color(QtGui.QPalette.Disabled, QtGui.QPalette.Text))
    self._normal_brush = QtGui.QBrush()      # invalid -> default palette text
    self._unlock_target = None               # cid the Delete->Unlock button acts on
    self._last_rename_revert = None          # (cid, pre-edit title) for revert

  # ---- data ----------------------------------------------------------------

  def set_conversations(self, metas, locked_ids=()):
    """Rebuild the list. ``locked_ids`` name conversations open in another
    phenix.chat process; they render DIMMED but stay selectable, so a click can
    re-check the lock (the chat window does this) and, if still held, turn the
    Delete button into Unlock. Locked rows are not renamable."""
    self._metas = list(metas)
    self._meta_by_id = {m.id: m for m in self._metas}
    self._locked_ids = set(locked_ids)
    self._list.blockSignals(True)
    self._list.clear()
    for m in self._metas:
      item = QtWidgets.QListWidgetItem(m.title or "Untitled")
      item.setData(QtCore.Qt.UserRole, m.id)
      item.setIcon(self._check_icon if m.id == self._active_id
                   else self._blank_icon)
      self._apply_lock_style(item, m, m.id in self._locked_ids)
      self._list.addItem(item)
    self._list.blockSignals(False)
    self.set_delete_button()                 # caller re-establishes the selection

  def _apply_lock_style(self, item, meta, locked):
    """Style one row for its lock state: locked -> dimmed, selectable, NOT
    editable; unlocked -> normal text, editable (so the rename triggers work).
    The _LOCKED_ROLE flag lets the delegate keep a locked row dimmed even when it
    is the selected row."""
    item.setData(_LOCKED_ROLE, locked)
    if locked:
      item.setToolTip("Open in another phenix.chat window (select it to unlock)")
      item.setForeground(self._dim_brush)
      item.setFlags((item.flags() | QtCore.Qt.ItemIsSelectable
                     | QtCore.Qt.ItemIsEnabled) & ~QtCore.Qt.ItemIsEditable)
    else:
      if meta is not None:
        item.setToolTip("%s - %s" % (meta.profile_name, meta.model))
      item.setForeground(self._normal_brush)
      item.setFlags(item.flags() | QtCore.Qt.ItemIsEditable)

  def _item_for(self, cid):
    """The QListWidgetItem whose conversation id is ``cid``, or None."""
    for i in range(self._list.count()):
      item = self._list.item(i)
      if item.data(QtCore.Qt.UserRole) == cid:
        return item
    return None

  def _set_row_locked(self, cid, locked):
    """Shared body of mark_locked / mark_unlocked: flip ``cid``'s locked styling
    and membership to ``locked``. No-op if the row isn't listed or is already in
    that state. When locking, close any in-place editor open on the row so its
    commit can't fire a rename for a now-locked conversation."""
    if (cid in self._locked_ids) == locked:
      return
    item = self._item_for(cid)
    if item is None:
      return
    if locked:
      self._locked_ids.add(cid)
      self._close_editor(item)
    else:
      self._locked_ids.discard(cid)
    meta = self._meta_by_id.get(cid)
    self._list.blockSignals(True)
    self._apply_lock_style(item, meta, locked)
    self._list.blockSignals(False)

  def mark_unlocked(self, cid):
    """Re-style ``cid`` as no longer locked (normal text, editable) and drop it
    from the locked set. No-op if it wasn't locked. Called by the chat window
    when a re-check finds the lock gone, or right after the user unlocks it."""
    self._set_row_locked(cid, False)

  def mark_locked(self, cid):
    """Re-style ``cid`` as locked (dimmed, not editable) and add it to the locked
    set. No-op if already locked. Called by the chat window when an on-select or
    poll re-check finds a conversation now open in another process, so it dims
    immediately and stays dimmed."""
    self._set_row_locked(cid, True)

  def _close_editor(self, item):
    """Close any in-place rename editor open on ``item`` (e.g. the lock-poll just
    dimmed the row mid-rename). closePersistentEditor also dismisses the
    transient editItem editor, discarding the pending edit -- otherwise its
    commit would still fire a rename for the now-locked conversation."""
    if self._list.state() == QtWidgets.QAbstractItemView.EditingState \
        and self._list.currentItem() is item:
      self._list.closePersistentEditor(item)

  def listed_ids(self):
    """The conv_ids currently shown as rows, in display order (used by the chat
    window's lock-poll to refresh dimming)."""
    return [self._list.item(i).data(QtCore.Qt.UserRole)
            for i in range(self._list.count())]

  def is_row_locked(self, cid):
    """True if row ``cid`` is currently rendered locked (dimmed) -- lets the poll
    sync the Delete/Unlock button without re-reading the marker it just read."""
    return cid in self._locked_ids

  def set_unlock_button(self, cid):
    """Turn the Delete button into Unlock, acting on ``cid``, and disable Rename
    (a locked conversation can't be renamed). Relabels ONLY on a real transition
    -- a repeated poll in the same state must not keep re-arming the debounce."""
    if self._unlock_target == cid:
      return
    self._unlock_target = cid
    self._del_btn.setText("Unlock")
    self._rename_btn.setEnabled(False)
    self._debounce_action_button()

  def set_delete_button(self):
    """Restore the Delete button (clears any Unlock target) and re-enable Rename.
    Relabels ONLY on a real transition (see set_unlock_button)."""
    if self._unlock_target is None:
      return
    self._unlock_target = None
    self._del_btn.setText("Delete")
    self._rename_btn.setEnabled(True)
    self._debounce_action_button()

  def _debounce_action_button(self):
    """Disable the Delete/Unlock button for _ACTION_DEBOUNCE_MS after its label
    flips, so a click aimed at the old label (e.g. the 2s lock-poll relabelled it
    under the cursor) drops harmlessly instead of firing the other action."""
    self._del_btn.setEnabled(False)
    self._action_debounce.start()

  def _enable_action_button(self):
    self._del_btn.setEnabled(True)

  def _on_del_button(self):
    """The shared Delete/Unlock button: Unlock emits unlock_requested for its
    target; Delete deletes the current selection."""
    if self._unlock_target is not None:
      self.unlock_requested.emit(self._unlock_target)
    else:
      self.click_delete()

  def title_for(self, cid):
    """The display title of row ``cid`` (its list text), or None if not listed.
    Lets the delete confirmation name the conversation."""
    item = self._item_for(cid)
    return item.text() if item is not None else None

  def select_index(self, i):
    """Select row ``i`` (used by the widget tests to drive a selection)."""
    if 0 <= i < self._list.count():
      self._list.setCurrentRow(i)

  def set_active(self, cid):
    """Mark the row whose id == ``cid`` as the active conversation (checkmark),
    clearing the marker from the others. Icon writes are wrapped in
    ``blockSignals`` -- ``setIcon`` is a ``DecorationRole`` data change that
    fires ``itemChanged``, which would otherwise emit a spurious
    ``rename_requested`` for an empty-title row (see the widget test).

    A no-op when ``cid`` is already the active id: ``set_conversations`` paints
    each row's icon from ``self._active_id`` and ``set_active`` is the only other
    icon writer, so the icons already reflect ``cid`` -- and the common caller
    (``_populate_sidebar`` on every turn-end rebuild, active unchanged) would
    otherwise re-run the O(N) icon pass to write identical icons."""
    if cid == self._active_id:
      return
    self._active_id = cid
    self._list.blockSignals(True)
    for i in range(self._list.count()):
      item = self._list.item(i)
      active = item.data(QtCore.Qt.UserRole) == cid
      item.setIcon(self._check_icon if active else self._blank_icon)
    self._list.blockSignals(False)

  def selected_id(self):
    item = self._list.currentItem()
    if item is None:
      return None
    return item.data(QtCore.Qt.UserRole)

  # ---- buttons -------------------------------------------------------------

  def click_new(self):
    self.new_requested.emit()

  def click_rename(self):
    """Start the in-place editor on the currently selected row.

    The actual rename is emitted from ``_on_item_changed`` once the
    user commits (Enter / focus-out).
    """
    item = self._list.currentItem()
    if item is None:
      return
    self._list.editItem(item)

  def click_delete(self):
    cid = self.selected_id()
    if cid:
      self.delete_requested.emit(cid)

  # ---- internal ------------------------------------------------------------

  def _on_row_changed(self, row):
    if row < 0:
      return
    item = self._list.item(row)
    if item is None:
      return
    self.selected.emit(item.data(QtCore.Qt.UserRole))

  def _on_item_double_clicked(self, item):
    """Double-click the ACTIVE conversation -> rename it (nothing to switch to);
    double-click any other UNLOCKED one -> activate (switch to) it. A locked row
    does nothing -- unlock it first via the Delete->Unlock button."""
    cid = item.data(QtCore.Qt.UserRole)
    if cid is None or cid in self._locked_ids:
      return
    if cid == self._active_id:
      self._list.editItem(item)
    else:
      self.activated.emit(cid)

  def select_id(self, cid):
    """Select the row whose id == ``cid`` (no signals; a switch is driven by
    ``activated``, not selection). A ``cid`` not present is a no-op."""
    for i in range(self._list.count()):
      if self._list.item(i).data(QtCore.Qt.UserRole) == cid:
        self._list.blockSignals(True)
        self._list.setCurrentRow(i)
        self._list.blockSignals(False)
        return

  def eventFilter(self, obj, event):
    """Enter/Return on the list activates a non-active selection (parallel to
    double-click) and is a no-op on the active one -- always consumed so it
    never falls through to QAbstractItemView, which on macOS maps Return to
    'start editing'. A NON-left double-click on the viewport is swallowed so only
    a left double-click activates/renames (itemDoubleClicked is button-agnostic,
    so a right/middle double-click would otherwise switch or open the editor)."""
    if obj is self._list.viewport() \
        and event.type() == QtCore.QEvent.MouseButtonDblClick \
        and event.button() != QtCore.Qt.LeftButton:
      return True                                    # ignore non-left dbl-click
    if obj is self._list and event.type() == QtCore.QEvent.KeyPress \
        and event.key() in (QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter):
      # Primary guard: while the editor is open, let it commit.
      if self._list.state() == QtWidgets.QAbstractItemView.EditingState:
        return False
      # Backup guard: swallow the one Return that some Qt bindings propagate to
      # the view after an edit commits (arms on closeEditor; auto-clears).
      if self._suppress_next_return:
        self._suppress_next_return = False
        return True
      cid = self.selected_id()
      if (cid is not None and cid != self._active_id
          and cid not in self._locked_ids):
        self.activated.emit(cid)                     # locked -> unlock first
      return True                                    # consume (activate or no-op)
    return super().eventFilter(obj, event)

  def _arm_return_suppression(self, *args):
    """Arm the post-rename Return guard, bounded to the current event-loop turn.
    closeEditor also fires on focus-out/Escape (no Return follows), so the flag
    must self-clear or a later genuine Return would be silently eaten."""
    self._suppress_next_return = True
    QtCore.QTimer.singleShot(0, self._clear_return_suppression)

  def _clear_return_suppression(self):
    self._suppress_next_return = False

  def _on_item_changed(self, item):
    """Handle a committed in-place rename of a list item.

    The list item's text changed -- either the user committed an
    in-place rename (Enter / focus-out) or ``set_conversations`` was
    called without blocking signals. The ``blockSignals`` wrapper
    around ``set_conversations`` means this slot only fires for real
    user commits; ``set_active`` and ``select_id`` rely on the same
    ``blockSignals`` invariant to keep their programmatic item changes
    from re-entering this slot.
    """
    cid = item.data(QtCore.Qt.UserRole)
    if not cid:
      return
    if cid in self._locked_ids:
      # A poll/select can lock a row mid-rename. _close_editor normally dismisses
      # the editor first; guard the commit too so a locked conversation is never
      # renamed -- revert the row to its stored title.
      self._list.blockSignals(True)
      item.setText(self._cached_title(cid) or "Untitled")
      self._list.blockSignals(False)
      return
    new_title = (item.text() or "").strip()
    old_title = self._cached_title(cid)
    if not new_title:
      # Reject empty: revert in place without re-firing this slot.
      self._list.blockSignals(True)
      item.setText(old_title or "Untitled")
      self._list.blockSignals(False)
      return
    if new_title == old_title:
      return
    # Remember the pre-edit title so a chat-window refusal (e.g. the row locked
    # mid-edit) can revert this row without a full rebuild -- see
    # revert_last_rename. Then optimistically update the cache so subsequent
    # compares see the new state without waiting for a round-trip refresh.
    self._last_rename_revert = (cid, old_title)
    m = self._meta_by_id.get(cid)
    if m is not None:
      m.title = new_title
    self.rename_requested.emit(cid, new_title)

  def revert_last_rename(self, cid):
    """Undo the most recent in-place rename of ``cid`` back to its pre-edit title.
    The chat window calls this when it refuses the rename (e.g. the row locked
    mid-edit) so the optimistic title doesn't linger even when a full rebuild is
    unavailable -- a transient list failure would skip one."""
    saved = getattr(self, "_last_rename_revert", None)
    if not saved or saved[0] != cid:
      return
    old_title = saved[1]
    item = self._item_for(cid)
    if item is None:
      return
    self._list.blockSignals(True)
    item.setText(old_title or "Untitled")
    self._list.blockSignals(False)
    m = self._meta_by_id.get(cid)
    if m is not None:
      m.title = old_title
    self._last_rename_revert = None

  def _cached_title(self, cid):
    m = self._meta_by_id.get(cid)
    return (m.title or "") if m is not None else ""
