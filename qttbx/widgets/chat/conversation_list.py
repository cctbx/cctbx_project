"""Sidebar listing of saved conversations."""

from qttbx.qt import QtCore, QtWidgets


class ConversationList(QtWidgets.QWidget):

  selected = QtCore.Signal(str)                    # conv_id
  new_requested = QtCore.Signal()
  delete_requested = QtCore.Signal(str)            # conv_id
  rename_requested = QtCore.Signal(str, str)       # conv_id, new_title

  def __init__(self, parent=None):
    super().__init__(parent)
    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(4, 4, 4, 4)
    self._list = QtWidgets.QListWidget(self)
    self._list.currentRowChanged.connect(self._on_row_changed)
    # Double-click or F2 starts an in-place rename editor; itemChanged
    # fires once the editor commits. The default trigger set on
    # QListWidget is NoEditTriggers, so we have to opt in explicitly
    # for the items' ItemIsEditable flag (set in set_conversations) to
    # take effect.
    self._list.setEditTriggers(
      QtWidgets.QAbstractItemView.DoubleClicked
      | QtWidgets.QAbstractItemView.EditKeyPressed)
    self._list.itemChanged.connect(self._on_item_changed)
    layout.addWidget(self._list, stretch=1)
    button_row = QtWidgets.QHBoxLayout()
    self._new_btn = QtWidgets.QPushButton("New", self)
    self._rename_btn = QtWidgets.QPushButton("Rename", self)
    self._del_btn = QtWidgets.QPushButton("Delete", self)
    self._new_btn.clicked.connect(self.click_new)
    self._rename_btn.clicked.connect(self.click_rename)
    self._del_btn.clicked.connect(self.click_delete)
    button_row.addWidget(self._new_btn)
    button_row.addWidget(self._rename_btn)
    button_row.addWidget(self._del_btn)
    button_row.addStretch(1)
    layout.addLayout(button_row)
    self._metas = []

  # ---- data ----------------------------------------------------------------

  def set_conversations(self, metas):
    self._metas = list(metas)
    self._list.blockSignals(True)
    self._list.clear()
    for m in self._metas:
      item = QtWidgets.QListWidgetItem(m.title or "Untitled")
      item.setToolTip("%s - %s" % (m.profile_name, m.model))
      item.setData(QtCore.Qt.UserRole, m.id)
      # ItemIsEditable lets the rename triggers (double-click / F2 /
      # the Rename button via editItem) open the in-place editor.
      item.setFlags(item.flags() | QtCore.Qt.ItemIsEditable)
      self._list.addItem(item)
    self._list.blockSignals(False)

  def select_index(self, i):
    if 0 <= i < self._list.count():
      self._list.setCurrentRow(i)

  def selected_id(self):
    item = self._list.currentItem()
    if item is None:
      return None
    return item.data(QtCore.Qt.UserRole)

  # ---- buttons -------------------------------------------------------------

  def click_new(self):
    self.new_requested.emit()

  def click_rename(self):
    """Start the in-place editor on the currently selected row. The
    actual rename is emitted from ``_on_item_changed`` once the user
    commits (Enter / focus-out)."""
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

  def _on_item_changed(self, item):
    """The list item's text changed -- either the user committed an
    in-place rename (Enter / focus-out) or set_conversations was
    called without blocking signals. The blockSignals wrapper around
    set_conversations means this slot only fires for real user
    commits."""
    cid = item.data(QtCore.Qt.UserRole)
    if not cid:
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
    # Update the cached meta so subsequent compares see the new state
    # without waiting for the chat window to round-trip a refresh.
    for m in self._metas:
      if m.id == cid:
        m.title = new_title
        break
    self.rename_requested.emit(cid, new_title)

  def _cached_title(self, cid):
    for m in self._metas:
      if m.id == cid:
        return m.title or ""
    return ""
