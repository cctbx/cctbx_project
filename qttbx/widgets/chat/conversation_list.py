"""Sidebar listing of saved conversations."""

from qttbx.qt import QtCore, QtWidgets


class ConversationList(QtWidgets.QWidget):

  selected = QtCore.Signal(str)                    # conv_id
  new_requested = QtCore.Signal()
  delete_requested = QtCore.Signal(str)            # conv_id

  def __init__(self, parent=None):
    super().__init__(parent)
    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(4, 4, 4, 4)
    self._list = QtWidgets.QListWidget(self)
    self._list.currentRowChanged.connect(self._on_row_changed)
    layout.addWidget(self._list, stretch=1)
    button_row = QtWidgets.QHBoxLayout()
    self._new_btn = QtWidgets.QPushButton("New", self)
    self._del_btn = QtWidgets.QPushButton("Delete", self)
    self._new_btn.clicked.connect(self.click_new)
    self._del_btn.clicked.connect(self.click_delete)
    button_row.addWidget(self._new_btn)
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
