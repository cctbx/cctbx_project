"""
Implements a line edit widget that maintains a history of entered commands
"""
from PySide2.QtWidgets import QLineEdit
from PySide2.QtCore import Qt


class HistoryLineEdit(QLineEdit):
  def __init__(self):
    super().__init__()
    self.history = []
    self.index = -1

  def keyPressEvent(self, event):
    if event.key() == Qt.Key_Up:
      self.index = min(self.index + 1, len(self.history) - 1)
      if self.index >= 0:
        self.setText(self.history[self.index])
    elif event.key() == Qt.Key_Down:
      self.index = max(self.index - 1, -1)
      if self.index == -1:
        self.clear()
      else:
        self.setText(self.history[self.index])
    else:
      super().keyPressEvent(event)

