from PySide2.QtCore import Qt
from PySide2.QtWidgets import  QLineEdit


class HistoryLineEdit(QLineEdit):
  def __init__(self):
    super(HistoryLineEdit, self).__init__()
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
      super(HistoryLineEdit, self).keyPressEvent(event)

