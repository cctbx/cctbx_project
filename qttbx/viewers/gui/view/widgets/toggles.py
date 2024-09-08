from PySide2.QtGui import QIcon
from PySide2.QtCore import Slot, Signal
from PySide2.QtWidgets import (
    QPushButton
)

class ToggleIconButton(QPushButton):
  def __init__(self, off_icon, on_icon, parent=None):
    super().__init__(parent)
    self.on_icon = on_icon
    self.off_icon = off_icon
    self.current_icon = self.off_icon
    self.setIcon(self.current_icon)
    self.clicked.connect(self.toggle_icon)

  

  @property
  def is_on(self):
    if self.current_icon == self.on_icon:
      return True
    else:
      return False

  @Slot()
  def toggle_icon(self):
    if self.current_icon == self.on_icon:
      self.current_icon = self.off_icon
    else:
      self.current_icon = self.on_icon
    self.setIcon(self.current_icon)

