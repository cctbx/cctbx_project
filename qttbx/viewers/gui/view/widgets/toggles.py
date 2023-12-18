from PySide2.QtWidgets import QPushButton, QLabel
from PySide2.QtGui import QIcon
from PySide2.QtCore import Slot, Signal

class ToggleIconButton(QPushButton):
  def __init__(self, on_path, off_path, parent=None):
    super().__init__(parent)
    self.on_icon = QIcon(str(on_path))
    self.off_icon = QIcon(str(off_path))
    self.current_icon = self.on_icon
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


class ToggleIconLabel(QLabel):
  # TODO: Make a widget to include spacers around label icons for centering
  stateChanged = Signal(bool)  # Define custom signal

  def __init__(self, checked_icon_path, unchecked_icon_path, *args, **kwargs):
    super(ToggleIconLabel, self).__init__(*args, **kwargs)
    self.checked_icon_path = checked_icon_path
    self.unchecked_icon_path = unchecked_icon_path
    self._is_checked = False
    self._is_destroyed = False
    self.update_icon()

  @property
  def is_checked(self):
    return self._is_checked

  @is_checked.setter
  def is_checked(self, value):
    if not self.is_destroyed:
      assert isinstance(value, bool), "is_checked must be boolean"
      self._is_checked = value
      self.update_icon()

  @property
  def is_destroyed(self):
    return self._is_destroyed

  @is_destroyed.setter
  def is_destroyed(self,value):
    self._is_destroyed = value

  def mousePressEvent(self, event):
    self.is_checked = not self.is_checked
    self.stateChanged.emit(self.is_checked)

  def update_icon(self):
    icon_path = self.checked_icon_path if self._is_checked else self.unchecked_icon_path
    self.setStyleSheet(f'image: url("{icon_path}")')
