from PySide2.QtGui import QIcon
from PySide2.QtCore import Slot, Signal
from PySide2.QtWidgets import (
    QLabel,
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





class ToggleIconLabel(QLabel):
  # TODO: Make a widget to include spacers around label icons for centering
  stateChanged = Signal(bool)  # Define custom signal

  def __init__(self, checked_icon, unchecked_icon, *args, **kwargs):
    super(ToggleIconLabel, self).__init__(*args, **kwargs)
    self.checked_icon = checked_icon
    self.unchecked_icon = unchecked_icon
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
    icon = self.checked_icon if self._is_checked else self.unchecked_icon
    # Convert QIcon to QPixmap
    pixmap = icon.pixmap(64, 64)  # Adjust size as necessary

    # Save the pixmap to a file (or you can use another in-memory method)
    temp_image_path = "temp_icon.png"
    pixmap.save(temp_image_path)

    # Use setStyleSheet with the saved image path
    self.setStyleSheet(f'image: url("{temp_image_path}")')
