from PySide2.QtGui import QIcon
from PySide2.QtCore import Slot, Signal
from PySide2.QtCore import QSize
from PySide2.QtWidgets import (
    QLabel,
    QPushButton,
)


class ToggleIconButton(QPushButton):
  def __init__(self, on_icon, off_icon, parent=None):
    super().__init__(parent)
    self.on_icon = on_icon
    self.off_icon = off_icon
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

  def __init__(self, checked_icon, unchecked_icon, *args, **kwargs):
    super().__init__(*args, **kwargs)
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
    self.set_icon(icon)

  def set_icon(self, icon):

    # Get the widget size and calculate the ideal icon size
    ideal_size = self.sizeHint()
    if not ideal_size.isValid() or ideal_size.isEmpty():
      ideal_size = QSize(64, 64)  # Default fallback size

    # Get the best-matching pixmap size from the QIcon
    pixmap = icon.pixmap(ideal_size)

    # Set the QPixmap on the QLabel
    self.setPixmap(pixmap)

  def resizeEvent(self, event):
    """Handle resizing of the widget and adjust icon size dynamically."""
    super(MyCustomWidget, self).resizeEvent(event)
    self.update_icon_size()

  def update_icon_size(self):
    """Update the icon size when the widget is resized."""
    icon = self.icon_label.pixmap()
    if icon:
      ideal_size = self.sizeHint()