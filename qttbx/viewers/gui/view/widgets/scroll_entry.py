from pathlib import Path

from PySide2.QtCore import Qt, QSize
from PySide2.QtWidgets import QWidget, QHBoxLayout
from PySide2.QtGui import QPalette, QPainter

from .toggles import ToggleIconLabel
from .editable_label import EditableLabel


class ScrollEntryView(QWidget):
  def __init__(self,parent=None,active_toggle=True):
    super().__init__(parent)
    # state
    self.parent_explicit = parent
    self._is_destroyed = False


    # Start Components
    self.layout = QHBoxLayout(self)
    self.setLayout(self.layout)


    # Styling
    self.palette = self.palette()
    # self.color = self.palette.color(QPalette.Background)
    # self.palette.setColor(QPalette.Background, self.color)
    self.setPalette(self.palette)
    self.setAutoFillBackground(True)
    self.setMinimumHeight(42)
    self.setMaximumHeight(42)
    self.layout.setContentsMargins(5, 5, 5, 5)  # left, top, right, bottom
    #self.layout.setAlignment(Qt.AlignVCenter)
    self._all_button_height = 32
    self._all_button_width = 48




    # Active toggle
    self.active_toggle = None
    if active_toggle:
      icon_path_checked = Path(__file__).parent / '../assets/icons/material/radio_checked.svg'
      icon_path_unchecked = Path(__file__).parent / '../assets/icons/material/radio_unchecked.svg'
      self.active_toggle = ToggleIconLabel(str(icon_path_checked),str(icon_path_unchecked))
      scale = 0.6
      self.active_toggle.setMaximumSize(QSize(self.height()*scale,self.height()*scale))
      self.active_toggle.setMinimumSize(QSize(self.height()*scale,self.height()*scale))
      self.active_toggle.setToolTip("Is active?")
      self.layout.addWidget(self.active_toggle)

    # Name Label Widget
    name = "entry name"
    tooltip = 'entry tooltip'
    self.visible_name = self._truncate_string(name,max_len=30).ljust(30)

    self.label_name = EditableLabel(self,self.visible_name)
    # self.label_name.setMinimumHeight(self.height()*0.3)
    # self.label_name.setMaximumHeight(self.height()*0.3)
    #self.label_name.resize(300, self.height()*0.6)  # width, height
    self.label_name.setToolTip(tooltip)
    self.layout.addWidget(self.label_name)


    # Add stretch to pushsubsequent widgets to the right
    self.layout.addStretch()



  @property
  def is_destroyed(self):
    return self._is_destroyed

  @is_destroyed.setter
  def is_destroyed(self,value):
    self._is_destroyed = value
    if self.active_toggle:
      self.active_toggle.is_destroyed = True




  def paintEvent(self, event):
    self.color = self.palette.color(QPalette.Background)
    self.palette.setColor(QPalette.Background, self.color)
    painter = QPainter(self)
    painter.setRenderHint(QPainter.Antialiasing)
    painter.setPen(Qt.NoPen)
    painter.setBrush(self.color)
    # Draw a rounded rectangle
    painter.drawRoundedRect(self.rect(), 7, 7)  # The last two arguments specify the x and y radius


  def _truncate_string(self,path, max_len=50):
    if len(path) > max_len:
      return path[:max_len // 2] + "..." + path[-max_len // 2:]
    else:
      return path
