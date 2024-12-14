from pathlib import Path
import qtawesome as qta
from PySide2.QtGui import QIcon
from PySide2.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton
)

from .scroll_list import ScrollableListView, ScrollEntryView
from .widgets.representation_select import RepresentationSelect
from .widgets.toggles import ToggleIconButton
from .widgets.checkbox import ConditionalCheckBox

class GeoCheckBox(ConditionalCheckBox):
  def __init__(self, label, parent=None):
    super().__init__(label, parent)

  def can_toggle(self):
    return False # just disable

class ModelLikeEntryView(ScrollEntryView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)

    # Create the vertical separator
    separator = QFrame()
    separator.setFrameShape(QFrame.VLine)
    separator.setFrameShadow(QFrame.Sunken)
    self.layout.addWidget(separator)


     # Visibility
    on_icon = qta.icon("mdi.eye")
    off_icon = qta.icon("mdi.eye-off")
    self.button_viz = ToggleIconButton(on_icon, off_icon, parent=self)
    self.button_viz.setToolTip("Toggle visibility")
    self.button_viz.setFixedSize(self._all_button_width,self._all_button_height)
    #button_color.setContentsMargins(0,0,0,0)
    self.layout.addWidget(self.button_viz)



    # Color theme widget # TODO: fix this
    # self.button_theme = ColorThemeButton(parent=self)
    # self.button_theme.button.setFixedSize(self._all_button_width,self._all_button_height)
    # self.layout.addWidget(self.button_theme)

    # # Representations
    # self.button_rep = RepresentationSelect(parent=self)
    # self.button_rep.button.setFixedSize(self._all_button_width,self._all_button_height)
    # self.layout.addWidget(self.button_rep)

    # Color picking widget
    self.button_color = QPushButton()
    icon = qta.icon("mdi.format-color-fill")
    self.button_color.setIcon(icon)
    self.button_color.setToolTip("Color fill")
    self.button_color.setFixedSize(self._all_button_width,self._all_button_height)
    #button_color.setContentsMargins(0,0,0,0)
    #button_color.setMaximumSize(QSize(maxs2,maxs2))
    self.layout.addWidget(self.button_color)

    # Create the second vertical separator
    separator2 = QFrame()
    separator2.setFrameShape(QFrame.VLine)
    separator2.setFrameShadow(QFrame.Sunken)
    self.layout.addWidget(separator2)

    # Close
    self.button_close = QPushButton()
    icon = qta.icon("mdi.close")
    self.button_close.setIcon(icon)
    self.button_close.setToolTip("Remove")
    self.button_close.setFixedSize(self._all_button_width,self._all_button_height)
    self.layout.addWidget(self.button_close)


    self._insert_index = 2 # a hint on where to insert widgets for subclasses. From back


    

class ModelEntryView(ModelLikeEntryView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)


class ModelListView(ScrollableListView):
  def __init__(self,parent=None,title="Models"):
    super().__init__(parent=parent)
    header_layout = QHBoxLayout()
    label = QLabel(title)
    current_font = label.font()
    current_font.setPointSize(16)
    current_font.setBold(False)
    label.setFont(current_font)

    self.load_button = QPushButton()
    load_icon = qta.icon("mdi.plus")
    self.load_button.setIcon(load_icon)
    self.load_button.setMaximumSize(50, 50)
    self.load_button.setContentsMargins(10, 10, 0, 0)
    header_layout.addWidget(label)
    header_layout.addWidget(self.load_button)

    self.layout.insertLayout(0, header_layout)
