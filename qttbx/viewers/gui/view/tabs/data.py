from PySide2.QtWidgets import QVBoxLayout

from ..widgets.tab import GUITab
from ..widgets.scroll_list import ScrollableListView
from PySide2.QtWidgets import (QHBoxLayout, QVBoxLayout, QLabel, QPushButton)
from PySide2.QtGui import QIcon
from pathlib import Path

class DataTabView(GUITab):
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    layout = QVBoxLayout()


    self.data_list_view = GenericDataListView(self)
    layout.addWidget(self.data_list_view)
    self.setLayout(layout)

    # self.model_list_view = ModelListView(self)
    # layout.addWidget(self.model_list_view)
    # self.setLayout(layout)

    # self.map_list_view = MapListView(self)
    # layout.addWidget(self.map_list_view)



class GenericDataListView(ScrollableListView):
  def __init__(self,parent=None,title=""):
    super().__init__(parent=parent)
    header_layout = QHBoxLayout()
    label = QLabel(title)
    current_font = label.font()
    current_font.setPointSize(16)
    current_font.setBold(False)
    label.setFont(current_font)

    self.load_button = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/plus.svg'
    load_icon = QIcon(str(icon_path))
    self.load_button.setIcon(load_icon)
    self.load_button.setMaximumSize(50, 50)
    self.load_button.setContentsMargins(10, 10, 0, 0)
    header_layout.addWidget(label)
    header_layout.addWidget(self.load_button)

    self.layout.insertLayout(0, header_layout)