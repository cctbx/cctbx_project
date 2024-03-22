from ..widgets.tab import GUITab
from PySide2.QtWidgets import (QHBoxLayout, QVBoxLayout, QLabel, QPushButton)
from PySide2.QtGui import QIcon

from ..tabs.selection import SelectionListView, SelectionEntryView
from pathlib import Path

class RestraintStagingEntryView(SelectionEntryView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)

class RestraintStagingListView(SelectionListView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    header_layout = QHBoxLayout()
    label = QLabel("Restraint edits")
    current_font = label.font()
    current_font.setPointSize(16)
    current_font.setBold(False)
    label.setFont(current_font)

    self.load_button = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/save.svg'
    load_icon = QIcon(str(icon_path))
    self.load_button.setIcon(load_icon)
    self.load_button.setMaximumSize(50, 50)
    self.load_button.setContentsMargins(10, 10, 0, 0)
    self.load_button.setToolTip("Save as restraint edits file")
    header_layout.addWidget(label)
    header_layout.addWidget(self.load_button)

    self.layout.insertLayout(0, header_layout)



class RestraintStagingTabView(GUITab):
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    layout = QVBoxLayout()

    self.list_view = RestraintStagingListView(self)
    layout.addWidget(self.list_view)
    self.setLayout(layout)
