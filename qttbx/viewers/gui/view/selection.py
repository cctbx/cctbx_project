from pathlib import Path

from PySide2.QtGui import QIcon
from PySide2.QtWidgets import (
    QPushButton,
    QVBoxLayout
)

from .widgets.tab import GUITab
from .scroll_list import ScrollableListView
from .models import ModelLikeEntryView


class SelectionEntryView(ModelLikeEntryView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)


    # Get info
    self.button_info = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/info.svg'
    icon = QIcon(str(icon_path))
    self.button_info.setIcon(icon)
    self.button_info.setToolTip("View selection info")
    self.button_info.setFixedSize(self._all_button_width,self._all_button_height)
    insert_index = max(self.layout.count() - 1, 0) # second to last
    self.layout.insertWidget(insert_index,self.button_info)


class SelectionListView(ScrollableListView):
  def __init__(self,parent=None):
        super().__init__(parent=parent)

class SelectionsTabView(GUITab):
  def __init__(self, parent=None):
    super().__init__(parent=parent)

    layout = QVBoxLayout()

    self.selection_list_view = SelectionListView(self)

    layout.addWidget(self.selection_list_view)
    self.setLayout(layout)

