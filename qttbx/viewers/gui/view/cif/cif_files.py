from pathlib import Path

from PySide2.QtGui import QIcon
from PySide2.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout
)

from ..widgets.tab import GUITab
from ..widgets.scroll_list import ScrollableListView
from ..widgets.scroll_entry import ScrollEntryView

class CifFileEntryView(ScrollEntryView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)


    # Close
    self.button_close = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/close.svg'
    icon = QIcon(str(icon_path))
    self.button_close.setIcon(icon)
    self.button_close.setToolTip("Remove")
    self.button_close.setFixedSize(self._all_button_width,self._all_button_height)
    self.layout.addWidget(self.button_close)


class CifFileListView(ScrollableListView):
  def __init__(self,parent=None,title="CIF Files"):
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


class CifFileTabView(GUITab):
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    layout = QVBoxLayout()

    self.list_view = CifFileListView(self)
    layout.addWidget(self.list_view)
    self.setLayout(layout)