from pathlib import Path

from PySide2.QtGui import QIcon
from PySide2.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QSpacerItem,
    QVBoxLayout
)

from ..widgets.tab import GUITab
from ..widgets.tab import GUITab
from ..widgets.scroll_list import ScrollableListView
from ..widgets.scroll_entry import ScrollEntryView


class RestraintFileEntryView(ScrollEntryView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)

    # Already have active toggle and label at this point
    #  will add open file button and filepath label
    self.remove_stretch()
    # Open in folder
    ################
    self.button_files = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/folder.svg'
    icon = QIcon(str(icon_path))
    self.button_files.setIcon(icon)
    self.button_files.setToolTip("Open containing folder")
    self.button_files.setMaximumWidth(self._all_button_width)
    self.layout.addWidget(self.button_files)

    # Full filepath label
    tooltip = 'Location'
    name = 'Restraint file from GeoStd + Monomer Library'
    self.visible_name = self._truncate_string(name,max_len=30).ljust(30)
    self.label_filepath = QLabel(self.visible_name)
    self.label_filepath.setToolTip(tooltip)
    self.layout.addWidget(self.label_filepath)

    # Spacer
    spacer = QSpacerItem(20, 20, QSizePolicy.Expanding, QSizePolicy.Minimum) # expand on x axis
    self.layout.addItem(spacer)

    # Close
    self.button_close = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/close.svg'
    icon = QIcon(str(icon_path))
    self.button_close.setIcon(icon)
    self.button_close.setToolTip("Remove")
    self.button_close.setFixedSize(self._all_button_width,self._all_button_height)
    self.layout.addWidget(self.button_close)

  def resizeEvent(self, event):
    # Resize child objects
    # Call the base class implementation
    super().resizeEvent(event)

    # Set the width of the child widget to 50% of the parent widget's width
    parent_width = self.width()
    child_width = int(parent_width * 0.05)
    self.label_name.setFixedWidth(child_width)


  def remove_stretch(self):
    # Fix for problems with widget alignment across entries
    # Iterate through the items in the horizontal layout
    for i in reversed(range(self.layout.count())):
      item = self.layout.itemAt(i)
      if isinstance(item, QSpacerItem):
        # Remove the spacer item
        self.layout.removeItem(item)
        break  # Assuming there's only one stretch to remove

class RestraintFileListView(ScrollableListView):
  def __init__(self,parent=None,title="Restraints"):
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


class RestraintFileTabView(GUITab):
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    layout = QVBoxLayout()

    self.list_view = RestraintFileListView(self)
    layout.addWidget(self.list_view)
    self.setLayout(layout)
