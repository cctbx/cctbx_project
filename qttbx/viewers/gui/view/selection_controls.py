"""
The view for presenting the selection controls
"""
from pathlib import Path

import qtawesome as qta
from PySide2.QtGui import QIcon
from PySide2.QtCore import Qt, QEvent
from PySide2.QtCore import QObject, Signal
from PySide2.QtGui import QIcon
from PySide2.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLayoutItem,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpacerItem,
    QVBoxLayout,
    QWidget
)

from qttbx.viewers.gui.view.widgets.history_line_edit import HistoryLineEdit

# Enable high DPI scaling
QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps)

  
class ScrollableHierarchyWidget(QWidget):
  selected_change = Signal(str) # the value that was changed to

  def __init__(self,parent=None):
    super().__init__(parent=parent)
    self.setFocusPolicy(Qt.StrongFocus)

    # Create the main layout
    self.layout = QVBoxLayout()

    # Create a QListWidget instance
    self.list_widget = QListWidget()
    self.list_widget.keyPressFilter = KeyPressFilter(self.list_widget)
    self.list_widget.installEventFilter(self.list_widget.keyPressFilter)

    # Create a scroll area and set the list widget as its widget
    self.scroll_area = QScrollArea()
    self.scroll_area.setWidgetResizable(True)
    self.scroll_area.setWidget(self.list_widget)

    # Add the scroll area to the main layout
    self.layout.addWidget(self.scroll_area)

    # Add a label to show the current selection
    self.label = QLabel("Selected: None")
    self.layout.addWidget(self.label)

    # Set the layout for the main window
    self.setLayout(self.layout)



  def update_list_widget(self, items):
    # Clear the current items
    self.list_widget.clear()
    # Add new items to the list widget
    for item in items:
      self.list_widget.addItem(QListWidgetItem(str(item)))
  



# End edits dialog
class SearchSelectDialog(QDialog):
  def __init__(self,title="Selection search",parent=None):
    super().__init__(parent)
    self.setWindowTitle(title)
    self.center_on_right_half()
    self.list_width = 150
    self.chain_scroller = ScrollableHierarchyWidget(self)
    self.comp_scroller = ScrollableHierarchyWidget(self)
    self.seq_scroller = ScrollableHierarchyWidget(self)
    self.atom_scroller = ScrollableHierarchyWidget(self)
    layout = QVBoxLayout(self)
    # navigate selection
    hierarchy_label = QLabel("Hierarchy filters:")
    layout.addWidget(hierarchy_label)
    hierarchy_layout  = QHBoxLayout()
    # Chains
    chain_layout = QVBoxLayout()
    chain_label = QLabel("    Chain")
    chain_list = self.chain_scroller
    chain_list.setFixedWidth(self.list_width)
    chain_layout.addWidget(chain_label)
    chain_layout.addWidget(chain_list)
    hierarchy_layout.addLayout(chain_layout)

    # Residue components
    comp_layout = QVBoxLayout()
    comp_label = QLabel("   Residue")
    comp_list = self.comp_scroller
    comp_list.setFixedWidth(self.list_width)

    comp_layout.addWidget(comp_label)
    comp_layout.addWidget(comp_list)
    hierarchy_layout.addLayout(comp_layout)


    # Residues
    seq_layout = QVBoxLayout()
    seq_label = QLabel("   Sequence")
    seq_list = self.seq_scroller
    seq_list.setFixedWidth(self.list_width)

    seq_layout.addWidget(seq_label)
    seq_layout.addWidget(seq_list)
    hierarchy_layout.addLayout(seq_layout)

    # Atoms
    atom_layout = QVBoxLayout()
    atom_label = QLabel("   Atom")
    atom_list = self.atom_scroller
    atom_list.setFixedWidth(self.list_width)

    atom_layout.addWidget(atom_label)
    atom_layout.addWidget(atom_list)
    hierarchy_layout.addLayout(atom_layout)

    layout.addLayout(hierarchy_layout)



    # String Selection
    selection_edit = parent.selection_edit
    
    search_label = QLabel("Phenix selection string:")
    selection_edit.setPlaceholderText("Enter selection string...")
    selection_edit.setToolTip("Enter a Phenix selection string")
    layout.addWidget(search_label)
    layout.addWidget(selection_edit)

  def center_on_right_half(self):
    screen_geometry = QApplication.desktop().screenGeometry()
    screen_width = screen_geometry.width()
    screen_height = screen_geometry.height()

    dialog_width = self.width()
    dialog_height = self.height()

    # Calculate the center position for the right half of the screen
    right_half_center_x = screen_width * 1 / 2 - dialog_width / 2
    center_y = screen_height / 2 - dialog_height / 2

    self.move(int(right_half_center_x), int(center_y))

  

class KeyPressFilter(QObject):
  space_pressed = Signal()
  def __init__(self, parent=None):
    super().__init__(parent)

  def eventFilter(self, obj, event):
    if event.type() == QEvent.KeyPress:
      if event.key() == Qt.Key_Space:
        self.space_pressed.emit()
        return True  # Event is handled
    return super().eventFilter(obj, event)
