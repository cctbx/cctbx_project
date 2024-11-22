from pathlib import Path
import qtawesome as qta
from PySide2.QtCore import Slot
from PySide2.QtGui import QIcon
from PySide2.QtWidgets import (
    QAction,
    QMenu,
    QPushButton,
    QVBoxLayout,
    QWidget
)

class RepresentationSelect(QWidget):
  def __init__(self,parent=None):
    super().__init__(parent)
    self.parent_explicit = parent
    self.button = QPushButton(self)
    icon = qta.icon("mdi.palette")
    self.button.setIcon(icon)
    self.button.setToolTip("Visual styling")
    self.menu = QMenu(self)

    self.options= { # Display name: molstar name
      'Ball & Stick':'ball-and-stick',
      'Cartoon':'cartoon'
    }
    self.actions = {}
    for key,value in self.options.items():
      action = QAction(key, self, checkable=False)
      self.menu.addAction(action)
      self.actions[value] = action

    self.button.clicked.connect(self.show_menu)

    layout = QVBoxLayout(self)
    layout.setContentsMargins(0,0,0,0)
    layout.addWidget(self.button)
    # maintain a registry of checked state
    self.selected_options = {value:False for key,value in self.options.items()}

    # initialize with all checked as false
    for key,value in self.selected_options.items():
      self.actions[key].setChecked(False)


  @Slot()
  def show_menu(self):
    global_pos = self.button.mapToGlobal(self.button.pos())
    self.menu.exec_(global_pos)
