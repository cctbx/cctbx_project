from PySide2.QtWidgets import QWidget, QVBoxLayout, QPushButton, QMenu, QAction
from PySide2.QtCore import Slot
from PySide2.QtGui import QIcon

from pathlib import Path

class RepresentationSelect(QWidget):
  def __init__(self,parent=None):
    super().__init__(parent)
    self.parent_explicit = parent
    self.button = QPushButton(self)
    icon_path = Path(__file__).parent / '../assets/icons/material/color_palette.svg'
    icon = QIcon(str(icon_path))
    self.button.setIcon(icon)
    self.button.setToolTip("Visual styling")
    self.menu = QMenu(self)

    self.options= { # Display name: molstar name
      'Ball & Stick':'ball-and-stick',
      'Cartoon':'cartoon'
    }
    self.actions = {}
    for key,value in self.options.items():
      action = QAction(key, self, checkable=True)
      self.menu.addAction(action)
      self.actions[value] = action 

    self.button.clicked.connect(self.show_menu)

    layout = QVBoxLayout(self)
    layout.setContentsMargins(0,0,0,0)
    layout.addWidget(self.button)

    self.selected_options = {value:False for key,value in self.options.items()}
    self.selected_options['ball-and-stick'] = True
    self.actions['ball-and-stick'].setChecked(True)


  @Slot()
  def show_menu(self):
    global_pos = self.button.mapToGlobal(self.button.pos())
    self.menu.exec_(global_pos)


    


    
