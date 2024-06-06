from PySide2.QtWidgets import QWidget, QVBoxLayout, QPushButton, QMenu, QAction
from PySide2.QtCore import Slot
from PySide2.QtGui import QIcon

from pathlib import Path

class ColorThemeButton(QWidget):
  def __init__(self,parent=None):
    super().__init__(parent)
    self.parent_explicit = parent
    self.button = QPushButton(self)
    icon_path = Path(__file__).parent / '../assets/icons/material/color_palette.svg'
    icon = QIcon(str(icon_path))
    self.button.setIcon(icon)
    self.button.setToolTip("Visual styling")
    self.menu = QMenu(self)

    # for more: mol-theme/color.js
    self.themes = { # Display name: molstar name
      'Uniform':'uniform',
      'By element':'element-symbol',
      'By Atom Id': 'atom-id',
      'By Chain Id': 'chain-id',
      'Occupancy' : 'occupancy',
      "Hydrophobicity": 'hydrophobicity'
    }
    self.actions = []
    for key,value in self.themes.items():
      action = QAction(key, self, checkable=True)
      self.menu.addAction(action)
      action.triggered.connect(self.option_selected)
      self.actions.append(action)

    self.button.clicked.connect(self.show_menu)

    layout = QVBoxLayout(self)
    layout.setContentsMargins(0,0,0,0)
    layout.addWidget(self.button)

  @property
  def state(self):
    return self.parent_explicit.state

  @property
  def ref(self):
    return self.parent_explicit.ref

  @Slot()
  def show_menu(self):
    global_pos = self.button.mapToGlobal(self.button.pos())
    self.menu.exec_(global_pos)

  @Slot()
  def option_selected(self):
    for action in self.actions:
      action.setChecked(False)

    action = self.sender()
    if action:
      key = action.text()
      self.themes[key]
      color = '#FFFF00'
      raise NotImplementedError
