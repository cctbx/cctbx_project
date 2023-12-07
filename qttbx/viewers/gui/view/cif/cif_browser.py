from pathlib import Path

import pandas as pd
from PySide2.QtWidgets import  QFrame, QLabel,QPushButton, QHBoxLayout,QVBoxLayout, QApplication, QMainWindow, QTabWidget, QTableView, QWidget, QComboBox
from PySide2.QtGui import QStandardItemModel, QStandardItem, QIcon


from ..widgets import  FastTableView, PandasTableModel
from ..widgets.tab import GUITab,GUITabWidget


class CifBrowserTabView(GUITab):
  """
  View cif structure
  """
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    layout = QVBoxLayout()
    self.layout = layout


    # Horizontal layout for comboboxes and save button
    combo_button_layout = QHBoxLayout()

    # Create a combobox for data level keys
    self.combobox_data = QComboBox()
    combo_button_layout.addWidget(self.combobox_data,4)

    # Create a combobox for block level keys
    self.combobox_block = QComboBox()
    combo_button_layout.addWidget(self.combobox_block,4)

    #Create the vertical separator
    separator = QFrame()
    separator.setFrameShape(QFrame.VLine)
    separator.setFrameShadow(QFrame.Sunken)
    combo_button_layout.addWidget(separator,1)

    # Save button
    self.save_button = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/save.svg'
    save_icon = QIcon(str(icon_path))
    self.save_button.setIcon(save_icon)
    # self.save_button.setMaximumSize(32, 48)
    # combo_button_layout.setContentsMargins(5, 5, 5, 5)  # left, top, right, bottom
    #self.layout.setAlignment(Qt.AlignVCenter)
    combo_button_layout.addWidget(self.save_button,1)

    # Add the horizontal layout to the main layout
    layout.addLayout(combo_button_layout)

    self.setLayout(layout)

    # add empty dataframe
    table = FastTableView()
    self.layout.addWidget(table)

