from pathlib import Path

import pandas as pd
from PySide2.QtWidgets import QLabel,QPushButton, QHBoxLayout,QVBoxLayout, QApplication, QMainWindow, QTabWidget, QTableView, QWidget, QComboBox
from PySide2.QtGui import QStandardItemModel, QStandardItem, QIcon


from ..widgets import  FastTableView, PandasTableModel
from ..widgets.tab import GUITab,GUITabWidget


class CifTabView(GUITab):
  """
  View cif structure
  """
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    layout = QVBoxLayout()
    self.layout = layout
    # header with buttons (rename to save)
    header_layout = QHBoxLayout()
    label = QLabel("")
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
    header_layout.addWidget(label)
    header_layout.addWidget(self.load_button)
    
    self.layout.insertLayout(0, header_layout)

    # Create a combobox for top-level keys
    self.combobox = QComboBox()
    layout.addWidget(self.combobox)
    self.setLayout(layout)


    # add empty dataframe
    table = FastTableView()
    self.layout.addWidget(table)


  