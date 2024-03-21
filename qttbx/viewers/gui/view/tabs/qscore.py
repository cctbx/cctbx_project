import numpy as np
from PySide2.QtWidgets import QWidget, QGridLayout, QTableView, QComboBox

from ..widgets import ClickableHistogramSeaborn
from ..widgets.tab import GUITab
from ..tabs.preprocess import PreprocessTab
from ..widgets import  PandasTableView

class QscoreTab(PreprocessTab):
  def __init__(self,parent=None):
    title = "Qscore"
    action = "Calculate Qscore"
    tool_tip = "Calculate Qscore to view results"
    super().__init__(parent=parent,title=title,action=action,tool_tip=tool_tip)

    # Create a combobox for atom/residue toggling
    self.combobox = QComboBox()
    self.combobox.setMaximumSize(200,48)
    self.combobox.addItems(["Q-Atoms","Q-Residues"])
    self.layout.addWidget(self.combobox)

    # Table
    self.table = PandasTableView()
    self.table.setSelectionBehavior(QTableView.SelectRows)

    self.layout.addWidget(self.table)


  def add_histogram(self):
    # Histogram
    hist_widget = QWidget()
    self.hist_widget =  hist_widget
    g_layout = QGridLayout()


    ch = ClickableHistogramSeaborn(self.table.model()._df["Qscore"].values)

    g_layout.addWidget(ch,0,0)

    hist_widget.setLayout(g_layout)

    self.layout.insertWidget(0,hist_widget)
