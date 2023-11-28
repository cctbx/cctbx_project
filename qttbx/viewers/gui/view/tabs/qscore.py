import numpy as np
from PySide2.QtWidgets import QWidget, QGridLayout, QTableView

from ..widgets import ClickableHistogramSeaborn
from ..widgets.tab import GUITab
from ..tabs.preprocess import PreprocessTab
from ..widgets import  FastTableView

class QscoreTab(PreprocessTab):
  def __init__(self,parent=None):
    title = "Qscore"
    action = "Calculate Qscore"
    tool_tip = "Calculate Qscore to view results"
    super().__init__(parent=parent,title=title,action=action,tool_tip=tool_tip)

  
    # Table
    self.table = FastTableView()
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
  