from PySide2.QtWidgets import (
    QComboBox,
    QGridLayout,
    QTableView,
    QWidget
)

from ..widgets import ClickableHistogramSeaborn
from ..tabs.preprocess import PreprocessTab
from ..table import PandasTableView

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
    self.table_view = PandasTableView()
    self.table_view.setSelectionBehavior(QTableView.SelectRows)

    self.layout.addWidget(self.table_view)


  def add_histogram(self):
    # Histogram
    hist_widget = QWidget()
    self.hist_widget =  hist_widget
    g_layout = QGridLayout()


    ch = ClickableHistogramSeaborn(self.table_view.model()._df["Qscore"].values)

    g_layout.addWidget(ch,0,0)

    hist_widget.setLayout(g_layout)

    self.layout.insertWidget(0,hist_widget)
