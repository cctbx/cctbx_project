from PySide2.QtWidgets import QTableView

from ..widgets import  FastTableView
from ..restraints.restraints import RestraintsTab

class BondsTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Bonds"):
    super().__init__(parent=parent,title=title)
    self.table = FastTableView()
    self.table.setSelectionBehavior(QTableView.SelectRows)

    self.layout.addWidget(self.table)
