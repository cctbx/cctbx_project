from PySide2.QtWidgets import QTableView

from .widgets import  PandasTableView
from .restraints.restraints import RestraintsTab

class BondsTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Bonds"):
    super().__init__(parent=parent,title=title)
    self.table = PandasTableView()
    self.table.setSelectionBehavior(QTableView.SelectRows)

    self.layout.addWidget(self.table)

class AnglesTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Angles"):
    super().__init__(parent=parent,title=title)
    self.table = PandasTableView()
    self.table.setSelectionBehavior(QTableView.SelectRows)

    self.layout.addWidget(self.table)

class DihedralsTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Dihedrals"):
    super().__init__(parent=parent,title=title)
    self.table = PandasTableView()
    self.table.setSelectionBehavior(QTableView.SelectRows)

    self.layout.addWidget(self.table)

class ChiralsTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Chirals"):
    super().__init__(parent=parent,title=title)
    self.table = PandasTableView()
    self.table.setSelectionBehavior(QTableView.SelectRows)

    self.layout.addWidget(self.table)

class PlanarityTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Planes"):
    super().__init__(parent=parent,title=title)
    self.table = PandasTableView()
    self.table.setSelectionBehavior(QTableView.SelectRows)

    self.layout.addWidget(self.table)

class NonbondedTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Nonbonded"):
    super().__init__(parent=parent,title=title)
    self.table = PandasTableView()
    self.table.setSelectionBehavior(QTableView.SelectRows)

    self.layout.addWidget(self.table)