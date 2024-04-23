from PySide2.QtWidgets import QTableView
from PySide2.QtCore import QObject, QAbstractTableModel,  Qt, QTimer, QPoint, Signal
from PySide2.QtWidgets import QApplication, QPushButton, QMenu, QMainWindow, QVBoxLayout, QWidget
from PySide2.QtCore import QObject, QAbstractTableModel,  Qt, QTimer, QPoint, Signal

from ..table import PandasTableView
from PySide2.QtWidgets import (QHBoxLayout, QVBoxLayout, QLabel, QPushButton)
from ..widgets.tab import GUITab


"""
View for the Pandas table tabs for geometry
"""


class RestraintsTableView(PandasTableView):
  addEdit = Signal(dict) # row dict corresponding to the restraint being edited

  def __init__(self, parent=None, default_col_width=75):
    super().__init__(parent=parent)

    self.setContextMenuPolicy(Qt.CustomContextMenu)
    self.customContextMenuRequested.connect(self.showContextMenu)

  def showContextMenu(self, position):
    index = self.indexAt(position)
    if not index.isValid():
        return

    menu = QMenu(self)
    #editAction = menu.addAction('Edit')
    #deleteAction = menu.addAction('Delete')
    toggleEditAction = menu.addAction("Create Edit")

    action = menu.exec_(self.viewport().mapToGlobal(position))

    # if action == deleteAction:
    #   self.deleteRow(index.row())
    if action == toggleEditAction:
      df = self.model().df
      # get named tuple for the row for edit
      specific_row = next(df.iloc[[index.row()]].itertuples(index=False), None)
      row_dict = specific_row._asdict()
      self.addEdit.emit(row_dict)



class RestraintsTab(GUITab):
  """
  Base class which presents the Dataframe. Generic functionality for dataframes
  exists in the PandasTableView class
  """
  def __init__(self,parent=None,title="Geometry"):
    super().__init__(parent=parent)
    self.title = title
    self.layout = QVBoxLayout()
    self.setLayout(self.layout)
    self.table_view = RestraintsTableView()
    self.layout.addWidget(self.table_view)


class BondsTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Bonds"):
    super().__init__(parent=parent,title=title)


class AnglesTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Angles"):
    super().__init__(parent=parent,title=title)


class DihedralsTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Dihedrals"):
    super().__init__(parent=parent,title=title)


class ChiralsTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Chirals"):
    super().__init__(parent=parent,title=title)


class PlanarityTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Planes"):
    super().__init__(parent=parent,title=title)


class NonbondedTableTabView(RestraintsTab):
  def __init__(self,parent=None,title="Nonbonded"):
    super().__init__(parent=parent,title=title)
