from PySide2.QtCore import Qt, Signal
from PySide2.QtWidgets import (
    QHBoxLayout,
    QMenu,
    QPushButton,
    QVBoxLayout
)

from ..widgets.tab import GUITab
from ..table import PandasTableView


class EditsTableView(PandasTableView):
  removeEdit = Signal(dict) # edits_ref, row dict corresponding to the restraint being edited

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
    removeAction = menu.addAction("Remove")

    action = menu.exec_(self.viewport().mapToGlobal(position))

    # if action == deleteAction:
    #   self.deleteRow(index.row())
    if action == removeAction:
      df = self.model().df
      # get named tuple for the row for edit
      specific_row = next(df.iloc[[index.row()]].itertuples(index=False), None)
      row_dict = specific_row._asdict()
      self.removeEdit.emit(row_dict)


class EditsTab(GUITab):
  """
  Base class which presents the Dataframe. Generic functionality for dataframes
  exists in the PandasTableView class
  """
  def __init__(self,parent=None,title="Geometry"):
    super().__init__(parent=parent)
    self.title = title
    self.layout = QVBoxLayout()
    self.setLayout(self.layout)



    # add controls at bottom to write edits file
    hlayout = QHBoxLayout()
    self.read_button = QPushButton("Read edits")
    self.read_button.setToolTip("Read a geometry edits phil file")
    self.read_button.setMaximumWidth(128)
    hlayout.addWidget(self.read_button)

    self.write_button = QPushButton("Write edits")
    # icon_path = Path(__file__).parent / '../assets/icons/material/save.svg'
    # load_icon = QIcon(str(icon_path))
    # self.write_button.setIcon(load_icon)
    self.write_button.setToolTip("Write edits file")
    self.write_button.setMaximumWidth(128)
    hlayout.addWidget(self.write_button)
    self.layout.addLayout(hlayout)


    # Edits table
    self.table_view = EditsTableView()
    self.layout.addWidget(self.table_view)


class BondsTableTabView(EditsTab):

  def __init__(self,parent=None,title="Bonds"):
    super().__init__(parent=parent,title=title)



class AnglesTableTabView(EditsTab):
  def __init__(self,parent=None,title="Angles"):
    super().__init__(parent=parent,title=title)


class DihedralsTableTabView(EditsTab):
  def __init__(self,parent=None,title="Dihedrals"):
    super().__init__(parent=parent,title=title)


class ChiralsTableTabView(EditsTab):
  def __init__(self,parent=None,title="Chirals"):
    super().__init__(parent=parent,title=title)


class PlanarityTableTabView(EditsTab):
  def __init__(self,parent=None,title="Planes"):
    super().__init__(parent=parent,title=title)


class NonbondedTableTabView(EditsTab):
  def __init__(self,parent=None,title="Nonbonded"):
    super().__init__(parent=parent,title=title)
