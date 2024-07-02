import qtawesome as qta
from PySide2.QtCore import Qt, Signal
from PySide2.QtWidgets import (
    QMenu,
    QLabel,
    QWidget,
    QSpacerItem,
    QSizePolicy,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton
)

from ..table import PandasTableView
from ..widgets.tab import GUITab
from ..widgets.filter import TableFilter


"""
View for the Pandas table tabs for geometry
"""
class GeometryTableView(PandasTableView):
  addEdit = Signal(dict) # row dict corresponding to the restraint being edited

  def __init__(self, parent=None):
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

class EditControls(QWidget):
  """
  Controls to edit entries
  """
  def __init__(self, parent=None):
    super().__init__(parent=parent)
    
    # Main vertical layout
    self.layout = QHBoxLayout()

    spacer = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
    self.layout.addSpacerItem(spacer)

    # Dropdown for filter selections
    edit_layout = QVBoxLayout()

    edit_layout.addWidget(QLabel("Modify/delete entry"))
    self.edit_button = QPushButton("Create Edit")
    icon = qta.icon("mdi.pencil")
    self.edit_button.setIcon(icon)
    edit_layout.addWidget(self.edit_button)
    self.layout.addLayout(edit_layout)

    self.setLayout(self.layout)



class GeometryTab(GUITab):
  """
  Base class which presents the Dataframe. Generic functionality for dataframes
  exists in the PandasTableView class
  """
  def __init__(self,parent=None,title="Geometry"):
    super().__init__(parent=parent)
    self.title = title
    self.layout = QVBoxLayout()
    self.setLayout(self.layout)

    # Create edit button
    self.edit_controls =  EditControls()


    # Filter controls
    self.filter = TableFilter()

    controls_layout = QHBoxLayout()
    controls_layout.addWidget(self.edit_controls)
    controls_layout.addWidget(self.filter)
    self.layout.addLayout(controls_layout)

    # Table
    self.table_view = GeometryTableView()
    self.layout.addWidget(self.table_view)


class BondsTableTabView(GeometryTab):
  def __init__(self,parent=None,title="Bonds"):
    super().__init__(parent=parent,title=title)


class AnglesTableTabView(GeometryTab):
  def __init__(self,parent=None,title="Angles"):
    super().__init__(parent=parent,title=title)


class DihedralsTableTabView(GeometryTab):
  def __init__(self,parent=None,title="Dihedrals"):
    super().__init__(parent=parent,title=title)


class ChiralsTableTabView(GeometryTab):
  def __init__(self,parent=None,title="Chirals"):
    super().__init__(parent=parent,title=title)


class PlanarityTableTabView(GeometryTab):
  def __init__(self,parent=None,title="Planes"):
    super().__init__(parent=parent,title=title)


class NonbondedTableTabView(GeometryTab):
  def __init__(self,parent=None,title="Nonbonded"):
    super().__init__(parent=parent,title=title)
