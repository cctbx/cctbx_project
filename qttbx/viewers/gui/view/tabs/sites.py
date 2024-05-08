from PySide2.QtWidgets import  QVBoxLayout, QTableView, QComboBox

from ..table import  PandasTableView
from ..widgets.tab import GUITab

class SitesTabView(GUITab):
  """
  Click on atoms, make changes
  """
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    layout = QVBoxLayout()
    # Atoms tab
    self.table_view = PandasTableView()
    #self.tab2_content.clicked.connect(self.on_atom_select)
    self.table_view.setSelectionBehavior(QTableView.SelectRows)
    layout.addWidget(self.table_view)
    self.setLayout(layout)


  def on_first_visit(self):
    #msg = QMessageBox.information(self,"Notice", "The atom sites tab is currently very slow for large structures.")
    pass
