from PySide2.QtWidgets import  QVBoxLayout, QTableView

from ..widgets import  FastTableView
from ..widgets.tab import GUITab

class SitesTabView(GUITab):
  """
  Click on atoms, make changes
  """
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    layout = QVBoxLayout()
    # Atoms tab
    self.table = FastTableView()
    #self.tab2_content.clicked.connect(self.on_atom_select)
    self.table.setSelectionBehavior(QTableView.SelectRows)
    
  
    layout.addWidget(self.table)
    self.setLayout(layout)



  