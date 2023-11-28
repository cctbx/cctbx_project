from PySide2.QtWidgets import QVBoxLayout

from ..widgets.tab import GUITab
from ..models import ModelListView
from ..maps import MapListView


class DataTabView(GUITab):
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    layout = QVBoxLayout()
    
    self.model_list_view = ModelListView(self)
    layout.addWidget(self.model_list_view)
    self.setLayout(layout)

    self.map_list_view = MapListView(self)
    layout.addWidget(self.map_list_view)



