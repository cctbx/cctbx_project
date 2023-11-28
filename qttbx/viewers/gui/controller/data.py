
from .scroll_entry import ScrollEntryController
from .scroll_list import ScrollableListController
from .models import ModelEntryController, ModelListController
from .maps import MapEntryController, MapListController
from ..view.models import ModelEntryView, ModelListView
from ..view.maps import MapEntryView, MapListView
from .controller import Controller

class DataTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.model_list_controller = ModelListController(parent=self,view=self.view.model_list_view)
    self.map_list_controller = MapListController(parent=self,view=self.view.map_list_view)
