from pathlib import Path

from PySide2.QtWidgets import QFileDialog

from ...view.geometry.geo_files import GeoFileEntryView
from ..scroll_entry import ScrollEntryController
from ..scroll_list import ScrollableListController
from ...state.geometry import Geometry
from ...state.ref import GeometryRef



class GeoFileEntryController(ScrollEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)


  def toggle_active_func(self,is_checked):
    assert False, "Deprecated"
    # geo = Geometry.from_geo_file(self.ref.data.filepath)
    # ref = GeometryRef(data=geo,model_ref= self.state.active_model_ref)
    # self.state.add_ref(ref)


class GeoFileListController(ScrollableListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    # # Load button
    # self.view.list_view.load_button.clicked.connect(self.showFileDialog)

    # update list
    self.state.signals.geofile_change.connect(self.update)



  def update(self):
    for ref in self.state.references_geo:
      if ref not in self.refs:
        if ref.show:
          entry_view = GeoFileEntryView()
          entry_controller = GeoFileEntryController(parent=self,view=entry_view,ref=ref)
          self.add_entry(entry_controller)

          # Make active (checked) and load data
          entry_controller.active = True

