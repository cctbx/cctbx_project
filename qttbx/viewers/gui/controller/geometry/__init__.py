from pathlib import Path

from PySide2.QtWidgets import QFileDialog

from ...view.geometry import GeoEntryView
from ..scroll_entry import ScrollEntryController
from ..scroll_list import ScrollableListController
from ...state.geometry import Geometry
from ...state.ref import GeometryRef



class GeoEntryController(ScrollEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)


  def toggle_active_func(self,is_checked):
    pass


class GeoListController(ScrollableListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    # Load button
    self.view.list_view.load_button.clicked.connect(self.showFileDialog)

    # update list
    self.state.signals.geo_change.connect(self.update)


  def showFileDialog(self):
    home_dir = Path.home()  # Cross-platform home directory
    fname = QFileDialog.getOpenFileName(self.view, 'Open file', str(home_dir))
    if fname[0]:
      filename = fname[0]
      filepath = Path(filename).absolute()

      # add to state
      data = Geometry(filepath=filepath)
      ref = GeometryRef(data=data,show=True)
      self.state.add_ref(ref)



  def update(self):
    for ref in self.state.references_geo:
      if ref not in self.refs:
        if ref.show:
          entry_view = GeoEntryView()
          entry_controller = GeoEntryController(parent=self,view=entry_view,ref=ref)
          self.add_entry(entry_controller)

          # Make active (checked) and load data
          entry_controller.active = True

