import os
from pathlib import Path

from PySide2.QtWidgets import QFileDialog, QColorDialog

from ..models import ModelEntryController, ModelListController
from ..maps import MapEntryController, MapListController
from ...view.models import ModelEntryView, ModelListView
from ...view.cif.cif_files import CifFileEntryView
from ...view.maps import MapEntryView, MapListView
from ..scroll_entry import ScrollEntryController
from ..scroll_list import ScrollableListController
from ..controller import Controller
from ...state.cif import CifFileData
from ...state.ref import CifFileRef



class CifFileEntryController(ScrollEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)



  def toggle_active_func(self,is_checked):
    pass
    # # TODO: here
    # if is_checked:
    #   self.state.active_ciffile_ref = self.ref
    # else:
    #   if self.state.active_ciffile_ref == self.ref:
    #     self.state.active_ciffile_ref = None
        

class CifFileListController(ScrollableListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    # Load button
    self.view.load_button.clicked.connect(self.showFileDialog)

    # update list
    self.state.signals.ciffile_change.connect(self.update)


  def showFileDialog(self):
    home_dir = os.path.expanduser("~")  # Cross-platform home directory
    fname = QFileDialog.getOpenFileName(self.view, 'Open file', home_dir)
    if fname[0]:
      filename = fname[0]
      filepath = str(Path(filename).absolute())
      #
      data = CifFileData(filepath=filepath)
      ref = CifFileRef(data=data)
      ref.label = data.filename
      self.state.add_ref(ref)



  def update(self):
    references_ciffiles = [value for key,value in self.state.references.items() if isinstance(value,CifFileRef)]
    for ref in references_ciffiles:
      if ref not in self.refs:
        if ref.show_in_list:
          entry_view = CifFileEntryView()
          entry_controller = CifFileEntryController(parent=self,view=entry_view,ref=ref)
          self.add_entry(entry_controller)
