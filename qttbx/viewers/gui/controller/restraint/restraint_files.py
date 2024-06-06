import os
from pathlib import Path
import platform
import subprocess

from PySide2.QtWidgets import QFileDialog, QColorDialog

from ..models import ModelEntryController, ModelListController
from ..maps import MapEntryController, MapListController
from ...view.models import ModelEntryView, ModelListView
from ...view.restraint.restraint_files import RestraintFileEntryView
from ...view.maps import MapEntryView, MapListView
from ..scroll_entry import ScrollEntryController
from ..scroll_list import ScrollableListController
from ..controller import Controller
from ...state.cif import CifFileData
from ...state.ref import RestraintRef
from ...state.restraint import Restraint


class RestraintFileEntryController(ScrollEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)


    # Files
    self.view.button_files.clicked.connect(lambda: self.open_file_explorer(str(Path(self.ref.data.filepath).parent)))

    # Filepath label
    self.view.label_filepath.setText(str(self.ref.data.filepath))
    

  def toggle_active_func(self,*args):
    self.state.signals.restraint_activated.emit(self.ref)
    # switch to browser sub tab
    self.parent.parent.view.setCurrentIndex(1)

  def open_file_explorer(self,path):
    if platform.system() == 'Windows':
      subprocess.run(['explorer', path])
    elif platform.system() == 'Darwin':
      subprocess.run(['open', path])
    elif platform.system() == 'Linux':
      subprocess.run(['xdg-open', path])

class RestraintFileListController(ScrollableListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    # Load button
    self.view.load_button.clicked.connect(self.showFileDialog)

    # update list
    #self.state.signals.ciffile_change.disconnect(self.update) # disconnect superclass
    self.state.signals.restraints_change.connect(self.update)


  def showFileDialog(self):
    home_dir = Path.home()  # Cross-platform home directory
    fname = QFileDialog.getOpenFileName(self.view, 'Open file', str(home_dir))
    if fname[0]:
      filename = fname[0]
      filepath = Path(filename).absolute()

      # add to state
      restraint = Restraint(comp_id=filepath.stem,filepath=filepath)
      ref = RestraintRef(data=restraint,show=True)
      self.state.add_ref(ref)



  def update(self):
    # Go through all (plural) restraint references in state
    #   and update the list with the singular form
    for ref in self.state.references_restraint:
      if ref not in self.refs:
        if ref.show:
          ref.label = ref.data.comp_id
          entry_view = RestraintFileEntryView()
          entry_controller = RestraintFileEntryController(parent=self,view=entry_view,ref=ref)
          self.add_entry(entry_controller)
