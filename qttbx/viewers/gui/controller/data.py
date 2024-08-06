from pathlib import Path

from PySide2.QtWidgets import QFileDialog

from .scroll_list import ScrollableListController
from .controller import Controller


class DataTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    #self.model_list_controller = ModelListController(parent=self,view=self.view.model_list_view)
    #self.map_list_controller = MapListController(parent=self,view=self.view.map_list_view)

    self.data_list_controller = GenericDataListController(parent=self,view=self.view.data_list_view)



class GenericDataListController(ScrollableListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self._first_load_completed = False
    self.openFileDialog = None
    # Load button
    self.view.load_button.clicked.connect(self.showFileDialog)

    # update list
    self.state.signals.model_change.connect(self.update)
    self.state.signals.update.connect(self.update) # generic update
    self.state.signals.geo_change.connect(self.update)


  def showFileDialog(self):
    str(Path.home())

    self.openFileDialog = QFileDialog(self.view)
    self.openFileDialog.setFileMode(QFileDialog.AnyFile)
    if self.openFileDialog.exec_():
        file_path = self.openFileDialog.selectedFiles()[0]
        if str(file_path).endswith(".geo"):
          self.log("Unable to open .geo file from here")
          # # add to state
          # data = Geometry(filepath=file_path)
          # ref = GeometryRef(data=data,show=True)
          # self.state.add_ref(ref)
        else:
          # All other files handled through datamanager
          filepath = str(Path(file_path).absolute())
          self.log(f"File selected: {filepath}")
          _ = self.state.data_manager.process_model_file(filepath)
          self.state._data_manager_changed()



  def check_only_model_is_active(self):
    # Upon first load, a single model should be made active
    if self._first_load_completed:
      return 
    else:
      if len(self.model_entries)==1:
        entry = self.model_entries[0]
        entry.active = True
        self._first_load_completed = True


  def update(self):
    # Main update of entries
    aspects = ["model"]
    for aspect in aspects:
      refs = getattr(self.state,f"references_{aspect}")
      for ref in refs:
        if ref not in self.refs:
          entry_view = ref.EntryViewClass()
          entry_controller = ref.EntryControllerClass(parent=self,view=entry_view,ref=ref)
          self.add_entry(entry_controller)
    
    # Other
    self.check_only_model_is_active()
