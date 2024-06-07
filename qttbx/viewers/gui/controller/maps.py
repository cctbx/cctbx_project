from pathlib import Path
import time

from PySide2.QtWidgets import (
    QFileDialog,
    QMessageBox
)

from ..view.maps import MapEntryView
from .scroll_entry import ScrollEntryController
from .scroll_list import ScrollableListController
from .widgets import ISOWidgetController, OpacityWidgetController


class MapEntryController(ScrollEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)

    self.iso_widget_controller = ISOWidgetController(parent=self,view=self.view.iso_widget,map_ref=ref)
    self.opacity_widget_controller = OpacityWidgetController(parent=self,view=self.view.opacity_widget,map_ref=ref)

    # Signals
    self.state.signals.map_change.connect(self.handle_map_change)


  def handle_map_change(self,ref):
    if ref is None or ref.id != self.ref.id:
      self.view.active_toggle.is_checked = False

  def toggle_active_func(self,is_checked):
    if is_checked:
      if self.state.active_model_ref is None:
        QMessageBox.information(self.view, 'Stopped', 'Cannot activate a map without an active model.')
        self.view.active_toggle.is_checked = False
      else:
        if self.state.active_map_ref is not None:
          if self.state.active_map_ref != self.ref:
            _active_model_ref = self.state.active_model_ref
            t_wait = 10 # time delay hack until async loading is completely sorted out
            t0 = time.time()
            msg = 'When switching maps, the viewer is fully cleared and reloaded with the active map and model.'
            self.state.signals.clear.emit(msg)
            #print("Now need to load model: ",_active_model_ref.id, "and map: ",self.ref.id)
            #print(self.state.active_model_ref, self.state.active_map_ref)
            self.state.active_model_ref = _active_model_ref
            self.state.active_model_ref.entry.view.active_toggle.is_checked = True
            t1 = time.time()
            t_diff = t1-t0
            if t_diff < t_wait:
              remaining_time = t_wait - t_diff
              time.sleep(remaining_time)

            self.state.active_map_ref = self.ref
            self.state.active_map_ref.entry.view.active_toggle.is_checked = True
        else:
          self.state.active_map_ref = self.ref


    else:
      if self.state.active_map_ref == self.ref:
        self.state.active_map_ref = None

class MapListController(ScrollableListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    # Load button
    self.view.load_button.clicked.connect(self.showFileDialog)

    # update list
    self.state.signals.map_change.connect(self.update)
    self.state.signals.update.connect(self.update)

  def showFileDialog(self):
    str(Path.home())

    self.openFileDialog = QFileDialog(self.view)
    self.openFileDialog.setFileMode(QFileDialog.AnyFile)
    if self.openFileDialog.exec_():
        file_path = self.openFileDialog.selectedFiles()[0]

        filepath = str(Path(file_path).absolute())
        print(f"File selected: {filepath}")
        _ = self.state.data_manager.process_real_map_file(filepath)
        self.state._data_manager_changed()

  def update(self):
    for ref in self.state.references_map:
      if ref not in self.refs:
        if ref.show:
          entry_view = MapEntryView()
          entry_controller = MapEntryController(parent=self,view=entry_view,ref=ref)
          self.add_entry(entry_controller)
          if ref == self.state.active_map_ref or len(self.state.references_map)==1:
            entry_controller.active = True # display radio button
