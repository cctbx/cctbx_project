from functools import partial
from dataclasses import replace
import platform
import subprocess
from pathlib import Path

from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
    QColorDialog,
    QFileDialog,
    QMenu,
    QPushButton
)

from ..view.models import ModelEntryView
from .scroll_entry import ScrollEntryController
from .scroll_list import ScrollableListController
from ..state.ref import SelectionRef, ModelRef, GeometryRef
from ..state.geometry import Geometry


class ModelLikeEntryController(ScrollEntryController):
  # An abstract base class
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)

    # Signals
    self.view.button_viz.clicked.connect(self.toggle_visibility)

    # Representation
    for key,action in self.view.button_rep.actions.items():
      action.triggered.connect(partial(self.representation_selected, action))


    # Color
    self.view.button_color.clicked.connect(self.show_color_dialog)

    # Folder
    file_based_ref = None
    if isinstance(self.ref,SelectionRef):
      if isinstance(self.ref.model_ref,ModelRef):
        file_based_ref = self.ref.model_ref
    else:
      file_based_ref = self.ref

    folder = Path.home()
    if hasattr(file_based_ref.data,"filepath"):
      if file_based_ref.data.filepath is not None:
        folder = Path(file_based_ref.data.filepath).parent

    if folder is None or not folder.exists():
      folder = Path.home()

    #self.view.button_files.clicked.connect(lambda: self.open_file_explorer(str(folder)))

  #   # Geometry
  #   self.view.button_restraints.setContextMenuPolicy(Qt.CustomContextMenu)
  #   self.view.button_restraints.customContextMenuRequested.connect(self.showContextMenu)
  #   self.view.button_restraints.mousePressEvent = self.buttonMousePressEvent  # Override mousePressEvent
    
  #   # set geo checkbox if applicable
  #   self.update_geo()

  # def update_geo(self,**args):
  #   ref = self.ref
  #   if hasattr(self.ref,"model_ref"):
  #     ref = self.ref.model_ref
  #   if ref.geometry_ref is not None:
  #     self.view.geo_checkbox.setChecked(True)

  def toggle_visibility(self,event):

      is_on= self.view.button_viz.is_on
      # For some reason this is inverted
      is_on = not is_on
      print("Toggling visibility, is on? : ",is_on)
      if is_on:
        self.state.signals.selection_hide.emit(self.ref)
      else:
        self.state.signals.selection_show.emit(self.ref)

  def representation_selected(self,action):

    #action = self.sender() # won't work why?
    if action:
      key = action.text()
      if key.startswith("Ball"):
        key = 'ball-and-stick'
      else:
        key = 'ribbon'
      self.state.signals.selection_rep_show.emit(self.ref,key)
      #self.parent.parent.parent.molstar.viewer.representation_query(self.ref.id_molstar,self.ref.query.to_json(),key)

  def open_file_explorer(self,path):
    if platform.system() == 'Windows':
      subprocess.run(['explorer', path])
    elif platform.system() == 'Darwin':
      subprocess.run(['open', path])
    elif platform.system() == 'Linux':
      subprocess.run(['xdg-open', path])

  def buttonMousePressEvent(self, event):
    if event.button() == Qt.LeftButton:
        self.showContextMenu(event.pos())
    else:
        QPushButton.mousePressEvent(self.view.button_restraints, event)  # Call the original mousePressEvent

  def showContextMenu(self, position):
    contextMenu = QMenu(self.view)

    # Add actions to the context menu
    action1 = contextMenu.addAction("Generate geometry")
    action2 = contextMenu.addAction("Open geometry file")

    # Connect actions to functions/slots
    action1.triggered.connect(self.process_and_make_restraints)
    action2.triggered.connect(self.load_geometry)

    # Show the context menu at the button's position
    contextMenu.exec_(self.view.button_restraints.mapToGlobal(position))

  def process_and_make_restraints(self):
    assert False, "Deprecated"
      # if self.ref.has_geometry:
      #   self.state.notify("Already have restraints loaded. New processing will replace existing restraints.")
      # try:
      #   restraints = Geometry.from_model_ref(self.ref)
      #   ref = GeometryRef(restraints,self.state.active_model_ref)
      #   self.state.add_ref(ref)
      # except:
      #   self.state.notify("Failed to process and make restraints.")
      #   raise

  def load_geometry(self):
      if self.ref.has_geometry:
        self.state.notify("Already have restraints loaded. Will now replace existing restraints.")
      self.open_geometry_file_dialog()

  def open_geometry_file_dialog(self):
    assert False, "Deprecated"
    # self.openFileDialog = QFileDialog(self.view)
    # self.openFileDialog.setFileMode(QFileDialog.AnyFile)
    # if self.openFileDialog.exec_():
    #     file_path = self.openFileDialog.selectedFiles()[0]

    #     filepath = str(Path(file_path).absolute())
    #     print(f"Geometry file selected: {filepath}")

    #     data = Geometry.from_geo_file(filepath)
    #     ref = GeometryRef(data,self.state.active_model_ref)
    #     self.state.add_ref(ref)

  def show_color_dialog(self):
    color_dialog = QColorDialog(self.view)
    button_pos = self.view.button_color.mapToGlobal(self.view.button_color.rect().topLeft())
    color_dialog.move(button_pos)
    color = color_dialog.getColor()

    if color.isValid():
      #print("User selected color:", color.name())

      style = replace(self.ref.style,color_theme='uniform',color=color.name())
      self.ref.style = style




class ModelEntryController(ModelLikeEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)

    self.state.signals.model_change.connect(self.handle_model_change)


  def handle_model_change(self,ref):
    if ref is None or ref.id != self.ref.id:
      self.view.active_toggle.is_checked = False

  def toggle_active_func(self,is_checked):
    if is_checked:
      self.state.active_model_ref = self.ref
    else:
      # Deactivate this as active model ref
      if self.state.active_model_ref == self.ref:
        self.state.active_model_ref = None

class ModelListController(ScrollableListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self._first_load_completed = False
    self.openFileDialog = None
    # Load button
    self.view.load_button.clicked.connect(self.showFileDialog)

    # update list
    self.state.signals.model_change.connect(self.update)
    self.state.signals.update.connect(self.update) # generic update


  def showFileDialog(self):
    str(Path.home())

    self.openFileDialog = QFileDialog(self.view)
    self.openFileDialog.setFileMode(QFileDialog.AnyFile)
    if self.openFileDialog.exec_():
        file_path = self.openFileDialog.selectedFiles()[0]

        filepath = str(Path(file_path).absolute())
        print(f"File selected: {filepath}")
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
    for ref in self.state.references_model:
      if ref not in self.refs:
        if ref.show:
          entry_view = ModelEntryView()
          entry_controller = ModelEntryController(parent=self,view=entry_view,ref=ref)
          self.add_entry(entry_controller)
          if ref == self.state.active_model_ref:
            entry_controller.active = True # display radio button
    self.check_only_model_is_active()