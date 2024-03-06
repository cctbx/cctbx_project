
from functools import partial
from dataclasses import replace
import platform
import subprocess
from pathlib import Path

from PySide2.QtWidgets import QFileDialog, QColorDialog

from ..view.models import ModelEntryView
from .scroll_entry import ScrollEntryController
from .scroll_list import ScrollableListController
from ..state.ref import SelectionRef, ModelRef

class ModelLikeEntryController(ScrollEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)

    # Signals
    self.view.button_viz.clicked.connect(self.toggle_visibility)
    #self.state.signals.update.connect(self.update)

    # Representation
    for key,action in self.view.button_rep.actions.items():
      action.triggered.connect(partial(self.representation_selected, action))
      for rep_name in ref.style.representation:
        # if representations are preselected in the style, enable them
        self.view.button_rep.selected_options[rep_name] = True
        self.view.button_rep.actions[rep_name].setChecked(True)

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

    self.view.button_files.clicked.connect(lambda: self.open_file_explorer(str(folder)))


  def toggle_visibility(self,event):
      value= self.view.button_viz.is_on
      #self.state.emitter.signal_viz_change.emit(self.ref.id,value)
      style = replace(self.ref.style,visible=(not self.ref.style.visible))
      self.ref.style = style

  def representation_selected(self,action):

    #action = self.sender() # won't work why?
    if action:
      key = action.text()
      # option = self.view.button_rep.options[key]
      # selected_options = self.view.button_rep.selected_options
      # current_state = selected_options[option]
      # action.setChecked(not current_state)
      # selected_options[option] = not selected_options[option]
      # reps = [key for key,value in selected_options.items() if value]
      # #print("Replacing data for style from ref: ",self.ref, self.ref.id)
      # #print("Old:")
      # #print(self.ref.style.to_json(indent=2))
      # style = replace(self.ref.style,representation=reps)
      # #print("New:")
      # #print(style.to_json(indent=2))
      # self.ref.style = style
      #self.state.emitter.signal_repr_change.emit(*emission)
      # DEBUG:
      if key.startswith("Ball"):
        key = 'ball-and-stick'
      else:
        key = 'ribbon'
      self.parent.parent.parent.molstar.viewer.representation_query(self.ref.id_molstar,self.ref.query.to_json(),key)

  def open_file_explorer(self,path):
    if platform.system() == 'Windows':
      subprocess.run(['explorer', path])
    elif platform.system() == 'Darwin':
      subprocess.run(['open', path])
    elif platform.system() == 'Linux':
      subprocess.run(['xdg-open', path])

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
    # TODO: here
    if is_checked:
      self.state.active_model_ref = self.ref
    else:
      if self.state.active_model_ref == self.ref:
        self.state.active_model_ref = None

class ModelListController(ScrollableListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.openFileDialog = None
    # Load button
    self.view.load_button.clicked.connect(self.showFileDialog)

    # update list
    self.state.signals.model_change.connect(self.update)
    self.state.signals.update.connect(self.update) # generic update


  def showFileDialog(self):
    self.state.is_updating = False
    home_dir = str(Path.home())

    self.openFileDialog = QFileDialog(self.view)
    self.openFileDialog.setFileMode(QFileDialog.AnyFile)
    if self.openFileDialog.exec_():
        file_path = self.openFileDialog.selectedFiles()[0]

        filepath = str(Path(file_path).absolute())
        print(f"File selected: {filepath}")
        _ = self.state.data_manager.process_model_file(filepath)
        self.state._data_manager_changed()



  def update(self):
    for ref in self.state.references_model:
      if ref not in self.refs:
        if ref.show_in_list:
          entry_view = ModelEntryView()
          entry_controller = ModelEntryController(parent=self,view=entry_view,ref=ref)
          self.add_entry(entry_controller)
          if ref == self.state.active_model_ref:
            entry_controller.active = True # display radio button
