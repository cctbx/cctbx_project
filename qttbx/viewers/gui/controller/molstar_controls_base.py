from pathlib import Path
import platform
import subprocess
import time

from PySide2.QtCore import Slot
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
    QFileDialog,
    QMenu,
    QPushButton
)

import iotbx
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str

from .controller import Controller
from .selection_controls import ScrollableHierarchyController, SearchSelectDialogController
from ..view.selection_controls import SearchSelectDialog
from ..model.ref import Ref, SelectionRef
from ..model.selection import Selection


class MolstarControlsController(Controller):
  """
  A controller which manages the minimum base level controls present in the molstar GUI tab
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.last_selection_exec = time.time()

    ## Signals
    
    # Enable return key to execute selection
    self.view.selection_edit.returnPressed.connect(self.execute_selection)
    self.state.signals.model_change.connect(self.active_model_change)
    self.view.button_search.clicked.connect(self.search_select_dialog)
    self.view.button_cancel.clicked.connect(self.viewer.deselect_all)

    # Hierarchy Selector
    self.search_select_dialog_controller = None
    self.picking_level = "atom"
  
  



  @property
  def viewer(self):
    return self.parent

  def clear_viewer(self):
    self.parent.clear_viewer("Clearing the viewer.")


  def active_model_change(self,model_ref):
    if model_ref:
      if model_ref.data:
        if model_ref.data.filename:
          label = model_ref.data.filename
          self.view.active_model_label.setText(label)


  @Slot()
  def select_active_selection(self):
    #self.highlight_persist(None,value=True) # set the highlight persist automatically when selecting
    if self.state.active_selection_ref is None:
      self.viewer.deselect_all()
    else:
      self.viewer.select_from_query(self.state.active_selection_ref.query)

  def clear_selection(self):
    self.state.active_selection_ref = None
    self.viewer.deselect_all()

  @Slot()
  def start_selecting(self,is_checked):
    # selection button was clicked. is_checked is NEW state
    if is_checked:
      self.viewer.toggle_selection_mode(True)
    else:
      self.viewer.toggle_selection_mode(False)

  @Slot()
  def execute_selection(self):
    """
    This is a selection from the text box
    """
    if time.time()-self.last_selection_exec<5:
      return
    text = self.view.selection_edit.text()
    if text:
      selection = Selection.from_selection_string(text,model=self.state.active_model)
      passed, fail_reason = selection.is_validated, selection.fail_reason
      
      if not passed:
        self.warning({"msg":f"Unable to select: {fail_reason}","uuid": Ref._generate_uuid()})
        self.view.selection_edit.clear()
        self.view.selection_edit.setPlaceholderText(f"Unsupported selection: {fail_reason}: {text}")
        return
      if text.startswith("select"):
        text = text[7:]
      elif text.startswith("sel "):
        text = text[4:]
      try:
        selection = Selection.from_selection_string(text,model=self.state.active_model)
        self.viewer.select_from_selection(selection)
        self.viewer.focus_selected()
        self.save_text_to_history()
        self.view.selection_edit.setPlaceholderText(text)
      except:
        #raise
        self.view.selection_edit.clear()
        self.view.selection_edit.setPlaceholderText(f"Unsupported selection: syntax not currently supported: {text}")
  

  def save_text_to_history(self):
    text = self.view.selection_edit.text()
    if text:
      self.view.selection_edit.history.insert(0, text)
      self.view.selection_edit.index = -1
      self.view.selection_edit.clear()


  def search_select_dialog(self):
    dialog = SearchSelectDialog(title="Selection search", parent=self.view)
    self.view.search_select_dialog = dialog
    self.search_select_dialog_controller = SearchSelectDialogController(parent=self,view=self.view.search_select_dialog)
    default_width = dialog.width()
    new_width = default_width * 6
    dialog.setMinimumWidth(new_width)
    dialog.exec_()
