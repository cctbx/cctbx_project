import subprocess
import platform
from functools import partial
from dataclasses import replace

from PySide2.QtCore import QObject
from PySide2.QtWidgets import QMessageBox

from .controller import Controller
from ..state.ref import ModelRef, MapRef


class ScrollEntryController(Controller, QObject):
  """
  Controls views/scroll_entry.py
  """
  def __init__(self,parent=None,view=None,ref=None,show=True):
    QObject.__init__(self) # QObject
    Controller.__init__(self,parent=parent,view=view)
    assert ref is not None, f"An Entry is a gui analog to a Ref. Must provide ref"
    self.ref = ref
    self.ref.entry = self
    self._is_active = ref.active
    self.show = True
    self._is_destroyed = False

    # Connect signals from view

    # Active
    self.view.active_toggle.stateChanged.connect(self._toggle_active_func)


    # Close
    self.view.button_close.clicked.connect(self.remove_entry)

    # Display name
    self.view.label_name.setText(self.ref.label)


  @property
  def is_active(self):
    return self._is_active

  @is_active.setter
  def is_active(self,value):
    if not self.view.is_destroyed:
      assert isinstance(value,bool), "Active must be boolean"
      self._is_active = value
      self.view.active_toggle.is_checked=value

  @property
  def is_destroyed(self):
    return self._is_destroyed

  @is_destroyed.setter
  def is_destroyed(self,value):
    self._is_destroyed = value
    self.view.is_destroyed = value

  @property
  def parent_list(self):
    return self.parent


  def remove_entry(self):
    print("removing entry for ref: ",self.ref.id)
    # TODO: move all  this to state
    if hasattr(self.ref,"show_in_list"):
      self.ref.show_in_list = False
    

    # remove from active
    if (self.ref == self.state.active_model_ref):
      self.state.active_model_entry = None
    if (self.ref == self.state.active_map_ref):
      self.state.active_map_entry = None
    if (self.ref == self.state.active_selection_ref):
      self.state.active_selection_entry = None

    # delete from data manager
    if isinstance(self.ref,ModelRef):
      self.state.data_manager.remove_model(self.ref.data.filepath)
    if isinstance(self.ref,MapRef):
      self.state.data_manager.remove_real_map(self.ref.data.filepath)


    # delete children
    
    for ref_id,ref in list(self.state.references.items()):
      print("Ref: ",ref,ref.id)
      if not isinstance(ref,ModelRef):
        if hasattr(ref,"model_ref") and ref.model_ref == self.ref:
          print(ref,ref.model_ref.id,self.ref.id)
          if ref.entry in ref.entry.parent_list.entries:
            ref.entry.parent_list.remove_entry(ref.entry)
    if self.ref.id in self.state.references:
      del self.state.references[self.ref.id]

    # harsh reset of viewer, the problem is that the viewer will send back old pairings if same model
    self.parent_list.remove_entry(self) # remove from gui
    self.state.signals.clear.emit("Currently, closing an entry triggers a full viewer reset.")



  def _toggle_active_func(self,is_checked):
    self.ref.active = is_checked
    if self.view._is_destroyed:
        return
    else:
      return self.toggle_active_func(is_checked)

  def toggle_active_func(self,is_checked):
    # implement for subclasses. Called when toggle is switched
    raise NotImplementedError

  


    