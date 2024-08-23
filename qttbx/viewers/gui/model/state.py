"""
The state object is a singleton container of references. 
It also functions as a "signal bus" for Qt Signals and Slots. This allows controllers to communicate with
each other in a way that is asynchronous and limited in scope.
"""
from pathlib import Path
import json
import os
from dataclasses import dataclass
from typing import Dict

from PySide2.QtCore import QObject, Signal
from PySide2.QtWidgets import QMessageBox

from iotbx.data_manager import DataManager

from .base import DataClassBase
from .molstar import Reference, Structure, Component, Representation
from .ref import *

from .selection import Selection
from .data import MolecularModelData
from .molstar import MolstarState



class StateSignals(QObject):
  """
  The 'signal bus' for the state object. The signals defined here are the preferred way
  for controllers to communicate with each other. 
  s
  Importantly, the State instance can notify all controllers about any changes to state.
  """
  # Messages
  warning = Signal(dict) # A warning message
  notification = Signal(dict) # A notification
  error = Signal(dict) # An error
  tab_change = Signal(str) 
  model_change = Signal(object) # model ref
  data_change = Signal() # The data manager has changed
  has_synced = Signal(bool)
  remove_ref = Signal(object) # ref

  # Selection signals
  select_all = Signal(bool)
  deselect_all = Signal(bool)
  selection_added = Signal(object) # selection ref *
  selection_activated = Signal(object)
  selection_deactivated = Signal(object)

  selection_hide = Signal(object)
  selection_show = Signal(object)
  selection_focus= Signal(object)
  selection_rep_show = Signal(object,str) # Ref, representation_name

  color_selection = Signal(object,object) # Ref, Color

class State:
  """
  The sole place where GUI state exists. Which is mostly just a collection of references to data (Refs).

  A typical program flow example would be a controller modifies the state singleton (single instance),
    say by adding a new SelectionRef. The state then broadcasts a signal that the state has changed.
    Any additional controller that cares about selection references can then respond to that signal with
    a connected method.

  Ideally the state object could be written to disk and then reloaded later to give the same GUI state.
  """
  @classmethod
  def from_empty(cls):
    dm = DataManager()
    return cls(dm)

  def __init__(self, data_manager,params=None,log=None,debug=False):
    self._log = log
    self._params = params
    self._active_model_ref = None
    self._active_selection_ref = None
    self._data_manager = data_manager
    self.debug = debug

    # Python GUI state
    self.references = {} # dictionary of all 'objects' tracked by the State

    # External molstar state
    self._molstar_state = MolstarState(references={},has_synced=False)


    # Signals
    self.signals = StateSignals()
    self.signals.data_change.connect(self._data_manager_changed)
    self.signals.remove_ref.connect(self.remove_ref)

  def init_from_datamanager(self):
    # models
    if hasattr(self.data_manager,"get_model_names"):
      for name in self.data_manager.get_model_names():
        model = self.data_manager.get_model(filename=name)
        self.add_ref_from_mmtbx_model(model,filename=name)
        # set the first model in dm as active
        if self.active_model_ref is None:
          self.active_model_ref = self.references_model[0]



  @property
  def params(self):
    return self._params

  def log(self,*args):
    if self.debug:
      print(*args)

  
  @property
  def mmcif_column_map(self):
    return self.parameters.core_map_to_mmcif

  @property
  def molstarState(self):
    return self._molstar_state

  @molstarState.setter
  def molstarState(self,molstarState):
    self._molstar_state = molstarState
    if isinstance(molstarState,MolstarState): # we got a valid phenix state
      if molstarState.has_synced:
        self.signals.has_synced.emit(True)

  @property
  def external_loaded(self):
    return {"molstar":[ref.identifiers["molstar"] for ref_id,ref in self.references.items()]}


  def update_from_remote_dict(self,remote_state_dict):
    for ref_id,ref in remote_state_dict["references"].items():
      if ref_id in self.references:
        local_ref = self.state.references[ref_id]
        local_ref.external_ids.update(ref["external_ids"])
        local_ref.style.update(**ref["style"])


  def start(self):
    # get data out of datamanager
    self.init_from_datamanager()
    # Run this after all controllers are initialized
    self.signals.model_change.emit(self.active_model_ref)


  def add_ref(self,ref,emit=True):
    """
    Add a reference object to the state. 
      If emit=True, a signal will be emitted about the change.
    """

    assert isinstance(ref,Ref), "Must add instance of Ref or subclass"
    self.references[ref.uuid] = ref

    if isinstance(ref,ModelRef):
      #self.active_model_ref = ref
      #self.signals.model_change.emit(ref)

      # Add 'all' selection ref
      selection = Selection.from_selection_string("all",model=ref.data.model)
      selection_ref = SelectionRef(selection,ref,show=True)
      self.add_ref(selection_ref)

    elif isinstance(ref,SelectionRef):
      if emit:
        self.signals.selection_added.emit(ref)
    
    else:
      raise ValueError(f"ref provided not among those expected: {ref}")



  def remove_ref(self,ref):
    assert isinstance(ref,Ref), f"Must add instance of Ref or subclass, not {type(ref)}"
    #self.remove_node_and_downstream(ref) #TODO: Deleting refs could lead to dependency issues
    if ref.entry:
      ref.entry.remove()
    del self.references[ref.uuid]
  

  def add_ref_from_phenix_selection_string(self,phenix_string):
    selection = Selection.from_selection_string(phenix_string)
    return self.add_ref_from_selection(selection)

  def add_ref_from_selection(self,selection,model_ref=None,show=True,make_active=False):
    if model_ref is None:
      model_ref = self.active_model_ref
    ref = SelectionRef(data=selection,model_ref=model,show=show)
    self.add_ref(ref)
    if make_active:
      self.active_selection_ref = ref
    return ref


  def add_ref_from_mmtbx_model(self,model,label=None,filename=None):
    if model.get_number_of_atoms()==0:
      self.log("Model with zero atoms, not added")
      return

    if filename is not None:
      filepath = Path(filename).absolute()

    # add the model
    data = MolecularModelData(filepath=filepath,model=model)
    model_ref = ModelRef(data=data,show=True)
    self.add_ref(model_ref)
    return model_ref


  @property
  def references_model(self):
    return [value for key,value in self.references.items() if isinstance(value,ModelRef)]

  @property
  def references_selection(self):
    return [value for key,value in self.references.items() if isinstance(value,SelectionRef)]

  @property
  def state(self):
    return self

  @property
  def data_manager(self):
    return self._data_manager

  @property
  def dm(self):
    # alias
    return self.data_manager

  def _data_manager_changed(self):
    # Call this to emit a signal the data_manager changed.
    #
    # Unless putting signals in data manager, this must
    # be called explicitly/manually if the data manager changes.
    #self.signals.data_manager_changed.emit()
    model_refs = [ref for ref in self.references_model]
    model_keys = [str(ref.data.filepath) for ref in model_refs]
    for filename in self.data_manager.get_model_names():
      if str(filename) not in model_keys:
        self.log(f"New file found in data manager: {filename} and not found in references: {model_keys}")
        model = self.data_manager.get_model(filename=filename)
        ref = self.add_ref_from_mmtbx_model(model,filename=filename)

        self.signals.model_change.emit(self._active_model_ref) # No change, just trigger update


  #####################################
  # Models / Mols
  #####################################

  @property
  def active_model_ref(self):
    return self._active_model_ref

  @active_model_ref.setter
  def active_model_ref(self,value):

    if value == self._active_model_ref:
      return # do nothing if already the active ref

    
    if value is None:
      self._active_model_ref = None
    else:
      assert isinstance(value,ModelRef), "Set active_model_ref with instance of Ref or subclass"
      assert value in self.references.values(), "Cannot set active ref before adding to state"
      self._active_model_ref = value


      # emit signals
      self.signals.model_change.emit(self.active_model_ref)
      
  @property
  def active_model(self):
    if self.active_model_ref is not None:
      return self.active_model_ref.model



  #####################################
  # Selections
  #####################################

  @property
  def active_selection_ref(self):
    return self._active_selection_ref

  @active_selection_ref.setter
  def active_selection_ref(self,value):
    """
    Assigning a SelectionRef to be 'active' emits a signal handled by the graphics controller

    """
    if value is None:
      self._active_selection_ref =None
      self.signals.deselect_all.emit(True)
    else:
      assert isinstance(value,SelectionRef), "Set active_model_ref with instance of Ref or subclass"
      assert value in self.references.values(), "Cannot set active ref before adding to state"
      self._active_selection_ref = value
      self.signals.selection_activated.emit(value)
      self.signals.selection_focus.emit(value)

  