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

from .refs import Ref, DataManagerRef, ModelRef, SelectionRef


class ActiveEmitDict(dict):
  def __init__(self,parent_state):
    self.parent_state = parent_state
    super().__init__()

  def __setitem__(self, key, value):
    assert isinstance(value,Ref)
    super().__setitem__(key, value)
    self.parent_state.signals.active_change.emit(value)


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
  active_change = Signal(object) # A ref instance
  has_synced = Signal(bool)
  add_ref = Signal(Ref)
  remove_ref = Signal(Ref)


  # Selection signals
  picking_level = Signal(str) # picking granularity ('atom', 'residue')
  selection_added = Signal(SelectionRef)
  selection_activated = Signal(object)
  select_all = Signal()
  select_none = Signal()

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
    self.active_selection = None
    self._data_manager = data_manager
    self.debug = debug

    # Signals
    self.signals = StateSignals()
    self.signals.remove_ref.connect(self.remove_ref)

    
    # Python GUI state
    self.refs = {} # dictionary of all 'objects' tracked by the State
    self.active_refs = ActiveEmitDict(self) # RefClass: Ref instance (max one active ref per class)

    # Build from datamanager
    dm_ref = DataManagerRef(data=self.data_manager)
    for ref in dm_ref.get_all_refs(): # includes dm_ref
      self.add_ref(ref)


    # External molstar state
    self._molstar_state = None




  @property
  def params(self):
    return self._params

  def log(self,*args):
    if self.debug:
      print(*args)

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
    return {"molstar":[ref.identifiers["molstar"] for ref_id,ref in self.refs.items()]}

  @property
  def active_model_ref(self):
    return self._active_model_ref

  @active_model_ref.setter
  def active_model_ref(self,value):
    assert isinstance(value,(Ref,type(None)))
    self._active_model_ref = value
    self.signals.active_change.emit(value)

  @property
  def active_model(self):
    return self.active_model_ref.model

  @property
  def active_selection_ref(self):
    return self._active_selection_ref

  @active_selection_ref.setter
  def active_selection_ref(self,value):
    assert isinstance(value,(Ref,type(None)))
    self._active_selection_ref = value
    self.signals.active_change.emit(value)
    if value is None:
      self.signals.select_none.emit()

  def add_ref(self,ref,emit=True):
    """
    Add a reference object to the state. 
      If emit=True, a signal will be emitted about the change.
    """
    self.refs[ref.uuid] = ref
    ref.adopt_state(self)

  def remove_ref(self,ref):
    assert isinstance(ref,Ref), f"Must add instance of Ref or subclass, not {type(ref)}"
    #self.remove_node_and_downstream(ref) #TODO: Deleting refs could lead to dependency issues
    if ref.entry:
      ref.entry.remove()
    del self.refs[ref.uuid]

  def get_refs(self,ref_class=None):
    if not ref_class:
      return list(self.refs.values())
    else:
      return [value for value in self.refs.values() if isinstance(value,ref_class)]

  def get_active_ref(self,ref_class=ModelRef):
    if ref_class in self.active_refs:
      return self.active_refs[ref_class]

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

  @property
  def references_selection(self):
    return [ref for ref_id,ref in self.refs.items() if isinstance(ref,SelectionRef)]