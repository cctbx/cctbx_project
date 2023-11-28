"""
The state object is a container of references. 
The origin of any signal related to changes of state.
"""
import copy
from collections import UserDict
import uuid
import json
from pathlib import Path

from PySide2.QtCore import QObject, QTimer, Signal, Slot

from .ref import Ref,ModelRef,MapRef,SelectionRef
from ...last.mol import MolDF
from ...last.selection_utils import SelectionQuery
from ...last.python_utils import DotDict

from iotbx.data_manager import DataManager

class StateSignalEmitter(QObject):
  """
  These are signals emitted from the state object (the Model in MVC). 
  It is the responsibility of the viewer to listen and implement these signals.
  """
  signal_dm_changed = Signal()
  signal_model_change = Signal()
  signal_map_change = Signal()
  signal_selection_change = Signal()
  signal_references_change = Signal() # the list of references changed, not just change in active
  signal_load_active_model = Signal()
  signal_load_active_map = Signal()

  # Style properties
  signal_iso_change = Signal(str,float) # (ref_id, iso value)
 
  signal_repr_change = Signal(str,list) # (ref_id, list of desired reprs)
  signal_viz_change = Signal(str,bool) # change visibility (ref_id, on/off)

class StateSignals(QObject):
  style_change = Signal(str)# json
  tab_change = Signal(str) # tab name TOOD: Move? Not really a state thing...
  color_change = Signal(Ref)
  restraints_change = Signal()

class State:

  @classmethod
  def from_empty(cls):
    dm = DataManager()
    return cls(dm)

  def __init__(self, data_manager):
    # Props are more complex data structures
    self._active_model_ref_id = None
    self._active_map_ref_id = None
    self._active_selection_ref_id = None
    self._data_manager = data_manager
    self._model = None
    self._map_manager = None
    self._iso = 0.5
    self.associations = {} # model: map associations
    self.references = {} # dictionary of all 'objects' tracked by the State
    self.params = DotDict()
    self.params.default_format = 'pdb'

    self.emitter= StateSignalEmitter()
    self.signals = StateSignals()

    # Initialize references

    # models
    for key in self.data_manager.get_model_names():
      self.add_model_ref_from_key(key)

    # maps
    for key in self.data_manager.get_real_map_names():
      self.add_map_ref_from_key(key)

  def add_ref(self,ref):
    assert isinstance(ref,Ref), "Must add instance of Ref or subclass"
    self.references[ref.id] = ref
    self.emitter.signal_references_change.emit()

  def remove_ref(self,ref):
    assert isinstance(ref,Ref), "Must add instance of Ref or subclass"
    del self.references[ref.id]
    self.emitter.signal_references_change.emit()

  def add_model_ref_from_key(self,model_key):
    ref = ModelRef.from_filename(self,model_key)
    self.add_ref(ref)
    
  def add_map_ref_from_key(self,map_key):
    ref = MapRef.from_filename(self,map_key)
    self.add_ref(ref)

  def add_selection_ref_from_query(self,query: SelectionQuery,show=False):
    ref = SelectionRef(self,query)
    ref._show_in_selections = show
    self.add_ref(ref)

  @property
  def references_model(self):
    return [value for key,value in self.references.items() if value.kind=="model"]
  @property
  def references_map(self):
    return [value for key,value in self.references.items() if value.kind=="map"]

  @property
  def references_selection(self):
    return [value for key,value in self.references.items() if value.kind=="selection"]

  @property
  def state(self):
    return self

  @property
  def data_manager(self):
    return self._data_manager

  @data_manager.setter
  def data_manager(self,value):
    self._data_manager = value
    self._data_manager_changed()

  @property
  def dm(self):
    # alias
    return self.data_manager

  def _data_manager_changed(self):
    # Call this to emit a signal the data_manager changed.
    # 
    # Unless putting signals in data manager, this must
    # be called if the data manager changes.
    #self.emitter.signal_dm_changed.emit()
    model_keys = [ref.key for ref in self.references_model]
    for key in self.data_manager.get_model_names():
      if key not in model_keys:
        self.add_model_ref_from_key(key)

    map_keys = [ref.key for ref in self.references_map]
    for key in self.data_manager.get_real_map_names():
      if key not in map_keys:
        self.add_map_ref_from_key(key)

#####################################
  # Models / Mols
  #####################################

  @property
  def active_model_ref(self):
    if self._active_model_ref_id is None:
      return None
    return self.references[self.active_model_ref_id]

  @active_model_ref.setter
  def active_model_ref(self,value):
    if value is None:
      self.active_model_ref_id = None
    else:
      assert isinstance(value,Ref), "Set active_model_ref with instance of Ref or subclass"
      if value not in self.references.values():
        self.add_ref(value)
      self.active_model_ref_id = value.id

  @property
  def active_model_ref_id(self):
    return self._active_model_ref_id

  @active_model_ref_id.setter
  def active_model_ref_id(self,value):
    self._active_model_ref_id = value
    print("Emitting model changed")
    self.emitter.signal_model_change.emit()
    self.signals.restraints_change.emit()


  @property
  def active_model(self):
    if self.active_model_ref is not None:
      return self.active_model_ref.model

  @property
  def model(self):
    # alias
    return self.active_model

  @property
  def active_mol(self):
    # the active mol
    if self.active_model_ref is not None:
      return self.active_model_ref.mol

  @property
  def mol(self):
    # alias
    return self.active_mol
  
  #####################################
  # Maps
  #####################################
  @property
  def active_map_ref(self):
    if self.active_map_ref_id is None:
      return None
    return self.references[self.active_map_ref_id]

  @active_map_ref.setter
  def active_map_ref(self,value):
    if value is None:
       self.active_map_ref_id = None
    else:
      assert isinstance(value,Ref), "Set active_model_ref with instance of Ref or subclass"
      if value not in self.references.values():
        self.add_ref(value)
      self.active_map_ref_id = value.id

  @property
  def active_map_ref_id(self):
    return self._active_map_ref_id

  @active_map_ref_id.setter
  def active_map_ref_id(self,value):
    print("Changing map key")
    self._active_map_ref_id = value
    self.emitter.signal_map_change.emit()

    # some logic to turn on model that is paired with map
    # if no models are on
    if self.active_model_ref_id is None:
      if self.active_model_ref_id in self.associations:
        model_key = self.associations[0]
        self.active_model_ref_id = model_key

  @property
  def active_map(self):
    if self.active_map_ref is not None:
      return self.active_map_ref.map_manager

  @property
  def map_manager(self):
    # alias
    return self.active_map

  #####################################
  # Selections
  #####################################

  @property
  def iso(self):
    return self._iso
  @iso.setter
  def iso(self,value):
    self._iso = value
    self.emitter.signal_iso_change.emit(value)

  @property
  def active_selection_ref(self):
    if self.active_selection_ref_id is None:
      return None
    return self.references[self.active_selection_ref_id]

  @active_selection_ref.setter
  def active_selection_ref(self,value):
    if value is None:
      self.active_selection_ref_id =None
    else:
      assert isinstance(value,Ref), "Set active_model_ref with instance of Ref or subclass"
      if value not in self.references.values():
        self.add_ref(value)
      self.active_selection_ref_id = value.id

  @property
  def active_selection_ref_id(self):
    return self._active_selection_ref_id

  @active_selection_ref_id.setter
  def active_selection_ref_id(self,value):
    self._active_selection_ref_id = value
    self.emitter.signal_selection_change.emit()

  @property
  def active_selection_json(self):
    if self.active_selection_ref is None:
      return None
    return self.active_selection_ref.query.to_json()

  @property
  def associations_by_id(self):
    return {key.id:value.id for key,value in self.associations.items()}

  def _guess_associations(self):
    # Try to guess map/model associations if not specified.
    # TODO: This is horrible
    for model_ref,map_ref in zip(self.references_model,self.references_map):
      self.associations[model_ref] = map_ref
      self.associations[map_ref] = [model_ref]
    
    # TODO: remove this, load maps no matter what
    for map_ref in self.references_map:
      if map_ref not in self.associations:
        model_ref = self.references_model[0]
        self.associations[model_ref] = map_ref
        self.associations[map_ref] = [model_ref]

  #####################################
  # Restraints
  #####################################