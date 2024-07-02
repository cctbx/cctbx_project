"""
The state object is a container of references.
The origin of any signal related to changes of state.
"""
from pathlib import Path
import json
import os
from dataclasses import dataclass
from typing import Dict

import networkx as nx
from PySide2.QtCore import QObject, Signal
from PySide2.QtWidgets import QMessageBox

from iotbx.data_manager import DataManager

from .base import DataClassBase
from .reference import Reference
from .structure import Structure
from .component import Component
from .representation import Representation
from .reference import Reference
from .ref import Ref,ModelRef,MapRef,SelectionRef, GeometryRef,  CifFileRef, EditsRef, RestraintRef, RestraintsRef, ResultsRef
from ...core.python_utils import DotDict
from ...core.selection import Selection
from .data import MolecularModelData, RealSpaceMapData
from .cif import CifFileData
from ...core.parameters import params


class StateSignals(QObject):
  style_change = Signal(object,str)# (ref,style_json)
  tab_change = Signal(str) # tab name TOOD: Move? Not really a state thing...
  color_change = Signal(Ref)
  model_change = Signal(object) # model ref
  map_change = Signal(object) # map ref
  data_change = Signal() # The data manager has changed
  geometry_change = Signal(object) # restraits ref


  ### Start new
  # Generic state changes
  new_active_ref = Signal(object) # Any ref type, now active

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

  # Geometry signals
  edits_added = Signal(object) # edits ref
  geometry_filter_from_restraint = Signal(object) # a FilterObj, apply to geometry filters

  # Restraints signals
  restraints_change = Signal(object) # Restraint ref
  restraint_activated = Signal(object)

  # Cif file signals
  #ciffile_change = Signal(object)


  #  End new
  references_change = Signal() # generic, TODO: refactor out
  results_change = Signal(object)
  ciffile_change = Signal(object) # cif file ref
  geo_change = Signal(object) # geo file ref
  repr_change = Signal(str,list) # (ref_id, list of desired reprs)
  viz_change = Signal(str,bool) # change visibility (ref_id, on/off)
  picking_level = Signal(str) # one of 'residue', 'atom'
  sync = Signal()
  clear = Signal(str) # reload all active objects, send the messagebox message
  #select = Signal(object) # select a ref object
  remove_ref = Signal(object) # ref
  update = Signal(object)
  stage_restraint = Signal(object,str) # send selection ref,type to be staged as restraint
  filter_update = Signal(object,bool) # emit a filter controller to apply, with debug flag
  #select_ref = Signal(object)
  



@dataclass(frozen=True)
class PhenixState(DataClassBase):
  references: Dict[str,Reference]

  @classmethod
  def from_dict(cls,state_dict):
    phenix_state =  cls(references = {})
    # update state from dict
    for ref_dict in state_dict["references"]:

      id_molstar = ref_dict["molstarKey"]
      id_viewer = ref_dict["phenixKey"]


      ref = Reference(id_molstar=id_molstar, id_viewer = id_viewer, structures = [])
      phenix_state.references[id_viewer] = ref
      # now modify within a ref
      for structure_dict in ref_dict['structures']:
        structure = Structure(
          phenixKey=structure_dict["phenixKey"],
          phenixReferenceKey=structure_dict['phenixReferenceKey'],
          data_id=structure_dict["data_id"],
          key=structure_dict["key"],
          components=[])
        ref.structures.append(structure)
        for component_dict in structure_dict['components']:
          component = Component(phenixKey=component_dict["phenixKey"],representations=[],key=component_dict["key"])
          structure.components.append(component)

          for representation_dict in component_dict["representations"]:
            representation = Representation(phenixKey=representation_dict["phenixKey"],name=representation_dict["name"])

            component.representations.append(representation)

    return phenix_state

class State:

  @classmethod
  def from_empty(cls):
    dm = DataManager()
    return cls(dm)

  def __init__(self, data_manager,log=None):
    self.log = log
    self._active_model_ref = None
    self._active_map_ref = None
    self._active_selection_ref = None
    #self._active_ciffile_ref = None
    self._data_manager = data_manager
    self._model = None
    self._map_manager = None
    self._has_synced = False
    self._phenix_state = PhenixState(references={})
    #self.is_updating = True

    # Graph
    self.G = nx.DiGraph()

      #self.associations = {} # model: map associations
    self.references = {} # dictionary of all 'objects' tracked by the State
    #self.external_loaded = defaultdict(list) # external name: [internal ref_ids]
    self.params = params

    # Signals
    self.signals = StateSignals()
    self.signals.data_change.connect(self._data_manager_changed)
    self.signals.remove_ref.connect(self.remove_ref)
    Ref.signals.ref_connect.connect(self._ref_connect_graph)
    Ref.signals.ref_connect.connect(self._ref_connect_signals)
    Ref.signals.style_change.connect(self.apply_style)

  def init_from_datamanager(self):
    # models
    for name in self.data_manager.get_model_names():
      model = self.data_manager.get_model(filename=name)
      self.add_ref_from_mmtbx_model(model,filename=name)
      # set the first model in dm as active
      if self.active_model_ref is None:
        self.active_model_ref = self.references_model[0]

    # maps
    for name in self.data_manager.get_real_map_names():
      map_manager = self.data_manager.get_real_map(filename=name)
      label = os.path.basename(name)
      self.add_ref_from_map_manager(map_manager=map_manager,filepath=name,label=label)

      # set the first map in dm as active
      if self.active_map_ref is None:
        self.active_map_ref = self.references_map[0]


  def apply_style(self,ref,style):
    self.signals.style_change.emit(ref,style.to_json())

  def _ref_connect_graph(self,ref_upstream,ref_downstream):
    """
    Two refs are being connected, update nx graph
    """
    if ref_upstream not in self.state.G:
      self.state.G.add_node(ref_upstream,id=ref_upstream.id)
    if ref_downstream not in self.state.G:
      self.state.G.add_node(ref_downstream,id=ref_downstream.id)
    self.state.G.add_edge(ref_upstream,ref_downstream)
    ref_upstream._is_top_level = self.G.in_degree(ref_upstream) == 0

  def _ref_connect_signals(self,ref_upstream,ref_downstream):
    """
    Two refs are being connected, determine which signals to emit
    """
    # adding restraints to model
    if (isinstance(ref_upstream,ModelRef) and
        isinstance(ref_downstream,GeometryRef)):
      self.signals.geometry_change.emit(ref_downstream)
    elif (isinstance(ref_upstream,ModelRef) and
        isinstance(ref_downstream,CifFileRef)):
      self.signals.ciffile_change.emit(ref_downstream)


  def remove_node_and_downstream(self, start_node):
    graph = self.G

    # Step 1: Find all nodes downstream of the start_node (including start_node)
    downstream_nodes = set(nx.dfs_preorder_nodes(graph, source=start_node))

    # Step 2: Identify nodes to remove - those that are not "top level"
    nodes_to_remove = set()
    for node in downstream_nodes:
      # Check if the node has incoming edges from nodes not in downstream_nodes
      has_external_incoming_edges = any(
        pred not in downstream_nodes for pred in graph.predecessors(node)
      )
      # If a node has no incoming edges from outside downstream_nodes, mark it for removal
      if not has_external_incoming_edges and not node._is_top_level:
        nodes_to_remove.add(node)

    # Collect edges that will be removed
    edges_to_remove = set()
    for node in nodes_to_remove:
      # Incoming edges
      edges_to_remove.update(graph.in_edges(node))
      # Outgoing edges
      edges_to_remove.update(graph.out_edges(node))

    # Perform checks on edges_to_remove before actually removing them
    for upstream_ref,downstream_ref in edges_to_remove:
      upstream_ref._clear_ref_connections(downstream_ref)

    # Optionally, based on checks, decide whether to proceed with node removal
    # For now, we'll proceed with removing the identified nodes
    graph.remove_nodes_from(nodes_to_remove)
    for node in nodes_to_remove:
      if node.id in self.references:
        del self.references[node.id]

  @property
  def mmcif_column_map(self):
    from ...core.selection_utils import core_keys_to_mmcif_keys_default
    return core_keys_to_mmcif_keys_default

  @property
  def phenixState(self):
    return self._phenix_state

  @phenixState.setter
  def phenixState(self,phenixState):
    self._phenix_state = phenixState
    if isinstance(phenixState,PhenixState): # we got a valid phenix state
      # Sync with remote
      self.has_synced = True
    #   for ref_id,ref in self.references.items():
    #     if ref_id in phenixState.references:
    #       ref.reference = phenixState.references[ref_id]
    #       ref.style = replace(ref.style,representation=ref.reference.representations)

  @property
  def external_loaded(self):
    return {"molstar":[ref.id_molstar for ref_id,ref in self.references.items()]}

  @property
  def has_synced(self):
    return self._has_synced

  @has_synced.setter
  def has_synced(self,value):
    self._has_synced = value


  def to_dict(self):
    d= {
      "class_name": self.__class__.__name__,
      "references": {ref_id: ref.to_dict() for ref_id,ref in self.references.items()},
      #"external_loaded": self.external_loaded,
      #"active_model_ref": self.active_model_ref.id if self.active_model_ref else None,
      #"active_map_ref": self.active_map_ref.id if self.active_map_ref else None,
      #"active_selection_ref": self.active_selection_ref.id if self.active_selection_ref else None
    }
    return d


  def to_json(self,indent=None):
    d = self.to_dict()
    return json.dumps(d,indent=indent)

  def show(self):
    print(self.to_json(indent=2))

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
    self.signals.map_change.emit(self.active_map_ref)

  def notify(self,message):
    # Run a notification popup that must be click closed
    msg = QMessageBox()
    msg.setWindowTitle("Notification")
    msg.setText(message)
    msg.setIcon(QMessageBox.Information)
    msg.setStandardButtons(QMessageBox.Ok)
    msg.exec_()

  def add_ref(self,ref,emit=True):
    """
    Add a reference object to the state. 
      If emit=True, a signal will be emitted about the change.
    """

    assert isinstance(ref,Ref), "Must add instance of Ref or subclass"
    self.references[ref.id] = ref

    if isinstance(ref,ModelRef):
      #self.active_model_ref = ref
      #self.signals.model_change.emit(ref)

      # Add 'all' selection ref
      selection = Selection.from_phenix_string("all")
      selection_ref = SelectionRef(selection,ref)
      self.add_ref(selection_ref)

      # Optionally add a cif file ref
      if ref.file_ref is not None:
        self.add_ref(ref.file_ref)

    elif isinstance(ref,MapRef):
      #self.active_map_ref = ref
      #self.signals.model_change.emit(ref)
      pass
    elif isinstance(ref,SelectionRef):
      if emit:
        self.signals.selection_added.emit(ref)
    elif isinstance(ref,GeometryRef):
      if emit:
        self.signals.geo_change.emit(ref)

    elif isinstance(ref,RestraintRef):
      # Singular
      if emit:
        self.signals.restraints_change.emit(ref)

    elif isinstance(ref,RestraintsRef):
      # Plural, add sub refs then emit
      for restraint in ref.data.restraints:
        sub_ref = RestraintRef(data=restraint,restraints_ref=ref,show=True)
        ref.restraints.append(sub_ref)
        self.add_ref(sub_ref,emit=False) # don't emit until the end

      # emit at the end
      self.signals.restraints_change.emit(ref)

    elif isinstance(ref,(ResultsRef)):
      if emit:
        self.signals.results_change.emit(ref)

    elif isinstance(ref,(CifFileRef)):
      if emit:
        self.signals.ciffile_change.emit(ref)

    elif isinstance(ref,EditsRef):
      if emit:
        self.signals.edits_added.emit(ref)

    else:
      raise ValueError(f"ref provided not among those expected: {ref}")
    # # update other controllers
    # if emit:
    #   self.signals.update.emit(ref)






  def remove_ref(self,ref):
    assert isinstance(ref,Ref), f"Must add instance of Ref or subclass, not {type(ref)}"
    #self.remove_node_and_downstream(ref) #TODO: Deleting refs could lead to dependency issues
    if ref.entry:
      ref.entry.remove()
    del self.references[ref.id]

  def add_ref_from_phenix_selection_string(self,phenix_string):
    selection = Selection.from_phenix_string(phenix_string)
    return self.add_ref_from_selection(selection)

  def add_ref_from_selection(self,selection,model_ref=None,show=True,make_active=False):
    if model_ref is None:
      model_ref = self.active_model_ref
    ref = SelectionRef(data=selection,model_ref=model_ref,show=show)
    self.add_ref(ref)
    if make_active:
      self.active_selection_ref = ref
    return ref

  def add_ref_from_model_file(self,filename=None,label=None,format=None):
    if label is None:
      label = filename
      #label = "model_"+str(id(filename))+str(time.time())
    dm = DataManager()
    dm.process_model_file(filename)
    model = dm.get_model()
    return self.add_ref_from_mmtbx_model(model,label=label,filename=None)

  def add_ref_from_model_string(self,model_string,label=None,format=None):
    #assert label is not None
      #label = "model_"+str(id(model_string))+str(time.time())
    dm = DataManager()
    dm.process_model_str(label,model_string)
    model = dm.get_model()

    return self.add_ref_from_mmtbx_model(model,label=label,filename=None)

  def add_ref_from_mmtbx_model(self,model,label=None,filename=None):
    if model.get_number_of_atoms()==0:
      print("Model with zero atoms, not added")
      return

    if filename is not None:
      filepath = Path(filename).absolute()


    # Add the file
    data = CifFileData(filepath=filepath)
    cif_ref = CifFileRef(data=data,show=True)
    self.add_ref(cif_ref)

    # add the model
    data = MolecularModelData(filepath=filepath,model=model)
    model_ref = ModelRef(data=data,ciffile_ref=cif_ref,show=True)
    self.add_ref(model_ref)



    return model_ref


  def add_ref_from_map_manager(self,filepath=None,map_manager=None,volume_id=None,model_id=None,label=None):
    data = RealSpaceMapData(filepath=filepath,map_manager=map_manager,label=label)
    ref = MapRef(data=data,model_ref=None)
    self.add_ref(ref)
    return ref


  def add_ref_from_map_file(self,filename=None,volume_id=None,model_id=None,label=None):
    # filepath = None
    # name = filename
    # if name is None:
    #   if map_manager.file_name not in ["",None]:
    #     filename = Path(map_manager.file_name)
    #     if filename.exists():
    #       name = filename
    if filename is not None:
      filepath = str(Path(filename).absolute())

    self.data_manager.process_real_map_file(filename=filename)
    map_manager = self.data_manager.get_real_map(filename=filename)

    data = RealSpaceMapData(filepath=filepath,map_manager=map_manager,label=label)
    ref = MapRef(data=data,model_ref=None)
    self.add_ref(ref)
    return ref

  @property
  def references_model(self):
    return [value for key,value in self.references.items() if isinstance(value,ModelRef)]
  @property
  def references_map(self):
    return [value for key,value in self.references.items() if isinstance(value,MapRef)]

  @property
  def references_selection(self):
    return [value for key,value in self.references.items() if isinstance(value,SelectionRef)]

  @property
  def references_ciffile(self):
    return [value for key,value in self.references.items() if isinstance(value,CifFileRef)]

  @property
  def references_edits(self):
    return [value for key,value in self.references.items() if isinstance(value,EditsRef)]

  @property # Singular
  def references_restraint(self):
    return [value for key,value in self.references.items() if isinstance(value,RestraintRef)]
  
  @property # Plural
  def references_restraints(self):
    return [value for key,value in self.references.items() if isinstance(value,RestraintsRef)]

  @property
  def references_geo(self):
    return [value for key,value in self.references.items() if isinstance(value,GeometryRef)]

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
    # be called explicitly/manually if the data manager changes.
    #self.signals.data_manager_changed.emit()
    model_refs = [ref for ref in self.references_model]
    model_keys = [str(ref.data.filepath) for ref in model_refs]
    for filename in self.data_manager.get_model_names():
      if str(filename) not in model_keys:
        print(f"New file found in data manager: {filename} and not found in references: {model_keys}")
        model = self.data_manager.get_model(filename=filename)
        ref = self.add_ref_from_mmtbx_model(model,filename=filename)

        self.signals.model_change.emit(self._active_model_ref) # No change, just trigger update


    map_refs = [ref for ref in self.references_map]
    map_keys = [ref.data.filepath for ref in map_refs]
    for filename in self.data_manager.get_real_map_names():
      if filename not in map_keys:
        print(f"New file found in data manager: {filename} and not found in references: {map_keys}")
        map_manager = self.data_manager.get_real_map(filename=filename)
        self.add_ref_from_map_manager(map_manager=map_manager,filepath=filename,label=os.path.basename(filename))
        self.signals.map_change.emit(self.active_map_ref) # No change, just trigger update
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
      assert isinstance(value,Ref), "Set active_model_ref with instance of Ref or subclass"
      assert value in self.references.values(), "Cannot set active ref before adding to state"
      self._active_model_ref = value

      # check if cif file for cif editor
      if self.active_model_ref.file_ref is not None:
        if isinstance(self.active_model_ref.file_ref,CifFileRef):
          self.signals.ciffile_change.emit(self.active_model_ref.file_ref)

      # emit signals
      self.signals.model_change.emit(self.active_model_ref)

    self.signals.references_change.emit()

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

  @property
  def sites(self):
    # alias
    return self.mol.sites
  #####################################
  # Maps
  #####################################
  @property
  def active_map_ref(self):
    return self._active_map_ref

  @active_map_ref.setter
  def active_map_ref(self,value):
    if value is None:
      self._active_map_ref = None
    else:
      assert isinstance(value,Ref), "Set active_model_ref with instance of Ref or subclass"
      assert value in self.references.values(), "Cannot set active ref before adding to state"
      self._active_map_ref = value

      # signals
      self.signals.map_change.emit(self.active_map_ref)
    self.signals.references_change.emit()


  @property
  def active_map(self):
    if self.active_map_ref is not None:
      return self.active_map_ref.map_manager

  # def _guess_associations(self):
  #   # Try to guess map/model associations if not specified.
  #   # TODO: This is horrible
  #   for model_ref,map_ref in zip(self.references_model,self.references_map):
  #     self.associations[model_ref] = map_ref
  #     self.associations[map_ref] = model_ref

  #   # TODO: remove this, load maps no matter what
  #   for map_ref in self.references_map:
  #     if map_ref not in self.associations:
  #       model_ref = self.references_model[0]
  #       self.associations[model_ref] = map_ref
  #       self.associations[map_ref] = model_ref

  #####################################
  # Selections
  #####################################


  @property
  def active_selection_ref(self):
    return self._active_selection_ref

  @active_selection_ref.setter
  def active_selection_ref(self,value):
    """
    Assigning a SelectionRef to be 'active' emits a signal.

    """
    if value is None:
      self._active_selection_ref =None
    else:
      assert isinstance(value,Ref), "Set active_model_ref with instance of Ref or subclass"
      assert value in self.references.values(), "Cannot set active ref before adding to state"
      self._active_selection_ref = value
      self.signals.selection_activated.emit(value)
      self.signals.selection_focus.emit(value)
    #self.signals.references_change.emit()

  #####################################
  # Geometry
  #####################################

  # @property
  # def active_restraint_ref(self):
  #   return self.active_model_ref.restraints


  #####################################
  # Cif Files
  #####################################

  # @property
  # def active_ciffile_ref(self):
  #   return self._active_ciffile_ref

  # @active_ciffile_ref.setter
  # def active_ciffile_ref(self,value):

  #   ref = value
  #   if value is None:
  #     self._active_ciffile_ref = None
  #   else:
  #     assert isinstance(value,CifFileRef), "Set active_model_ref with instance of Ref or subclass"
  #     assert value in self.references.values(), "Cannot set active ref before adding to state"
  #     self._active_ciffile_ref = value

  #     # control active flag
  #     ref.active = True
  #     for key,other_ref in self.references.items():
  #       if isinstance(ref,CifFileRef):
  #         if other_ref != ref and other_ref.active:
  #           other_ref.active = False

  #     # Emit signals
  #     self.signals.ciffile_change.emit(self.active_ciffile_ref)
  #   self.signals.references_change.emit()
