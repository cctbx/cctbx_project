"""
A Ref is a top level container for data. It provides:
 1. Compositions of basic data types
 2. 'Instances' of data, where identical data appears in multiple objects
 3. An unique identifier to keep track of data
"""

from __future__ import annotations # backwards compat string literal types
from pathlib import Path
import uuid
import hashlib
import json
from dataclasses import replace
from typing import Optional

import pandas as pd
from PySide2.QtCore import QObject, Signal

from mmtbx.geometry_restraints.geo_file_parsing import add_i_seq_columns_from_id_str

from .style import Style
from ...core.selection import Selection
from .data import MolecularModelData, RealSpaceMapData
from .cif import CifFileData
from .base import DataClassBase, ObjectFrame
from .geometry import Geometry
from ...core.mol import MolDataFrame

class RefSignals(QObject):
  """
  Signals that Ref objects emit upstream to the State singleton
  """
  ref_connect = Signal(object,object) # upstream ref, downstream ref
  style_change = Signal(object,object) # Ref, style

class Ref:
  """
  The fundamental object, manages the relationship between data and defines
  the data present in the State singleton
  """
  # Signals
  signals = RefSignals()
  _class_label_name = ""
  EntryViewClass = None
  EntryControllerClass = None

  def __init__(self,data: DataClassBase,style: Optional[Style] = None, show: bool = False):
    self._data = data
    self._id = self._generate_uuid()
    self._external_ids = {}
    self._file_ref = None
    self._label = None
    self._show = show
    self.entry = None # set later the entry controller object
    self.results = {} # program_name: result_ref
    self._active = False
    self._query = None

    self._other_ref_connections = {} # Keep track of connections to other refs
    self._is_top_level = False # whether a ref was created without upstream nodes
    self._reference = None


    # Style
    assert isinstance(style,(Style,type(None)))
    if style is None:
      style = Style.from_default()
    self._style = style
    self.style_history = []
 

  @property
  def EntryControllerClass(self):
    raise NotImplementedError

  @property
  def EntryViewClass(self):
    raise NotImplementedError

  def _clear_ref_connections(self, value):
    for name,ref in self._other_ref_connections.items():
      if value == ref:
        setattr(self,name,None)

  @property
  def file_ref(self):
    return self._file_ref

  @file_ref.setter
  def file_ref(self,value):
    self._file_ref = value



  @property
  def active(self):
    return self._active

  @active.setter
  def active(self,value):
    # self.log("#"*80)
    # self.log(f"setting ref: {self} to active: ",value)
    # self.log("#"*80)

    # toggle entry as active
    if self.entry:
      self.entry.active = value

    if self.file_ref:
      self.file_ref.active = value

    self._active = value


  # def _set_active_ref(self,ref):
  #   if ref == self:
  #     self.active = True
  #   else:
  #     self.active = False

  def __setattr__(self, name, value):
    # Emit a signal if connecting two refs
    if isinstance(value, Ref):
      self.signals.ref_connect.emit(self,value)
      self._other_ref_connections[name] = value
    super().__setattr__(name, value)

  @property
  def id(self):
    return self._id

  @property
  def data(self):
    return self._data

  @property
  def external_ids(self):
    return self._external_ids

  @property
  def reference(self):
    return self._reference
    # if self.state.phenixState.references:
    #   if self.id in self.state.phenixState.references:
    #     return self.state.phenixState.references[self.id]
  @reference.setter
  def reference(self,value):
    self._reference = value

  @property
  def id_molstar(self):
    if self.reference:
      return self.reference.id_molstar


  @property
  def state(self):
    assert False, "Ref should not have direct access to state"
    # assert self._state is not None, "State was never set for this ref"
    # return self._state

  @state.setter
  def state(self,value):
    assert False, "Ref should not have direct access to state"
    # self._state = value
    # self._connections()

  @property
  def style(self):
    return self._style

  @style.setter
  def style(self,value):
    assert isinstance(value,Style), "Set with Style class"
    # sync
    #if self.state.phenixState.references and self.id in self.state.phenixState.references:
      #reference = self.state.phenixState.references[self.id]
    if self.reference:
      value = replace(value,representation=self.reference.representations)


    self.signals.style_change.emit(self,value)
    self.style_history.append(self.style)
    self._style = value

  @property
  def last_visible_style(self):
    for style in reversed(self.style_history):
      if style.visible:
        return style

  @property
  def show(self):
    return self._show

  @ show.setter
  def show(self,value):
    self._show = value

  @property
  def label(self):
    if self._label is None:
      if hasattr(self.data,"label") and self.data.label is not None:
        self._label = self.data.label

      elif hasattr(self.data,"filename") and self.data.filename is not None:
        try:
          # if key is a path, get the stem
          self._label = self._truncate_string(Path(self.data.filename).name)
        except:
          self._label = self._truncate_string(self.data.filename)
      else:
        self._label = self._class_label_name

    return self._label

  @label.setter
  def label(self,value):
    self._label = value

  @property
  def query(self):
    return self._query


  def to_dict(self):
    d = {
      "id":self.id,
      "class_name" : self.__class__.__name__,
      "external_ids": self.external_ids,
      "data":self.data.to_dict(),
      "style":self.style.to_dict(),
    }
    return d

  def to_json(self,indent=None):
    return json.dumps(self.to_dict(),indent=indent)

  def _truncate_string(self,path, max_len=20):
    if len(path) > max_len:
      return path[:max_len // 2] + "..." + path[-max_len // 2:]
    else:
      return path


  @staticmethod
  def _generate_uuid(length: int=24):
    # Generate a UUID
    full_uuid = str(uuid.uuid4())

    # Hash the UUID
    hashed_uuid = hashlib.sha1(full_uuid.encode()).hexdigest()

    # Truncate to the desired length
    short_uuid = hashed_uuid[:length]
    return short_uuid


class CifFileRef(Ref):
  _class_label_name = 'file'
  def __init__(self,data: CifFileData,show: bool = False):
    super().__init__(data=data,show=show)

  @property
  def label(self):
    if self.data.filename is not None:
      return self.data.filename
    return super().label


class GeometryRef(Ref):
  # TODO: Move most of this to the actual table not the ref, adding columns here
  #       modifies the restraints data which is bad!
  #  1. you cannot have just one chain, res, seq value per restraint... (inter-residue atoms)
  #  2. Do need to go through i_seqs? Could go straight from id_str labels to columns?
  #  3. Can we rely on i_seqs obtained here (directly from mmtbx model) to match i_seqs in sites?
  #  4. Where is the appropriate place for all/any of this to happen?
  #       - The data object constructor, but then cannot instantiate without mmtbx model and/or sites
  #       - The ref object (current version), but then misusing the Ref class? Maybe not.
  #       - The controller for the restraints table, only do when needed. But then misusing controller?
  #
  #   Most of this stuff requires sites or model anyway. The clean way would be to get i_seqs, then extract
  #     all core attrs from sites into the geometry table with appropriate suffixes.
  #   Besides the experimental features, this whole thing requires an mmtbx model as the starting point.
  #  Flow 
  #   1. get model
  #   2. fast model => sites dataframe
  #   3. fast sites dataframe => downstream dataframes
  #
  #  Locations:
  #    1. Little things like fixing names and incrementing indices can go in the data class
  #    2. Big things like converting labels <=> i_seqs or building out downstream dataframes go here (In ref)
  #    3. Its ok to modify the geometry data at the beginning, do it then. 
  #    Ultimately want this functionality:
  #      selection = get_sel_from_row(i,atom_sites)
  #      row = get_row_from_sel(sel,atom_sites)
  #       
  #      Or: 
  #        selection = sites.sel_obj_from_iseq(i)
  #        sites_sel = sites.select_from_selection(selection)
  #    
  _class_label_name = "geometry"
  def __init__(self,
      data: Geometry, 
      model_ref: 'ModelRef'):
    assert model_ref is not None
    super().__init__(data=data)
    self.model_ref = model_ref
    self.has_synchronized_with_model = False
    self.synchronize_with_model() # add individual i_seq columns from model
    self.has_incremented_column_suffixes = False
    self.increment_column_suffixes() # Make consistent per-atom suffixes (i_seq_1, i_seq_2, etc)
    self.add_iseqs_column() # Add a plural i_seqs column (column of lists)
    
    self.add_other_columns(prefix="asym_id")
    self.add_other_columns(prefix="comp_id")
    self.add_other_columns(prefix="atom_id")
    self.add_other_columns(prefix="alt_id")
    self.add_other_columns(prefix="residue_class")

    self.add_compositional_labels()

    # Set model restraint_ref value (do last, will emit signal that assumes a complete object)
    self.model_ref.geometry_ref = self

  @property
  def EntryControllerClass(self):
    from ..controller.geometry import GeoEntryController
    return GeoEntryController

  @property
  def EntryViewClass(self):
      from ..view.geometry import GeoEntryView
      return GeoEntryView

  @property
  def dfs(self):
    return self.data.dataframes

  def synchronize_with_model(self):
    # reset i_seq values from model if possible
    if self.data.has_id_str:
      for name,df in self.data.dataframes.items():
        # TODO: bug with c-beta parsing
        if name in ["c-beta","cbeta"] and 'i_seq_0' in df and df["i_seq_0"][0]=="dihedral":
          df.drop(columns=["i_seq_0"],inplace=True)

      _ = add_i_seq_columns_from_id_str(self.data.dataframes,self.model_ref.model)
    self.has_synchronized_with_model = True

  def increment_column_suffixes(self):
    for name, df in self.data.dataframes.items():
      for prefix in ["id_str", "i_seq"]:
        # Track columns to rename
        columns_to_rename = {}
        for col in df.columns:
          if col.startswith(prefix) and col.split("_")[-1].isdigit():
            parts = col.split("_")
            if parts[-1] == "0":  # Check for zero indexing
              zero_indexed = True
            else:
              zero_indexed = False
            if zero_indexed or parts[-1].isdigit():
              # Increment the suffix digit
              new_digit = int(parts[-1]) + 1
              # Reconstruct the column name with the new digit
              new_col = "_".join(parts[:-1] + [str(new_digit)])
              columns_to_rename[col] = new_col

        # Rename columns after determining all necessary changes
        df.rename(columns=columns_to_rename, inplace=True)
    self.has_incremented_column_suffixes = True

  def add_iseqs_column(self):
    # add a new column that is a list of the i_seqs
    for name,df in self.dfs.items():
      if df is not None:
        i_seq_cols = [col for col in df.columns if "i_seq" in col and col != 'i_seqs']
        if len(i_seq_cols)>0:
          df['i_seqs'] = df.apply(lambda row: [row[col] for col in i_seq_cols], axis=1)


  def add_other_columns(self,prefix):
    # prefix is a column name assumed to be in sites
    # It addes columns prefix_1, prefix_2, etc for each atom
    for name,df in self.dfs.items():
      if df is not None and self.model_ref is not None:
        i_seq_cols = [col for col in df.columns if "i_seq" in col and col != 'i_seqs']
        if len(i_seq_cols)>0:
          i_seq_suffixes = [col.replace('i_seq_','') for col in i_seq_cols]
          name_cols = [f"{prefix}_{suffix}" for suffix in i_seq_suffixes]
          for i_seq_col,name_col in zip(i_seq_cols,name_cols):
     

            df[name_col] = pd.NA
            df.reset_index(drop=True, inplace=True)
            not_na = df[i_seq_col].notna()

            indices = df.loc[not_na, i_seq_col].to_numpy()  # Extract as numpy array for direct access
            vals = self.model_ref.mol.sites[prefix].iloc[indices].to_numpy()  # Extract values as numpy array

            # Directly assign values to df where not_na is True, bypassing index-based reindexing
            df.loc[not_na, name_col] = vals

  def add_compositional_labels(self):
    for name,df in self.dfs.items():
      if df is not None and self.model_ref is not None:
        i_seq_cols = [col for col in df.columns if "i_seq" in col and col != 'i_seqs']
        if len(i_seq_cols)>0:
          i_seq_suffixes = [col.replace('i_seq_','') for col in i_seq_cols]
          name_cols = [f"Label_{suffix}" for suffix in i_seq_suffixes]
          for i_seq_col,name_col in zip(i_seq_cols,name_cols):

            df[name_col] = pd.NA
            df.reset_index(drop=True, inplace=True)
            not_na = df[i_seq_col].notna()

            i_seqs = df.loc[not_na, i_seq_col].to_numpy()  # Extract as numpy array for direct access
            sites = self.model_ref.mol.sites
            labels_compositional = sites.select_from_i_seqs(i_seqs).to_labels_compositional()

            # Directly assign values to df where not_na is True, bypassing index-based reindexing
            df.loc[not_na, name_col] = labels_compositional
  
     
class EditsRef(Ref):
  # A collection of edits of a single type
  _class_label_name = "edits"
  def __init__(self,data: ObjectFrame,geometry_ref: GeometryRef):
    super().__init__(data=data)
    self.geometry_ref = geometry_ref

class BondEditsRef(EditsRef):
  _class_label_name = "bond_edits"

class AngleEditsRef(EditsRef):
  _class_label_name = "angle_edits"

class DihedralEditsRef(EditsRef):
  _class_label_name = "dihedral_edits"

class ChiralEditsRef(EditsRef):
  _class_label_name = "chiral_edits"

class PlaneEditsRef(EditsRef):
  _class_label_name = "plane_edits"

class StyleRef(Ref):
  # A collection of edits of a single type
  _class_label_name = "style"
  def __init__(self,data: Style):
    super().__init__(data=data)

  @classmethod
  def from_default(cls):
    style = Style.from_default()
    return cls(data=style)



class RestraintsRef(Ref):
  """
  The plural reference to a collection of restraint refs
  """
  _class_label_name = "restraints"
  def __init__(self,data: Restraints, model_ref: ModelRef,show=True):
    super().__init__(data=data,show=show)
    self.model_ref = model_ref
    self.model_ref.restraints_ref = self
    self.restraints = [] # RestraintRef list


class RestraintRef(CifFileRef):
  """
  A singular reference to a restraint file
  """
  _class_label_name = "restraint"
  def __init__(self,data: Restraint, restraints_ref: Optional[RestraintsRef] = None,show=True):
    super().__init__(data=data,show=show)
    self.restraints_ref = restraints_ref


class ModelRef(Ref):
  _class_label_name = "model"
  def __init__(self,
    data: MolecularModelData,
    ciffile_ref: Optional[CifFileRef]=None,
    geometry_ref: Optional[GeometryRef] = None, 
    restraints_ref: Optional[RestraintsRef] = None, 
    show=True):
    super().__init__(data=data,show=show)
    self._mol = None
    self._geometry_ref = geometry_ref
    self._restraints_ref = restraints_ref
    self._ciffile_ref = ciffile_ref


  @property
  def EntryControllerClass(self):
    from ..controller.models import ModelEntryController
    return ModelEntryController

  @property
  def EntryViewClass(self):
      from ..view.models import ModelEntryView
      return ModelEntryView


  @classmethod
  def from_filename(cls,filename: str):
    model_data = MolecularModelData(filename=filename)
    return cls(model_data)

  @property
  def model(self):
    return self.data.model


  @property
  def ciffile_ref(self):
    return self._ciffile_ref

  @property
  def geometry_ref(self):
    return self._geometry_ref

  @geometry_ref.setter
  def geometry_ref(self,value):
    assert isinstance(value,(GeometryRef,type(None)))
    self._geometry_ref = value


  @property
  def has_geometry(self):
    return self.geometry_ref is not None


  @property
  def restraints_ref(self):
    return self._restraints_ref

  @restraints_ref.setter
  def restraints_ref(self,value):
    assert isinstance(value,(RestraintsRef,type(None)))
    self._restraints_ref = value

  @property
  def has_restraints(self):
    return self.restraints_ref is not None

  @property
  def mol(self):
    if self._mol is None:
      mol = MolDataFrame.from_mmtbx_model(self.model)
      self._mol = mol
    return self._mol

  # @property
  # def label(self):
  #   if self._label is None:
  #     if self.data.filename is not None:
  #       try:
  #         # if key is a path, get the stem
  #         self._label = self._truncate_string(Path(self.data.filename).name)
  #       except:
  #         self._label = self._truncate_string(self.data.filename)
  #     else:
  #       self._label ="model"

  #   return self._label

  @property
  def model_ref(self):
    return self
  @property
  def query(self):
    query = SelectionQuery.from_model_ref(self)
    query.params.refId = self.id
    return query

  # def _connections(self):
  #   super()._connections()
  #   self.state.signals.model_change.connect(self._set_active_ref)


class MapRef(Ref):
  _class_label_name = 'map'
  def __init__(self,data: RealSpaceMapData, model_ref: Optional[ModelRef] = None,style: Optional[Style] = None):
    super().__init__(data=data,style=style)
    self._model_ref = model_ref


  @property
  def map_manager(self):
    return self.data.map_manager

  @property
  def model_ref(self):
    return self._model_ref


  @model_ref.setter
  def model_ref(self,value):
    self._model_ref = value

  # @property
  # def label(self):
  #   if self._label is None:
  #     if self.data.filename is not None:
  #       try:
  #         # if key is a path, get the stem
  #         self._label = self._truncate_string(Path(self.data.filename).name)
  #       except:
  #         self._label = self._truncate_string(self.data.filename)
  #     else:
  #       self._label = "map"

  #   return self._label

  @property
  def query(self):
    return None

class SelectionRef(Ref):
  _class_label_name = 'selection'
  def __init__(self,data: Selection,model_ref: ModelRef,  show: Optional[bool] = True):
    assert show is not None, "Be explicit about whether to show selection ref in list of selections or not"
    assert ModelRef is not None, "Selection Ref cannot exist without reference model"
    assert data is not None, "Provide a Selection query as the data"
    super().__init__(data=data, style=model_ref.style, show=show)
    #self._debug_data = data
    # # TODO: Setting model_ref before super init is messy, can this be fixed?
    # self._model_ref = model_ref
    # self._query_sel = data
    # self._query_sel.params.refId=model_ref.id

    self._model_ref = model_ref
    self._sites_sel = None
    self._number_of_atoms = None



  @property
  def selection(self):
    return self.data

  @property
  def model_ref(self):
    return self._model_ref


  def to_json(self,indent=2):
    raise NotImplementedError
    return json.dumps(d,indent=indent)
  @property
  def sites_sel(self):
    if not self._sites_sel:
      self._sites_sel = self.model_ref.mol.select_from_selection(self.selection)
    return self._sites_sel
  @property
  def number_of_atoms(self):
    if not self._number_of_atoms:
      self._number_of_atoms = len(self.sites_sel)
    return self._number_of_atoms

# class GeometryRef(Ref):
#   _class_label_name = 'restraint'
#   def __init__(self,data: Geometry, model_ref: ModelRef):
#     super().__init__(data=data,style=model_ref.style)
#     self.model_ref = model_ref
#     self.supported_restraints = ['bond','angle']
#     for field in fields(self.data.__class__):
#       field_name = field.name
#       field_value = getattr(self.data, field_name)
#       if field_name in self.supported_restraints and field_value is not None:
#           field_value = [GeometryRef(data=data,model_ref=model_ref) for data in field_value]
#       setattr(self,field_name,field_value)

#   @property
#   def dfs(self):
#     dfs = {}
#     for field in fields(self.data.__class__):
#       field_name = field.name
#       field_value = getattr(self, field_name)
#       if field_name in self.supported_restraints and field_value is not None:
#         df = pd.DataFrame([dataclass_instance.data.__dict__ for dataclass_instance in field_value])
#         df = df.round(2)
#         # Drop some columns for now
#         columns_to_drop = ['labels', 'slack', 'rt_mx',"sym_op_i"]
#         columns_to_drop = [col for col in columns_to_drop if col in df.columns]
#         df.drop(columns=columns_to_drop, inplace=True)
#         # Sort columns
#         sort_order = ['residual',"delta",'ideal','model','i_seq']
#         matched_columns = []
#         for string in sort_order:
#             matched_columns += [col for col in df.columns if string in col]

#         # append the remaining columns
#         remaining_columns = [col for col in df.columns if col not in matched_columns]
#         new_column_order = matched_columns + remaining_columns
#         df = df[new_column_order]

#         # sort rows
#         df = df.sort_values(by='residual',ascending=False)

#         dfs[field_name] = df
#     return dfs


#   @classmethod
#   def from_model_ref(cls,model_ref):
#     data = Geometry.from_model_ref(model_ref)
#     return cls(data=data,model_ref=model_ref)



# class QscoreRef(ModelResultsRef):
#     def __init__(self,data: QscoreResult,model_ref: ModelRef, selection_ref: Optional[SelectionRef] = None):
#       super().__init__(data=data,model_ref=model_ref,selection_ref=selection_ref)
#       self.model_ref.results["qscore"] = self



  # def _connections(self):
  #   self.state.signals.ciffile_change.connect(self._set_active_ref)

  # @property
  # def label(self):
  #   if self._label is None:
  #     if self.data.filename is not None:
  #       try:
  #         # if key is a path, get the stem
  #         self._label = self._truncate_string(Path(self.data.filename).name)
  #       except:
  #         self._label = self._truncate_string(self.data.filename)
  #     else:
  #       self._label = "file"

  #   return self._label

class ResultsRef(Ref):
  def __init__(self,data: Result):
    super().__init__(data=data)

    for model_ref in self.data.model_refs:
      model_ref.results[self.data.program_name] = self
    for map_ref in self.data.map_refs:
      map_ref.results[self.data.program_name] = self