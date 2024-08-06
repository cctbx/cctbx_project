from dataclasses import dataclass, fields, field
from typing import List, Optional
import tempfile
from pathlib import Path

import pandas as pd
import numpy as np
from mmtbx.geometry_restraints.geo_file_parsing import add_i_seq_columns_from_id_str
from mmtbx.geometry_restraints.geo_file_parsing import (
  parse_geo_file
)

from .base import DataClassBase, ObjectFrame

geometry_names = [
      "bond",
      "angle",
      "dihedral",
      "chirality",
      "cbeta",
      "plane",
      "nonbonded",
    ]

@dataclass(frozen=True)
class Geo(DataClassBase):
  #labels: Optional[List[str]] = None
  i_seqs: Optional[List[int]] = None
  #id_strs: Optional[List[str]] = None

@dataclass(frozen=True)
class NonBondedGeometry(Geo):
  name = 'nonbonded'

@dataclass(frozen=True)
class BondGeometry(Geo):
  name = 'bond'
  ideal: Optional[float] = 0.0
  sigma: Optional[float] = 1.0
  weight: Optional[float] = 0.0
  slack: Optional[float] = 0.0
  sym_op_j: Optional[object] = None
  rt_mx: Optional[object] = None



@dataclass(frozen=True)
class AngleGeometry(Geo):
  name = "angle"
  ideal: Optional[float] = 0.0
  sigma: Optional[float] = 0.0
  weight: Optional[float] = 0.0

@dataclass(frozen=True)
class DihedralGeometry(Geo):
  name = "dihedral"
  ideal: Optional[float] = 0.0
  sigma: Optional[float] = 0.0
  weight: Optional[float] = 0.0
  harmonic: Optional[float] = 0.0

  @property
  def sinusoidal(self):
    return self.harmonic

@dataclass(frozen=True)
class ChiralGeometry(Geo):
  name = "chirality"
  ideal: Optional[float] = 0.0
  sigma: Optional[float] = 0.0
  weight: Optional[float] = 0.0
  both_signs: Optional[bool] = False

@dataclass(frozen=True)
class PlaneGeometry(Geo):
  name = "plane"
  default_weight = 2500.0
  weights: Optional[List[float]] = field(default_factory=lambda: [0.0])

@dataclass(frozen=True)
class Geometry(DataClassBase):
  # Container for geometry dataframes
  filepath: Optional[Path] = None
  bond: Optional[pd.DataFrame] = None
  angle: Optional[pd.DataFrame] = None
  dihedral: Optional[pd.DataFrame] = None
  chirality: Optional[pd.DataFrame] = None
  cbeta: Optional[pd.DataFrame] = None
  plane: Optional[pd.DataFrame] = None
  nonbonded: Optional[pd.DataFrame] = None
  #manager: Optional[object] = None # GRM



  def __post_init__(self):
    # Rename for consistency
    for f in fields(self):
      value = getattr(self, f.name)
      if isinstance(value,pd.DataFrame) and "restraint_type" in value:
        value["geometry_type"] = value["restraint_type"]
        value.drop(columns=["geometry_type"],inplace=True)

  def select_from_selection(self,df_name,selection,sites,exclusive=False):
    sel_sites = sites.select_from_selection(selection)
    i_seqs = list(sel_sites.index)
    return self.select_from_iseqs(df_name,i_seqs,exclusive=exclusive)

  def select_from_string(self,df_name,sel_str,sites,exclusive=False):
    sel_sites = sites.select_from_phenix_string(sel_str)
    i_seqs = list(sel_sites.index)
    return self.select_from_iseqs(df_name,i_seqs,exclusive=exclusive)

  def _select_from_iseqs(self,df,i_seqs,exclusive=False):
    i_seq_cols = [col for col in df.columns if col.startswith("i_seq_")]
    i_seq_df = df[i_seq_cols]
    i_seqs = np.unique(i_seqs)
    mask = np.isin(i_seq_df,i_seqs)
    if not exclusive:
      # Inclusive, either component i_seqs in selection
      mask = mask.any(axis=1)
    else:
      # Exclusive, both component i_seqs must be in selection
      mask = mask.all(axis=1)
    return df[mask]

  def select_from_iseqs(self,df_name,i_seqs,exclusive=False):
    return self._select_from_iseqs(getattr(self,df_name),i_seqs,exclusive=exclusive)

  @staticmethod
  def add_i_seqs_from_model(self,model):
    # add i_seqs
    d = {}
    for field in fields(self):
      key = field.name
      value = getattr(self, key)
      if isinstance(value,pd.DataFrame):
        d[key] = value
    add_i_seq_columns_from_id_str(d,model)


    # collapse to list
    for key,value in d.items():
      object.__setattr__(self,key,ObjectFrame.collapse_cols(value,"i_seq"))

  @property
  def filename(self):
    if self.filepath is not None:
      return self.filepath.name

  @property
  def dataframes(self):
    d = {
      name:getattr(self,name) for name in self.geometry_names
    }
    return d

  @property
  def geometry_names(self):
    return geometry_names

  # Functions to determine if id_str or i_seq is present. This is hard to
  # determine because how to handle mixed presence? For now, if any are present
  # these functions return True
  @property
  def has_id_str(self):
    ok = []
    for name in self.geometry_names:
      value = getattr(self,name)
      if isinstance(value,pd.DataFrame):
        for col in value.columns:
          if "id_str" in col:
            is_ok = value[col].notna().any()
            ok.append(is_ok)

    if not any(ok) or len(ok)==0:
      return False
    else:
      return True

  @property
  def has_i_seq(self):
    # TODO: a bug in cbeta parsing gives an i_seq column with 'dihedral'
    results = []
    for name in self.geometry_names:
      result = False
      value = getattr(self,name)
      if isinstance(value,pd.DataFrame):
        for col in value.columns:
          if "i_seq" in col:
            result = value[col].notna().any()
      results.append(result)

    if all(results) and len(results)>0:
      return True
    else:
      return False

  @staticmethod
  def guess_and_convert_to_float(df, threshold=1.0):
    """
    Iterates through each column of a DataFrame and converts it to float if appropriate.
    df (pd.DataFrame): The input DataFrame.
    threshold (float): The minimum proportion of valid numeric values needed to convert the column to float.

    """
    for idx, col in enumerate(df.columns):
      try:
        # Attempt to convert the column to numeric values
        converted_col = pd.to_numeric(df[col], errors='coerce')
        # Determine the ratio of numeric values to the total number of values
        non_na_ratio = converted_col.notna().sum() / len(converted_col)
        if non_na_ratio > threshold:  # If ratio of numeric values is above the threshold, convert to float
          df[col] = converted_col
      except:
        pass
    return df

  @classmethod
  def from_geo_file(cls,geo_file_path):
    dfs = parse_geo_file(geo_file_path,return_format="df")
    if dfs is None or all([v is None for k,v in dfs.items()]):
      assert False, "Unable to parse geo file"

    # TODO: Converting here to float guesses, should be defined clearly by field name
    #dfs = {name:cls.guess_and_convert_to_float(df) for name,df in dfs.items()}


    # rename c beta for dataclasses
    if "c-beta" in dfs:
      dfs["cbeta"] = dfs.pop("c-beta")
    kwargs = {k:v for k,v in dfs.items() if k in {field.name for field in fields(cls)}}
    kwargs.update({"filepath":Path(geo_file_path)})
    return cls(**kwargs)

  @classmethod
  def from_mmtbx_model(cls,model,process=True,pdb_interpretation_params=None):
    # Using a .geo file intermediate for geometry to be identical across file/grm for now
    if process:
      model.process(make_restraints=True,pdb_interpretation_params=pdb_interpretation_params)
    grm = model.get_restraints_manager()
    # TODO: replace tempfile with tested direct-from-grm functionality
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
      grm.write_geo_file(sites_cart=model.get_sites_cart(),
                          site_labels = model.get_xray_structure().scatterers().extract_labels(),
                          file_descriptor=temp_file)
      temp_file.flush()  # Ensure all data is written to disk
      obj = cls.from_geo_file(temp_file.name)
      cls.add_i_seqs_from_model(obj,model)
      return obj

  # Add new geometry edits
  def add_plane_geometry(self,geometry_data):
    df = self.plane
    new_row = pd.DataFrame([pd.Series([pd.NA] * len(df.columns), index=df.columns)], index=[0])
    df = pd.concat([new_row, df]).reset_index(drop=True)
    for prefix in ["i_seq","weight"]:
      for i in range(1,len(geometry_data.i_seqs)+1):
        col = f"{prefix}_{i}"
        if col in df.columns:
          value = getattr(geometry_data,prefix+"s")[i-1]
          df.loc[0,col] = value

    if "action" not in df:
      df["action"] = pd.NA
    df.loc[0,"action"] = "add"
    self.plane = df


  def add_bond_geometry(self,geometry_data):
    n_atoms = 2
    df = self.bond
    new_row = pd.DataFrame([pd.Series([pd.NA] * len(df.columns), index=df.columns)], index=[0])
    df = pd.concat([new_row, df]).reset_index(drop=True)
    i_seq_cols = [f"i_seq_{i}" for i in range(1,n_atoms+1)]
    i_seq_values = list(geometry_data.i_seqs)
    cols = [field.name for field in fields(geometry_data) if field.name in df.columns]
    values = [getattr(geometry_data,col) for col in cols]
    cols = i_seq_cols+cols
    values = i_seq_values+values
    for col,value in zip(cols,values):
      # if isinstance(value,list):
      #   value = [value]
      df.at[0,col] = value

    if "action" not in df:
      df["action"] = pd.NA
    df.loc[0,"action"] = "add"
    self.bond = df
  @staticmethod
  def add_iseqs_column(df):
    # add a new column that is a list of the i_seqs
    i_seq_cols = [col for col in df.columns if "i_seq" in col and col != 'i_seqs']
    if len(i_seq_cols)>0:
      df['i_seqs'] = df.apply(lambda row: [row[col] for col in i_seq_cols], axis=1)