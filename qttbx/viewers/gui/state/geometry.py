from dataclasses import dataclass, fields
from typing import List, Optional
from io import StringIO
import tempfile
from pathlib import Path

import pandas as pd


from mmtbx.geometry_restraints.geo_file_parsing import (
  parse_geo_file,
  add_i_seq_columns_from_id_str
)
from .base import DataClassBase
from ...core.cctbx_utils import get_restraint_dfs_from_model


@dataclass(frozen=True)
class GeoFileData(DataClassBase):
  filepath: Optional[Path] = None
    
  @property
  def filename(self):
    if self.filepath is not None:
      return self.filepath.name

  


@dataclass(frozen=True)
class Geo(DataClassBase):
  #labels: Optional[List[str]] = None
  i_seqs: Optional[List[int]] = None
  #id_strs: Optional[List[str]] = None


@dataclass(frozen=True)
class BondGeometry(Geo):
  ideal: Optional[float] = None
  sigma: Optional[float] = None
  weight: Optional[float] = None
  slack: Optional[float] = None
  sym_op_j: Optional[object] = None
  rt_mx: Optional[object] = None

@dataclass(frozen=True)
class AngleGeometry(Geo):
  ideal: Optional[float] = None
  sigma: Optional[float] = None
  weight: Optional[float] = None

@dataclass(frozen=True)
class DihedralGeometry(Geo):
  ideal: Optional[float] = None
  sigma: Optional[float] = None
  weight: Optional[float] = None
  harmonic: Optional[float] = None

  @property
  def sinusoidal(self):
    return self.harmonic

@dataclass(frozen=True)
class ChiralGeometry(Geo):
  ideal: Optional[float] = None
  sigma: Optional[float] = None
  weight: Optional[float] = None
  both_signs: Optional[bool] = None

@dataclass(frozen=True)
class PlaneGeometry(Geo):
  default_weight = 2500.0
  weights: Optional[List[float]] = None


@dataclass(frozen=True)
class Geometry(DataClassBase):
  file: Optional[str] = None
  bond: Optional[pd.DataFrame] = None
  angle: Optional[pd.DataFrame] = None
  dihedral: Optional[pd.DataFrame] = None
  chirality: Optional[pd.DataFrame] = None
  cbeta: Optional[pd.DataFrame] = None
  plane: Optional[pd.DataFrame] = None
  nonbonded: Optional[pd.DataFrame] = None
  #manager: Optional[object] = None # GRM


  @property
  def dataframes(self):
    d = {
      name:getattr(self,name) for name in self.restraint_names
    }
    return d

  @property
  def restraint_names(self):
    return [
      "bond",
      "angle",
      "dihedral",
      "chirality",
      "cbeta",
      "plane",
      "nonbonded",
    ]

  # Functions to determine if id_str or i_seq is present. This is hard to
  # determine because how to handle mixed presence? For now, if any are present
  # these functions return True
  @property
  def has_id_str(self):
    ok = []
    for name in self.restraint_names:
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
    for name in self.restraint_names:
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


  @classmethod
  def from_geo_file(cls,geo_file_path):
    dfs = parse_geo_file(geo_file_path,return_format="df")
    if dfs is None or all([v is None for k,v in dfs.items()]):
      assert False, "Unable to parse geo file"
    # rename c beta for dataclasses
    if "c-beta" in dfs:
      dfs["cbeta"] = dfs.pop("c-beta")
    kwargs = {k:v for k,v in dfs.items() if k in {field.name for field in fields(cls)}}
    kwargs.update({"file":geo_file_path})
    return cls(**kwargs)

  @classmethod
  def from_model_ref(cls,model_ref):
    # from a grm
    model = model_ref.model
    model.process(make_restraints=True)
    grm = model.get_restraints_manager()
    # TODO: replace tempfile with tested direct-from-grm functionality
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
      grm.write_geo_file(sites_cart=model.get_sites_cart(),
                          site_labels = model.get_xray_structure().scatterers().extract_labels(),
                          file_descriptor=temp_file)
      temp_file.flush()  # Ensure all data is written to disk
      return cls.from_geo_file(temp_file.name)

  # Add new restraint edits
  def add_plane_restraint(self,restraint_data):
    df = self.plane
    new_row = pd.DataFrame([pd.Series([pd.NA] * len(df.columns), index=df.columns)], index=[0])
    df = pd.concat([new_row, df]).reset_index(drop=True)
    for prefix in ["i_seq","weight"]:
      for i in range(1,len(restraint_data.i_seqs)+1):
        col = f"{prefix}_{i}"
        if col in df.columns:
          value = getattr(restraint_data,prefix+"s")[i-1]
          df.loc[0,col] = value

    if "action" not in df:
      df["action"] = pd.NA
    df.loc[0,"action"] = "add"
    self.plane = df


  def add_bond_restraint(self,restraint_data):
    n_atoms = 2
    df = self.bond
    new_row = pd.DataFrame([pd.Series([pd.NA] * len(df.columns), index=df.columns)], index=[0])
    df = pd.concat([new_row, df]).reset_index(drop=True)
    i_seq_cols = [f"i_seq_{i}" for i in range(1,n_atoms+1)]
    i_seq_values = list(restraint_data.i_seqs)
    cols = [field.name for field in fields(restraint_data) if field.name in df.columns]
    values = [getattr(restraint_data,col) for col in cols]
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