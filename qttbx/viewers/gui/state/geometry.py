from dataclasses import dataclass, fields
from typing import List, Optional
import tempfile
from pathlib import Path

import pandas as pd


from mmtbx.geometry_restraints.geo_file_parsing import (
  parse_geo_file
)
from .base import DataClassBase

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
  ideal: Optional[float] = None
  sigma: Optional[float] = None
  weight: Optional[float] = None
  slack: Optional[float] = None
  sym_op_j: Optional[object] = None
  rt_mx: Optional[object] = None

@dataclass(frozen=True)
class AngleGeometry(Geo):
  name = "angle"
  ideal: Optional[float] = None
  sigma: Optional[float] = None
  weight: Optional[float] = None

@dataclass(frozen=True)
class DihedralGeometry(Geo):
  name = "dihedral"
  ideal: Optional[float] = None
  sigma: Optional[float] = None
  weight: Optional[float] = None
  harmonic: Optional[float] = None

  @property
  def sinusoidal(self):
    return self.harmonic

@dataclass(frozen=True)
class ChiralGeometry(Geo):
  name = "chirality"
  ideal: Optional[float] = None
  sigma: Optional[float] = None
  weight: Optional[float] = None
  both_signs: Optional[bool] = None

@dataclass(frozen=True)
class PlaneGeometry(Geo):
  name = "plane"
  default_weight = 2500.0
  weights: Optional[List[float]] = None


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


  @classmethod
  def from_geo_file(cls,geo_file_path):
    dfs = parse_geo_file(geo_file_path,return_format="df")
    if dfs is None or all([v is None for k,v in dfs.items()]):
      assert False, "Unable to parse geo file"
    # rename c beta for dataclasses
    if "c-beta" in dfs:
      dfs["cbeta"] = dfs.pop("c-beta")
    kwargs = {k:v for k,v in dfs.items() if k in {field.name for field in fields(cls)}}
    kwargs.update({"filepath":Path(geo_file_path)})
    return cls(**kwargs)

  @classmethod
  def from_mmtbx_model(cls,model,process=True):
    # from a grm
    if process:
      model.process(make_restraints=True)
    grm = model.get_restraints_manager()
    # TODO: replace tempfile with tested direct-from-grm functionality
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
      grm.write_geo_file(sites_cart=model.get_sites_cart(),
                          site_labels = model.get_xray_structure().scatterers().extract_labels(),
                          file_descriptor=temp_file)
      temp_file.flush()  # Ensure all data is written to disk
      return cls.from_geo_file(temp_file.name)

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