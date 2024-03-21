from dataclasses import dataclass, fields
from typing import List, Optional
from io import StringIO
import tempfile

import pandas as pd


from mmtbx.geometry_restraints.geo_file_parsing import (
  parse_geo_file,
  add_i_seq_columns_from_id_str
)
from .base import DataClassBase
from ...last.cctbx_utils import get_restraint_dfs_from_model

@dataclass(frozen=False)
class Restraint(DataClassBase):
  labels: Optional[List[str]] = None
  i_seqs: Optional[List[int]] = None
  id_strs: Optional[List[str]] = None



@dataclass(frozen=False)
class BondRestraint(Restraint):
  ideal: Optional[float] = None
  model: Optional[float] = None
  sigma: Optional[float] = None
  delta: Optional[float] = None
  residual: Optional[float] = None
  weight: Optional[float] = None
  slack: Optional[float] = None
  sym_op_j: Optional[object] = None
  rt_mx: Optional[object] = None

@dataclass(frozen=False)
class AngleRestraint(Restraint):
  ideal: Optional[float] = None
  model: Optional[float] = None
  delta: Optional[float] = None
  residual: Optional[float] = None
  sigma: Optional[float] = None
  weight: Optional[float] = None


@dataclass(frozen=False)
class Restraints(DataClassBase):
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


  # @classmethod
  # def from_geo_file(cls,geo_file_path):
  #   dfs = parse_geo_file(geo_file_path,return_format="df")
  #   instances_dict = {}
  #   for key,restraint_class in cls.components.items():
  #     df = dfs[key]
  #     instances = []
  #     for row in df.itertuples():
  #       d = {
  #       "i_seqs" :[v for k,v in row._asdict().items() if "i_seq" in k],
  #       "id_strs":[v for k,v in row._asdict().items() if "id_str" in k],
  #       }
  #       d.update({k: v for k, v in row._asdict().items() if k in {field.name for field in fields(restraint_class)}})
  #       instances.append(restraint_class(**d))

  #     instances_dict[key] = instances
  #   return cls(filename=None,**instances_dict,manager=None)