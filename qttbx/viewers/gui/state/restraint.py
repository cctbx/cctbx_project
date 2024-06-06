from dataclasses import dataclass, fields
from pathlib import Path
from typing import List, Optional

from elbow.command_line.where_is_that_cif_file import run as where_is_that_cif_file

from .base import DataClassBase
from .cif import CifFileData


@dataclass(frozen=True)
class Restraint(CifFileData):
  """
  A single chemical component and restraint file
  """
  comp_id: str
  # Also inherits filepath: str


@dataclass(frozen=True)
class Restraints(DataClassBase):
  """
  A collection of restraints
  """
  restraints: List[Restraint]

  @classmethod
  def from_sites(cls,sites):
    comp_ids = sites["comp_id"].unique()
    restraints = []
    for comp_id in comp_ids:
      file = where_is_that_cif_file(comp_id.upper(),file_name_only=True,verbose=0)
      restraint = Restraint(comp_id=comp_id,filepath=file)
      restraints.append(restraint)
    return cls(restraints=restraints)
    
