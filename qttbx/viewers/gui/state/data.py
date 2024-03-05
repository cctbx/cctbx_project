
from dataclasses import dataclass
from typing import Optional, Dict
from pathlib import Path

import pandas as pd

from .base import DataClassBase
from .cif import CifFileData


@dataclass
class MolecularModelData(DataClassBase):
  label: Optional[str] = None
  filepath: Optional[str] = None
  filename: Optional[str] = None
  model: Optional[object] = None
  #cif_data: Optional[CifFileData] = None

  def __post_init__(self):
    super().__post_init__()
    if self.filepath is not None and self.filename is None:
        self.filename = Path(self.filepath).name
    # if self.filename is not None:
    #   if ".cif" in Path(self.filename).suffixes or ".mmcif" in Path(self.filename).suffixes:
    #     cif_data = CifFileData(filepath = self.filepath)
    #     self.cif_data = cif_data

    # if self.filename is None:
    #    self.filename = "model_"+str(id(self))



@dataclass
class RealSpaceMapData(DataClassBase):
  label: Optional[str] = None
  filepath: Optional[str] = None
  filename: Optional[str] = None
  map_manager: Optional[object] = None

  def __post_init__(self):
    super().__post_init__()
    if self.filepath and not self.filename:
        self.filename = Path(self.filepath).name
