
from dataclasses import dataclass
from typing import Optional, Dict
from pathlib import Path

import pandas as pd

from .base import DataClassBase
from .cif import CifFileData


@dataclass(frozen=True)
class MolecularModelData(DataClassBase):
  filepath: Optional[Path] = None
  model: Optional[object] = None

  @property
  def filename(self):
     if self.filepath is not None:
      return self.filepath.name



@dataclass(frozen=True)
class RealSpaceMapData(DataClassBase):
  label: Optional[str] = None
  filepath: Optional[str] = None
  map_manager: Optional[object] = None

  @property
  def filename(self):
     return self.filepath.name
