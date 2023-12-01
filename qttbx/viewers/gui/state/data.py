
from dataclasses import dataclass
from typing import Optional
from pathlib import Path

from .base import DataClassBase


@dataclass
class MolecularModelData(DataClassBase):
  filepath: Optional[str] = None
  filename: Optional[str] = None
  model: Optional[object] = None
  
  def __post_init__(self):
    super().__post_init__()
    if self.filepath and not self.filename:
        self.filename = Path(self.filepath).name

   

@dataclass
class RealSpaceMapData(DataClassBase):
  filepath: Optional[str] = None
  filename: Optional[str] = None
  map_manager: Optional[object] = None

  def __post_init__(self):
    super().__post_init__()
    if self.filepath and not self.filename:
        self.filename = Path(self.filepath).name
    
  