"""
Thing wrappers around external data structures. ie: molecular model or density map
"""
from dataclasses import dataclass
from typing import Optional
from pathlib import Path

from mmtbx.model import manager as ModelManager

from .base import DataClassBase


@dataclass(frozen=True)
class MolecularModelData(DataClassBase):
  filepath: Optional[Path] = None
  model: Optional[ModelManager] = None

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
