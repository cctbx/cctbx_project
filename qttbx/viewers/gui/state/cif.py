
from dataclasses import dataclass
from typing import Optional, Dict
from pathlib import Path
from ...core.cif_io import CifInput

from .base import DataClassBase


@dataclass(frozen=True)
class CifFileData(DataClassBase):
  filepath: Optional[Path] = None

  @property
  def filename(self):
    if self.filepath is not None:
      return self.filepath.name

  @property
  def cif_input(self):
    cif_input = CifInput(self.filepath)
    return cif_input

  @property
  def dataframes(self):
    return self.cif_input
