

from dataclasses import dataclass
from pathlib import Path


from ...core.cif_io import CifInput
from .base import DataClassBase


@dataclass(frozen=True)
class CifFileData(DataClassBase):
  filepath: Path

  def __post_init__(self):
    if not isinstance(self.filepath,Path):
      # Workaround for frozen=True
      object.__setattr__(self, 'filepath', Path(self.filepath))

  @property
  def filename(self):
    return self.filepath.name

  @property
  def cif_input(self):
    cif_input = CifInput(self.filepath)
    return cif_input

  @property
  def dataframes(self):
    return self.cif_input
