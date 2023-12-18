
from dataclasses import dataclass
from typing import Optional, Dict
from pathlib import Path

from ...last.cif_io import read_cif_file
from .base import DataClassBase


@dataclass
class CifFileData(DataClassBase):
  filepath: Optional[str] = None
  filename: Optional[str] = None

  def __post_init__(self):
    super().__post_init__()
    if self.filepath and not self.filename:
        self.filename = Path(self.filepath).name

  @property
  def dataframes(self):
    if not hasattr(self,"_pandas"):
      self._pandas = read_cif_file(self.filepath,method="iotbx",return_as="pandas")
      #self._pandas = read_cif_file(self.filepath,method="cifpd",return_as="pandas")
    return self._pandas
