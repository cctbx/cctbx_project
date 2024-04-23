
from dataclasses import dataclass
from typing import Optional, Dict, List
from pathlib import Path

from ..state.base import ObjectRow

@dataclass(frozen=True)
class EditData(ObjectRow):
  pass


@dataclass(frozen=True)
class BondEdit(EditData):
  i_seqs: List[int]
  ideal_old: float
  ideal_new: float
  sigma_old: float
  sigma_new: float


