from dataclasses import dataclass
from typing import List, Dict
from .base import DataClassBase

@dataclass(frozen=True)
class Representation(DataClassBase):
  phenixKey: str
  name: str
