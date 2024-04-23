from dataclasses import dataclass
from typing import List, Dict
from .base import DataClassBase

@dataclass(frozen=True)
class Component(DataClassBase):
  key: str
  representations: List[str]