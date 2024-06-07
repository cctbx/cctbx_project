
from dataclasses import dataclass
from typing import List

from .base import DataClassBase

@dataclass(frozen=True)
class Component(DataClassBase):
  phenixKey: str
  key: str
  representations: List[str]