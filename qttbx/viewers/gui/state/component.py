from dataclasses import dataclass
from typing import List, Dict
from .base import DataClassBase

@dataclass
class Component(DataClassBase):
  key: str
  representations: List[str]