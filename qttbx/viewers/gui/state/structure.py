from dataclasses import dataclass
from typing import List, Dict, Optional
from .component import Component
from .base import DataClassBase

@dataclass
class Structure(DataClassBase):
  components: List[Component]

  @classmethod
  def from_default(cls):
    return cls(
      components = {}
    )