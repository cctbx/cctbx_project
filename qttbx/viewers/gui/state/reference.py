from dataclasses import dataclass
from typing import List, Optional
from .structure import Structure
from .base import DataClassBase

@dataclass
class Reference(DataClassBase):
  id_viewer: Optional[str]
  id_molstar: Optional[str]
  structures: List[Structure]

  @classmethod
  def from_default(cls):
    return cls(
      id_viewer = None,
      id_molstar = None,
      structure = Structure.from_default()
    )

  @property
  def representations(self):
    output = []
    for structure in self.structures:
      for component in structure.components:
        for repr_name in component.representations:
          output.append(repr_name)
    return output