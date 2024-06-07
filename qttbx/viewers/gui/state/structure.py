from dataclasses import dataclass
from typing import List

from .component import Component
from .base import DataClassBase

@dataclass(frozen=True)
class Structure(DataClassBase):
  phenixReferenceKey: str
  phenixKey: str
  data_id: str
  key: str
  components: List[Component]
