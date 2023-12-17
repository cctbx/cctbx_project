from dataclasses import dataclass
from typing import Optional, List

from .base import DataClassBase
@dataclass
class Style(DataClassBase):
  ref_id: str
  iso: float
  color_theme: str
  opacity: float
  representation: List[str] # one of: ['ball-and-stick', 'cartoon']
  visible: bool
  color: Optional[str] = None
  query: object = None


  @classmethod
  def from_default(cls,ref_id):
    return cls(
      ref_id = ref_id,
      #query = SelectionQuery.from_all(ref_id=ref_id),
      iso=0.5,
      color=None,
      color_theme='uniform',
      opacity=1.0,
      representation=['ball-and-stick'],
      visible=True
    )

