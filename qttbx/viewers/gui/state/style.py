from dataclasses import dataclass
from typing import Optional, List

from .base import DataClassBase
@dataclass(frozen=True)
class Style(DataClassBase):
  iso: float
  color_theme: str
  opacity: float
  representation: List[str] # one of: ['ball-and-stick', 'cartoon']
  visible: bool
  color: Optional[str] = None
  query: object = None


  @classmethod
  def from_default(cls):
    return cls(
      #query = SelectionQuery.from_all(ref_id=ref_id),
      iso=0.5,
      color=None,
      color_theme='uniform',
      opacity=1.0,
      representation=['ball-and-stick'],
      visible=True
    )

