
from dataclasses import dataclass

from .base import DataClassBase

@dataclass(frozen=True)
class Representation(DataClassBase):
  phenixKey: str
  name: str
