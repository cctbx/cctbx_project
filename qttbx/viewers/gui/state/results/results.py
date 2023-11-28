from dataclasses import dataclass

from ..base import DataClassBase

@dataclass(frozen=True)
class Result(DataClassBase):
  program_name: str

