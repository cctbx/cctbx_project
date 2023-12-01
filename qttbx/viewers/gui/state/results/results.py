from dataclasses import dataclass

from ..base import DataClassBase

@dataclass
class Result(DataClassBase):
  program_name: str

