from dataclasses import dataclass
from typing import List

from libtbx import group_args

from ..base import DataClassBase
from ..ref import Ref

@dataclass(frozen=True)
class Result(DataClassBase):
  program_name: str
  results: group_args
  params: group_args # or libtbx phil scope extract
  model_refs: List[Ref]
  map_refs: List[Ref]

