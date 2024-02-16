from typing import Optional

from ..ref import Ref, ModelRef, SelectionRef
from .results import Result

class ResultsRef(Ref):
  def __init__(self,data: Result):
    super().__init__(data=data)

    for model_ref in self.data.model_refs:
      model_ref.results[self.data.program_name] = self
    for map_ref in self.data.map_refs:
      map_ref.results[self.data.program_name] = self

# class ModelResultsRef(ResultsRef):
#   def __init__(self,data: Result, model_ref: ModelRef, selection_ref: Optional[SelectionRef] = None):
#     super().__init__(data=data)
#     self.model_ref = model_ref
#     self.selection_ref = selection_ref
