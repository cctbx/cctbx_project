
from dataclasses import dataclass
from typing import List, Optional

from ..state.base import ObjectRow

@dataclass(frozen=True)
class EditData(ObjectRow):
  action: str # "add", "mod"
  i_seqs: List[int]
  labels_compositional: List[str]
  sel_strings: List[str]  
    
  def to_edits_string(self):
    raise NotImplementedError

@dataclass(frozen=True)
class BondEdit(EditData):
  ideal_old: float
  ideal_new: float
  sigma_new: float = 0.001
  sigma_old: float = None


  def to_edits_string(self):
    s = f"""
bond {{
  action: *{self.action}
  atom_selection_1 = {self.sel_strings[0]}
  atom_selection_2 = {self.sel_strings[1]}
  symmetry_operation = {None}
  distance_ideal = {self.ideal_new}
  sigma = {self.sigma_new}
  slack = {None}
  }}
    """
    return s

@dataclass(frozen=True)
class AngleEdit(EditData):
  i_seqs: List[int]
  ideal_old: float
  ideal_new: float
  sigma_old: float
  sigma_new: float = 3.0
  
  def to_edits_string(self):
    s = f"""
angle {{
  action: *{self.action}
  atom_selection_1 = {self.sel_strings[0]}
  atom_selection_2 = {self.sel_strings[1]}
  atom_selection_3 = {self.sel_strings[2]}
  symmetry_operation = {None}
  distance_ideal = {self.ideal_new}
  sigma = {self.ideal_new}
  slack = {None}
  }}
    """
    return s

@dataclass(frozen=True)
class DihedralEdit(EditData):
  i_seqs: List[int]
  ideal_old: float
  ideal_new: float
  sigma_old: float
  sigma_new: float
  weight_old: Optional[float] = None
  weight_new: Optional[float] = None
  harmonic_old: Optional[float] = None
  harmonic_new: Optional[float] = None

  def to_edits_string(self):
    raise NotImplementedError

