
from dataclasses import dataclass, fields
from typing import List, Optional

import iotbx
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str

from ..state.base import ObjectRow



@dataclass(frozen=True)
class EditData(ObjectRow):
  action: str # "add", ...
  
  
  def as_phil_obj(self):
    if hasattr(self,"name"):
      gm_phil = iotbx.phil.parse(
        input_string     = grand_master_phil_str,
        process_includes = True)
      phil_obj = getattr(gm_phil.extract().geometry_restraints.edits,self.name)[0]
      for key in dir(phil_obj):
        if not key.startswith("_"):
          if hasattr(self,key):
            value = getattr(self,key)
            setattr(phil_obj,key,value)
      return phil_obj

  def add_labels_from_sites(self,sites):
    for f in fields(self):
      value = getattr(self, f.name)
      if "atom_selection_" in f.name:
        idx = f.name.replace("atom_selection_","")
        sites_sel = sites.select_from_phenix_string(value)
        i_seqs= sites_sel.index.values
        label = sites_sel.to_labels_compositional()[0]
        object.__setattr__(self,f"atom_label_{idx}",label)
        #object.__setattr__(self,f"i_seqs",self.i_seqs+i_seqs)

  def to_edits_string(self):
    raise NotImplementedError

@dataclass(frozen=True)
class BondEdit(EditData):
  atom_selection_1: str
  atom_selection_2: str
  distance_ideal: float
  sigma: float
  symmetry_operation: Optional[object] =  None
  slack: Optional[float] = None
  limit: Optional[float] = -1.0
  top_out: Optional[bool] = False
  # New labels
  atom_label_1: Optional[str] = None
  atom_label_2: Optional[str] = None
  #i_seqs: Optional[List[int]] = field(default_factory=lambda: [])



@dataclass(frozen=True)
class AngleEdit(EditData):
  atom_selection_1: str
  atom_selection_2: str
  atom_selection_3: str
  angle_ideal: float
  sigma: float
  symmetry_operation: Optional[object] =  None

  # New labels
  atom_label_1: Optional[str] = None
  atom_label_2: Optional[str] = None
  atom_label_3: Optional[str] = None

    

@dataclass(frozen=True)
class DihedralEdit(EditData):
  atom_selection_1: str
  atom_selection_2: str
  atom_selection_3: str
  atom_selection_4: str

  angle_ideal: float
  sigma: float
  weight:  Optional[float] = None
  harmonic:  Optional[float] = None
  symmetry_operation: Optional[object] =  None

  # New labels
  atom_label_1: Optional[str] = None
  atom_label_2: Optional[str] = None
  atom_label_3: Optional[str] = None
  atom_label_4: Optional[str] = None

      

@dataclass(frozen=True)
class ChiralEdit(EditData):
  atom_selection_1: str
  atom_selection_2: str
  atom_selection_3: str
  ideal: float
  sigma: float


@dataclass(frozen=True)
class PlaneEdit(EditData):
  i_seqs: List[int]


  def to_edits_string(self):
    raise NotImplementedError
