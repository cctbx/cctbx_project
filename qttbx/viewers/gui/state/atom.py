from typing import Type, Optional
from dataclasses import dataclass

from .base import DataClassBase


@dataclass(frozen=True)
class Atom(DataClassBase):
  name = 'atom'

  # Required fields
  id: str
  i_seq: int
  label_comp_id: str
  label_atom_id: str
  type_symbol: str   
  group_PDB: str
  label_alt_id: str

  # At least one asym_id
  auth_asym_id: Optional[str] = None
  label_asym_id: Optional[str] = None

  # At least one seq_id
  auth_seq_id: Optional[int] = None
  label_seq_id: Optional[int] = None


  # Not required fields
  pdbx_PDB_model_num: Optional[str] = None
  pdbx_PDB_ins_code: Optional[str] = None
  Cartn_x: Optional[float] = None
  Cartn_y: Optional[float] = None
  Cartn_z: Optional[float] = None
  B_iso_or_equiv: Optional[float] = None
  occupancy: Optional[float] = None
  pdbx_formal_charge: Optional[int] = None
  label_entity_id: Optional[str] = None

  def __post_init__(self):
    assert self.auth_asym_id or self.label_asym_id
    assert self.auth_seq_id or self.label_seq_id