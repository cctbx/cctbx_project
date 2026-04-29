from __future__ import absolute_import, division, print_function
from iotbx.pdb.mmcif import cif_input

def tst_bad_symmetry_1():
  """ NMR structure with question mark in unit cell parameters.
  Expected outcome - as if absent completely.
  """
  cif_txt_1 = """\
data_1JZC
_cell.entry_id           1JZC
_cell.length_a           ?
_cell.length_b           ?
_cell.length_c           ?
_cell.angle_alpha        ?
_cell.angle_beta         ?
_cell.angle_gamma        ?
_cell.Z_PDB              1
_cell.pdbx_unique_axis   ?
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1   O "O5'"  . G A 1 1  ? -4.337  -11.083 -6.294 1.00 2.13 ? 1  G A "O5'"  1
ATOM 2   C "C5'"  . G A 1 1  ? -5.629  -11.603 -6.617 1.00 2.29 ? 1  G A "C5'"  1
ATOM 3   C "C4'"  . G A 1 1  ? -6.489  -11.784 -5.370 1.00 2.23 ? 1  G A "C4'"  1
"""
  c_inp = cif_input(source_info=None, lines=cif_txt_1)
  cs = c_inp.crystal_symmetry()
  assert cs.unit_cell() == None
  assert cs.space_group() == None
  assert cs.is_empty()

def tst_bad_symmetry_2():
  """ Same as #1, but no parameters at all.
  Expected outcome - as if absent completely.
  """
  cif_txt_1 = """\
data_1JZC
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1   O "O5'"  . G A 1 1  ? -4.337  -11.083 -6.294 1.00 2.13 ? 1  G A "O5'"  1
ATOM 2   C "C5'"  . G A 1 1  ? -5.629  -11.603 -6.617 1.00 2.29 ? 1  G A "C5'"  1
ATOM 3   C "C4'"  . G A 1 1  ? -6.489  -11.784 -5.370 1.00 2.23 ? 1  G A "C4'"  1
"""
  c_inp = cif_input(source_info=None, lines=cif_txt_1)
  cs = c_inp.crystal_symmetry()
  assert cs.unit_cell() == None
  assert cs.space_group() == None
  assert cs.is_empty()


if __name__ == '__main__':
  tst_bad_symmetry_1()
  tst_bad_symmetry_2()
