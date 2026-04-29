from __future__ import absolute_import, division, print_function
import iotbx.pdb

cif_str = """\
data_1YJP
_cell.entry_id           1YJP
_cell.length_a           21.937
_cell.length_b           4.866
_cell.length_c           23.477
_cell.angle_alpha        90.00
_cell.angle_beta         107.08
_cell.angle_gamma        90.00
_cell.Z_PDB              2
_cell.pdbx_unique_axis   ?
_pdbx_database_status.recvd_initial_deposition_date   2005-01-15
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
ATOM   1  N N   . GLY A 1 1 ? -9.009  4.612  6.102  1.00 16.77 ? 1  GLY A N   1
ATOM   2  C CA  . GLY A 1 1 ? -9.052  4.207  4.651  1.00 16.57 ? 1  GLY A CA  1
ATOM   3  C C   . GLY C 1 1 ? -8.015  3.140  4.419  1.00 16.16 ? 1  GLY CC C   1
ATOM   4  O O   . GLY D 1 1 ? -7.523  2.521  5.381  1.00 16.78 ? 1  GLY DDDD O   1
"""


def exercise_extract_header_year(prefix="iotbx_tst_mmcif_segids"):
  cif_inp = iotbx.pdb.input(lines=cif_str, source_info=None)
  # print(cif_inp.deposition_date())
  # print(cif_inp.extract_header_year())
  assert cif_inp.deposition_date() == "2005-01-15"
  assert cif_inp.extract_header_year() == 2005

  # now the same but with ?

  cif_str_2 = cif_str.replace('2005-01-15', '?')
  cif_inp_2 = iotbx.pdb.input(lines=cif_str_2, source_info=None)
  # print(cif_inp_2.deposition_date())
  # print(cif_inp_2.extract_header_year())
  assert cif_inp_2.deposition_date() == '?'
  assert cif_inp_2.extract_header_year() == None

def exercise_label_to_auth_asym_id_dictionary():
  cif_inp = iotbx.pdb.input(lines=cif_str, source_info=None)
  res = cif_inp.label_to_auth_asym_id_dictionary()
  # print(res)
  assert res == {'A': 'A', 'C': 'CC', 'D': 'DDDD'}, res

if __name__ == "__main__":
  exercise_extract_header_year()
  exercise_label_to_auth_asym_id_dictionary()
