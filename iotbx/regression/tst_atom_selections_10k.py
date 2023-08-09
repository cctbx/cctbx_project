from __future__ import absolute_import, division, print_function

import iotbx.pdb

def tst_1():
  h = iotbx.pdb.input(source_info=None, lines="""\
data_1YJP
#
_cell.entry_id           1YJP
_cell.length_a           21.937
_cell.length_b           4.866
_cell.length_c           23.477
_cell.angle_alpha        90.00
_cell.angle_beta         107.08
_cell.angle_gamma        90.00
_cell.Z_PDB              2
_cell.pdbx_unique_axis   ?
#
_symmetry.entry_id                         1YJP
_symmetry.space_group_name_H-M             'P 1 21 1'
_symmetry.pdbx_full_space_group_name_H-M   ?
_symmetry.cell_setting                     ?
_symmetry.Int_Tables_number                4
_symmetry.space_group_name_Hall            ?
#
_atom_sites.entry_id                    1YJP
_atom_sites.fract_transf_matrix[1][1]   0.045585
_atom_sites.fract_transf_matrix[1][2]   0.000000
_atom_sites.fract_transf_matrix[1][3]   0.014006
_atom_sites.fract_transf_matrix[2][1]   0.000000
_atom_sites.fract_transf_matrix[2][2]   0.205508
_atom_sites.fract_transf_matrix[2][3]   0.000000
_atom_sites.fract_transf_matrix[3][1]   0.000000
_atom_sites.fract_transf_matrix[3][2]   0.000000
_atom_sites.fract_transf_matrix[3][3]   0.044560
_atom_sites.fract_transf_vector[1]      0.00000
_atom_sites.fract_transf_vector[2]      0.00000
_atom_sites.fract_transf_vector[3]      0.00000
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
ATOM   1  N N   . GLY A 1 1 ? -9.009  4.612  6.102  1.00 16.77 ? 1  GLY A N   1
ATOM   2  C CA  . GLY A 1 1 ? -9.052  4.207  4.651  1.00 16.57 ? 1  GLY A CA  1
ATOM   3  C C   . GLY A 1 1 ? -8.015  3.140  4.419  1.00 16.16 ? 1  GLY A C   1
ATOM   4  O O   . GLY A 1 1 ? -7.523  2.521  5.381  1.00 16.78 ? 1  GLY A O   1
HETATM 60 O O   . HOH B 2 . ? -6.471  5.227  7.124  1.00 22.62 ? 8  HOH A O   1
HETATM 61 O O   . HOH B 2 . ? 10.431  1.858  3.216  1.00 19.71 ? 9  HOH A O   1
HETATM 62 O O   . HOH B 2 . ? -11.286 1.756  -1.468 1.00 17.08 ? 10 HOH A O   1
HETATM 63 O O   . HOH B 2 . ? 11.808  4.179  9.970  1.00 23.99 ? 10011 HOH A O   1
HETATM 64 O O   . HOH B 2 . ? 13.605  1.327  9.198  1.00 26.17 ? 10012 HOH A O   1
HETATM 65 O O   . HOH B 2 . ? -2.749  3.429  10.024 1.00 39.15 ? 10013 HOH A O   1
HETATM 66 O O   . HOH B 2 . ? -1.500  0.682  10.967 1.00 43.49 ? 10014 HOH A O   1
#    """).construct_hierarchy()

  # for a in h.atoms():
  #   print (a.id_str(), a.parent().parent().resseq,  a.parent().parent().resseq_as_int())

  asc = h.atom_selection_cache()
  # single residues:
  #
  assert list(asc.selection("resseq 10014").iselection()) == [10]
  assert list(asc.selection("resseq A00E").iselection()) == [10]
  # note resid is different:
  assert list(asc.selection("resid 10014").iselection()) == []
  assert list(asc.selection("resid A00E").iselection()) == [10]

  # Range of residues:
  #
  assert list(asc.selection("resseq 10011:10014").iselection()) == [7,8,9,10]
  assert list(asc.selection("resseq 10:10014").iselection()) == [6,7,8,9,10]
  assert list(asc.selection("resseq 10:10012").iselection()) == [6,7,8]
  # With hy36
  assert list(asc.selection("resseq 10011:A00E").iselection()) == [7,8,9,10]
  assert list(asc.selection("resseq 10:A00E").iselection()) == [6,7,8,9,10]

if (__name__ == "__main__"):
  tst_1()
