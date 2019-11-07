from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from libtbx.test_utils import show_diff


test_cif_str = """\
data_default
_cell.angle_beta                  107.080
_cell.angle_gamma                 90.000
_cell.length_b                    4.866
_cell.length_c                    23.477
_cell.angle_alpha                 90.000
_cell.volume                      2395.534
_cell.length_a                    21.937
_space_group.crystal_system       monoclinic
_space_group.name_H-M_alt         'P 1 21 1'
_space_group.IT_number            4
_space_group.name_Hall            ' P 2yb'
_symmetry.space_group_name_H-M    'P 1 21 1'
_symmetry.Int_Tables_number       4
_symmetry.space_group_name_Hall   ' P 2yb'

_atom_site.group_PDB ATOM
_atom_site.id 1
_atom_site.label_atom_id N
_atom_site.label_alt_id A
_atom_site.label_comp_id GLY
_atom_site.auth_asym_id A
_atom_site.auth_seq_id 1
_atom_site.pdbx_PDB_ins_code I
_atom_site.Cartn_x   -9.13600
_atom_site.Cartn_y  4.55200
_atom_site.Cartn_z  5.82700
_atom_site.occupancy  1.000
_atom_site.B_iso_or_equiv  18.49000
_atom_site.type_symbol  N
_atom_site.pdbx_formal_charge  ?
_atom_site.label_asym_id  A
_atom_site.label_entity_id  ?
_atom_site.label_seq_id  1
_atom_site.pdbx_PDB_model_num  1

_atom_site_anisotrop.id  1
_atom_site_anisotrop.pdbx_auth_atom_id  N
_atom_site_anisotrop.pdbx_label_alt_id  A
_atom_site_anisotrop.pdbx_auth_comp_id  GLY
_atom_site_anisotrop.pdbx_auth_asym_id  A
_atom_site_anisotrop.pdbx_auth_seq_id  1
_atom_site_anisotrop.pdbx_PDB_ins_code  I
_atom_site_anisotrop.U[1][1]  0.29030
_atom_site_anisotrop.U[2][2]  0.13360
_atom_site_anisotrop.U[3][3]  0.27880
_atom_site_anisotrop.U[1][2]  0.06870
_atom_site_anisotrop.U[1][3]  0.10560
_atom_site_anisotrop.U[2][3]  0.05900
"""

pdb_answer = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N  AGLY A   1I     -9.136   4.552   5.827  1.00 18.49           N
ANISOU    1  N  AGLY A   1I    2903   1336   2788    687   1056    590       N
TER
END
"""

def exercise(prefix="tst_cif_as_pdb_1atom"):
  """ Note that cif contains only 1 atom and information is not wrapped into
  loops (loop_) as usual
  """
  cif_inp = iotbx.pdb.input(lines=test_cif_str, source_info=None)
  model = mmtbx.model.manager(model_input=cif_inp)
  pdb_txt = model.model_as_pdb()
  assert not show_diff(pdb_txt, pdb_answer)

if __name__ == '__main__':
  exercise()
