from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out

# ------------------------------------------------------------------------------

def run():
  """
  Test that MTRIX matrix is not applied if it is referred to as
  "given" and it is not under loop_.
  """
  inp = iotbx.pdb.input(lines=mmcif_str.split("\n"), source_info=None)
  #inp = iotbx.pdb.input('5aj2_test.cif')
  m = mmtbx.model.manager(model_input = inp,log = null_out())
  p = m.get_hierarchy()
  #p.overall_counts().show()
  chains = list(p.overall_counts().chain_ids.keys())
  chains.sort()
  answer = ['A', 'B', 'C']
  assert (chains == answer), '%s %s' % (chains, answer)

mmcif_str = '''
data_5AJ2
#
_cell.entry_id           5AJ2
_cell.length_a           1.000
_cell.length_b           1.000
_cell.length_c           1.000
_cell.angle_alpha        90.00
_cell.angle_beta         90.00
_cell.angle_gamma        90.00
_cell.Z_PDB              1
_cell.pdbx_unique_axis   ?
#
_symmetry.entry_id                         5AJ2
_symmetry.space_group_name_H-M             'P 1'
_symmetry.pdbx_full_space_group_name_H-M   ?
_symmetry.cell_setting                     ?
_symmetry.Int_Tables_number                1
#
_struct_ncs_oper.id             1
_struct_ncs_oper.code           given
_struct_ncs_oper.details        ?
_struct_ncs_oper.matrix[1][1]   0.858065
_struct_ncs_oper.matrix[1][2]   0.000000
_struct_ncs_oper.matrix[1][3]   0.513541
_struct_ncs_oper.matrix[2][1]   0.000000
_struct_ncs_oper.matrix[2][2]   1.000000
_struct_ncs_oper.matrix[2][3]   0.000000
_struct_ncs_oper.matrix[3][1]   -0.513541
_struct_ncs_oper.matrix[3][2]   0.000000
_struct_ncs_oper.matrix[3][3]   0.858065
_struct_ncs_oper.vector[1]      -131.70160
_struct_ncs_oper.vector[2]      5.57000
_struct_ncs_oper.vector[3]      173.48700
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
ATOM 1   C CA . VAL A 1 93  ? 265.910 287.538 371.741 1.00 65.63  ? 93   VAL A CA 1
ATOM 2   C CA . THR A 1 94  ? 263.468 290.283 372.749 1.00 75.26  ? 94   THR A CA 1
ATOM 3   C CA . GLU A 1 95  ? 263.384 291.733 376.277 1.00 36.70  ? 95   GLU A CA 1
ATOM 4   C CA . GLU A 1 96  ? 259.609 291.836 376.486 1.00 59.29  ? 96   GLU A CA 1
ATOM 5   C CA . ASP A 1 97  ? 259.526 288.058 375.892 1.00 55.35  ? 97   ASP A CA 1
ATOM 6   C CA . LEU A 1 98  ? 261.809 287.789 378.893 1.00 19.26  ? 98   LEU A CA 1
ATOM 7   C CA . ASN A 1 99  ? 259.447 289.832 381.026 1.00 24.16  ? 99   ASN A CA 1
ATOM 8   C CA . VAL A 1 100 ? 256.439 287.714 380.038 1.00 10.60  ? 100  VAL A CA 1
ATOM 264 C CA . ALA B 2 1   ? 253.610 250.832 414.360 1.00 26.29  ? 356  ALA B CA 1
ATOM 265 C CA . HIS B 2 2   ? 257.043 250.568 415.989 1.00 21.06  ? 357  HIS B CA 1
ATOM 266 C CA . THR B 2 3   ? 259.058 250.909 412.773 1.00 5.95   ? 358  THR B CA 1
ATOM 267 C CA . GLN B 2 4   ? 258.465 252.946 409.610 1.00 4.25   ? 359  GLN B CA 1
ATOM 268 C CA . THR B 2 5   ? 258.221 249.640 407.755 1.00 5.45   ? 360  THR B CA 1
ATOM 269 C CA . MET B 2 6   ? 255.270 248.643 409.910 1.00 34.55  ? 361  MET B CA 1
ATOM 270 C CA . LEU B 2 7   ? 253.772 252.095 409.294 1.00 5.13   ? 362  LEU B CA 1
ATOM 271 C CA . PHE B 2 8   ? 254.068 251.554 405.553 1.00 6.06   ? 363  PHE B CA 1
ATOM 272 C CA . GLN B 2 9   ? 252.690 248.038 406.042 1.00 19.38  ? 364  GLN B CA 1
ATOM 273 C CA . THR B 2 10  ? 249.682 249.259 408.001 1.00 6.76   ? 365  THR B CA 1
ATOM 489 C CA . GLN C 3 1   ? 291.310 230.862 413.727 1.00 13.26  ? 580  GLN C CA 1
ATOM 490 C CA . GLY C 3 2   ? 294.400 231.823 415.719 1.00 4.03   ? 581  GLY C CA 1
ATOM 491 C CA . LYS C 3 3   ? 294.558 235.261 414.129 1.00 4.01   ? 582  LYS C CA 1
ATOM 492 C CA . SER C 3 4   ? 294.108 238.876 415.238 1.00 10.15  ? 583  SER C CA 1
ATOM 493 C CA . LEU C 3 5   ? 291.813 241.792 414.393 1.00 6.63   ? 584  LEU C CA 1
ATOM 494 C CA . TYR C 3 6   ? 292.428 245.561 414.505 1.00 4.96   ? 585  TYR C CA 1
ATOM 495 C CA . ILE C 3 7   ? 289.616 248.010 415.243 1.00 4.50   ? 586  ILE C CA 1
ATOM 496 C CA . ASN C 3 8   ? 289.741 251.810 415.279 1.00 6.27   ? 587  ASN C CA 1
ATOM 497 C CA . SER C 3 9   ? 286.657 253.085 417.096 1.00 26.29  ? 588  SER C CA 1
ATOM 498 C CA . GLU C 3 10  ? 286.522 256.305 415.074 1.00 19.32  ? 589  GLU C CA 1
#
'''

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
