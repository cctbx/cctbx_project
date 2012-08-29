from __future__ import division
from cStringIO import StringIO

from libtbx.test_utils import approx_equal
from libtbx.test_utils import Exception_expected
from libtbx.test_utils import show_diff

import iotbx.cif
from iotbx.pdb.mmcif import pdb_hierarchy_builder


def exercise_pdb_hierachy_builder():
  input_1ab1 = """\
data_1AB1
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
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   17  C CA   A THR A 1 2 13.800 11.354 5.944  0.70 2.98  2  THR A CA   1
ATOM   18  C CA   B THR A 1 2 13.872 10.918 5.772  0.30 5.01  2  THR A CA   1
ATOM   19  C C    . THR A 1 2 14.118 10.664 7.274  1.00 2.92  2  THR A C    1
ATOM   20  O O    . THR A 1 2 15.013 9.787  7.376  1.00 3.70  2  THR A O    1
ATOM   21  C CB   A THR A 1 2 12.784 10.610 5.099  0.70 4.35  2  THR A CB   1
ATOM   22  C CB   B THR A 1 2 12.938 9.905  5.108  0.30 1.07  2  THR A CB   1
ATOM   23  O OG1  A THR A 1 2 13.364 9.360  4.756  0.70 5.60  2  THR A OG1  1
ATOM   24  O OG1  B THR A 1 2 12.559 10.450 3.917  0.30 5.43  2  THR A OG1  1
ATOM   25  C CG2  A THR A 1 2 12.451 11.323 3.797  0.70 7.65  2  THR A CG2  1
ATOM   26  C CG2  B THR A 1 2 11.642 9.729  5.862  0.30 2.71  2  THR A CG2  1
ATOM   27  H H    . THR A 1 2 15.500 10.728 5.006  1.00 1.74  2  THR A H    1
ATOM   28  H HA   . THR A 1 2 13.437 12.241 6.128  1.00 4.09  2  THR A HA   1
ATOM   29  H HB   A THR A 1 2 11.947 10.457 5.635  0.70 1.11  2  THR A HB   1
ATOM   30  H HB   B THR A 1 2 13.426 9.030  5.019  0.30 4.20  2  THR A HB   1
ATOM   31  H HG21 A THR A 1 2 13.241 11.856 3.510  0.70 1.00  2  THR A HG21 1
ATOM   32  H HG21 B THR A 1 2 11.196 10.615 5.970  0.30 2.01  2  THR A HG21 1
ATOM   33  H HG22 A THR A 1 2 12.211 10.692 3.092  0.70 5.54  2  THR A HG22 1
ATOM   34  H HG22 B THR A 1 2 11.027 9.120  5.405  0.30 10.62 2  THR A HG22 1
ATOM   35  H HG23 A THR A 1 2 11.693 11.979 3.928  0.70 9.26  2  THR A HG23 1
ATOM   36  H HG23 B THR A 1 2 11.801 9.381  6.801  0.30 8.15  2  THR A HG23 1
# ...trunacted...
HETATM 679 C C1   A EOH B 2 . 16.083 0.909  12.713 0.60 15.95 66 EOH A C1   1
HETATM 680 C C1   B EOH B 2 . 15.565 1.042  13.226 0.40 14.21 66 EOH A C1   1
HETATM 681 C C2   A EOH B 2 . 15.130 -0.269 13.016 0.60 14.85 66 EOH A C2   1
HETATM 682 C C2   B EOH B 2 . 14.830 -0.303 13.058 0.40 15.17 66 EOH A C2   1
HETATM 683 O O    A EOH B 2 . 15.378 1.984  12.104 0.60 12.07 66 EOH A O    1
HETATM 684 O O    B EOH B 2 . 14.811 2.078  12.602 0.40 5.53  66 EOH A O    1
"""
  cif_model = iotbx.cif.reader(input_string=input_1ab1).model()
  builder = pdb_hierarchy_builder(cif_model["1AB1"])
  #assert builder.crystal_symmetry is None
  hierarchy = builder.hierarchy
  s = StringIO()
  hierarchy.show(out=s)
  assert not show_diff(s.getvalue(), """\
model id="" #chains=2
  chain id="A" #residue_groups=1  ### WARNING: duplicate chain id ###
    resid="   2 " #atom_groups=3
      altloc="" resname="THR" #atoms=4
        " C  "
        " O  "
        " H  "
        " HA "
      altloc="A" resname="THR" #atoms=8
        " CA "
        " CB "
        " OG1"
        " CG2"
        " HB "
        "HG21"
        "HG22"
        "HG23"
      altloc="B" resname="THR" #atoms=8
        " CA "
        " CB "
        " OG1"
        " CG2"
        " HB "
        "HG21"
        "HG22"
        "HG23"
  chain id="A" #residue_groups=1  ### WARNING: duplicate chain id ###
    resid="  66 " #atom_groups=2
      altloc="A" resname="EOH" #atoms=3
        " C1 "
        " C2 "
        " O  "
      altloc="B" resname="EOH" #atoms=3
        " C1 "
        " C2 "
        " O  "
""")

  input_missing_mandatory_items = """\
data_1AB1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM   17  C A THR A 1 2 0.70 2.98
ATOM   18  C B THR A 1 2 0.30 5.01
ATOM   19  C . THR A 1 2 1.00 2.92
ATOM   20  O . THR A 1 2 1.00 3.70
"""
  cif_model = iotbx.cif.reader(input_string=input_missing_mandatory_items).model()
  try: pdb_hierarchy_builder(cif_model["1AB1"])
  except AssertionError, e: pass
  else:
    # TODO: raise a better error here
    raise Exception_expected

  input_4edr = """\
data_4EDR
_cell.length_a           150.582
_cell.length_b           150.582
_cell.length_c           38.633
_cell.angle_alpha        90.000
_cell.angle_beta         90.000
_cell.angle_gamma        120.000
#
_symmetry.space_group_name_H-M             'P 61'
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
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1    N  N     . SER A 1 1 21.138  -69.073 17.360  1.00 23.68 108 SER A N     1
ATOM   2    C  CA    . SER A 1 1 22.164  -68.793 18.358  1.00 22.98 108 SER A CA    1
ATOM   3    C  C     . SER A 1 1 23.173  -67.799 17.805  1.00 21.13 108 SER A C     1
ATOM   4    O  O     . SER A 1 1 23.251  -67.594 16.595  1.00 19.34 108 SER A O     1
ATOM   5    C  CB    . SER A 1 1 22.882  -70.080 18.766  1.00 22.68 108 SER A CB    1
ATOM   6    O  OG    . SER A 1 1 23.683  -70.569 17.703  1.00 24.00 108 SER A OG    1
HETATM 2650 MN MN    . MN  F 4 . 9.296   -44.783 -6.320  1.00 44.18 505 MN  A MN    1
#
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.type_symbol
_atom_site_anisotrop.pdbx_label_atom_id
_atom_site_anisotrop.pdbx_label_alt_id
_atom_site_anisotrop.pdbx_label_comp_id
_atom_site_anisotrop.pdbx_label_asym_id
_atom_site_anisotrop.pdbx_label_seq_id
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
1    N N   . SER A 1   0.4097 0.2916 0.1984 0.1130  0.0328  0.0375
2    C CA  . SER A 1   0.4035 0.2848 0.1847 0.1242  0.0297  0.0347
3    C C   . SER A 1   0.3744 0.2637 0.1648 0.1278  0.0215  0.0258
4    O O   . SER A 1   0.3494 0.2393 0.1463 0.1236  0.0185  0.0236
5    C CB  . SER A 1   0.4085 0.2738 0.1795 0.1267  0.0311  0.0373
6    O OG  . SER A 1   0.4276 0.2843 0.1998 0.1252  0.0273  0.0344
"""

  cif_model = iotbx.cif.reader(input_string=input_4edr).model()
  cif_block = cif_model["4EDR"]
  builder = pdb_hierarchy_builder(cif_block)
  assert approx_equal(builder.crystal_symmetry.unit_cell().parameters(),
                      (150.582, 150.582, 38.633, 90, 90, 120))
  assert builder.crystal_symmetry.space_group_info().symbol_and_number() == \
         'P 61 (No. 169)'
  hierarchy = builder.hierarchy
  residue_group = hierarchy.residue_groups().next()
  assert residue_group.resseq == ' 108'
  assert residue_group.resseq_as_int() == 108
  atoms = hierarchy.atoms()
  assert atoms[0].serial == '1'
  assert atoms[0].name == ' N  '
  assert atoms[0].b == 23.68
  assert atoms[0].uij_is_defined()
  assert approx_equal(
    atoms[0].uij, (0.4097, 0.2916, 0.1984, 0.113, 0.0328, 0.0375))
  assert atoms[1].serial == '2'
  assert atoms[1].name == ' CA '
  assert atoms[1].b == 22.98
  assert atoms[1].uij_is_defined()
  assert approx_equal(
    atoms[1].uij, (0.4035, 0.2848, 0.1847, 0.1242, 0.0297, 0.0347))
  assert not atoms[6].uij_is_defined()
  assert atoms[6].serial == '2650'
  assert atoms[6].name == 'MN  '
  assert atoms[6].b == 44.18
  assert approx_equal(atoms[6].uij, (-1, -1, -1, -1, -1, -1))
  s = StringIO()
  hierarchy.show(out=s)
  assert not show_diff(s.getvalue(), """\
model id="" #chains=2
  chain id="A" #residue_groups=1  ### WARNING: duplicate chain id ###
    resid=" 108 " #atom_groups=1
      altloc="" resname="SER" #atoms=6
        " N  "
        " CA "
        " C  "
        " O  "
        " CB "
        " OG "
  chain id="A" #residue_groups=1  ### WARNING: duplicate chain id ###
    resid=" 505 " #atom_groups=1
      altloc="" resname="MN" #atoms=1
        "MN  "
""")
  #
  input_1ezu = """\
data_1EZU
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
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   3422 N  N   . GLY C 2 164 ? -25.713 -9.765  41.937  1.00 38.56 584 GLY C N   1
ATOM   3423 C  CA  . GLY C 2 164 ? -27.027 -9.485  41.364  1.00 38.07 584 GLY C CA  1
ATOM   3424 C  C   . GLY C 2 164 ? -27.645 -10.541 40.476  1.00 35.83 584 GLY C C   1
ATOM   3425 O  O   . GLY C 2 164 ? -27.554 -11.745 40.752  1.00 35.99 584 GLY C O   1
ATOM   3426 N  N   . PHE C 2 165 A -28.265 -10.093 39.385  1.00 34.13 584 PHE C N   1
ATOM   3427 C  CA  . PHE C 2 165 A -28.937 -11.016 38.471  1.00 36.09 584 PHE C CA  1
ATOM   3428 C  C   . PHE C 2 165 A -28.697 -10.685 37.018  1.00 36.74 584 PHE C C   1
ATOM   3429 O  O   . PHE C 2 165 A -28.853 -9.536  36.627  1.00 39.19 584 PHE C O   1
ATOM   3430 C  CB  . PHE C 2 165 A -30.444 -10.982 38.748  1.00 35.84 584 PHE C CB  1
ATOM   3431 C  CG  . PHE C 2 165 A -30.802 -11.298 40.174  1.00 39.04 584 PHE C CG  1
ATOM   3432 C  CD1 . PHE C 2 165 A -30.914 -10.282 41.119  1.00 39.32 584 PHE C CD1 1
ATOM   3433 C  CD2 . PHE C 2 165 A -30.983 -12.614 40.579  1.00 36.86 584 PHE C CD2 1
ATOM   3434 C  CE1 . PHE C 2 165 A -31.203 -10.572 42.458  1.00 35.08 584 PHE C CE1 1
ATOM   3435 C  CE2 . PHE C 2 165 A -31.271 -12.924 41.904  1.00 32.57 584 PHE C CE2 1
ATOM   3436 C  CZ  . PHE C 2 165 A -31.379 -11.898 42.850  1.00 32.97 584 PHE C CZ  1
"""
  cif_model = iotbx.cif.reader(input_string=input_1ezu).model()
  cif_block = cif_model["1EZU"]
  builder = pdb_hierarchy_builder(cif_block)
  hierarchy = builder.hierarchy
  residue_groups = list(hierarchy.residue_groups())
  assert residue_groups[0].icode == " "
  assert residue_groups[1].icode == "A"
  s = StringIO()
  hierarchy.show(out=s)
  assert not show_diff(s.getvalue(), """\
model id="" #chains=1
  chain id="C" #residue_groups=2
    resid=" 584 " #atom_groups=1
      altloc="" resname="GLY" #atoms=4
        " N  "
        " CA "
        " C  "
        " O  "
    resid=" 584A" #atom_groups=1
      altloc="" resname="PHE" #atoms=11
        " N  "
        " CA "
        " C  "
        " O  "
        " CB "
        " CG "
        " CD1"
        " CD2"
        " CE1"
        " CE2"
        " CZ "
""")

  input_charges = """\
data_charges
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
_atom_site.Cartn_x_esd
_atom_site.Cartn_y_esd
_atom_site.Cartn_z_esd
_atom_site.occupancy_esd
_atom_site.B_iso_or_equiv_esd
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 2932 CA CA  . CA  J 2 .   ? 221.154 27.397 60.094 1.00 49.40  ? ? ? ? ? 2 1123 CA  C CA  1
HETATM 2996 O  O1S . MES D 4 .   ? 47.470 -6.157  23.319 1.00 13.84 ? ? ? ? ? -1 1653 MES F O1S 1
HETATM 2997 O  O2S . MES D 4 .   ? 47.327 -5.296  20.939 1.00 15.26 ? ? ? ? ? ?  1653 MES F O2S 1
"""
  cif_model = iotbx.cif.reader(input_string=input_charges).model()
  cif_block = cif_model["charges"]
  builder = pdb_hierarchy_builder(cif_block)
  hierarchy = builder.hierarchy
  atoms = hierarchy.atoms()
  assert atoms[0].charge == "2+"
  assert atoms[1].charge == "1-"
  assert atoms[2].charge == ""
  s = StringIO()
  hierarchy.show(out=s)
  assert not show_diff(s.getvalue(), """\
model id="" #chains=2
  chain id="C" #residue_groups=1
    resid="1123 " #atom_groups=1
      altloc="" resname="CA" #atoms=1
        "CA  "
  chain id="F" #residue_groups=1
    resid="1653 " #atom_groups=1
      altloc="" resname="MES" #atoms=2
        " O1S"
        " O2S"
""")

def run():
  exercise_pdb_hierachy_builder()


if __name__ == '__main__':
  run()
  print "OK"
