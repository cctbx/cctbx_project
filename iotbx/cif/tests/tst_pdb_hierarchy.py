from cStringIO import StringIO

from libtbx.test_utils import approx_equal
from libtbx.test_utils import Exception_expected
from libtbx.test_utils import show_diff

import iotbx.cif
from iotbx.cif.builders import pdb_hierarchy_builder


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
model id="1" #chains=1
  chain id="A" #residue_groups=2
    resid="   2 " #atom_groups=3
      altloc="A" resname="THR" #atoms=8
        "CA"
        "CB"
        "OG1"
        "CG2"
        "HB"
        "HG21"
        "HG22"
        "HG23"
      altloc="B" resname="THR" #atoms=8
        "CA"
        "CB"
        "OG1"
        "CG2"
        "HB"
        "HG21"
        "HG22"
        "HG23"
      altloc="" resname="THR" #atoms=4
        "C"
        "O"
        "H"
        "HA"
    resid="  66 " #atom_groups=2
      altloc="A" resname="EOH" #atoms=3
        "C1"
        "C2"
        "O"
      altloc="B" resname="EOH" #atoms=3
        "C1"
        "C2"
        "O"
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
  atoms = hierarchy.atoms()
  assert atoms[0].serial == '1'
  assert atoms[0].name == 'N'
  assert atoms[0].b == 23.68
  assert atoms[0].uij_is_defined()
  assert approx_equal(
    atoms[0].uij, (0.4097, 0.2916, 0.1984, 0.113, 0.0328, 0.0375))
  assert atoms[1].serial == '2'
  assert atoms[1].name == 'CA'
  assert atoms[1].b == 22.98
  assert atoms[1].uij_is_defined()
  assert approx_equal(
    atoms[1].uij, (0.4035, 0.2848, 0.1847, 0.1242, 0.0297, 0.0347))
  assert not atoms[6].uij_is_defined()
  assert atoms[6].serial == '2650'
  assert atoms[6].name == 'MN'
  assert atoms[6].b == 44.18
  assert approx_equal(atoms[6].uij, (-1, -1, -1, -1, -1, -1))
  s = StringIO()
  hierarchy.show(out=s)
  assert not show_diff(s.getvalue(), """\
model id="1" #chains=1
  chain id="A" #residue_groups=2
    resid=" 108 " #atom_groups=1
      altloc="" resname="SER" #atoms=6
        "N"
        "CA"
        "C"
        "O"
        "CB"
        "OG"
    resid=" 505 " #atom_groups=1
      altloc="" resname="MN" #atoms=1
        "MN"
""")


def run():
  exercise_pdb_hierachy_builder()


if __name__ == '__main__':
  run()
  print "OK"
