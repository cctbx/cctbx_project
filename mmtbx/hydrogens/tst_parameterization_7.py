from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
import iotbx.cif
from libtbx.utils import null_out
from six.moves import zip

def prepare_inputs(pdb_str, cif_str):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  cif_object = iotbx.cif.reader(input_string = cif_str).model()
  # bla.cif does not exist, but cif_objects needs a filename in first position
  # of the tuple
  cif_objects = [('bla.cif', cif_object)]

  model = mmtbx.model.manager(
    model_input=pdb_inp,
    restraint_objects = cif_objects,
    log = null_out())
  model.process(make_restraints=True)
  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()

  return model

#-----------------------------------------------------------------------------
# For this ligand, three H cannot be parameterized
# NH2 group is planar, but there are no dihedral restraints
# no dihedral restraint for H11
#-----------------------------------------------------------------------------

def exercise1(pdb_str, cif_str):
  model = prepare_inputs(pdb_str, cif_str)
  riding_h_manager = model.get_riding_h_manager()
  atoms = model.get_hierarchy().atoms()

  h_para = riding_h_manager.h_parameterization

  diagnostics = riding_h_manager.diagnostics(
    sites_cart = model.get_sites_cart(),
    threshold  = 0.05)
  h_distances   = diagnostics.h_distances
  type_list     = diagnostics.type_list

# number of H atoms
  number_h = model.get_hd_selection().count(True)
  number_h_para = len(h_para) - h_para.count(None)

  assert (number_h_para == number_h-3), 'Not all H atoms are parameterized'

  for ih in h_distances:
    # One atom is expected to be moved
    if (ih == 16):
      continue
    labels = atoms[ih].fetch_labels()
    assert (h_distances[ih] < 0.1), \
      'distance too large: %s  atom: %s (%s) residue: %s ' \
      % (h_para[ih].htype, atoms[ih].name, ih, labels.resseq.strip())

  for type1, type2 in zip(type_list, type_list_known1):
    assert (type1 == type2)
    #print "'%s'," % type1,

def exercise2(pdb_str, cif_str):
  model = prepare_inputs(pdb_str, cif_str)
  riding_h_manager = model.get_riding_h_manager()
  atoms = model.get_hierarchy().atoms()

  h_para = riding_h_manager.h_parameterization

  diagnostics = riding_h_manager.diagnostics(
    sites_cart = model.get_sites_cart(),
    threshold  = 0.05)
  h_distances   = diagnostics.h_distances
  type_list     = diagnostics.type_list

# number of H atoms
  number_h = model.get_hd_selection().count(True)
  number_h_para = len(h_para) - h_para.count(None)

  assert (number_h_para == 0), 'Not all H atoms are parameterized'

pdb_str1 = """\
REMARK iotbx.pdb.box_around_molecule --buffer-layer=5 "02.updated.nomin..pdb"
REMARK Date 2017-07-10 Time 12:11:01 PDT -0700 (1499713861.46 s)
CRYST1   15.938   15.073   11.550  90.00  90.00  90.00 P 1
SCALE1      0.062743  0.000000  0.000000        0.00000
SCALE2      0.000000  0.066344  0.000000        0.00000
SCALE3      0.000000  0.000000  0.086580        0.00000
HETATM    1  C1  3HA   401       7.842   9.143   5.222  1.00 83.17           C
HETATM    2  C2  3HA   401       7.240   7.911   5.459  1.00 81.93           C
HETATM    3  C3  3HA   401       8.027   6.770   5.589  1.00 82.00           C
HETATM    4  C4  3HA   401       9.411   6.864   5.486  1.00 82.52           C
HETATM    5  C5  3HA   401      10.012   8.097   5.252  1.00 83.82           C
HETATM    6  C6  3HA   401       9.226   9.238   5.118  1.00 84.24           C
HETATM    7  C7  3HA   401       5.719   7.834   5.619  1.00 81.09           C
HETATM    8  N10 3HA   401       7.463   5.583   5.793  1.00 82.31           N
HETATM    9  O11 3HA   401      10.183   5.752   5.619  1.00 82.09           O
HETATM   10  O8  3HA   401       5.197   6.868   6.174  1.00 81.05           O
HETATM   11  O9  3HA   401       5.000   8.734   5.190  1.00 78.80           O
HETATM   12  H1  3HA   401       7.315   9.910   5.134  1.00 83.17           H
HETATM   13  H11 3HA   401      10.719   5.854   6.276  1.00 82.09           H
HETATM   14  H5  3HA   401      10.938   8.151   5.134  1.00 83.82           H
HETATM   15  H6  3HA   401       9.630  10.073   5.000  1.00 84.24           H
HETATM   16 H101 3HA   401       7.572   5.175   6.550  1.00 82.31           H
HETATM   17 H102 3HA   401       7.435   5.000   5.152  1.00 82.31           H
TER
END
"""


cif_str1 = """
#
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
3HA 3HA "Unknown                  " ligand 17 11 .
#
data_comp_3HA
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
3HA        O8      O   OC    -1 .         -2.5304   -0.3544    2.0294
3HA        C7      C   C      0 .         -1.2760   -0.4611    2.0390
3HA        O9      O   O      0 .         -0.6571   -0.5146    3.1342
3HA        C2      C   CR6    0 .         -0.5023   -0.4939    0.7224
3HA        C1      C   CR16   0 .         -0.2734   -1.7018    0.0834
3HA        C6      C   CR16   0 .          0.4107   -1.7296   -1.1210
3HA        C5      C   CR16   0 .          0.8659   -0.5495   -1.6865
3HA        C4      C   CR6    0 .          0.6370    0.6583   -1.0475
3HA        O11     O   OH1    0 .          1.1017    1.8508   -1.6167
3HA        C3      C   CR6    0 .         -0.0471    0.6861    0.1569
3HA        N10     N   NH2    0 .         -0.3039    1.9536    0.8168
3HA        H1      H   HCR6   0 .         -0.6294   -2.6246    0.5256
3HA        H6      H   HCR6   0 .          0.5887   -2.6740   -1.6213
3HA        H5      H   HCR6   0 .          1.4009   -0.5713   -2.6283
3HA        H11     H   HOH1   0 .          2.0033    1.7449   -1.8770
3HA        H101    H   HNH2   0 .          0.0644    2.1253    1.7332
3HA        H102    H   HNH2   0 .         -0.8528    2.6556    0.3576
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
3HA  O8      C7     deloc         1.259 0.020     1.259
3HA  C7      O9     deloc         1.259 0.020     1.259
3HA  C7      C2     single        1.527 0.020     1.527
3HA  C2      C1     aromatic      1.386 0.020     1.386
3HA  C2      C3     aromatic      1.385 0.020     1.385
3HA  C1      C6     aromatic      1.385 0.020     1.385
3HA  C1      H1     single        0.930 0.020     1.080
3HA  C6      C5     aromatic      1.385 0.020     1.385
3HA  C6      H6     single        0.930 0.020     1.080
3HA  C5      C4     aromatic      1.385 0.020     1.385
3HA  C5      H5     single        0.930 0.020     1.080
3HA  C4      O11    single        1.401 0.020     1.401
3HA  C4      C3     aromatic      1.385 0.020     1.385
3HA  O11     H11    single        0.850 0.020     0.980
3HA  C3      N10    single        1.452 0.020     1.452
3HA  N10     H101   single        0.860 0.020     1.020
3HA  N10     H102   single        0.860 0.020     1.020
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
3HA  C2      C7      O9           120.00 3.000
3HA  C2      C7      O8           119.99 3.000
3HA  O9      C7      O8           120.00 3.000
3HA  C3      C2      C1           120.00 3.000
3HA  C3      C2      C7           120.00 3.000
3HA  C1      C2      C7           120.00 3.000
3HA  H1      C1      C6           120.00 3.000
3HA  H1      C1      C2           120.00 3.000
3HA  C6      C1      C2           120.00 3.000
3HA  H6      C6      C5           120.00 3.000
3HA  H6      C6      C1           120.00 3.000
3HA  C5      C6      C1           120.00 3.000
3HA  H5      C5      C4           120.00 3.000
3HA  H5      C5      C6           120.00 3.000
3HA  C4      C5      C6           120.00 3.000
3HA  C3      C4      O11          120.00 3.000
3HA  C3      C4      C5           120.00 3.000
3HA  O11     C4      C5           120.00 3.000
3HA  H11     O11     C4           109.47 3.000
3HA  N10     C3      C4           120.00 3.000
3HA  N10     C3      C2           120.00 3.000
3HA  C4      C3      C2           120.00 3.000
3HA  H102    N10     H101         120.00 3.000
3HA  H102    N10     C3           120.00 3.000
3HA  H101    N10     C3           120.00 3.000
#
loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
3HA CONST_01      C5      C6      C1      C2             0.00   0.0 0
3HA CONST_02      C5      C4      C3      C2             0.00   0.0 0
3HA CONST_03      C4      C3      C2      C1            -0.00   0.0 0
3HA CONST_04      C4      C5      C6      C1            -0.00   0.0 0
3HA CONST_05      C3      C2      C1      C6             0.00   0.0 0
3HA CONST_06      C3      C4      C5      C6             0.00   0.0 0
3HA CONST_07      C6      C1      C2      C7           179.02   0.0 0
3HA CONST_08      C4      C3      C2      C7          -179.02   0.0 0
3HA CONST_09      O11     C4      C3      C2          -179.76   0.0 0
3HA CONST_10      N10     C3      C2      C1           179.11   0.0 0
3HA CONST_11      O11     C4      C5      C6           179.76   0.0 0
3HA CONST_12      N10     C3      C4      C5          -179.11   0.0 0
3HA CONST_13      H6      C6      C1      C2          -179.93   0.0 0
3HA CONST_14      H5      C5      C6      C1           180.00   0.0 0
3HA CONST_15      H1      C1      C6      C5          -180.00   0.0 0
3HA CONST_16      H101    N10     C3      C2            61.62   0.0 0
3HA CONST_17      H102    N10     C3      C2          -118.12   0.0 0
3HA Var_01        C1      C2      C7      O8           -88.86  30.0 2
#
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
3HA plan-1  C7     0.020
3HA plan-1  C2     0.020
3HA plan-1  C1     0.020
3HA plan-1  C6     0.020
3HA plan-1  C5     0.020
3HA plan-1  C4     0.020
3HA plan-1  O11    0.020
3HA plan-1  C3     0.020
3HA plan-1  N10    0.020
3HA plan-1  H1     0.020
3HA plan-1  H6     0.020
3HA plan-1  H5     0.020
3HA plan-2  C3     0.020
3HA plan-2  N10    0.020
3HA plan-2  H101   0.020
3HA plan-2  H102   0.020
3HA plan-3  O8     0.020
3HA plan-3  C7     0.020
3HA plan-3  O9     0.020
3HA plan-3  C2     0.020
"""

#type_list_known1 = ['flat_2neigbs', 'alg1b', 'flat_2neigbs', 'flat_2neigbs',
#  'alg1a', 'alg1a']
type_list_known1 = ['flat_2neigbs', 'flat_2neigbs', 'flat_2neigbs']

# Zundelion ( [H2O -- H -- OH2]+)

pdb_str2 = """
CRYST1   12.247   13.433   12.428  90.00  90.00  90.00 P 1
SCALE1      0.081653  0.000000  0.000000        0.00000
SCALE2      0.000000  0.074444  0.000000        0.00000
SCALE3      0.000000  0.000000  0.080463        0.00000
HETATM    1  O   ZDL A   1       0.414   0.852  -0.974  1.00 30.00           O
HETATM    2  H1  ZDL A   1      -0.135   1.630  -0.847  1.00 30.00           H
HETATM    3  H2  ZDL A   1       1.338   1.093  -0.875  1.00 30.00           H
HETATM    4  O1  ZDL A   1      -0.233  -1.018   0.869  1.00 30.00           O
HETATM    5  H3  ZDL A   1      -0.566  -1.803   0.427  1.00 30.00           H
HETATM    6  H4  ZDL A   1      -0.909  -0.670   1.454  1.00 30.00           H
HETATM    7  H5  ZDL A   1       0.091  -0.084  -0.054  1.00 30.00           H
"""

cif_str2 = """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
ZDL        ZDL 'zundel ion                  ' ligand 7 6 .
#
data_comp_ZDL
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
ZDL         O      O   O      0    .       1.8400   38.7010   -3.2140
ZDL         H1     H   H      0    .       1.3350   39.4980   -3.3000
ZDL         H2     H   H      0    .       2.7520   38.9490   -3.2430
ZDL         O1     O   O      0    .       1.0410   36.8380   -1.5510
ZDL         H3     H   H      0    .       0.8620   35.9420   -1.8460
ZDL         H4     H   H      0    .       0.5460   37.0160   -0.7480
ZDL         H5     H   H     +1    .       1.2760   37.8090   -2.1650
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
ZDL   O       H1    single        0.960 0.020     0.960
ZDL   O       H2    single        0.960 0.020     0.960
ZDL   O       H5    single        1.350 0.350     1.350
ZDL   O1      H5    single        1.350 0.350     1.350
ZDL   O1      H4    single        0.960 0.020     0.960
ZDL   O1      H3    single        0.960 0.020     0.960
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
ZDL   H2      O       H1          109.47 3.000
ZDL   H3      O1      H4          109.47 3.000
ZDL   H1      O       H5          109.49 6.000
ZDL   H2      O       H5          109.49 6.000
ZDL   H3      O1      H5          109.49 6.000
ZDL   H4      O1      H5          109.49 6.000
ZDL   O       H5      O1          180.00 9.000
"""

if (__name__ == "__main__"):
  t0 = time.time()
  exercise1(pdb_str1, cif_str1)
  exercise2(pdb_str2, cif_str2)
  print("OK. Time: %8.3f"%(time.time()-t0))
