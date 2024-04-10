from __future__ import absolute_import, division, print_function
import time, sys
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex

def compare_XH_bond_length_to_ideal(model):
  geometry = model.get_restraints_manager().geometry
  atoms = model.get_hierarchy().atoms()
  sites_cart = model.get_sites_cart()
  bond_proxies_simple, asu = \
    geometry.get_all_bond_proxies(sites_cart = sites_cart)
  hd_selection = model.get_hd_selection()
  for bp in bond_proxies_simple:
    i, j = bp.i_seqs
    if hd_selection[i] or hd_selection[j]:
      assert approx_equal(atoms[i].distance(atoms[j]),
                          bp.distance_ideal,
                          eps=0.001)
def get_dist(s1, s2):
  return flex.sqrt((s1 - s2).dot())

def tst_0():
  """
  Check going back and forth does not accumulate errors.
  Also, confirm: the change cannot be undone.
  """
  pdb_inp = iotbx.pdb.input(lines=pdb_str1, source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  sites_cart_o = model.get_sites_cart()
  model.set_hydrogen_bond_length(
    use_neutron_distances=True, show=False, log=null_out)
  sites_cart_n = model.get_sites_cart()
  model.set_hydrogen_bond_length(
    use_neutron_distances=False, show=False, log=null_out)
  sites_cart_x = model.get_sites_cart()
  dist_mean = flex.mean(get_dist(sites_cart_n, sites_cart_x))
  assert dist_mean > 0.07
  #
  model.set_hydrogen_bond_length(
    use_neutron_distances=True, show=False, log=null_out)
  sites_cart_n1 = model.get_sites_cart()
  dist_mean = flex.mean(get_dist(sites_cart_n, sites_cart_n1))
  assert approx_equal(dist_mean, 0)
  #
  model.set_hydrogen_bond_length(
    use_neutron_distances=False, show=False, log=null_out)
  sites_cart_x1 = model.get_sites_cart()
  dist_mean = flex.mean(get_dist(sites_cart_x, sites_cart_x1))
  assert approx_equal(dist_mean, 0)
  # This is irreversible
  for sites_cart in [sites_cart_n, sites_cart_n1, sites_cart_x, sites_cart_x1]:
    dist_mean = flex.mean(get_dist(sites_cart, sites_cart_o))
    assert dist_mean > 0.005

#-------------------------------------------------------------------------------

def tst_1():
  '''
  input model: bond lengths close to neutron but not quite
  input grm: neutron bond lengths
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str1.split("\n"), source_info=None)
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(pdb_interpretation_params=params,
    make_restraints=True)
  # request neutron bond lengths
  model.set_hydrogen_bond_length(use_neutron_distances=True,
                                 show=False,
                                 log=sys.stdout)
  compare_XH_bond_length_to_ideal(model = model)
  #STOP()
  # request X-ray bond lengths
  pdb_inp = iotbx.pdb.input(lines=pdb_str1.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(pdb_interpretation_params=params, make_restraints=False)
  model.set_hydrogen_bond_length(use_neutron_distances=False,
                                 show=False,
                                 log=sys.stdout)
  compare_XH_bond_length_to_ideal(model = model)
#-------------------------------------------------------------------------------

def tst_2():
  '''
  input model: bond lengths close to neutron but not quite
  input grm: X-ray bond lengths
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str1.split("\n"), source_info=None)
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = False
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(pdb_interpretation_params = params)
  #
  # Check if selected individual bond lengths are correct
  atoms = model.get_hierarchy().atoms()
  assert (atoms[4].id_str() == 'pdb=" CB  ARG A  74 "')
  assert (atoms[18].id_str() == 'pdb=" HB2 ARG A  74 "')
  assert approx_equal(atoms[4].distance(atoms[18]), 1.08, eps = 0.005)
  assert (atoms[10].id_str() == 'pdb=" NH2 ARG A  74 "')
  assert (atoms[15].id_str() == 'pdb="DH21 ARG A  74 "')
  assert approx_equal(atoms[10].distance(atoms[15]), 0.99, eps = 0.005)
  # request neutron bond lengths
  model.set_hydrogen_bond_length(use_neutron_distances=True,
                                 show=False,
                                 log=sys.stdout)
  # Check if bond lengths are equal to ideal
  compare_XH_bond_length_to_ideal(model = model)
  # Check if selected individual bond lengths are correct
  atoms = model.get_hierarchy().atoms()
  assert approx_equal(atoms[4].distance(atoms[18]), 1.09, eps = 0.005)
  assert approx_equal(atoms[10].distance(atoms[15]), 1.02, eps = 0.005)
  # request X-ray bond lengths
  model.set_hydrogen_bond_length(use_neutron_distances=False,
                                 show=False,
                                 log=sys.stdout)
  compare_XH_bond_length_to_ideal(model = model)
  # Check if selected individual bond lengths are correct
  atoms = model.get_hierarchy().atoms()
  assert approx_equal(atoms[4].distance(atoms[18]), 0.97, eps = 0.005)
  assert approx_equal(atoms[10].distance(atoms[15]), 0.86, eps = 0.005)

#-------------------------------------------------------------------------------

def tst_3():
  '''
  Test if the modification works also when cif_objects are supplied
  (meaning there is a ligand cif file)
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str2.split("\n"), source_info=None)
  cif_object = iotbx.cif.reader(input_string = cif_str2).model()
  # bla.cif does not exist, but cif_objects needs a filename in first position
  # of the tuple
  cif_objects = [('bla.cif', cif_object)]
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = False
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    restraint_objects = cif_objects,
    log         = null_out())
  model.process(pdb_interpretation_params = params)
  # request neutron bond lengths
  model.set_hydrogen_bond_length(use_neutron_distances=True,
                                 show=False,
                                 log=sys.stdout)
  compare_XH_bond_length_to_ideal(model = model)
  # request X-ray bond lengths
  model.set_hydrogen_bond_length(use_neutron_distances=False,
                                 show=False,
                                 log=sys.stdout)
  compare_XH_bond_length_to_ideal(model = model)

#-------------------------------------------------------------------------------
pdb_str1 = """
CRYST1   17.692   14.832   14.146  90.00  90.00  90.00 P 1
SCALE1      0.056523  0.000000  0.000000        0.00000
SCALE2      0.000000  0.067422  0.000000        0.00000
SCALE3      0.000000  0.000000  0.070691        0.00000
ATOM      1  N   ARG A  74       5.000   7.873   6.729  1.00 26.20           N
ATOM      2  CA  ARG A  74       6.196   8.078   7.547  1.00 25.89           C
ATOM      3  C   ARG A  74       6.284   9.534   7.986  1.00 24.08           C
ATOM      4  O   ARG A  74       6.576   9.832   9.146  1.00 23.39           O
ATOM      5  CB  ARG A  74       7.477   7.738   6.773  1.00 28.75           C
ATOM      6  CG  ARG A  74       7.718   6.268   6.516  1.00 33.03           C
ATOM      7  CD  ARG A  74       9.095   6.031   5.871  1.00 36.08           C
ATOM      8  NE  ARG A  74      10.208   6.328   6.778  1.00 37.90           N
ATOM      9  CZ  ARG A  74      11.051   7.351   6.639  1.00 38.74           C
ATOM     10  NH1 ARG A  74      10.896   8.234   5.654  1.00 38.21           N
ATOM     11  NH2 ARG A  74      12.056   7.493   7.500  1.00 38.01           N
ATOM     12  D   ARG A  74       5.091   7.539   5.816  0.81 26.33           D
ATOM     13  DE  ARG A  74      10.345   5.711   7.532  0.24 37.60           D
ATOM     14 DH11 ARG A  74      10.140   8.133   5.004  0.17 38.08           D
ATOM     15 DH12 ARG A  74      11.533   9.006   5.559  1.00 36.47           D
ATOM     16 DH21 ARG A  74      12.692   8.250   7.395  0.05 38.56           D
ATOM     17 DH22 ARG A  74      12.178   6.838   8.249  0.77 39.36           D
ATOM     18  HA  ARG A  74       6.132   7.447   8.424  1.00 26.28           H
ATOM     19  HB2 ARG A  74       7.441   8.242   5.818  1.00 29.24           H
ATOM     20  HB3 ARG A  74       8.320   8.122   7.334  1.00 28.50           H
ATOM     21  HG2 ARG A  74       7.669   5.728   7.450  1.00 33.30           H
ATOM     22  HG3 ARG A  74       6.953   5.908   5.849  1.00 32.25           H
ATOM     23  HD2 ARG A  74       9.163   5.000   5.580  1.00 36.76           H
ATOM     24  HD3 ARG A  74       9.184   6.659   5.000  1.00 36.58           H
"""
# i_seqs for X-H bonds
#0 11 pdb=" N   ARG A  74 " pdb=" D   ARG A  74 "
#1 17 pdb=" CA  ARG A  74 " pdb=" HA  ARG A  74 "
#4 18 pdb=" CB  ARG A  74 " pdb=" HB2 ARG A  74 "
#4 19 pdb=" CB  ARG A  74 " pdb=" HB3 ARG A  74 "
#5 20 pdb=" CG  ARG A  74 " pdb=" HG2 ARG A  74 "
#5 21 pdb=" CG  ARG A  74 " pdb=" HG3 ARG A  74 "
#6 22 pdb=" CD  ARG A  74 " pdb=" HD2 ARG A  74 "
#6 23 pdb=" CD  ARG A  74 " pdb=" HD3 ARG A  74 "
#7 12 pdb=" NE  ARG A  74 " pdb=" DE  ARG A  74 "
#9 13 pdb=" NH1 ARG A  74 " pdb="DH11 ARG A  74 "
#9 14 pdb=" NH1 ARG A  74 " pdb="DH12 ARG A  74 "
#10 15 pdb=" NH2 ARG A  74 " pdb="DH21 ARG A  74 "
#10 16 pdb=" NH2 ARG A  74 " pdb="DH22 ARG A  74 "


pdb_str2 = """
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

cif_str2 = """
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

if (__name__ == "__main__"):
  t0 = time.time()
  tst_0()
  tst_1()
  tst_2()
  tst_3()
  print("OK. Time: %8.3f"%(time.time()-t0))
