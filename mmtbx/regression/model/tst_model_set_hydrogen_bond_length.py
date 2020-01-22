from __future__ import absolute_import, division, print_function
import time, sys
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal



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
    pdb_interpretation_params = params,
    log         = null_out())
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
    pdb_interpretation_params = params,
    log         = null_out())
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
  pass

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

if (__name__ == "__main__"):
  t0 = time.time()
  tst_1()
  tst_2()
#  tst_3()
  print("OK. Time: %8.3f"%(time.time()-t0))
