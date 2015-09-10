from __future__ import division
from cctbx.array_family import flex
from mmtbx.monomer_library import server, pdb_interpretation
from mmtbx.geometry_restraints import c_beta
import iotbx

pdb_str_1 = """\
CRYST1   26.960   29.455   29.841  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   TYR A  20      10.702  10.331   8.954  1.00 23.81           N
ATOM      2  CA  TYR A  20      11.170   9.547   7.817  1.00 35.69           C
ATOM      3  C   TYR A  20      10.082   8.601   7.319  1.00 28.94           C
ATOM      4  O   TYR A  20       9.756   7.616   7.982  1.00 25.26           O
ATOM      5  CB  TYR A  20      12.426   8.758   8.192  1.00 33.37           C
ATOM      6  CG  TYR A  20      13.016   7.958   7.052  1.00 31.48           C
ATOM      7  CD1 TYR A  20      13.847   8.558   6.115  1.00 36.91           C
ATOM      8  CD2 TYR A  20      12.747   6.602   6.916  1.00 21.36           C
ATOM      9  CE1 TYR A  20      14.389   7.832   5.071  1.00 36.56           C
ATOM     10  CE2 TYR A  20      13.285   5.868   5.875  1.00 29.56           C
ATOM     11  CZ  TYR A  20      14.105   6.487   4.957  1.00 35.08           C
ATOM     12  OH  TYR A  20      14.644   5.760   3.920  1.00 38.73           O
ATOM     13  N   ARG A  21       9.530   8.915   6.149  1.00 38.95           N
ATOM     14  CA  ARG A  21       8.473   8.118   5.525  1.00 38.77           C
ATOM     15  C   ARG A  21       7.250   7.963   6.428  1.00 27.69           C
ATOM     16  O   ARG A  21       7.134   6.992   7.176  1.00 22.82           O
ATOM     17  CB  ARG A  21       9.004   6.747   5.093  1.00 20.00           C
ATOM     18  CG  ARG A  21       8.015   5.920   4.285  1.00 20.00           C
ATOM     19  CD  ARG A  21       8.608   4.577   3.893  1.00 20.00           C
ATOM     20  NE  ARG A  21       7.671   3.771   3.116  1.00 20.00           N
ATOM     21  CZ  ARG A  21       7.939   2.556   2.649  1.00 20.00           C
ATOM     22  NH1 ARG A  21       9.121   2.001   2.879  1.00 20.00           N
ATOM     23  NH2 ARG A  21       7.025   1.895   1.951  1.00 20.00           N
ATOM     24  N   GLY A  22       6.340   8.929   6.351  1.00 24.85           N
ATOM     25  CA  GLY A  22       5.132   8.903   7.154  1.00 29.53           C
ATOM     26  C   GLY A  22       5.373   9.358   8.580  1.00 33.22           C
ATOM     27  O   GLY A  22       5.196  10.531   8.906  1.00 30.06           O
"""

pdb_str_2 = """\
CRYST1   72.072   33.173   34.033  90.00  90.00  90.00 P 1
ATOM    533  N   GLU B   3       7.910  11.239  20.396  1.00 20.00           N
ATOM    534  CA  GLU B   3       8.310  11.284  21.798  1.00 20.00           C
ATOM    534  CA  GLU B   3       8.310  11.284  21.798  1.00 20.00           C
ATOM    535  C   GLU B   3       9.344  10.190  21.979  1.00 20.00           C
ATOM    536  O   GLU B   3      10.197  10.267  22.867  1.00 20.00           O
ATOM    537  CB  GLU B   3       7.115  11.041  22.731  1.00 20.00           C
ATOM    538  HA  GLU B   3       8.761  12.248  22.034  1.00 20.00           H
ATOM    539  H  AGLU B   3       7.474  10.360  20.122  0.50 20.00           H
ATOM    540  D  BGLU B   3       7.474  10.360  20.122  0.50 20.00           D
"""

pdb_str_3 = """\
CRYST1   72.072   33.173   34.033  90.00  90.00  90.00 P 1
ATOM    533  N   GLU B   3       7.910  11.239  20.396  1.00 20.00           N
ATOM    534  CA AGLU B   3       8.310  11.284  21.798  1.00 20.00           C
ATOM    534  CA BGLU B   3       8.310  11.284  21.798  1.00 20.00           C
ATOM    535  C   GLU B   3       9.344  10.190  21.979  1.00 20.00           C
ATOM    536  O   GLU B   3      10.197  10.267  22.867  1.00 20.00           O
ATOM    537  CB  GLU B   3       7.115  11.041  22.731  1.00 20.00           C
"""

pdb_str_4 = """\
CRYST1   72.072   33.173   34.033  90.00  90.00  90.00 P 1
ATOM    533  N   GLU B   3       7.910  11.239  20.396  1.00 20.00           N
ATOM    534  CA AGLU B   3       8.310  11.284  21.798  1.00 20.00           C
ATOM    535  C   GLU B   3       9.344  10.190  21.979  1.00 20.00           C
ATOM    536  O   GLU B   3      10.197  10.267  22.867  1.00 20.00           O
ATOM    537  CB  GLU B   3       7.115  11.041  22.731  1.00 20.00           C
"""

def exercise_1():
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv              = server.server(),
    ener_lib                 = server.ener_lib(),
    raw_records              = flex.std_string(pdb_str_1.splitlines()),
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)
  grm = processed_pdb_file.geometry_restraints_manager()
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  # c-beta restraints are added by default!!!
  assert len(grm.get_c_beta_torsion_proxies()) == 4

  #test global selection and removing c-beta restraints
  tst_boolsel = pdb_hierarchy.atom_selection_cache().selection("resname TYR")
  tst_iselection = tst_boolsel.iselection()
  #test global selection
  grm2 = grm.select(iselection=tst_iselection)
  assert len(grm2.get_c_beta_torsion_proxies()) == 2
  grm2 = grm.select(selection=tst_boolsel)
  assert len(grm2.get_c_beta_torsion_proxies()) == 2
  #remove a selection
  grm.remove_c_beta_torsion_restraints_in_place(selection=tst_iselection)
  assert len(grm.get_c_beta_torsion_proxies()) == 2
  #add a selection
  grm.remove_c_beta_torsion_restraints_in_place()
  assert len(grm.get_c_beta_torsion_proxies()) == 0
  c_beta_torsion_proxies = c_beta.get_c_beta_torsion_proxies(
      pdb_hierarchy,
      selection=tst_iselection,
      sigma=2.5)
  assert len(c_beta_torsion_proxies) == 2

def exercise_2():
  """
  Testing with ACs
  """
  pdb_h = iotbx.pdb.input(
      source_info=None,
      lines=pdb_str_2).construct_hierarchy()
  c_beta_restrs = c_beta.get_c_beta_torsion_proxies(pdb_h)
  assert len(c_beta_restrs) == 2

  pdb_h = iotbx.pdb.input(
      source_info=None,
      lines=pdb_str_3).construct_hierarchy()
  c_beta_restrs = c_beta.get_c_beta_torsion_proxies(pdb_h)
  assert len(c_beta_restrs) == 4

  pdb_h = iotbx.pdb.input(
      source_info=None,
      lines=pdb_str_4).construct_hierarchy()
  c_beta_restrs = c_beta.get_c_beta_torsion_proxies(pdb_h)
  assert len(c_beta_restrs) == 2

if (__name__ == "__main__") :
  exercise_1()
  exercise_2()
  print "OK"
