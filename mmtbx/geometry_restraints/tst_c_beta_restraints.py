from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from mmtbx.geometry_restraints import c_beta
import mmtbx.model
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

pdb_str_5 = """\
ATOM      1  N   DAL A   1      -1.425  -0.261  -0.246  1.00 20.00      A    N
ATOM      2  CA  DAL A   1       0.026  -0.261  -0.246  1.00 20.00      A    C
ATOM      3  CB  DAL A   1       0.535  -0.261   1.194  1.00 20.00      A    C
ATOM      4  C   DAL A   1       0.535   0.986  -0.965  1.00 20.00      A    C
ATOM      5  O   DAL A   1       1.425   0.883  -1.851  1.00 20.00      A    O
ATOM      6  OXT DAL A   1       0.066   2.117  -0.674  1.00 20.00      A    O-1
"""

pdb_str_6 = """
ATOM      1  N   ALA A   1      -1.425  -0.261  -0.246  1.00 20.00      A    N
ATOM      2  CA  ALA A   1       0.026  -0.261  -0.246  1.00 20.00      A    C
ATOM      3  CB  ALA A   1       0.535  -0.261   1.194  1.00 20.00      A    C
ATOM      4  C   ALA A   1       0.535   0.986  -0.965  1.00 20.00      A    C
ATOM      5  O   ALA A   1       1.425   0.883  -1.851  1.00 20.00      A    O
ATOM      6  OXT ALA A   1       0.066   2.117  -0.674  1.00 20.00      A    O-1
"""

pdb_str_7 = """
CRYST1  405.833  405.833  236.667  90.00  90.00  90.00 P 1
ATOM  10021  N   ALA E1273     302.855 197.287 175.926  1.00 96.14           N
ATOM  10022  CA  ALA E1273     303.086 195.852 176.039  1.00 96.14           C
ATOM  10023  C   ALA E1273     303.598 195.351 174.696  1.00 96.14           C
ATOM  10024  O   ALA E1273     303.769 196.123 173.749  1.00 96.14           O
ATOM  10025  CB  ALA E1273     301.819 195.120 176.456  1.00174.39           C
ATOM  10026  N   HIS E1274     303.840 194.045 174.612  1.00 96.33           N
ATOM  10027  CA  HIS E1274     304.277 193.409 173.371  1.00 96.33           C
ATOM  10028  C   HIS E1274     303.842 191.944 173.370  1.00 96.33           C
ATOM  10029  O   HIS E1274     303.670 191.358 174.436  1.00 96.33           O
ATOM  10030  CB  HIS E1274     305.794 193.523 173.216  1.00174.39           C
ATOM  10031  CG  HIS E1274     306.463 192.245 172.818  1.00174.39           C
ATOM  10032  ND1 HIS E1274     306.808 191.271 173.730  1.00174.39           N
ATOM  10033  CD2 HIS E1274     306.848 191.779 171.606  1.00174.39           C
ATOM  10034  CE1 HIS E1274     307.379 190.261 173.098  1.00174.39           C
ATOM  10035  NE2 HIS E1274     307.416 190.544 171.809  1.00174.39           N
ATOM  10036  N   ARG E1275     303.671 191.339 172.195  1.00 96.53           N
ATOM  10037  CA  ARG E1275     303.301 189.931 172.174  1.00 96.53           C
ATOM  10038  C   ARG E1275     304.315 189.194 173.039  1.00 96.53           C
ATOM  10039  O   ARG E1275     305.308 188.665 172.530  1.00 96.53           O
ATOM  10040  CB  ARG E1275     303.279 189.375 170.749  1.00174.39           C
ATOM  10041  N   UNK E1276     304.062 189.142 174.347  1.00 95.20           N
ATOM  10042  CA  UNK E1276     305.103 188.914 175.332  1.00 95.20           C
ATOM  10043  C   UNK E1276     305.569 187.487 175.535  1.00 95.20           C
ATOM  10044  O   UNK E1276     305.814 187.071 176.671  1.00 95.20           O
ATOM  10045  CB  UNK E1276     304.640 189.472 176.669  1.00161.31           C
ATOM  10046  N   UNK E1277     305.690 186.727 174.455  1.00 94.99           N
ATOM  10047  CA  UNK E1277     306.151 185.349 174.528  1.00 94.99           C
ATOM  10048  C   UNK E1277     306.888 185.034 173.233  1.00 94.99           C
ATOM  10049  O   UNK E1277     307.382 185.940 172.554  1.00 94.99           O
ATOM  10050  CB  UNK E1277     304.991 184.389 174.752  1.00161.31           C
TER
END
"""

def exercise_1():
  pdb_inp = iotbx.pdb.input(lines=flex.std_string(pdb_str_1.splitlines()), source_info=None)
  model = mmtbx.model.manager(model_input = pdb_inp)
  model.process(make_restraints=True)
  grm = model.get_restraints_manager().geometry
  pdb_hierarchy = model.get_hierarchy()
  sites_cart = model.get_sites_cart()
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
  c_beta_restrs, c_beta_skip = c_beta.get_c_beta_torsion_proxies(pdb_h)
  assert len(c_beta_restrs) == 2

  pdb_h = iotbx.pdb.input(
      source_info=None,
      lines=pdb_str_3).construct_hierarchy()
  c_beta_restrs, c_beta_skip = c_beta.get_c_beta_torsion_proxies(pdb_h)
  assert len(c_beta_restrs) == 4

  pdb_h = iotbx.pdb.input(
      source_info=None,
      lines=pdb_str_4).construct_hierarchy()
  c_beta_restrs, c_beta_skip = c_beta.get_c_beta_torsion_proxies(pdb_h)
  assert len(c_beta_restrs) == 2

def exercise_3():
  """
  Testing d-peptide
  """
  pdb_h = iotbx.pdb.input(
      source_info=None,
      lines=pdb_str_5).construct_hierarchy()
  c_beta_restrs, c_beta_skip = c_beta.get_c_beta_torsion_proxies(pdb_h)
  assert len(c_beta_restrs) == 0
  assert len(c_beta_skip.get("d-peptide", []))==1
  assert len(c_beta_skip.get("-ve", []))==0

  pdb_h = iotbx.pdb.input(
      source_info=None,
      lines=pdb_str_6).construct_hierarchy()
  c_beta_restrs, c_beta_skip = c_beta.get_c_beta_torsion_proxies(pdb_h)
  assert len(c_beta_restrs) == 0
  assert len(c_beta_skip.get("d-peptide", []))==0
  assert len(c_beta_skip.get("-ve", []))==1

def exercise_4():
  """
  Testing UNK to be equivalent to ALA
  """
  pdb_h = iotbx.pdb.input(
      source_info=None,
      lines=pdb_str_7).construct_hierarchy()
  c_beta_restrs, c_beta_skip = c_beta.get_c_beta_torsion_proxies(pdb_h)
  assert c_beta_restrs.size()==10
  assert c_beta_skip == {}
  #
  pdb_h = iotbx.pdb.input(
      source_info=None,
      lines=pdb_str_7.replace("UNK","ALA")).construct_hierarchy()
  c_beta_restrs, c_beta_skip = c_beta.get_c_beta_torsion_proxies(pdb_h)
  assert c_beta_restrs.size()==10
  assert c_beta_skip == {}

if (__name__ == "__main__"):
  exercise_1()
  exercise_2()
  exercise_3()
  exercise_4()
  print("OK")
