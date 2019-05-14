from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
import os.path
import time

pdb_lines = """\
CRYST1  194.138  101.751   80.651  90.00 105.74  90.00 C 1 2 1       4
ATOM      1  N   MET A   1     -63.482  31.658  20.489  1.00 88.61           N
ATOM      2  CA  MET A   1     -64.594  31.042  19.724  1.00 87.28           C
ATOM      3  C   MET A   1     -65.689  30.614  20.691  1.00 85.59           C
ATOM      4  O   MET A   1     -65.419  30.412  21.876  1.00 85.68           O
ATOM      5  CB  MET A   1     -64.075  29.831  18.973  1.00 89.26           C
ATOM      6  CG  MET A   1     -62.902  30.150  18.065  1.00 92.21           C
ATOM      7  SD  MET A   1     -62.278  28.659  17.271  1.00 96.77           S
ATOM      8  CE  MET A   1     -63.798  27.825  16.920  1.00 94.61           C
ATOM      9  N   THR A   2     -66.924  30.488  20.211  1.00 83.63           N
ATOM     10  CA  THR A   2     -68.002  30.078  21.106  1.00 83.02           C
ATOM     11  C   THR A   2     -69.152  29.283  20.508  1.00 81.29           C
ATOM     12  O   THR A   2     -69.485  29.420  19.332  1.00 79.88           O
ATOM     13  CB  THR A   2     -68.642  31.264  21.829  1.00 84.34           C
ATOM     14  OG1 THR A   2     -69.678  31.807  21.006  1.00 86.38           O
ATOM     15  CG2 THR A   2     -67.602  32.328  22.151  1.00 86.86           C
ATOM     16  N   MET A   3     -69.768  28.490  21.387  1.00 79.99           N
ATOM     17  CA  MET A   3     -70.888  27.620  21.082  1.00 77.22           C
ATOM     18  C   MET A   3     -70.424  26.519  20.160  1.00 76.93           C
ATOM     19  O   MET A   3     -69.748  25.598  20.606  1.00 78.22           O
ATOM     20  CB  MET A   3     -72.019  28.431  20.476  1.00 76.32           C
ATOM     21  CG  MET A   3     -72.439  29.543  21.406  1.00 75.48           C
ATOM     22  SD  MET A   3     -72.756  28.964  23.098  1.00 75.99           S
ATOM     23  CE  MET A   3     -71.183  29.165  23.909  1.00 73.62           C
ATOM     24  N   ASP A   4     -70.753  26.588  18.880  1.00 76.31           N
ATOM     25  CA  ASP A   4     -70.290  25.529  18.001  1.00 75.92           C
ATOM     26  C   ASP A   4     -69.444  26.061  16.882  1.00 75.06           C
ATOM     27  O   ASP A   4     -69.938  26.667  15.939  1.00 75.09           O
ATOM     28  CB  ASP A   4     -71.451  24.727  17.426  1.00 76.22           C
ATOM     29  CG  ASP A   4     -70.984  23.449  16.763  1.00 76.96           C
ATOM     30  OD1 ASP A   4     -70.169  23.538  15.816  1.00 76.92           O
ATOM     31  OD2 ASP A   4     -71.428  22.362  17.205  1.00 77.39           O
ATOM     32  N   PHE A   5     -68.156  25.798  17.003  1.00 74.26           N
ATOM     33  CA  PHE A   5     -67.167  26.224  16.037  1.00 74.26           C
ATOM     34  C   PHE A   5     -67.515  25.771  14.623  1.00 72.79           C
ATOM     35  O   PHE A   5     -67.028  26.337  13.644  1.00 72.51           O
ATOM     36  CB  PHE A   5     -65.831  25.636  16.442  1.00 77.83           C
ATOM     37  CG  PHE A   5     -65.602  25.613  17.931  1.00 82.00           C
ATOM     38  CD1 PHE A   5     -64.713  24.690  18.492  1.00 84.76           C
ATOM     39  CD2 PHE A   5     -66.249  26.521  18.775  1.00 84.12           C
ATOM     40  CE1 PHE A   5     -64.473  24.672  19.871  1.00 85.57           C
ATOM     41  CE2 PHE A   5     -66.017  26.516  20.156  1.00 85.19           C
ATOM     42  CZ  PHE A   5     -65.127  25.589  20.708  1.00 86.02           C
HETATM 8902  C1  MLI A1001     -35.107  -8.183  19.606  1.00 48.70           C
HETATM 8903  C2  MLI A1001     -34.404  -7.934  20.966  1.00 47.76           C
HETATM 8904  C3  MLI A1001     -35.406  -7.002  18.730  1.00 51.13           C
HETATM 8905  O6  MLI A1001     -33.181  -7.738  21.141  1.00 48.11           O
HETATM 8906  O7  MLI A1001     -35.151  -7.962  21.978  1.00 44.44           O
HETATM 8907  O8  MLI A1001     -34.422  -6.274  18.392  1.00 48.51           O
HETATM 8908  O9  MLI A1001     -36.601  -6.959  18.453  1.00 53.02           O
"""

cif_lines = """\
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
MLI MLI "Unknown                  " ligand 9 7 .
#
data_comp_MLI
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
MLI        C1      C   CH2    0 .         -0.3524   -0.0000   -0.4983
MLI        C2      C   C      0 .          1.1747   -0.0000   -0.4983
MLI        C3      C   C      0 .         -0.8614   -0.0000    0.9415
MLI        O6      O   O      0 .          1.8042    1.0799   -0.3478
MLI        O7      O   OC    -1 .          1.8042   -1.0799   -0.6489
MLI        O8      O   O      0 .         -1.2112    1.0802    1.4855
MLI        O9      O   OC    -1 .         -0.9313   -1.0802    1.5844
MLI        H11     H   HCH2   0 .         -0.7135    0.8845   -1.0090
MLI        H12     H   HCH2   0 .         -0.7135   -0.8845   -1.0090
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
MLI  C1      C2     single        1.527 0.020     1.527
MLI  C1      C3     single        1.527 0.020     1.527
MLI  C1      H11    single        0.970 0.020     1.090
MLI  C1      H12    single        0.970 0.020     1.090
MLI  C2      O6     deloc         1.259 0.020     1.259
MLI  C2      O7     deloc         1.259 0.020     1.259
MLI  C3      O8     deloc         1.259 0.020     1.259
MLI  C3      O9     deloc         1.259 0.020     1.259
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
MLI  H12     C1      H11          109.47 3.000
MLI  H12     C1      C3           109.47 3.000
MLI  H11     C1      C3           109.47 3.000
MLI  H12     C1      C2           109.47 3.000
MLI  H11     C1      C2           109.47 3.000
MLI  C3      C1      C2           109.47 3.000
MLI  O7      C2      O6           120.00 3.000
MLI  O7      C2      C1           120.00 3.000
MLI  O6      C2      C1           120.00 3.000
MLI  O9      C3      O8           120.00 3.000
MLI  O9      C3      C1           120.00 3.000
MLI  O8      C3      C1           120.00 3.000
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
MLI Var_01        O8      C3      C1      C2            97.82  30.0 3
MLI Var_02        O6      C2      C1      C3           -82.06  30.0 3
#
loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
MLI plan-1  C1     0.020
MLI plan-1  C2     0.020
MLI plan-1  O6     0.020
MLI plan-1  O7     0.020
MLI plan-2  C1     0.020
MLI plan-2  C3     0.020
MLI plan-2  O8     0.020
MLI plan-2  O9     0.020
"""

def exercise_01(prefix="tst_mi_ligands_test_01"):
  """
  Simple run to a completion with reference map. no SS annotations.
  """
  pdb_file = open("%s_start.pdb" % prefix, "w")
  pdb_file.write(pdb_lines)
  pdb_file.close()
  cif_file = open("%s_start.cif" % prefix, "w")
  cif_file.write(cif_lines)
  cif_file.close()
  cmd = " ".join([
      "phenix.model_idealization",
      "%s_start.pdb" % prefix,
      "%s_start.cif" % prefix,
      "use_map_for_reference=True",
      "number_of_refinement_cycles=1",
      "run_minimization_first=False",
      "loop_idealization.number_of_ccd_trials=1",
      "n_macro=1",
      "debug=True",
      ">%s.log" % prefix])
  print(cmd)
  assert not easy_run.call(cmd)
  res_log = open("%s.log" % prefix, "r")
  log_lines = res_log.readlines()
  for l in [
      # "Secondary structure substitution step will be skipped\n",
      "  Minimizing...\n",
      "Using map as reference\n",
      # "Ramachandran outliers:      0.00      0.00      0.00      0.00      0.00\n",
      "All done.\n"]:
    assert l in log_lines, "'%s' not in log file." % l
  res_log.close()
  # assert os.path.isfile("%s_start.pdb_idealized.pdb" % prefix)

if (__name__ == "__main__"):
  t0 = time.time()
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping: probe not configured")
  else:
    exercise_01()
  print("Time: %.2f" % (time.time() - t0))
  print("OK")
