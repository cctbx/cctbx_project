from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
import os.path
import time

cryst_str = """\n
CRYST1   34.917   22.246   44.017  90.00  90.00  90.00 P 1
SCALE1      0.028639  0.000000  0.000000        0.00000
SCALE2      0.000000  0.044952  0.000000        0.00000
SCALE3      0.000000  0.000000  0.022718        0.00000
"""
# taken from phenix_regression/refinement/ncs/tst_ncs_0.py
pdb_str = """\
ATOM      1  N   ALA A   1      27.344  16.348  30.784  1.00 10.00           N
ATOM      2  CA  ALA A   1      26.429  15.281  31.335  1.00 10.00           C
ATOM      3  C   ALA A   1      26.610  14.025  30.603  1.00 10.00           C
ATOM      4  O   ALA A   1      26.479  13.979  29.356  1.00 10.00           O
ATOM      5  CB  ALA A   1      24.874  15.800  31.300  1.00 10.00           C
ATOM      1  N   ALA A   2      26.812  12.925  31.345  1.00 10.00           N
ATOM      2  CA  ALA A   2      27.084  11.577  30.797  1.00 10.00           C
ATOM      3  C   ALA A   2      25.856  10.737  30.707  1.00 10.00           C
ATOM      4  O   ALA A   2      25.741   9.860  29.891  1.00 10.00           O
ATOM      5  CB  ALA A   2      28.151  10.950  31.721  1.00 10.00           C
ATOM      1  N   ALA A   3      25.009  10.973  31.714  1.00 10.00           N
ATOM      2  CA  ALA A   3      23.621  10.543  31.560  1.00 10.00           C
ATOM      3  C   ALA A   3      23.023  11.008  30.214  1.00 10.00           C
ATOM      4  O   ALA A   3      22.786  10.233  29.249  1.00 10.00           O
ATOM      5  CB  ALA A   3      22.760  11.040  32.654  1.00 10.00           C
ATOM      1  N   ALA A   4      22.798  12.304  30.175  1.00 10.00           N
ATOM      2  CA  ALA A   4      22.329  13.084  28.981  1.00 10.00           C
ATOM      3  C   ALA A   4      23.116  12.816  27.721  1.00 10.00           C
ATOM      4  O   ALA A   4      22.533  12.805  26.670  1.00 10.00           O
ATOM      5  CB  ALA A   4      22.372  14.607  29.318  1.00 10.00           C
ATOM      1  N   ALA A   5      24.448  12.622  27.823  1.00 10.00           N
ATOM      2  CA  ALA A   5      25.228  12.407  26.573  1.00 10.00           C
ATOM      3  C   ALA A   5      25.222  10.947  26.143  1.00 10.00           C
ATOM      4  O   ALA A   5      25.386  10.664  24.983  1.00 10.00           O
ATOM      5  CB  ALA A   5      26.634  12.906  26.746  1.00 10.00           C
ATOM      1  N   ALA A   6      24.976  10.048  27.071  1.00 10.00           N
ATOM      2  CA  ALA A   6      24.857   8.614  26.805  1.00 10.00           C
ATOM      3  C   ALA A   6      23.537   8.349  26.054  1.00 10.00           C
ATOM      4  O   ALA A   6      23.439   7.570  25.057  1.00 10.00           O
ATOM      5  CB  ALA A   6      24.874   7.845  28.114  1.00 10.00           C
ATOM      1  N   ALA A   7      22.542   9.039  26.580  1.00 10.00           N
ATOM      2  CA  ALA A   7      21.228   8.903  25.942  1.00 10.00           C
ATOM      3  C   ALA A   7      21.329   9.698  24.628  1.00 10.00           C
ATOM      4  O   ALA A   7      20.707   9.383  23.632  1.00 10.00           O
ATOM      5  CB  ALA A   7      20.146   9.465  26.862  1.00 10.00           C
ATOM      1  N   ALA A   8      22.181  10.696  24.613  1.00 10.00           N
ATOM      2  CA  ALA A   8      22.526  11.372  23.378  1.00 10.00           C
ATOM      3  C   ALA A   8      23.351  10.555  22.448  1.00 10.00           C
ATOM      4  O   ALA A   8      23.618  10.883  21.252  1.00 10.00           O
ATOM      5  CB  ALA A   8      23.168  12.697  23.693  1.00 10.00           C
ATOM      1  N   ALA A   9      23.864   9.423  22.961  1.00 10.00           N
ATOM      2  CA  ALA A   9      24.785   8.541  22.264  1.00 10.00           C
ATOM      3  C   ALA A   9      24.057   7.451  21.484  1.00 10.00           C
ATOM      4  O   ALA A   9      24.127   7.381  20.257  1.00 10.00           O
ATOM      5  CB  ALA A   9      25.815   7.975  23.249  1.00 10.00           C
ATOM      1  N   ALA A  10      23.518   6.548  22.264  1.00 10.00           N
ATOM      2  CA  ALA A  10      22.629   5.525  21.690  1.00 10.00           C
ATOM      3  C   ALA A  10      21.549   6.308  21.009  1.00 10.00           C
ATOM      4  O   ALA A  10      21.114   5.933  19.930  1.00 10.00           O
ATOM      5  CB  ALA A  10      22.057   4.714  22.784  1.00 10.00           C
ATOM      1  N   ALA A  11      21.120   7.452  21.541  1.00 10.00           N
ATOM      2  CA  ALA A  11      20.186   8.260  20.874  1.00 10.00           C
ATOM      3  C   ALA A  11      20.978   9.215  19.937  1.00 10.00           C
ATOM      4  O   ALA A  11      20.386  10.177  19.507  1.00 10.00           O
ATOM      5  CB  ALA A  11      19.295   9.031  21.867  1.00 10.00           C
ATOM      1  N   ALA A  12      22.222   8.932  19.598  1.00 10.00           N
ATOM      2  CA  ALA A  12      22.896   9.709  18.563  1.00 10.00           C
ATOM      3  C   ALA A  12      22.924   8.925  17.308  1.00 10.00           C
ATOM      4  O   ALA A  12      22.982   9.445  16.193  1.00 10.00           O
ATOM      5  CB  ALA A  12      24.294  10.138  18.994  1.00 10.00           C
ATOM      1  N   ALA A  13      22.951   7.633  17.508  1.00 10.00           N
ATOM      2  CA  ALA A  13      22.709   6.629  16.554  1.00 10.00           C
ATOM      3  C   ALA A  13      21.275   6.673  16.206  1.00 10.00           C
ATOM      4  O   ALA A  13      20.870   6.521  15.092  1.00 10.00           O
ATOM      5  CB  ALA A  13      23.077   5.254  17.025  1.00 10.00           C
ATOM      1  N   ALA A  14      20.471   6.929  17.226  1.00 10.00           N
ATOM      2  CA  ALA A  14      19.039   6.992  17.025  1.00 10.00           C
ATOM      3  C   ALA A  14      18.676   8.380  16.528  1.00 10.00           C
ATOM      4  O   ALA A  14      17.748   8.556  15.761  1.00 10.00           O
ATOM      5  CB  ALA A  14      18.240   6.715  18.272  1.00 10.00           C
ATOM      1  N   ALA A  15      19.381   9.390  17.055  1.00 10.00           N
ATOM      2  CA  ALA A  15      19.204  10.743  16.669  1.00 10.00           C
ATOM      3  C   ALA A  15      19.407  10.807  15.174  1.00 10.00           C
ATOM      4  O   ALA A  15      18.402  10.987  14.424  1.00 10.00           O
ATOM      5  CB  ALA A  15      20.190  11.665  17.493  1.00 10.00           C
ATOM      1  N   ALA A  16      20.702  10.653  14.831  1.00 10.00           N
ATOM      2  CA  ALA A  16      21.206  10.546  13.480  1.00 10.00           C
ATOM      3  C   ALA A  16      20.484   9.612  12.585  1.00 10.00           C
ATOM      4  O   ALA A  16      20.380   9.918  11.386  1.00 10.00           O
ATOM      5  CB  ALA A  16      22.631  10.174  13.475  1.00 10.00           C
ATOM      1  N   ALA A  17      20.064   8.475  13.175  1.00 10.00           N
ATOM      2  CA  ALA A  17      19.355   7.473  12.426  1.00 10.00           C
ATOM      3  C   ALA A  17      17.924   7.807  12.064  1.00 10.00           C
ATOM      4  O   ALA A  17      17.535   7.721  10.871  1.00 10.00           O
ATOM      5  CB  ALA A  17      19.359   6.123  13.216  1.00 10.00           C
ATOM      1  N   ALA A  18      17.152   8.115  13.031  1.00 10.00           N
ATOM      2  CA  ALA A  18      15.835   8.594  12.861  1.00 10.00           C
ATOM      3  C   ALA A  18      15.811   9.835  11.861  1.00 10.00           C
ATOM      4  O   ALA A  18      15.020   9.889  10.868  1.00 10.00           O
ATOM      5  CB  ALA A  18      15.272   8.918  14.234  1.00 10.00           C
ATOM      1  N   ALA A  19      16.661  10.845  12.100  1.00 10.00           N
ATOM      2  CA  ALA A  19      16.435  12.061  11.275  1.00 10.00           C
ATOM      3  C   ALA A  19      17.004  11.815   9.833  1.00 10.00           C
ATOM      4  O   ALA A  19      16.334  12.117   8.857  1.00 10.00           O
ATOM      5  CB  ALA A  19      17.059  13.242  11.866  1.00 10.00           C
ATOM      1  N   ALA A  20      18.191  11.200   9.841  1.00 10.00           N
ATOM      2  CA  ALA A  20      19.091  11.247   8.697  1.00 10.00           C
ATOM      3  C   ALA A  20      19.549   9.835   8.231  1.00 10.00           C
ATOM      4  O   ALA A  20      20.670   9.692   7.663  1.00 10.00           O
ATOM      5  CB  ALA A  20      20.326  12.105   9.035  1.00 10.00           C
ATOM      1  N   ALA A  21      18.654   8.850   8.523  1.00 10.00           N
ATOM      2  CA  ALA A  21      18.827   7.437   8.168  1.00 10.00           C
ATOM      3  C   ALA A  21      17.565   6.607   8.282  1.00 10.00           C
ATOM      4  O   ALA A  21      16.485   6.992   7.820  1.00 10.00           O
ATOM      5  CB  ALA A  21      19.888   6.838   8.983  1.00 10.00           C
"""

def exercise_01(prefix="tst_mi_ext_map_test_01"):
  """
  Simple run to a completion with reference map. no SS annotations.
  """
  pdb_file = open("%s_start.pdb" % prefix, "w")
  pdb_file.write(cryst_str)
  pdb_file.write(pdb_str)
  pdb_file.close()
  # making map
  fmodel_cmd = " ".join([
      "phenix.fmodel",
      "%s_start.pdb" % prefix,
      "high_resolution=5",
      # "generate_fake_p1_symmetry=True",
      ])
  print(fmodel_cmd)
  assert not easy_run.call(fmodel_cmd)
  mtz2map_cmd = " ".join([
    "phenix.mtz2map",
    "%s_start.pdb.mtz" % prefix,
    "include_fmodel=True"])
  print(mtz2map_cmd)
  assert not easy_run.call(mtz2map_cmd)
  # STOP()
  cmd = " ".join([
      "phenix.model_idealization",
      "%s_start.pdb" % prefix,
      "%s_start.pdb_fmodel.ccp4" % prefix,
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
      "Preparing user map...\n",
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
