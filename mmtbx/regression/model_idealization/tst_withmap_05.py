from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
import os.path
import time

# taken from phenix_regression/refinement/ncs/tst_ncs_0.py
pdb_str = """\
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
ATOM   1395  N   GLY A  95      55.545  77.050  62.972  1.00 67.65           N
ATOM   1396  CA  GLY A  95      54.815  76.493  64.067  1.00 66.50           C
ATOM   1397  C   GLY A  95      55.509  76.948  65.306  1.00 65.94           C
ATOM   1398  O   GLY A  95      56.728  77.057  65.338  1.00 66.27           O
ATOM   1399  N   LYS A  96      54.733  77.220  66.367  1.00 65.59           N
ATOM   1400  CA  LYS A  96      55.263  77.542  67.659  1.00 67.02           C
ATOM   1401  C   LYS A  96      55.940  76.308  68.162  1.00 65.59           C
ATOM   1402  O   LYS A  96      56.897  76.361  68.935  1.00 65.24           O
ATOM   1403  CB  LYS A  96      54.152  77.941  68.646  1.00 69.25           C
ATOM   1404  CG  LYS A  96      53.150  78.931  68.042  1.00 68.89           C
ATOM   1405  CD  LYS A  96      51.947  79.224  68.942  1.00 67.92           C
ATOM   1406  CE  LYS A  96      50.602  79.121  68.220  1.00 71.86           C
ATOM   1407  NZ  LYS A  96      50.769  78.385  66.946  1.00 70.05           N
ATOM   1408  N   GLU A  97      55.411  75.151  67.727  1.00 63.83           N
ATOM   1409  CA  GLU A  97      55.703  73.848  68.246  1.00 62.40           C
ATOM   1410  C   GLU A  97      57.109  73.451  67.909  1.00 56.41           C
ATOM   1411  O   GLU A  97      57.815  74.145  67.178  1.00 54.82           O
ATOM   1412  CB  GLU A  97      54.732  72.787  67.705  1.00 69.38           C
ATOM   1413  CG  GLU A  97      53.280  73.274  67.765  1.00 77.38           C
ATOM   1414  CD  GLU A  97      52.416  72.412  66.855  1.00 83.20           C
ATOM   1415  OE1 GLU A  97      52.933  71.941  65.807  1.00 84.20           O
ATOM   1416  OE2 GLU A  97      51.218  72.220  67.198  1.00 84.08           O
ATOM   1417  N   ASP A  98      57.560  72.337  68.527  1.00 54.83           N
ATOM   1418  CA  ASP A  98      58.945  71.963  68.620  1.00 52.52           C
ATOM   1419  C   ASP A  98      59.391  71.254  67.375  1.00 45.70           C
ATOM   1420  O   ASP A  98      58.695  70.397  66.837  1.00 45.06           O
ATOM   1421  CB  ASP A  98      59.220  71.021  69.805  1.00 61.00           C
ATOM   1422  CG  ASP A  98      58.309  71.431  70.957  1.00 64.25           C
ATOM   1423  OD1 ASP A  98      57.123  71.006  70.957  1.00 65.51           O
ATOM   1424  OD2 ASP A  98      58.789  72.174  71.855  1.00 63.23           O
ATOM   1425  N   ALA A  99      60.588  71.627  66.874  1.00 40.59           N
ATOM   1426  CA  ALA A  99      61.188  71.036  65.709  1.00 34.76           C
ATOM   1427  C   ALA A  99      61.566  69.603  65.956  1.00 31.73           C
ATOM   1428  O   ALA A  99      61.242  68.714  65.169  1.00 24.61           O
ATOM   1429  CB  ALA A  99      62.480  71.767  65.318  1.00 35.50           C
ATOM   1430  N   ALA A 100      62.263  69.356  67.081  1.00 29.14           N
ATOM   1431  CA  ALA A 100      62.905  68.103  67.363  1.00 26.25           C
ATOM   1432  C   ALA A 100      63.894  67.769  66.291  1.00 24.22           C
ATOM   1433  O   ALA A 100      63.646  67.922  65.096  1.00 25.43           O
ATOM   1434  CB  ALA A 100      61.940  66.914  67.494  1.00 30.57           C
ATOM   1435  N   ASN A 101      65.061  67.266  66.721  1.00 24.03           N
ATOM   1436  CA  ASN A 101      66.043  66.814  65.791  1.00 25.16           C
ATOM   1437  C   ASN A 101      65.760  65.366  65.582  1.00 22.14           C
ATOM   1438  O   ASN A 101      66.599  64.515  65.877  1.00 22.24           O
ATOM   1439  CB  ASN A 101      67.447  66.846  66.398  1.00 30.60           C
ATOM   1440  CG  ASN A 101      68.318  67.948  65.831  1.00 33.80           C
ATOM   1441  OD1 ASN A 101      69.536  67.907  65.989  1.00 37.78           O
ATOM   1442  ND2 ASN A 101      67.719  68.969  65.169  1.00 34.81           N
ATOM   1443  N   ASN A 102      64.577  65.034  65.045  1.00 18.33           N
ATOM   1444  CA  ASN A 102      64.423  63.668  64.678  1.00 15.90           C
ATOM   1445  C   ASN A 102      63.381  63.583  63.621  1.00 14.00           C
ATOM   1446  O   ASN A 102      62.490  64.426  63.515  1.00 12.14           O
ATOM   1447  CB  ASN A 102      64.070  62.725  65.846  1.00 18.46           C
ATOM   1448  CG  ASN A 102      62.726  63.100  66.443  1.00 21.95           C
ATOM   1449  OD1 ASN A 102      62.629  63.981  67.295  1.00 27.01           O
ATOM   1450  ND2 ASN A 102      61.656  62.391  65.994  1.00 27.14           N
ATOM   1451  N   TYR A 103      63.516  62.534  62.799  1.00 10.81           N
ATOM   1452  CA  TYR A 103      62.708  62.269  61.655  1.00  6.91           C
ATOM   1453  C   TYR A 103      61.309  62.060  62.142  1.00  2.00           C
ATOM   1454  O   TYR A 103      60.349  62.519  61.526  1.00  2.00           O
ATOM   1455  CB  TYR A 103      63.217  60.995  60.944  1.00  3.80           C
ATOM   1456  CG  TYR A 103      62.318  60.589  59.826  1.00  9.36           C
ATOM   1457  CD1 TYR A 103      62.305  61.252  58.626  1.00  9.82           C
ATOM   1458  CD2 TYR A 103      61.478  59.505  59.974  1.00  8.69           C
ATOM   1459  CE1 TYR A 103      61.465  60.854  57.606  1.00  8.11           C
ATOM   1460  CE2 TYR A 103      60.651  59.094  58.966  1.00 13.62           C
ATOM   1461  CZ  TYR A 103      60.638  59.771  57.765  1.00 13.62           C
ATOM   1462  OH  TYR A 103      59.782  59.351  56.728  1.00 17.18           O
ATOM   1463  N   ALA A 104      61.153  61.334  63.263  1.00  4.07           N
ATOM   1464  CA  ALA A 104      59.857  60.824  63.598  1.00  6.65           C
ATOM   1465  C   ALA A 104      58.853  61.920  63.825  1.00  7.19           C
ATOM   1466  O   ALA A 104      57.788  61.896  63.211  1.00  6.00           O
ATOM   1467  CB  ALA A 104      59.885  59.966  64.870  1.00 11.87           C
ATOM   1468  N   ARG A 105      59.131  62.918  64.696  1.00  7.49           N
ATOM   1469  CA  ARG A 105      58.072  63.870  64.909  1.00 10.39           C
ATOM   1470  C   ARG A 105      57.876  64.604  63.636  1.00 11.29           C
ATOM   1471  O   ARG A 105      56.754  64.797  63.170  1.00 11.85           O
ATOM   1472  CB  ARG A 105      58.328  64.974  65.957  1.00 11.69           C
ATOM   1473  CG  ARG A 105      59.547  64.795  66.856  1.00 12.68           C
ATOM   1474  CD  ARG A 105      59.356  65.416  68.245  1.00 13.62           C
ATOM   1475  NE  ARG A 105      58.864  66.817  68.080  1.00 14.65           N
ATOM   1476  CZ  ARG A 105      57.573  67.139  68.396  1.00 16.77           C
ATOM   1477  NH1 ARG A 105      56.696  66.163  68.773  1.00 18.53           N
ATOM   1478  NH2 ARG A 105      57.160  68.439  68.341  1.00 18.02           N
END
"""

def exercise_04(prefix="tst_mi_map_test_05"):
  """
  Actually will try to CCD.
  Should include actual check for success.
  """
  # without cryst
  pdb_file = open("%s_start.pdb" % prefix, "w")
  pdb_file.write(pdb_str)
  pdb_file.close()
  cmd = " ".join([
      "phenix.model_idealization",
      "%s_start.pdb" % prefix,
      "use_map_for_reference=True",
      "run_minimization_first=False",
      "run_minimization_last=False",
      "loop_idealization.minimize_whole=False",
      "loop_idealization.number_of_ccd_trials=1",
      "number_of_refinement_cycles=1",
      "n_macro=1",
      "debug=True",
      # ">%s.log" % prefix,
      ])
  print(cmd)
  assert not easy_run.call(cmd)
  # assert os.path.isfile("%s_start.pdb_all_idealized.pdb" % prefix)
  # res_log = open("%s.log" % prefix, "r")
  # log_lines = res_log.readlines()
  # # NCS constraints with map are not implemented yet
  # for l in [
  #     # "Using ncs\n",
  #     "Using map as reference\n",
  #     "  Minimizing... (NCS)\n",
  #     # "Ramachandran outliers:      0.00      0.00      0.00      0.00      0.00\n",
  #     "All done.\n"]:
  #   assert l in log_lines, "'%s' not in log file." % l
  # res_log.close()

if (__name__ == "__main__"):
  t0 = time.time()
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping: probe not configured")
  else:
    exercise_04()
  print("Time: %.2f" % (time.time() - t0))
  print("OK")
