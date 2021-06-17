from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
import os.path
import time

# taken from 3zee multiplied
pdb_str = """\
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
ATOM    114  N   VALA1  13      28.866   7.875-103.684  1.00  0.00           N
ATOM    115  CA  VALA1  13      27.966   7.696-104.717  1.00  0.00           C
ATOM    116  C   VALA1  13      28.923   7.301-105.774  1.00  0.00           C
ATOM    117  O   VALA1  13      29.996   7.851-105.765  1.00  0.00           O
ATOM    118  CB  VALA1  13      27.357   9.044-105.115  1.00  0.00           C
ATOM    119  CG1 VALA1  13      28.173  10.400-104.901  1.00  0.00           C
ATOM    120  CG2 VALA1  13      26.772   9.026-106.573  1.00  0.00           C
ATOM    121  N   VALA1  14      28.597   6.282-106.581  1.00  0.00           N
ATOM    122  CA  VALA1  14      29.479   5.817-107.556  1.00  0.00           C
ATOM    123  C   VALA1  14      28.731   5.902-108.820  1.00  0.00           C
ATOM    124  O   VALA1  14      27.735   5.291-108.935  1.00  0.00           O
ATOM    125  CB  VALA1  14      29.816   4.426-107.266  1.00  0.00           C
ATOM    126  CG1 VALA1  14      30.981   4.043-108.197  1.00  0.00           C
ATOM    127  CG2 VALA1  14      30.218   4.344-105.790  1.00  0.00           C
ATOM    128  N   PROA1  15      29.286   6.585-109.794  1.00  0.00           N
ATOM    129  CA  PROA1  15      28.744   6.599-111.104  1.00  0.00           C
ATOM    130  C   PROA1  15      29.165   5.280-111.724  1.00  0.00           C
ATOM    131  O   PROA1  15      30.277   4.833-111.335  1.00  0.00           O
ATOM    132  CB  PROA1  15      29.324   7.870-111.770  1.00  0.00           C
ATOM    133  CG  PROA1  15      30.683   7.956-111.088  1.00  0.00           C
ATOM    134  CD  PROA1  15      30.341   7.551-109.692  1.00  0.00           C
ATOM    135  N   CYSA1  16      28.318   4.535-112.451  1.00  0.00           N
ATOM    136  CA  CYSA1  16      28.758   3.175-112.680  1.00  0.00           C
ATOM    137  C   CYSA1  16      28.739   2.854-114.087  1.00  0.00           C
ATOM    138  O   CYSA1  16      28.894   1.759-114.604  1.00  0.00           O
ATOM    139  CB  CYSA1  16      27.905   2.251-111.855  1.00  0.00           C
ATOM    140  SG  CYSA1  16      28.255   2.522-110.159  1.00  0.00           S
ATOM    141  N   GLYA1  17      28.299   3.834-114.840  1.00  0.00           N
ATOM    142  CA  GLYA1  17      27.972   3.877-116.194  1.00  0.00           C
ATOM    143  C   GLYA1  17      27.015   2.770-116.530  1.00  0.00           C
ATOM    144  O   GLYA1  17      25.791   3.024-116.400  1.00  0.00           O
ATOM    145  N   ASPA1  18      27.485   1.673-117.088  1.00  0.00           N
ATOM    146  CA  ASPA1  18      26.581   0.717-117.651  1.00  0.00           C
ATOM    147  C   ASPA1  18      26.603  -0.454-116.677  1.00  0.00           C
ATOM    148  O   ASPA1  18      25.726  -1.309-116.715  1.00  0.00           O
ATOM    149  CB  ASPA1  18      27.058   0.127-118.978  1.00  0.00           C
ATOM    150  CG  ASPA1  18      25.941  -0.447-119.867  1.00  0.00           C
ATOM    151  OD1 ASPA1  18      25.283  -1.471-119.437  1.00  0.00           O
ATOM    152  OD2 ASPA1  18      25.763   0.055-120.988  1.00  0.00           O
ATOM    153  N   GLYA1  19      27.598  -0.551-115.748  1.00  0.00           N
ATOM    154  CA  GLYA1  19      27.850  -1.623-114.809  1.00  0.00           C
ATOM    155  C   GLYA1  19      28.483  -2.857-115.392  1.00  0.00           C
ATOM    156  O   GLYA1  19      28.602  -3.897-114.714  1.00  0.00           O
ATOM    157  N   ARGA1  20      28.884  -2.776-116.690  1.00  0.00           N
ATOM    158  CA  ARGA1  20      29.503  -3.853-117.351  1.00  0.00           C
ATOM    159  C   ARGA1  20      30.956  -3.624-117.090  1.00  0.00           C
ATOM    160  O   ARGA1  20      31.604  -2.729-117.659  1.00  0.00           O
ATOM    161  CB  ARGA1  20      29.135  -3.792-118.837  1.00  0.00           C
ATOM    162  CG  ARGA1  20      27.632  -4.042-119.004  1.00  0.00           C
ATOM    163  CD  ARGA1  20      27.125  -5.547-118.721  1.00  0.00           C
ATOM    164  NE  ARGA1  20      27.534  -6.358-119.907  1.00  0.00           N
ATOM    165  CZ  ARGA1  20      28.739  -6.931-119.948  1.00  0.00           C
ATOM    166  NH1 ARGA1  20      29.387  -7.363-118.805  1.00  0.00           N
ATOM    167  NH2 ARGA1  20      29.293  -7.238-121.145  1.00  0.00           N
ATOM    168  N   META1  21      31.444  -4.500-116.136  1.00  0.00           N
ATOM    169  CA  META1  21      32.768  -4.515-115.553  1.00  0.00           C
ATOM    170  C   META1  21      32.830  -5.593-114.505  1.00  0.00           C
ATOM    171  O   META1  21      31.843  -5.906-113.868  1.00  0.00           O
ATOM    172  CB  META1  21      32.981  -3.224-114.734  1.00  0.00           C
ATOM    173  CG  META1  21      34.438  -2.959-114.287  1.00  0.00           C
ATOM    174  SD  META1  21      34.506  -2.151-112.763  1.00  0.00           S
ATOM    175  CE  META1  21      34.190  -0.475-113.397  1.00  0.00           C
ATOM    176  N   LYSA1  22      34.049  -6.099-114.200  1.00  0.00           N
ATOM    177  CA  LYSA1  22      34.284  -6.987-113.056  1.00  0.00           C
ATOM    178  C   LYSA1  22      34.369  -6.049-111.900  1.00  0.00           C
ATOM    179  O   LYSA1  22      34.435  -4.865-112.055  1.00  0.00           O
ATOM    180  CB  LYSA1  22      35.583  -7.836-113.206  1.00  0.00           C
ATOM    181  CG  LYSA1  22      35.330  -9.254-113.722  1.00  0.00           C
ATOM    182  CD  LYSA1  22      36.457 -10.351-113.280  1.00  0.00           C
ATOM    183  CE  LYSA1  22      35.864 -11.774-113.255  1.00  0.00           C
ATOM    184  NZ  LYSA1  22      36.511 -12.735-112.250  1.00  0.00           N
TER
ATOM    784  N   VALA2  13      26.276 -14.311-100.152  1.00  0.00           N
ATOM    785  CA  VALA2  13      25.503 -13.817-101.185  1.00  0.00           C
ATOM    786  C   VALA2  13      25.920 -14.765-102.242  1.00  0.00           C
ATOM    787  O   VALA2  13      27.074 -15.112-102.233  1.00  0.00           O
ATOM    788  CB  VALA2  13      25.997 -12.423-101.583  1.00  0.00           C
ATOM    789  CG1 VALA2  13      27.525 -12.011-101.369  1.00  0.00           C
ATOM    790  CG2 VALA2  13      25.563 -12.031-103.041  1.00  0.00           C
ATOM    791  N   VALA2  14      24.979 -15.275-103.049  1.00  0.00           N
ATOM    792  CA  VALA2  14      25.293 -16.221-104.024  1.00  0.00           C
ATOM    793  C   VALA2  14      24.813 -15.641-105.288  1.00  0.00           C
ATOM    794  O   VALA2  14      23.671 -15.392-105.403  1.00  0.00           O
ATOM    795  CB  VALA2  14      24.573 -17.458-103.734  1.00  0.00           C
ATOM    796  CG1 VALA2  14      25.148 -18.540-104.665  1.00  0.00           C
ATOM    797  CG2 VALA2  14      24.806 -17.795-102.258  1.00  0.00           C
ATOM    798  N   PROA2  15      25.685 -15.533-106.262  1.00  0.00           N
ATOM    799  CA  PROA2  15      25.304 -15.148-107.572  1.00  0.00           C
ATOM    800  C   PROA2  15      24.695 -16.390-108.192  1.00  0.00           C
ATOM    801  O   PROA2  15      25.187 -17.483-107.803  1.00  0.00           O
ATOM    802  CB  PROA2  15      26.603 -14.632-108.238  1.00  0.00           C
ATOM    803  CG  PROA2  15      27.643 -15.511-107.556  1.00  0.00           C
ATOM    804  CD  PROA2  15      27.116 -15.567-106.160  1.00  0.00           C
ATOM    805  N   CYSA2  16      23.568 -16.341-108.919  1.00  0.00           N
ATOM    806  CA  CYSA2  16      22.944 -17.627-109.148  1.00  0.00           C
ATOM    807  C   CYSA2  16      22.707 -17.845-110.555  1.00  0.00           C
ATOM    808  O   CYSA2  16      22.061 -18.743-111.072  1.00  0.00           O
ATOM    809  CB  CYSA2  16      21.688 -17.702-108.323  1.00  0.00           C
ATOM    810  SG  CYSA2  16      22.128 -17.749-106.627  1.00  0.00           S
ATOM    811  N   GLYA2  17      23.069 -16.834-111.308  1.00  0.00           N
ATOM    812  CA  GLYA2  17      22.862 -16.576-112.662  1.00  0.00           C
ATOM    813  C   GLYA2  17      21.406 -16.712-112.998  1.00  0.00           C
ATOM    814  O   GLYA2  17      20.698 -15.681-112.868  1.00  0.00           O
ATOM    815  N   ASPA2  18      20.985 -17.829-113.556  1.00  0.00           N
ATOM    816  CA  ASPA2  18      19.670 -17.892-114.119  1.00  0.00           C
ATOM    817  C   ASPA2  18      18.875 -18.753-113.145  1.00  0.00           C
ATOM    818  O   ASPA2  18      17.651 -18.762-113.183  1.00  0.00           O
ATOM    819  CB  ASPA2  18      19.605 -18.648-115.446  1.00  0.00           C
ATOM    820  CG  ASPA2  18      18.403 -18.289-116.335  1.00  0.00           C
ATOM    821  OD1 ASPA2  18      17.219 -18.571-115.905  1.00  0.00           O
ATOM    822  OD2 ASPA2  18      18.622 -17.803-117.456  1.00  0.00           O
ATOM    823  N   GLYA2  19      19.526 -19.511-112.216  1.00  0.00           N
ATOM    824  CA  GLYA2  19      18.965 -20.459-111.277  1.00  0.00           C
ATOM    825  C   GLYA2  19      18.567 -21.788-111.860  1.00  0.00           C
ATOM    826  O   GLYA2  19      17.933 -22.621-111.182  1.00  0.00           O
ATOM    827  N   ARGA2  20      18.912 -22.007-113.158  1.00  0.00           N
ATOM    828  CA  ARGA2  20      18.614 -23.213-113.819  1.00  0.00           C
ATOM    829  C   ARGA2  20      19.820 -24.054-113.558  1.00  0.00           C
ATOM    830  O   ARGA2  20      20.907 -23.857-114.127  1.00  0.00           O
ATOM    831  CB  ARGA2  20      18.390 -22.914-115.305  1.00  0.00           C
ATOM    832  CG  ARGA2  20      17.132 -22.053-115.472  1.00  0.00           C
ATOM    833  CD  ARGA2  20      15.724 -22.788-115.189  1.00  0.00           C
ATOM    834  NE  ARGA2  20      15.458 -23.656-116.375  1.00  0.00           N
ATOM    835  CZ  ARGA2  20      15.930 -24.904-116.416  1.00  0.00           C
ATOM    836  NH1 ARGA2  20      16.099 -25.664-115.273  1.00  0.00           N
ATOM    837  NH2 ARGA2  20      16.117 -25.509-117.613  1.00  0.00           N
ATOM    838  N   META2  21      19.565 -25.024-112.604  1.00  0.00           N
ATOM    839  CA  META2  21      20.509 -25.952-112.021  1.00  0.00           C
ATOM    840  C   META2  21      19.807 -26.772-110.973  1.00  0.00           C
ATOM    841  O   META2  21      18.879 -26.314-110.336  1.00  0.00           O
ATOM    842  CB  META2  21      21.558 -25.168-111.202  1.00  0.00           C
ATOM    843  CG  META2  21      22.792 -25.985-110.755  1.00  0.00           C
ATOM    844  SD  META2  21      23.400 -25.450-109.231  1.00  0.00           S
ATOM    845  CE  META2  21      24.334 -24.022-109.865  1.00  0.00           C
ATOM    846  N   LYSA2  22      20.337 -27.981-110.668  1.00  0.00           N
ATOM    847  CA  LYSA2  22      19.891 -28.784-109.524  1.00  0.00           C
ATOM    848  C   LYSA2  22      20.602 -28.167-108.368  1.00  0.00           C
ATOM    849  O   LYSA2  22      21.470 -27.358-108.523  1.00  0.00           O
ATOM    850  CB  LYSA2  22      20.240 -30.297-109.674  1.00  0.00           C
ATOM    851  CG  LYSA2  22      19.076 -31.144-110.190  1.00  0.00           C
ATOM    852  CD  LYSA2  22      19.129 -32.716-109.748  1.00  0.00           C
ATOM    853  CE  LYSA2  22      17.715 -33.332-109.723  1.00  0.00           C
ATOM    854  NZ  LYSA2  22      17.517 -34.473-108.718  1.00  0.00           N
TER
END
"""

def exercise_04(prefix="tst_mi_map_test_06"):
  """
  Actually will try to CCD. Same as 05, but with NCS copies.
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
  res_log = open("%s_start.pdb.log" % prefix, "r")
  log_lines = res_log.readlines()
  res_log.close()
  # NCS constraints with map are not implemented yet
  # should be there
  for l in [
      "Using ncs\n",
      "Idealizing chain A1\n",
      # 'Working on pdbres="VALA1  13 " OUTLIER\n',
      # "Using map as reference\n",
      # "  Minimizing... (NCS)\n",
      # "Ramachandran outliers:      0.00      0.00      0.00      0.00      0.00\n",
      # "Rotamer outliers     :     12.50     12.50     12.50     12.50      0.00      0.00\n",
      "All done.\n"]:
    assert l in log_lines, "'%s' not in log file." % l
  # should not be there
  for l in [
      "Idealizing chain A2\n",
      'Working on pdbres="VALA2  13 " OUTLIER\n']:
    assert l not in log_lines, "'%s' should not be in log file." % l

if (__name__ == "__main__"):
  t0 = time.time()
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping: probe not configured")
  else:
    exercise_04()
  print("Time: %.2f" % (time.time() - t0))
  print("OK")
