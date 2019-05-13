from __future__ import division
from __future__ import print_function
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
import time
import mmtbx.nci.hbond
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

pdb_str_00 = """
ATOM      1  N   GLY A   1      -5.606  -2.251 -12.878  1.00  0.00           N
ATOM      2  CA  GLY A   1      -5.850  -1.194 -13.852  1.00  0.00           C
ATOM      3  C   GLY A   1      -5.186  -1.524 -15.184  1.00  0.00           C
ATOM      4  O   GLY A   1      -5.744  -1.260 -16.249  1.00  0.00           O
ATOM      0  H1  GLY A   1      -6.104  -2.109 -12.154  1.00  0.00           H   new
ATOM      0  H2  GLY A   1      -5.819  -3.038 -13.234  1.00  0.00           H   new
ATOM      0  H3  GLY A   1      -4.746  -2.252 -12.650  1.00  0.00           H   new
ATOM      0  HA2 GLY A   1      -6.805  -1.081 -13.981  1.00  0.00           H   new
ATOM      0  HA3 GLY A   1      -5.507  -0.352 -13.515  1.00  0.00           H   new
ATOM      1  N   GLY A   2      -3.992  -2.102 -15.115  1.00  0.00           N
ATOM      2  CA  GLY A   2      -3.261  -2.499 -16.313  1.00  0.00           C
ATOM      3  C   GLY A   2      -3.961  -3.660 -17.011  1.00  0.00           C
ATOM      4  O   GLY A   2      -4.016  -3.716 -18.240  1.00  0.00           O
ATOM      0  H   GLY A   2      -3.585  -2.274 -14.378  1.00  0.00           H   new
ATOM      0  HA2 GLY A   2      -3.191  -1.745 -16.920  1.00  0.00           H   new
ATOM      0  HA3 GLY A   2      -2.356  -2.756 -16.075  1.00  0.00           H   new
ATOM      1  N   GLY A   3      -4.492  -4.585 -16.219  1.00  0.00           N
ATOM      2  CA  GLY A   3      -5.216  -5.731 -16.755  1.00  0.00           C
ATOM      3  C   GLY A   3      -6.531  -5.289 -17.389  1.00  0.00           C
ATOM      4  O   GLY A   3      -6.939  -5.814 -18.425  1.00  0.00           O
ATOM      0  H   GLY A   3      -4.443  -4.566 -15.361  1.00  0.00           H   new
ATOM      0  HA2 GLY A   3      -4.669  -6.185 -17.416  1.00  0.00           H   new
ATOM      0  HA3 GLY A   3      -5.392  -6.369 -16.046  1.00  0.00           H   new
ATOM      1  N   GLY A   4      -7.189  -4.323 -16.758  1.00  0.00           N
ATOM      2  CA  GLY A   4      -8.442  -3.785 -17.273  1.00  0.00           C
ATOM      3  C   GLY A   4      -8.205  -3.003 -18.561  1.00  0.00           C
ATOM      4  O   GLY A   4      -9.007  -3.065 -19.492  1.00  0.00           O
ATOM      0  H   GLY A   4      -6.924  -3.963 -16.024  1.00  0.00           H   new
ATOM      0  HA2 GLY A   4      -9.066  -4.509 -17.439  1.00  0.00           H   new
ATOM      0  HA3 GLY A   4      -8.848  -3.207 -16.608  1.00  0.00           H   new
ATOM      1  N   GLY A   5      -7.099  -2.269 -18.604  1.00  0.00           N
ATOM      2  CA  GLY A   5      -6.735  -1.498 -19.787  1.00  0.00           C
ATOM      3  C   GLY A   5      -6.358  -2.423 -20.939  1.00  0.00           C
ATOM      4  O   GLY A   5      -6.687  -2.157 -22.094  1.00  0.00           O
ATOM      0  H   GLY A   5      -6.541  -2.204 -17.953  1.00  0.00           H   new
ATOM      0  HA2 GLY A   5      -7.477  -0.932 -20.051  1.00  0.00           H   new
ATOM      0  HA3 GLY A   5      -5.990  -0.912 -19.580  1.00  0.00           H   new
ATOM      1  N   GLY A   6      -5.665  -3.509 -20.614  1.00  0.00           N
ATOM      2  CA  GLY A   6      -5.268  -4.493 -21.614  1.00  0.00           C
ATOM      3  C   GLY A   6      -6.485  -5.236 -22.153  1.00  0.00           C
ATOM      4  O   GLY A   6      -6.565  -5.533 -23.345  1.00  0.00           O
ATOM      0  H   GLY A   6      -5.413  -3.695 -19.813  1.00  0.00           H   new
ATOM      0  HA2 GLY A   6      -4.804  -4.051 -22.343  1.00  0.00           H   new
ATOM      0  HA3 GLY A   6      -4.645  -5.125 -21.223  1.00  0.00           H   new
ATOM      1  N   GLY A   7      -7.430  -5.532 -21.267  1.00  0.00           N
ATOM      2  CA  GLY A   7      -8.660  -6.212 -21.655  1.00  0.00           C
ATOM      3  C   GLY A   7      -9.529  -5.303 -22.518  1.00  0.00           C
ATOM      4  O   GLY A   7     -10.158  -5.756 -23.474  1.00  0.00           O
ATOM      0  H   GLY A   7      -7.377  -5.346 -20.429  1.00  0.00           H   new
ATOM      0  HA2 GLY A   7      -8.447  -7.022 -22.144  1.00  0.00           H   new
ATOM      0  HA3 GLY A   7      -9.151  -6.479 -20.862  1.00  0.00           H   new
ATOM      1  N   GLY A   8      -9.559  -4.021 -22.172  1.00  0.00           N
ATOM      2  CA  GLY A   8     -10.324  -3.039 -22.930  1.00  0.00           C
ATOM      3  C   GLY A   8      -9.706  -2.819 -24.306  1.00  0.00           C
ATOM      4  O   GLY A   8     -10.416  -2.660 -25.299  1.00  0.00           O
ATOM      0  H   GLY A   8      -9.139  -3.697 -21.495  1.00  0.00           H   new
ATOM      0  HA2 GLY A   8     -11.241  -3.341 -23.027  1.00  0.00           H   new
ATOM      0  HA3 GLY A   8     -10.352  -2.199 -22.445  1.00  0.00           H   new
ATOM      1  N   GLY A   9      -8.378  -2.810 -24.356  1.00  0.00           N
ATOM      2  CA  GLY A   9      -7.658  -2.641 -25.613  1.00  0.00           C
ATOM      3  C   GLY A   9      -7.843  -3.861 -26.508  1.00  0.00           C
ATOM      4  O   GLY A   9      -7.980  -3.734 -27.725  1.00  0.00           O
ATOM      0  H   GLY A   9      -7.872  -2.901 -23.667  1.00  0.00           H   new
ATOM      0  HA2 GLY A   9      -7.978  -1.847 -26.070  1.00  0.00           H   new
ATOM      0  HA3 GLY A   9      -6.714  -2.505 -25.435  1.00  0.00           H   new
"""

pdb_str_01 = """
CRYST1   44.527   44.527   62.306  90.00  90.00 120.00 P 3 2 1
ATOM    484  CB  ASP A 729      24.227   3.198  28.416  1.00 40.30           C
ATOM    485  CG  ASP A 729      23.109   3.028  27.407  1.00 47.94           C
ATOM    486  OD1 ASP A 729      21.948   3.391  27.749  1.00 53.20           O
ATOM    487  OD2 ASP A 729      23.414   2.492  26.298  1.00 59.90           O
ATOM    658  CZ  ARG A 740      20.407  -3.170  38.924  1.00 70.50           C
ATOM    660  NH2 ARG A 740      20.911  -2.860  37.739  1.00 64.32           N
ATOM    672 HH21 ARG A 740      21.872  -2.629  37.580  1.00 64.32           H
ATOM    673 HH22 ARG A 740      20.350  -2.943  36.924  1.00 64.32           H
TER
END
"""

pdb_str_02 = """
CRYST1   54.464   54.312   58.062  90.00  90.00  90.00 P 1
ATOM      8  NE  ARG A  96      37.134  48.090   7.933  1.00 51.90           N
ATOM      9  CZ  ARG A  96      36.368  47.222   8.585  1.00 55.75           C
ATOM     10  NH1 ARG A  96      35.225  46.792   8.057  1.00 55.83           N
ATOM     11  NH2 ARG A  96      36.755  46.779   9.765  1.00 46.04           N
ATOM     20  HE  ARG A  96      37.867  48.340   8.306  1.00 51.90           H
ATOM     21 HH11 ARG A  96      34.974  47.076   7.285  1.00 55.82           H
ATOM     22 HH12 ARG A  96      34.737  46.230   8.488  1.00 55.82           H
ATOM     23 HH21 ARG A  96      37.496  47.052  10.105  1.00 46.04           H
ATOM     24 HH22 ARG A  96      36.266  46.217  10.195  1.00 46.04           H
ATOM     44  CB  ASP A 118      34.771  43.337  12.171  1.00 50.91           C
ATOM     45  CG  ASP A 118      34.959  44.392  11.096  1.00 59.09           C
ATOM     46  OD1 ASP A 118      34.559  44.164   9.925  1.00 60.01           O
ATOM     47  OD2 ASP A 118      35.520  45.449  11.443  1.00 59.21           O
ATOM     50  HB2 ASP A 118      34.866  43.750  13.043  1.00 50.91           H
ATOM     51  HB3 ASP A 118      35.474  42.674  12.093  1.00 50.91           H
TER
END
"""

pdb_str_03 = """
ATOM      1  N   GLY A   1      -5.423  -6.897  -8.991  1.00  0.00           N
ATOM      2  CA  GLY A   1      -3.976  -6.554  -9.489  1.00  0.00           C
ATOM      3  C   GLY A   1      -3.642  -5.220 -10.016  1.00  0.00           C
ATOM      4  O   GLY A   1      -2.895  -4.472  -9.450  1.00  0.00           O
ATOM      5  H1  GLY A   1      -6.043  -6.804  -9.631  1.00  0.00           H
ATOM      6  H2  GLY A   1      -5.535  -7.752  -8.722  1.00  0.00           H
ATOM      7  H3  GLY A   1      -5.708  -6.370  -8.334  1.00  0.00           H
ATOM      8  HA2 GLY A   1      -3.306  -6.803  -8.861  1.00  0.00           H
ATOM      9  HA3 GLY A   1      -3.777  -7.117 -10.209  1.00  0.00           H
ATOM     10  N   GLY A   2      -4.208  -4.919 -11.204  1.00  0.00           N
ATOM     11  CA  GLY A   2      -3.776  -3.707 -11.889  1.00  0.00           C
ATOM     12  C   GLY A   2      -4.936  -2.864 -12.308  1.00  0.00           C
ATOM     13  O   GLY A   2      -5.832  -2.544 -11.494  1.00  0.00           O
ATOM     14  H   GLY A   2      -4.692  -5.484 -11.621  1.00  0.00           H
ATOM     15  HA2 GLY A   2      -3.186  -3.174 -11.317  1.00  0.00           H
ATOM     16  HA3 GLY A   2      -3.316  -3.955 -12.670  1.00  0.00           H
ATOM     17  N   GLY A   3      -4.859  -2.354 -13.528  1.00  0.00           N
ATOM     18  CA  GLY A   3      -5.899  -1.438 -13.954  1.00  0.00           C
ATOM     19  C   GLY A   3      -6.867  -2.105 -14.900  1.00  0.00           C
ATOM     20  O   GLY A   3      -7.641  -3.024 -14.476  1.00  0.00           O
ATOM     21  H   GLY A   3      -4.224  -2.449 -14.045  1.00  0.00           H
ATOM     22  HA2 GLY A   3      -6.449  -1.066 -13.206  1.00  0.00           H
ATOM     23  HA3 GLY A   3      -5.590  -0.635 -14.420  1.00  0.00           H
ATOM     24  N   GLY A   4      -6.766  -1.789 -16.194  1.00  0.00           N
ATOM     25  CA  GLY A   4      -7.861  -2.087 -17.072  1.00  0.00           C
ATOM     26  C   GLY A   4      -7.699  -1.579 -18.516  1.00  0.00           C
ATOM     27  O   GLY A   4      -7.751  -0.367 -18.770  1.00  0.00           O
ATOM     28  H   GLY A   4      -6.135  -1.331 -16.478  1.00  0.00           H
ATOM     29  HA2 GLY A   4      -7.943  -3.006 -16.969  1.00  0.00           H
ATOM     30  HA3 GLY A   4      -8.656  -1.754 -16.709  1.00  0.00           H
ATOM     31  N   GLY A   5      -7.809  -2.518 -19.502  1.00  0.00           N
ATOM     32  CA  GLY A   5      -7.740  -2.130 -20.882  1.00  0.00           C
ATOM     33  C   GLY A   5      -6.478  -2.653 -21.539  1.00  0.00           C
ATOM     34  O   GLY A   5      -5.802  -3.533 -20.997  1.00  0.00           O
ATOM     35  H   GLY A   5      -7.838  -3.379 -19.455  1.00  0.00           H
ATOM     36  HA2 GLY A   5      -8.468  -2.578 -21.312  1.00  0.00           H
ATOM     37  HA3 GLY A   5      -7.866  -1.212 -20.932  1.00  0.00           H
ATOM     38  N   GLY A   6      -6.117  -2.123 -22.731  1.00  0.00           N
ATOM     39  CA  GLY A   6      -4.987  -2.658 -23.491  1.00  0.00           C
ATOM     40  C   GLY A   6      -5.160  -2.640 -24.995  1.00  0.00           C
ATOM     41  O   GLY A   6      -4.908  -1.653 -25.672  1.00  0.00           O
ATOM     42  H   GLY A   6      -6.592  -1.532 -23.131  1.00  0.00           H
ATOM     43  HA2 GLY A   6      -4.142  -2.201 -23.305  1.00  0.00           H
ATOM     44  HA3 GLY A   6      -4.889  -3.560 -23.232  1.00  0.00           H
ATOM     45  N   GLY A   7      -5.638  -3.711 -25.585  1.00  0.00           N
ATOM     46  CA  GLY A   7      -5.960  -3.699 -27.015  1.00  0.00           C
ATOM     47  C   GLY A   7      -7.032  -4.738 -27.283  1.00  0.00           C
ATOM     48  O   GLY A   7      -7.311  -5.619 -26.444  1.00  0.00           O
ATOM     49  H   GLY A   7      -5.810  -4.476 -25.186  1.00  0.00           H
ATOM     50  HA2 GLY A   7      -6.227  -2.814 -27.234  1.00  0.00           H
ATOM     51  HA3 GLY A   7      -5.157  -3.875 -27.536  1.00  0.00           H
ATOM     52  N   GLY A   8      -7.771  -4.520 -28.396  1.00  0.00           N
ATOM     53  CA  GLY A   8      -8.795  -5.411 -28.892  1.00  0.00           C
ATOM     54  C   GLY A   8     -10.113  -5.074 -28.202  1.00  0.00           C
ATOM     55  O   GLY A   8     -10.343  -3.940 -27.807  1.00  0.00           O
ATOM     56  H   GLY A   8      -7.622  -3.806 -28.902  1.00  0.00           H
ATOM     57  HA2 GLY A   8      -8.940  -5.392 -29.823  1.00  0.00           H
ATOM     58  HA3 GLY A   8      -8.673  -6.343 -28.574  1.00  0.00           H
ATOM     59  N   GLY A   9     -11.025  -6.011 -28.143  1.00  0.00           N
ATOM     60  CA  GLY A   9     -12.304  -5.594 -27.670  1.00  0.00           C
ATOM     61  C   GLY A   9     -13.026  -4.589 -28.491  1.00  0.00           C
ATOM     62  O   GLY A   9     -13.941  -3.837 -28.020  1.00  0.00           O
ATOM     63  H   GLY A   9     -10.890  -6.885 -28.344  1.00  0.00           H
ATOM     64  HA2 GLY A   9     -12.901  -6.331 -27.544  1.00  0.00           H
ATOM     65  HA3 GLY A   9     -12.107  -5.113 -26.821  1.00  0.00           H
TER
END
"""

pdb_str_04 = """
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
ATOM     86  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     87  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     88  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     89  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     90  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     91  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     92  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     93  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     94  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     95  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     96  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     97  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     98  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
ATOM     99  H   TYR A   7       8.109   0.980   6.070  1.00 14.70           H
ATOM    100  HA  TYR A   7       8.852   2.975   7.692  1.00 15.18           H
ATOM    101  HB2 TYR A   7       9.327   0.215   7.987  1.00 15.35           H
ATOM    102  HB3 TYR A   7       9.645   1.296   9.108  1.00 15.35           H
ATOM    103  HD1 TYR A   7       7.058  -0.530   7.635  1.00 15.68           H
ATOM    104  HD2 TYR A   7       7.794   2.369  10.304  1.00 14.80           H
ATOM    105  HE1 TYR A   7       4.899  -0.713   8.417  1.00 13.46           H
ATOM    106  HE2 TYR A   7       5.612   2.197  11.109  1.00 14.33           H
ATOM    107  HH  TYR A   7       3.618   1.146  10.902  1.00 14.39           H
TER
HETATM  108  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM  109  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM  110  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM  111  O   HOH A  11      11.808   4.179   9.970  1.00 23.99           O
HETATM  112  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
HETATM  113  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM  114  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
END
"""

def core(pdb_str, pair_proxies = None):
  pdb_inp = iotbx.pdb.input(source_info = None, lines = pdb_str)
  model = mmtbx.model.manager(
    model_input   = pdb_inp,
    process_input = True,
    log           = null_out())
  return mmtbx.nci.hbond.find(model=model, pair_proxies=pair_proxies)

def exercise_00():
  r = core(pdb_str=pdb_str_00)
  #r.show()
  assert len(r.result) == 5
  d_HA = flex.double()
  d_AD = flex.double()
  a_DHA = flex.double()
  for r in r.result:
    d_HA .append(r.d_HA )
    d_AD .append(r.d_AD )
    a_DHA.append(r.a_DHA)
  assert approx_equal(flex.mean(d_HA), 2.104, 1.e-2)
  assert approx_equal(flex.mean(d_AD), 2.899, 1.e-2)
  assert approx_equal(flex.mean(a_DHA), 153.2, 1.e-1)

def exercise_01():
  r = core(pdb_str=pdb_str_01)
  #r.show()
  assert len(r.result) == 2
  r0 = r.result[0]
  r1 = r.result[1]
  for r_ in r.result:
    assert approx_equal(r_.d_HA, 2.206, 1.e-2)
    assert approx_equal(r_.d_AD, 3.065, 1.e-2)
    assert approx_equal(r_.a_DHA, 142.98, 1.e-1)
    assert str(r_.symop) == "x-y,-y,-z+1"
  assert [r0.i, r0.j] == [3,6], [r0.i, r0.j]
  assert [r1.i, r1.j] == [6,3], [r1.i, r1.j]

def exercise_02():
  r = core(pdb_str=pdb_str_02)
  #r.show()
  #r.as_pymol()
  assert len(r.result) == 2
  r0 = r.result[0]
  r1 = r.result[1]
  assert approx_equal(r0.d_HA, 2.523, 1.e-2)
  assert approx_equal(r0.d_AD, 3.292, 1.e-2)
  assert approx_equal(r0.a_DHA, 149.378, 1.e-1)
  assert str(r0.symop) == "x,y,z"
  assert approx_equal(r1.d_HA, 1.644, 1.e-2)
  assert approx_equal(r1.d_AD, 2.472, 1.e-2)
  assert approx_equal(r1.a_DHA, 160.475, 1.e-1)
  assert str(r1.symop) == "x,y,z"

def exercise_03():
  r1 = core(pdb_str=pdb_str_00)
  assert len(r1.result) == 5
  r2 = core(pdb_str=pdb_str_03, pair_proxies = r1.pair_proxies)
  assert len(r2.result) == 5, len(r2.result)
  r2.as_pymol()
  for a,b in zip(r1.result, r2.result):
    assert a.i == b.i
    assert a.j == b.j

def exercise_04():
  """
  Case when A is not bonded to anything (based on 1yjp).
  """
  r = core(pdb_str=pdb_str_04)

if __name__ == '__main__':
  t0 = time.time()
  exercise_00()
  exercise_01()
  exercise_02()
  exercise_03()
  exercise_04()
  print("OK. Time: %6.3f"%(time.time()-t0))
