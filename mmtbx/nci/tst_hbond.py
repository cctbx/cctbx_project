from __future__ import division
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

def core(pdb_str):
  pdb_inp = iotbx.pdb.input(source_info = None, lines = pdb_str)
  model = mmtbx.model.manager(
    model_input   = pdb_inp,
    process_input = True,
    log           = null_out())
  return mmtbx.nci.hbond.find(model=model)

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
  r = core(pdb_str=pdb_str_03)
  r.show()
  r.as_pymol()

if __name__ == '__main__':
  t0 = time.time()
  exercise_00()
  exercise_01()
  exercise_02()
  print "OK. Time: %6.3f"%(time.time()-t0)
