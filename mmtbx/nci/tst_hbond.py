from __future__ import division
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
import time
import mmtbx.nci.hbond

pdb_str_00 = """
ATOM      1  N   GLY A   1      -5.606  -2.251 -12.878  1.00  0.00           N
ATOM      2  CA  GLY A   1      -5.850  -1.194 -13.852  1.00  0.00           C
ATOM      3  C   GLY A   1      -5.186  -1.524 -15.184  1.00  0.00           C
ATOM      4  O   GLY A   1      -5.744  -1.260 -16.249  1.00  0.00           O
ATOM      1  N   GLY A   2      -3.992  -2.102 -15.115  1.00  0.00           N
ATOM      2  CA  GLY A   2      -3.261  -2.499 -16.313  1.00  0.00           C
ATOM      3  C   GLY A   2      -3.961  -3.660 -17.011  1.00  0.00           C
ATOM      4  O   GLY A   2      -4.016  -3.716 -18.240  1.00  0.00           O
ATOM      1  N   GLY A   3      -4.492  -4.585 -16.219  1.00  0.00           N
ATOM      2  CA  GLY A   3      -5.216  -5.731 -16.755  1.00  0.00           C
ATOM      3  C   GLY A   3      -6.531  -5.289 -17.389  1.00  0.00           C
ATOM      4  O   GLY A   3      -6.939  -5.814 -18.425  1.00  0.00           O
ATOM      1  N   GLY A   4      -7.189  -4.323 -16.758  1.00  0.00           N
ATOM      2  CA  GLY A   4      -8.442  -3.785 -17.273  1.00  0.00           C
ATOM      3  C   GLY A   4      -8.205  -3.003 -18.561  1.00  0.00           C
ATOM      4  O   GLY A   4      -9.007  -3.065 -19.492  1.00  0.00           O
ATOM      1  N   GLY A   5      -7.099  -2.269 -18.604  1.00  0.00           N
ATOM      2  CA  GLY A   5      -6.735  -1.498 -19.787  1.00  0.00           C
ATOM      3  C   GLY A   5      -6.358  -2.423 -20.939  1.00  0.00           C
ATOM      4  O   GLY A   5      -6.687  -2.157 -22.094  1.00  0.00           O
ATOM      1  N   GLY A   6      -5.665  -3.509 -20.614  1.00  0.00           N
ATOM      2  CA  GLY A   6      -5.268  -4.493 -21.614  1.00  0.00           C
ATOM      3  C   GLY A   6      -6.485  -5.236 -22.153  1.00  0.00           C
ATOM      4  O   GLY A   6      -6.565  -5.533 -23.345  1.00  0.00           O
ATOM      1  N   GLY A   7      -7.430  -5.532 -21.267  1.00  0.00           N
ATOM      2  CA  GLY A   7      -8.660  -6.212 -21.655  1.00  0.00           C
ATOM      3  C   GLY A   7      -9.529  -5.303 -22.518  1.00  0.00           C
ATOM      4  O   GLY A   7     -10.158  -5.756 -23.474  1.00  0.00           O
ATOM      1  N   GLY A   8      -9.559  -4.021 -22.172  1.00  0.00           N
ATOM      2  CA  GLY A   8     -10.324  -3.039 -22.930  1.00  0.00           C
ATOM      3  C   GLY A   8      -9.706  -2.819 -24.306  1.00  0.00           C
ATOM      4  O   GLY A   8     -10.416  -2.660 -25.299  1.00  0.00           O
ATOM      1  N   GLY A   9      -8.378  -2.810 -24.356  1.00  0.00           N
ATOM      2  CA  GLY A   9      -7.658  -2.641 -25.613  1.00  0.00           C
ATOM      3  C   GLY A   9      -7.843  -3.861 -26.508  1.00  0.00           C
ATOM      4  O   GLY A   9      -7.980  -3.734 -27.725  1.00  0.00           O
TER
"""

pdb_str_01 = """
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

def core(pdb_str, suffix):
  pdb_inp = iotbx.pdb.input(source_info = None, lines = pdb_str)
  model = mmtbx.model.manager(
    model_input   = pdb_inp,
    process_input = True,
    log           = null_out())
  get_h_bonds = mmtbx.nci.hbond.get_hydrogen_bonds(model=model)
  get_h_bonds.write_restrains_file(pdb_file_name="exercise_%s.eff"%suffix)
  return get_h_bonds.get_hydrogen_bonds_pairs()

def exercise_00():
  r1 = core(pdb_str=pdb_str_00, suffix=1)
  r2 = core(pdb_str=pdb_str_01, suffix=2)
  assert len(r1) == len(r2)
  for r in r1:
    print r
  print
  for r in r2:
    print r

if __name__ == '__main__':
  t0 = time.time()
  exercise_00()
  print "OK. Time: %6.3f"%(time.time()-t0)
