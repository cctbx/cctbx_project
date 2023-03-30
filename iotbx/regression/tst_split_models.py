from __future__ import absolute_import, division, print_function

import os
import iotbx.pdb
from libtbx import easy_run
from libtbx.test_utils import approx_equal


def exercise():
  if (os.path.isfile("multi_model_1.pdb")):
    os.remove("multi_model_1.pdb")
  if (os.path.isfile("multi_model_2.pdb")):
    os.remove("multi_model_2.pdb")
  f = open("multi_model.pdb", "w")
  f.write("""\
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      1.000000  0.000000  0.000000        0.00000
SCALE2      0.000000  1.000000  0.000000        0.00000
SCALE3      0.000000  0.000000  1.000000        0.00000
MODEL        1
ATOM      1  N   GLY P  -1     -22.866  -2.627  15.217  1.00  0.00           N
ATOM      2  CA  GLY P  -1     -22.714  -3.068  16.621  1.00  0.00           C
ATOM      3  C   GLY P  -1     -21.276  -3.457  16.936  1.00  0.00           C
ATOM      4  O   GLY P  -1     -20.538  -3.887  16.047  1.00  0.00           O
ATOM      5  H1  GLY P  -1     -22.583  -3.364  14.590  1.00  0.00           H
ATOM      6  H2  GLY P  -1     -22.293  -1.817  15.040  1.00  0.00           H
ATOM      7  H3  GLY P  -1     -23.828  -2.392  15.027  1.00  0.00           H
ATOM      8  HA2 GLY P  -1     -23.016  -2.261  17.288  1.00  0.00           H
ATOM      9  HA3 GLY P  -1     -23.352  -3.933  16.803  1.00  0.00           H
ENDMDL
MODEL        2
ATOM      1  N   GLY P  -1       7.889 -28.444  15.733  1.00  0.00           N
ATOM      2  CA  GLY P  -1       7.573 -27.021  15.481  1.00  0.00           C
ATOM      3  C   GLY P  -1       6.614 -26.458  16.523  1.00  0.00           C
ATOM      4  O   GLY P  -1       5.922 -27.210  17.215  1.00  0.00           O
ATOM      5  H1  GLY P  -1       8.513 -28.794  15.023  1.00  0.00           H
ATOM      6  H2  GLY P  -1       7.043 -28.994  15.722  1.00  0.00           H
ATOM      7  H3  GLY P  -1       8.329 -28.551  16.633  1.00  0.00           H
ATOM      8  HA2 GLY P  -1       8.492 -26.437  15.501  1.00  0.00           H
ATOM      9  HA3 GLY P  -1       7.112 -26.918  14.499  1.00  0.00           H
ENDMDL
END""")
  f.close()
  out = easy_run.fully_buffered("iotbx.pdb.split_models multi_model.pdb")
  assert (len(out.stdout_lines) == 2) and (len(out.stderr_lines) == 0)
  assert (os.path.isfile("multi_model_1.pdb") and
          os.path.isfile("multi_model_2.pdb"))
  sites = []
  for file in ["multi_model_1.pdb", "multi_model_2.pdb"] :
    pdb_in = iotbx.pdb.input(file)
    h = pdb_in.construct_hierarchy()
    assert (len(h.models()) == 1)
    sites.append(h.atoms().extract_xyz())
  assert approx_equal(sites[0].rms_difference(sites[1]), 38.9216638)

if (__name__ == "__main__"):
  exercise()
  print("OK")
