from __future__ import absolute_import, division, print_function
import mmtbx
import iotbx.pdb
import mmtbx.model
from libtbx import easy_run
from libtbx.utils import null_out
import mmtbx.atomic_environment_vectors as aev

pdb_str = """
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
ATOM      1  N   GLY A  10      -7.846  -5.040 -25.897  1.00  0.00           N
ATOM      2  CA  GLY A  10      -8.046  -6.284 -26.631  1.00  0.00           C
ATOM      3  C   GLY A  10      -9.473  -6.375 -27.160  1.00  0.00           C
ATOM      4  O   GLY A  10      -9.704  -6.850 -28.272  1.00  0.00           O
"""

def run():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(make_restraints=True)
  a = aev.AEV(model = model)
  print('forward AEVs')
  print(a.BAEVs)
  print('backward AEVs')
  print(a.EAEVs)
  gv = a.get_values()
  for key in gv:
    outl = key
    for cc in gv[key]:
      outl += ' %0.3f' % cc
    print(outl)
  # print('backward AEVs')
  # print(a.EAEVs)
  # print('middle AEVs')
  # print(a.MAEVs)
  b = aev.compare(a)
  print(b)
  # assert b['GLY  A   1']['E'] is None
  print(b['GLY  A   5'])
  assert b['GLY  A   5']['E'] > 0.99
  assert b['GLY  A   5']['M'] > 0.99
  assert b['GLY  A   5']['B'] > 0.99
  # assert b['GLY  A  10']['B'] is None
  if 0:
    recs = aev.format_HELIX_records_from_AEV(b, 0.9)
    assert len(recs) == 1
    r="HELIX    1   1 GLY  A   2  GLY  A   9                                      8"
    assert r == recs[0]
  #
  if 0:
    with open("development.aev.pdb", "w") as fo:
      fo.write(pdb_str)
    assert not easy_run.call("mmtbx.development.aev development.aev.pdb 0.9")

if __name__ == '__main__':
  run()
