from __future__ import absolute_import, division, print_function
import libtbx.load_env
import iotbx.pdb
import mmtbx.model.statistics
from StringIO import StringIO
from libtbx.test_utils import show_diff
import os

def test_1():
  regression_pdb = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/pdb/2qxs.pdb",
      test=os.path.isfile)
  pdb_inp = iotbx.pdb.input(file_name=regression_pdb)
  h = pdb_inp.construct_hierarchy()
  stats = mmtbx.model.statistics.geometry(pdb_hierarchy=h)
  out = StringIO()
  stats.show(log=out)
  assert not show_diff(out.getvalue(), """
GEOMETRY RESTRAINTS LIBRARY: NONE
DEVIATIONS FROM IDEAL VALUES.
  BOND      :  0.000   0.000      0
  ANGLE     :  0.000   0.000      0
  CHIRALITY :  0.000   0.000      0
  PLANARITY :  0.000   0.000      0
  DIHEDRAL  :  0.000   0.000      0
  MIN NONBONDED DISTANCE : 0.000

MOLPROBITY STATISTICS.
  ALL-ATOM CLASHSCORE : 7.78
  RAMACHANDRAN PLOT:
    OUTLIERS :  0.68 %
    ALLOWED  :  0.90 %
    FAVORED  : 98.42 %
  ROTAMER OUTLIERS :  4.35 %
  CBETA DEVIATIONS :  0.87 %
  PEPTIDE PLANE:
    CIS-PROLINE     : 7.14 %
    CIS-GENERAL     : 0.00 %
    TWISTED PROLINE : 0.00 %
    TWISTED GENERAL : 0.42 %
""")

if __name__ == '__main__':
  test_1()
