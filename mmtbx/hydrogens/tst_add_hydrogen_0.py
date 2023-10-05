from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from mmtbx.hydrogens import reduce_hydrogen
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal

from mmtbx.hydrogens.tst_add_hydrogen import compare_models

# ------------------------------------------------------------------------------

def run():
  test_000()
  # test_001()
  # test_002()
  # test_003()
  # test_004()
  # test_005()
  # test_006()
  # test_007()
  # test_008()

# ------------------------------------------------------------------------------

def test_000():
  '''
    Added a torsion to NH2_CTERM
  '''
  compare_models(pdb_str = pdb_str_000)

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------

pdb_str_000 = """
CRYST1   56.445   72.085  123.593  90.00  90.00  90.00 P 21 21 21
SCALE1      0.017716  0.000000  0.000000        0.00000
SCALE2      0.000000  0.013873  0.000000        0.00000
SCALE3      0.000000  0.000000  0.008091        0.00000
HETATM    1  N   DLE E  12     -46.710 -11.701  15.468  1.00 34.07           N
HETATM    2  CA  DLE E  12     -47.200 -11.850  14.105  1.00 39.64           C
HETATM    3  C   DLE E  12     -48.679 -11.482  14.009  1.00 41.13           C
HETATM    4  O   DLE E  12     -49.335 -11.798  13.021  1.00 42.80           O
HETATM    5  CB  DLE E  12     -46.390 -10.976  13.144  1.00 43.98           C
HETATM    6  CG  DLE E  12     -44.950 -11.415  12.889  1.00 46.65           C
HETATM    7  CD1 DLE E  12     -44.919 -12.765  12.197  1.00 46.88           C
HETATM    8  CD2 DLE E  12     -44.198 -10.363  12.085  1.00 48.25           C
HETATM   10  HA  DLE E  12     -47.098 -12.777  13.837  1.00 39.64           H
HETATM   11  HB2 DLE E  12     -46.356 -10.078  13.509  1.00 43.98           H
HETATM   12  HB3 DLE E  12     -46.844 -10.969  12.287  1.00 43.98           H
HETATM   13  HG  DLE E  12     -44.493 -11.511  13.739  1.00 46.65           H
HETATM   14 HD11 DLE E  12     -43.996 -13.021  12.046  1.00 46.88           H
HETATM   15 HD12 DLE E  12     -45.386 -12.696  11.349  1.00 46.88           H
HETATM   16 HD13 DLE E  12     -45.356 -13.420  12.763  1.00 46.88           H
HETATM   17 HD21 DLE E  12     -44.646 -10.236  11.234  1.00 48.25           H
HETATM   18 HD22 DLE E  12     -43.289 -10.668  11.938  1.00 48.25           H
HETATM   19 HD23 DLE E  12     -44.192  -9.530  12.582  1.00 48.25           H
HETATM   20  N   NH2 E 202     -49.193 -10.809  15.033  1.00 40.90           N
HETATM   21  HN1 NH2 E 202     -48.680 -10.596  15.728  1.00 40.90           H
HETATM   22  HN2 NH2 E 202     -50.051 -10.572  15.024  1.00 40.90           H
END
"""

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
