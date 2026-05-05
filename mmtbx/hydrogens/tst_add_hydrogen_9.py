from __future__ import absolute_import, division, print_function
import time
from mmtbx.hydrogens.tst_add_hydrogen_1 import compare_models

def run():
  test_002()
  test_001()

# ------------------------------------------------------------------------------

def test_001():
  '''
    CD1 is missing --> H, HD1, HE1 can't be parameterized --> should not be added.
  '''
  compare_models(pdb_str      = pdb_str_001,
                 not_contains = ' HE1')

# ------------------------------------------------------------------------------

def test_002():
  '''
    EKB ligand (geostd cif): the four primary-amine H atoms HN1A/B/E/F have
    only CONST (period=0, ESD=0) torsions as third-neighbor source. They
    are placed correctly via the CONST-proxy connectivity path; pre-fix
    they were silently dropped.
  '''
  compare_models(pdb_str  = pdb_str_002,
                 contains = 'HN1A')

# ------------------------------------------------------------------------------

pdb_str_001 = """
REMARK CD1 is missing --> HD1, HE1 cannot not be placed; H is not placed
CRYST1   17.955   13.272   13.095  90.00  90.00  90.00 P 1
ATOM      1  N   TYR A 139      10.241   7.920   5.000  1.00 10.00           N
ATOM      2  CA  TYR A 139      10.853   7.555   6.271  1.00 10.00           C
ATOM      3  C   TYR A 139      12.362   7.771   6.227  1.00 10.00           C
ATOM      4  O   TYR A 139      12.955   8.272   7.181  1.00 10.00           O
ATOM      5  CB  TYR A 139      10.540   6.098   6.617  1.00 10.00           C
ATOM      6  CG  TYR A 139       9.063   5.805   6.749  1.00 10.00           C
ATOM      7  CD2 TYR A 139       8.414   5.943   7.969  1.00 10.00           C
ATOM      8  CE1 TYR A 139       6.966   5.122   5.770  1.00 10.00           C
ATOM      9  CE2 TYR A 139       7.064   5.676   8.095  1.00 10.00           C
ATOM     10  CZ  TYR A 139       6.345   5.266   6.993  1.00 10.00           C
ATOM     11  OH  TYR A 139       5.000   5.000   7.113  1.00 10.00           O
ATOM     12  H   TYR A 139       9.382   7.879   5.001  1.00 10.00           H
ATOM     13  HA  TYR A 139      10.487   8.115   6.973  1.00 10.00           H
ATOM     14  HB2 TYR A 139      10.961   5.881   7.464  1.00 10.00           H
ATOM     15  HB3 TYR A 139      10.893   5.529   5.916  1.00 10.00           H
ATOM     16  HD2 TYR A 139       8.896   6.220   8.714  1.00 10.00           H
ATOM     17  HE2 TYR A 139       6.643   5.772   8.919  1.00 10.00           H
ATOM     18  HH  TYR A 139       4.752   5.127   7.905  1.00 10.00           H
TER
"""

pdb_str_002 = """
REMARK EKB ligand: HN1A/B/E/F driven by CONST torsions (period=0, ESD=0).
REMARK H atom positions are reduce_hydrogen output (idealized) - not the
REMARK QM-optimized positions in the geostd cif - so the compare_models
REMARK roundtrip is self-consistent within its 0.01 A tolerance.
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
HETATM    1  N1  EKB A   1      -3.949   1.249  -1.562  1.00 20.00           N
HETATM    2  C2  EKB A   1      -5.127   0.815  -1.114  1.00 20.00           C
HETATM    3  N3  EKB A   1      -5.332   0.048  -0.036  1.00 20.00           N
HETATM    4  C4  EKB A   1      -4.257  -0.317   0.650  1.00 20.00           C
HETATM    5  C5  EKB A   1      -2.971   0.072   0.283  1.00 20.00           C
HETATM    6  C6  EKB A   1      -2.872   0.886  -0.878  1.00 20.00           C
HETATM    7  C1A EKB A   1      -4.316  -2.688   1.438  1.00 20.00           C
HETATM    8  C1B EKB A   1       3.560   3.153   1.638  1.00 20.00           C
HETATM    9  C1C EKB A   1       3.010  -3.103  -1.785  1.00 20.00           C
HETATM   10  C1D EKB A   1       6.031   0.189  -0.851  1.00 20.00           C
HETATM   11  N1E EKB A   1      -6.214   1.187  -1.806  1.00 20.00           N
HETATM   12  N1F EKB A   1      -1.681   1.304  -1.322  1.00 20.00           N
HETATM   13  C1G EKB A   1      -0.768  -0.601   1.537  1.00 20.00           C
HETATM   14  C1H EKB A   1      -1.803  -0.310   0.988  1.00 20.00           C
HETATM   15  C1I EKB A   1       2.264   0.674   1.426  1.00 20.00           C
HETATM   16  C1J EKB A   1       2.082  -1.438   0.274  1.00 20.00           C
HETATM   17  C1K EKB A   1      -4.478  -1.218   1.827  1.00 20.00           C
HETATM   18  C1L EKB A   1       0.503  -0.953   2.161  1.00 20.00           C
HETATM   19  O1O EKB A   1       3.959   2.221   0.652  1.00 20.00           O
HETATM   20  O1P EKB A   1       3.608  -1.840  -1.563  1.00 20.00           O
HETATM   21  O1Q EKB A   1       4.727   0.559  -1.277  1.00 20.00           O
HETATM   22  C1S EKB A   1       1.661  -0.565   1.268  1.00 20.00           C
HETATM   23  C1V EKB A   1       3.307   1.048   0.580  1.00 20.00           C
HETATM   24  C1W EKB A   1       3.123  -1.068  -0.575  1.00 20.00           C
HETATM   25  C1Y EKB A   1       3.732   0.181  -0.430  1.00 20.00           C
HETATM   26  H1A EKB A   1      -4.978  -2.938   0.774  1.00 20.00           H
HETATM   27  H1B EKB A   1       3.660   2.799   2.535  1.00 20.00           H
HETATM   28  H1C EKB A   1       3.093  -3.684  -1.012  1.00 20.00           H
HETATM   29  H1D EKB A   1       6.127  -0.774  -0.782  1.00 20.00           H
HETATM   30  H1I EKB A   1       1.962   1.242   2.098  1.00 20.00           H
HETATM   31  H1J EKB A   1       1.662  -2.263   0.186  1.00 20.00           H
HETATM   32  H1K EKB A   1      -3.844  -1.003   2.528  1.00 20.00           H
HETATM   33  H1L EKB A   1       0.584  -0.508   3.019  1.00 20.00           H
HETATM   34 H1AA EKB A   1      -4.429  -3.259   2.214  1.00 20.00           H
HETATM   35 H1AB EKB A   1      -3.434  -2.846   1.067  1.00 20.00           H
HETATM   36 H1BA EKB A   1       2.630   3.410   1.541  1.00 20.00           H
HETATM   37 H1BB EKB A   1       4.082   3.970   1.599  1.00 20.00           H
HETATM   38 H1CA EKB A   1       3.415  -3.569  -2.533  1.00 20.00           H
HETATM   39 H1CB EKB A   1       2.063  -3.026  -1.979  1.00 20.00           H
HETATM   40 H1DA EKB A   1       6.244   0.560   0.020  1.00 20.00           H
HETATM   41 H1DB EKB A   1       6.713   0.500  -1.467  1.00 20.00           H
HETATM   42 H1KA EKB A   1      -5.370  -1.080   2.183  1.00 20.00           H
HETATM   43 H1LA EKB A   1       0.529  -1.906   2.342  1.00 20.00           H
HETATM   44 HN1A EKB A   1      -6.997   0.933  -1.557  1.00 20.00           H
HETATM   45 HN1B EKB A   1      -1.632   1.791  -2.029  1.00 20.00           H
HETATM   46 HN1E EKB A   1      -6.137   1.683  -2.504  1.00 20.00           H
HETATM   47 HN1F EKB A   1      -0.964   1.085  -0.901  1.00 20.00           H
END
"""

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
