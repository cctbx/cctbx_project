from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.hydrogens.tst_add_hydrogen_1 import compare_models

def run():
  test_001()
  test_002()
  test_003()

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
    N-terminal residue with two altloc conformers --> both conformers receive
    NH3 propeller (H1, H2, H3); neither retains the standard backbone 'H'.
    Regression test for shadowed nested loop in place_n_terminal_propeller
    that only removed the placeholder 'H' from the first altloc, leaving the
    other altloc with [H, H2, H3] -- a combination pdb_interpretation cannot
    parameterize, so H2 and H3 were silently dropped.
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_002.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model = model.select(~model.get_hd_selection())
  add_h = reduce_hydrogen.place_hydrogens(model=model)
  add_h.run()
  model_h = add_h.get_model()

  altloc_h_names = {}
  rg = model_h.get_hierarchy().models()[0].chains()[0].residue_groups()[0]
  assert rg.resseq_as_int() == 1
  for ag in rg.atom_groups():
    if ag.altloc == '': continue
    altloc_h_names[ag.altloc] = set(
      a.name.strip() for a in ag.atoms() if a.element.strip() in ('H', 'D'))

  assert set(altloc_h_names) == {'A', 'B'}, \
    "expected altlocs A and B, got %s" % sorted(altloc_h_names)
  for alt in ('A', 'B'):
    for h in ('H1', 'H2', 'H3'):
      assert h in altloc_h_names[alt], \
        "altloc %s SER 1 is missing %s; got %s" % (
          alt, h, sorted(altloc_h_names[alt]))
    assert 'H' not in altloc_h_names[alt], \
      "altloc %s SER 1 retained plain backbone 'H' instead of NH3 propeller" % alt

# ------------------------------------------------------------------------------

def test_003():
  '''
    A nucleotide with a missing dihedral.
  '''
  compare_models(pdb_str = pdb_str_003)

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
REMARK N-terminal SER in double conformation (altlocs A/B) -- both should get NH3
CRYST1   46.893   46.017   95.993  90.00  98.40  90.00 P 1 21 1
ATOM      1  C   SER A   1     -12.619  34.880  11.188  1.00 18.10           C
ATOM      2  O   SER A   1     -12.707  34.745   9.948  1.00 22.52           O
ATOM      3  N  ASER A   1     -14.822  36.058  11.168  0.55 24.23           N
ATOM      4  CA ASER A   1     -13.798  35.463  11.990  0.55 23.39           C
ATOM      5  CB ASER A   1     -14.561  34.377  12.750  0.55 24.81           C
ATOM      6  OG ASER A   1     -15.918  34.652  13.022  0.55 34.66           O
ATOM      7  N  BSER A   1     -14.115  36.726  11.048  0.46 26.00           N
ATOM      8  CA BSER A   1     -13.579  35.763  11.997  0.46 23.88           C
ATOM      9  CB BSER A   1     -14.690  34.985  12.696  0.46 23.97           C
ATOM     10  OG BSER A   1     -15.372  34.116  11.815  0.46 29.41           O
ATOM     11  N   TYR A   2     -11.649  34.403  11.891  1.00 10.98           N
ATOM     12  CA  TYR A   2     -10.600  33.615  11.277  1.00  9.36           C
ATOM     13  C   TYR A   2     -11.075  32.217  10.979  1.00  8.66           C
ATOM     14  O   TYR A   2     -11.965  31.658  11.655  1.00  9.72           O
ATOM     15  CB  TYR A   2      -9.378  33.538  12.232  1.00  8.96           C
ATOM     16  CG  TYR A   2      -8.734  34.820  12.606  1.00  8.86           C
ATOM     17  CD1 TYR A   2      -8.761  35.972  11.813  1.00  9.47           C
ATOM     18  CD2 TYR A   2      -7.981  34.882  13.783  1.00  9.74           C
ATOM     19  CE1 TYR A   2      -8.113  37.147  12.185  1.00  9.54           C
ATOM     20  CE2 TYR A   2      -7.335  36.044  14.145  1.00  9.97           C
ATOM     21  CZ  TYR A   2      -7.406  37.182  13.362  1.00  8.39           C
ATOM     22  OH  TYR A   2      -6.760  38.306  13.784  1.00  9.22           O
ATOM     23  N   THR A   3     -10.473  31.593   9.961  1.00  8.54           N
ATOM     24  CA  THR A   3     -10.697  30.229   9.657  1.00  8.82           C
ATOM     25  C   THR A   3      -9.356  29.515   9.445  1.00  7.80           C
ATOM     26  O   THR A   3      -8.339  30.115   9.091  1.00  8.51           O
ATOM     27  CB  THR A   3     -11.537  30.085   8.363  1.00 11.47           C
ATOM     28  OG1 THR A   3     -10.781  30.625   7.267  1.00 13.80           O
ATOM     29  CG2 THR A   3     -12.878  30.751   8.531  1.00 16.31           C
TER
"""

pdb_str_003 = '''
CRYST1   60.683   61.851   76.893  90.00  90.00  90.00 P 21 21 21
SCALE1      0.016479  0.000000  0.000000        0.00000
SCALE2      0.000000  0.016168  0.000000        0.00000
SCALE3      0.000000  0.000000  0.013005        0.00000
HETATM    1  C1' ADP A1311      -7.459  14.326  10.821  0.60 14.02           C
HETATM    2  C2  ADP A1311      -5.449  12.545   7.224  0.60 15.85           C
HETATM    3  C2' ADP A1311      -7.594  15.611  11.604  0.60 14.04           C
HETATM    4  C3' ADP A1311      -8.097  15.126  12.927  0.60 13.75           C
HETATM    5  C4  ADP A1311      -6.903  13.885   8.453  0.60 14.26           C
HETATM    6  C4' ADP A1311      -8.964  13.965  12.513  0.60 14.88           C
HETATM    7  C5  ADP A1311      -7.414  14.430   7.198  0.60 14.79           C
HETATM    8  C5' ADP A1311     -10.403  14.410  12.351  0.60 12.89           C
HETATM    9  C6  ADP A1311      -6.877  13.954   5.942  0.60 15.37           C
HETATM   10  C8  ADP A1311      -8.458  15.356   8.818  0.60 13.69           C
HETATM   11  N1  ADP A1311      -5.907  13.025   6.029  0.60 16.08           N
HETATM   12  N3  ADP A1311      -5.925  12.949   8.425  0.60 15.93           N
HETATM   13  N6  ADP A1311      -7.358  14.471   4.786  0.60 16.26           N
HETATM   14  N7  ADP A1311      -8.348  15.322   7.484  0.60 15.22           N
HETATM   15  N9  ADP A1311      -7.603  14.521   9.383  0.60 14.10           N
HETATM   16  O1A ADP A1311     -12.593  15.064  10.083  0.60  6.78           O
HETATM   17  O1B ADP A1311     -11.480  19.808  11.509  0.30 13.50           O
HETATM   18  O2' ADP A1311      -6.345  16.239  11.737  0.60 12.16           O
HETATM   19  O2A ADP A1311     -12.829  15.834  12.325  0.60  8.51           O
HETATM   20  O2B ADP A1311     -10.576  17.778  12.692  0.30 17.96           O
HETATM   21  O3' ADP A1311      -7.034  14.659  13.737  0.60 15.01           O
HETATM   22  O3A ADP A1311     -11.739  17.443  10.617  0.60 12.19           O
HETATM   23  O3B ADP A1311     -13.078  18.233  12.535  0.20 17.22           O
HETATM   24  O4' ADP A1311      -8.522  13.493  11.245  0.60 15.01           O
HETATM   25  O5' ADP A1311     -10.557  15.534  11.511  0.60 12.05           O
HETATM   26  PA  ADP A1311     -12.004  15.966  11.098  0.60 10.28           P
HETATM   27  PB  ADP A1311     -11.709  18.390  11.914  0.30 15.89           P
HETATM   28  H2  ADP A1311      -4.662  11.801   7.209  0.60 15.85           H
HETATM   29  H1' ADP A1311      -6.457  13.916  11.009  0.60 14.02           H
HETATM   30  H2' ADP A1311      -8.248  16.358  11.133  0.60 14.04           H
HETATM   31  H3' ADP A1311      -8.605  15.901  13.517  0.60 13.75           H
HETATM   32  H4' ADP A1311      -8.894  13.181  13.280  0.60 14.88           H
HETATM   33  H8  ADP A1311      -9.154  15.983   9.361  0.60 13.69           H
HETATM   34 H5'1 ADP A1311     -10.983  13.581  11.940  0.60 12.89           H
HETATM   35 H5'2 ADP A1311     -10.810  14.648  13.336  0.60 12.89           H
HETATM   36 HN61 ADP A1311      -6.985  14.160   3.900  0.60 16.26           H
HETATM   37 HN62 ADP A1311      -8.089  15.167   4.810  0.60 16.26           H
HETATM   38 HO2' ADP A1311      -6.432  17.010  12.315  0.60 12.16           H
HETATM   39 HO3' ADP A1311      -6.479  15.405  14.002  0.60 15.01           H
'''

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
