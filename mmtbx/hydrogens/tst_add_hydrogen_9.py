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

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
