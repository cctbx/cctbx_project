from __future__ import absolute_import, division, print_function
import time, os
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2 as reduce2
from libtbx.utils import null_out
from mmtbx.hydrogens.tst_add_hydrogen_1 import compare_models

# ------------------------------------------------------------------------------

def run():
  test_000()
  test_001()

# ------------------------------------------------------------------------------

def test_000():
  '''
  Check that output spacegroup is correct
  '''
  model_fn = "tst_reduce2.pdb"
  with open(model_fn, "w") as f:
    f.write(pdb_str_000)

  result = run_program(program_class=reduce2.Program,args=[model_fn,
    'overwrite=True'], logger = null_out())

  m = result.model

  assert(m.crystal_symmetry().unit_cell().parameters() == \
    (19.465, 21.432, 29.523, 90.0, 90.0, 90.0))
  assert(str(m.crystal_symmetry().space_group_info()) == 'P 21 21 21')

  os.remove('tst_reduce2H.txt')
  os.remove('tst_reduce2H.pdb')
  os.remove('tst_reduce2.pdb')

# ------------------------------------------------------------------------------

def test_001():
  '''
    5b1a : SAC is residue number 1. This molecule has a N nitrogen atom but
    only a single H atom (not a propeller) because it is not at the end of the
    molecule.
  '''
  compare_models(pdb_str = pdb_str_001)

# ------------------------------------------------------------------------------

def test_002():
  '''
    1eko : AYA is residue number 1. This molecule has a N nitrogen atom but
    only a single H atom (not a propeller) because it is not at the end of the
    molecule.
  '''
  compare_models(pdb_str = pdb_str_002)

# ------------------------------------------------------------------------------

pdb_str_000 = '''
REMARK based on PDB model 3njw
CRYST1   19.465   21.432   29.523  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   GLY A   1       6.011  23.726   5.538  1.00  4.36           N
ATOM      2  CA  GLY A   1       7.279  24.418   5.504  1.00  4.79           C
ATOM      3  C   GLY A   1       8.370  23.751   6.291  1.00  4.66           C
ATOM      4  O   GLY A   1       9.449  24.344   6.484  1.00  5.69           O
ATOM      5  H1  GLY A   1       5.585  23.650   6.316  1.00  4.36           H
ATOM      8  HA2 GLY A   1       7.158  25.311   5.862  1.00  4.79           H
ATOM      9  HA3 GLY A   1       7.575  24.484   4.583  1.00  4.79           H
ATOM     10  N   ASP A   9       3.861  20.962   2.856  1.00  4.41           N
ATOM     11  CA  ASP A   9       3.324  22.217   3.387  1.00  4.26           C
ATOM     12  C   ASP A   9       3.240  23.233   2.275  1.00  4.25           C
ATOM     13  O   ASP A   9       3.655  23.005   1.142  1.00  4.88           O
ATOM     14  CB  ASP A   9       4.061  22.659   4.637  1.00  3.84           C
ATOM     15  CG  ASP A   9       5.454  23.195   4.444  1.00  3.99           C
ATOM     16  OD1 ASP A   9       6.011  23.127   3.350  1.00  4.95           O
ATOM     17  H   ASP A   9       3.281  20.481   2.442  1.00  4.41           H
ATOM     18  HA  ASP A   9       2.416  22.112   3.712  1.00  4.26           H
ATOM     19  HB2 ASP A   9       3.543  23.362   5.059  1.00  3.84           H
ATOM     20  HB3 ASP A   9       4.131  21.896   5.231  1.00  3.84           H
END
'''

pdb_str_001 = '''
REMARK based on 5b1a
CRYST1   18.410   15.860   17.467  90.00  90.00  90.00 P 1
SCALE1      0.054318  0.000000  0.000000        0.00000
SCALE2      0.000000  0.063052  0.000000        0.00000
SCALE3      0.000000  0.000000  0.057251        0.00000
HETATM    1  N   SAC I   1       7.347   7.485  10.838  1.00 72.85           N
HETATM    2  CA  SAC I   1       8.265   8.536  10.273  1.00 65.93           C
HETATM    3  C   SAC I   1       9.533   7.934   9.714  1.00 61.51           C
HETATM    4  O   SAC I   1       9.919   6.721   9.909  1.00 54.09           O
HETATM    5  CB  SAC I   1       8.769   9.656  11.243  1.00 75.40           C
HETATM    6  OG  SAC I   1       7.699   9.995  12.090  1.00 74.06           O
HETATM    7  C1A SAC I   1       6.609   6.651  10.059  1.00 74.63           C
HETATM    8  C2A SAC I   1       5.724   5.626  10.777  1.00 82.80           C
HETATM    9  OAC SAC I   1       6.705   6.749   8.717  1.00 70.21           O
HETATM   10  H   SAC I   1       7.552   7.191  11.620  1.00 72.85           H
HETATM   11  HA  SAC I   1       7.694   8.955   9.610  1.00 65.93           H
HETATM   12  HB2 SAC I   1       9.532   9.332  11.747  1.00 75.40           H
HETATM   13  HB3 SAC I   1       9.071  10.421  10.729  1.00 75.40           H
HETATM   14  HG  SAC I   1       7.429   9.283  12.467  1.00 74.06           H
HETATM   15 H2A1 SAC I   1       6.226   5.100  11.419  1.00 82.80           H
HETATM   16 H2A2 SAC I   1       5.317   5.000  10.157  1.00 82.80           H
HETATM   17 H2A3 SAC I   1       5.000   6.050  11.265  1.00 82.80           H
ATOM     18  N   THR I   2      10.131   8.757   8.935  1.00 51.98           N
ATOM     19  CA  THR I   2      11.426   8.484   8.402  1.00 54.97           C
ATOM     20  C   THR I   2      12.126   9.761   8.554  1.00 59.30           C
ATOM     21  O   THR I   2      11.599  10.860   8.274  1.00 47.76           O
ATOM     22  CB  THR I   2      11.514   7.889   6.907  1.00 58.76           C
ATOM     23  OG1 THR I   2      12.896   7.485   6.503  1.00 58.81           O
ATOM     24  CG2 THR I   2      11.011   8.856   5.886  1.00 63.25           C
ATOM     25  H   THR I   2       9.803   9.512   8.685  1.00 51.98           H
ATOM     26  HA  THR I   2      11.858   7.749   8.865  1.00 54.97           H
ATOM     27  HB  THR I   2      10.953   7.098   6.919  1.00 58.76           H
ATOM     28  HG1 THR I   2      13.410   8.149   6.535  1.00 58.81           H
ATOM     29 HG21 THR I   2      10.083   9.074   6.065  1.00 63.25           H
ATOM     30 HG22 THR I   2      11.078   8.467   5.000  1.00 63.25           H
ATOM     31 HG23 THR I   2      11.537   9.671   5.913  1.00 63.25           H
TER
END
'''

pdb_str_002 = '''
CRYST1   16.414   13.704   17.728  90.00  90.00  90.00 P 1
SCALE1      0.060924  0.000000  0.000000        0.00000
SCALE2      0.000000  0.072971  0.000000        0.00000
SCALE3      0.000000  0.000000  0.056408        0.00000
HETATM    1  N   AYA A   1       9.296   7.617  11.004  1.00 79.54           N
HETATM    2  CA  AYA A   1       7.999   6.974  10.799  1.00 77.17           C
HETATM    3  C   AYA A   1       7.913   6.349   9.406  1.00 75.41           C
HETATM    4  O   AYA A   1       8.493   5.293   9.127  1.00 75.33           O
HETATM    5  CB  AYA A   1       6.879   7.999  10.982  1.00 77.63           C
HETATM    6  CM  AYA A   1      10.269   5.465  11.808  1.00 79.53           C
HETATM    7  CT  AYA A   1      10.339   6.982  11.537  1.00 80.57           C
HETATM    8  OT  AYA A   1      11.414   7.576  11.665  1.00 78.24           O
HETATM    9  H   AYA A   1       9.252   8.471  11.098  1.00 79.54           H
HETATM   10  HA  AYA A   1       7.886   6.268  11.455  1.00 77.17           H
HETATM   11  HB1 AYA A   1       6.007   7.584  10.896  1.00 77.63           H
HETATM   12  HB2 AYA A   1       6.937   8.704  10.318  1.00 77.63           H
HETATM   13  HB3 AYA A   1       6.923   8.414  11.858  1.00 77.63           H
HETATM   14  HM1 AYA A   1       9.365   5.148  11.655  1.00 79.53           H
HETATM   15  HM2 AYA A   1      10.520   5.287  12.728  1.00 79.53           H
HETATM   16  HM3 AYA A   1      10.877   5.000  11.213  1.00 79.53           H
ATOM     17  N   SER A   2       7.170   6.983   8.504  1.00 72.40           N
ATOM     18  CA  SER A   2       7.016   6.502   7.133  1.00 67.69           C
ATOM     19  C   SER A   2       7.368   7.668   6.208  1.00 63.08           C
ATOM     20  O   SER A   2       7.553   7.516   5.000  1.00 61.07           O
ATOM     21  CB  SER A   2       5.572   6.031   6.895  1.00 68.39           C
ATOM     22  OG  SER A   2       5.010   6.600   5.725  1.00 68.99           O
ATOM     23  H   SER A   2       6.736   7.708   8.663  1.00 72.40           H
ATOM     24  HA  SER A   2       7.588   5.741   6.946  1.00 67.69           H
ATOM     25  HB2 SER A   2       5.031   6.291   7.657  1.00 68.39           H
ATOM     26  HB3 SER A   2       5.572   5.066   6.802  1.00 68.39           H
ATOM     27  HG  SER A   2       5.000   7.437   5.789  1.00 68.99           H
TER
END'''

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
