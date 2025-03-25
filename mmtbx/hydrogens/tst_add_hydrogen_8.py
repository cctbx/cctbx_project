from __future__ import absolute_import, division, print_function
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2 as reduce2
from libtbx.utils import null_out
import time, os

# ------------------------------------------------------------------------------

def run():
  test_000()

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

# ------------------------------------------------------------------------------

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
