from __future__ import absolute_import, division, print_function
import time, os
import mmtbx.model
import iotbx.pdb
from libtbx.utils import Sorry
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2 as reduce2

# ------------------------------------------------------------------------------

def run():
  test_000()

# ------------------------------------------------------------------------------

def test_000():
  '''
    Run reduce on a single atom model --> no H will be placed
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_000.split("\n"), source_info=None)
  # initial model
  model = mmtbx.model.manager(model_input = pdb_inp,
                              log         = null_out())
  fn = 'single_atom.cif'
  with open(fn, 'w') as fp:
    fp.write(model.model_as_mmcif())

  args=["overwrite=True",
        "%s" % fn]

  sorry_msg = ''
  try:
    result = run_program(program_class=reduce2.Program,
                         args=args,
                         logger = null_out())
  except Sorry as e:
    sorry_msg = str(e)

  assert('Is this a single atom model?' in sorry_msg)

  os.remove(fn)

pdb_str_000 = '''
CRYST1  140.400  109.000  320.000  90.00  90.00  90.00 P 21 21 2     4
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.007123  0.000000  0.000000        0.00000
SCALE2      0.000000  0.009174  0.000000        0.00000
SCALE3      0.000000  0.000000  0.003125        0.00000
ATOM      1  CA  LYS A  36      38.866  60.745 -26.642  1.00  0.00           C
ATOM      2  CA  ARG A  37      40.620  57.631 -25.579  1.00  0.00           C
ATOM      3  CA  LYS A  38      38.386  54.647 -26.046  1.00  0.00           C
ATOM      4  CA  GLU A  39      40.092  53.260 -23.023  1.00  0.00           C
ATOM      5  CA  LYS A  40      38.483  55.954 -20.912  1.00  0.00           C
ATOM      6  CA  LEU A  41      35.218  56.618 -22.675  1.00  0.00           C
ATOM      7  CA  GLU A  42      34.406  52.949 -23.237  1.00  0.00           C
ATOM      8  CA  ASN A  43      36.549  51.336 -20.567  1.00  0.00           C
ATOM      9  CA  MET A  44      35.713  53.288 -17.402  1.00  0.00           C
ATOM     10  CA  LYS A  45      33.715  50.178 -16.389  1.00  0.00           C
ATOM     11  CA  LYS A  46      36.615  48.063 -15.140  1.00  0.00           C
ATOM     12  CA  GLU A  47      35.801  49.624 -11.771  1.00  0.00           C
'''

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
