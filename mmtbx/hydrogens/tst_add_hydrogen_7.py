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
  test_001()
  test_002()
  test_003()

# ------------------------------------------------------------------------------

def test_000():
  '''
    Run reduce on a CA-only model --> stops before H placement
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

  assert('The model contains only CA atoms' in sorry_msg), sorry_msg

  os.remove(fn)

# ------------------------------------------------------------------------------

def sorry_for(pdb_str):
  '''Run reduce2 on pdb_str (written as mmCIF); return the Sorry message.'''
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
  fn = 'trace.cif'
  with open(fn, 'w') as fp:
    fp.write(model.model_as_mmcif())
  sorry_msg = ''
  try:
    run_program(program_class=reduce2.Program,
                args=["overwrite=True", fn], logger=null_out())
  except Sorry as e:
    sorry_msg = str(e)
  os.remove(fn)
  return sorry_msg

def test_001():
  '''
    Run reduce on a P-only model --> stops before H placement
    (it used to place nothing and end with "Is this a single atom model?")
  '''
  sorry_msg = sorry_for(pdb_str_001)
  assert('The model contains only P atoms' in sorry_msg), sorry_msg

def test_002():
  '''
    Run reduce on a CA + P model (protein-RNA/DNA trace) --> stops before H
    placement
  '''
  sorry_msg = sorry_for(pdb_str_002)
  assert('The model contains only CA and P atoms' in sorry_msg), sorry_msg

def test_003():
  '''
    Waters and ions do not hide a trace: CA-only and P-only models with HOH
    and MG/ZN still stop before H placement
  '''
  extra = '''
HETATM  901  O   HOH C   1       5.000   5.000   5.000  1.00 20.00           O
HETATM  902  O   HOH C   2       8.000   5.000   5.000  1.00 20.00           O
HETATM  903 MG    MG D   1      11.000   5.000   5.000  1.00 20.00          MG
HETATM  904 ZN    ZN D   2      14.000   5.000   5.000  1.00 20.00          ZN
'''
  sorry_msg = sorry_for(pdb_str_000 + extra)
  assert('The model contains only CA atoms' in sorry_msg), sorry_msg
  sorry_msg = sorry_for(pdb_str_001 + extra)
  assert('The model contains only P atoms' in sorry_msg), sorry_msg

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

pdb_str_001 = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM     64  P    DG A   3      -0.558   7.198   8.645  1.00  0.00           P
ATOM     97  P    DG A   4      -2.235  11.036   3.182  1.00  0.00           P
ATOM    130  P    DT A   5      -2.469  12.378  -2.597  1.00  0.00           P
ATOM    162  P    DT A   6       4.702  14.634  -3.585  1.00  0.00           P
ATOM    194  P    DA A   7       7.706  10.079   1.228  1.00  0.00           P
ATOM    259  P    DG A   9       6.964   4.072  -7.633  1.00  0.00           P
'''

pdb_str_002 = '''
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM      2  CA  PRO A   0     -12.809  22.034  19.642  1.00 14.11           C
ATOM      9  CA  GLY A   1     -11.125  18.632  19.764  1.00  8.35           C
ATOM     13  CA  VAL A   2      -8.060  20.002  18.002  1.00  6.86           C
ATOM     20  CA  GLN A   3      -6.818  23.398  16.847  1.00  7.58           C
ATOM     29  CA  VAL A   4      -4.007  23.489  14.304  1.00  7.48           C
ATOM     36  CA AGLU A   5      -1.702  26.529  14.261  0.50  7.49           C
ATOM     64  P    DG B   3      -0.558   7.198   8.645  1.00  0.00           P
ATOM     97  P    DG B   4      -2.235  11.036   3.182  1.00  0.00           P
ATOM    130  P    DT B   5      -2.469  12.378  -2.597  1.00  0.00           P
ATOM    162  P    DT B   6       4.702  14.634  -3.585  1.00  0.00           P
ATOM    194  P    DA B   7       7.706  10.079   1.228  1.00  0.00           P
ATOM    259  P    DG B   9       6.964   4.072  -7.633  1.00  0.00           P
'''

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
