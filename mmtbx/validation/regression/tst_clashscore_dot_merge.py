from __future__ import absolute_import, division, print_function
"""
clashscore and clashscore2 must find the same clashes with condensed_probe=False (the
default, which reads probe's per-dot output) and with condensed_probe=True.

The per-dot parsers merge consecutive dot lines of the same atom pair, and the merged group
keeps the type of its first line. If a pair's 'so' (small overlap) lines come right before
its 'bo' (bad overlap) lines, the pair must still be a clash (cctbx_project#413).
"""
import os
import time
import libtbx.load_env
from iotbx import pdb
from iotbx.cli_parser import CCTBXParser
from libtbx.utils import null_out, multi_out
from mmtbx.validation import clashscore

# LYS 83 - LEU 84 of chain A of PDB entry 104L with the H that clashscore adds to the whole
# model, rotated 90 degrees about z: (x, y, z) -> (y, -x, z). In probe's per-dot output, the
# dots of A 84 LEU HD23 -> HA are 8 'so' lines followed by 2 'bo' lines (gap -0.486).
pdb_str = """
ATOM    652  N   LYS A  83      21.442   0.177  36.860  1.00 45.97           N
ATOM    653  CA  LYS A  83      22.314  -0.548  37.708  1.00 44.04           C
ATOM    654  C   LYS A  83      23.457   0.308  38.189  1.00 55.28           C
ATOM    655  O   LYS A  83      24.210  -0.194  39.028  1.00 64.07           O
ATOM    656  CB  LYS A  83      22.803  -1.760  37.002  1.00 43.91           C
ATOM    657  CG  LYS A  83      22.168  -1.927  35.615  1.00 27.47           C
ATOM    658  CD  LYS A  83      22.853  -3.060  34.852  1.00100.00           C
ATOM    659  CE  LYS A  83      23.917  -3.789  35.700  1.00100.00           C
ATOM    660  NZ  LYS A  83      24.003  -5.261  35.516  1.00 57.70           N
ATOM      0  H   LYS A  83      21.758   0.454  36.110  1.00 45.97           H
ATOM      0  HA  LYS A  83      21.821  -0.824  38.497  1.00 44.04           H
ATOM      0  HB2 LYS A  83      23.767  -1.710  36.909  1.00 43.91           H
ATOM      0  HB3 LYS A  83      22.610  -2.544  37.540  1.00 43.91           H
ATOM      0  HG2 LYS A  83      21.221  -2.116  35.706  1.00 27.47           H
ATOM      0  HG3 LYS A  83      22.245  -1.099  35.115  1.00 27.47           H
ATOM      0  HD2 LYS A  83      22.184  -3.699  34.559  1.00100.00           H
ATOM      0  HD3 LYS A  83      23.271  -2.701  34.054  1.00100.00           H
ATOM      0  HE2 LYS A  83      24.784  -3.405  35.498  1.00100.00           H
ATOM      0  HE3 LYS A  83      23.738  -3.608  36.636  1.00100.00           H
ATOM      0  HZ1 LYS A  83      24.640  -5.590  36.043  1.00 57.70           H
ATOM      0  HZ2 LYS A  83      23.223  -5.633  35.727  1.00 57.70           H
ATOM      0  HZ3 LYS A  83      24.198  -5.443  34.667  1.00 57.70           H
ATOM    661  N   LEU A  84      23.576   1.580  37.674  1.00 49.40           N
ATOM    662  CA  LEU A  84      24.628   2.488  38.126  1.00 41.55           C
ATOM    663  C   LEU A  84      24.175   3.267  39.309  1.00 46.70           C
ATOM    664  O   LEU A  84      24.850   3.214  40.389  1.00 49.33           O
ATOM    665  CB  LEU A  84      25.296   3.389  37.063  1.00 31.16           C
ATOM    666  CG  LEU A  84      26.154   2.468  36.183  1.00 60.76           C
ATOM    667  CD1 LEU A  84      27.141   3.211  35.331  1.00 14.11           C
ATOM    668  CD2 LEU A  84      26.931   1.423  36.978  1.00 17.79           C
ATOM      0  H   LEU A  84      23.055   1.908  37.074  1.00 49.40           H
ATOM      0  HA  LEU A  84      25.353   1.893  38.373  1.00 41.55           H
ATOM      0  HB2 LEU A  84      24.627   3.849  36.532  1.00 31.16           H
ATOM      0  HB3 LEU A  84      25.843   4.071  37.484  1.00 31.16           H
ATOM      0  HG  LEU A  84      25.500   2.028  35.617  1.00 60.76           H
ATOM      0 HD11 LEU A  84      27.651   2.579  34.800  1.00 14.11           H
ATOM      0 HD12 LEU A  84      26.667   3.819  34.742  1.00 14.11           H
ATOM      0 HD13 LEU A  84      27.744   3.716  35.899  1.00 14.11           H
ATOM      0 HD21 LEU A  84      27.451   0.875  36.370  1.00 17.79           H
ATOM      0 HD22 LEU A  84      27.526   1.867  37.602  1.00 17.79           H
ATOM      0 HD23 LEU A  84      26.310   0.862  37.468  1.00 17.79           H
END
"""

expected = [' A  84  LEU  HA   A  84  LEU HD23 :0.486']


def exercise_clashscore():
  found = {}
  for cp in [True, False]:
    h = pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
    cs = clashscore.clashscore(pdb_hierarchy=h, keep_hydrogens=True, fast=False,
                               condensed_probe=cp, out=null_out())
    found[cp] = [c.format_old() for c in cs.results]
  assert found[False] == found[True] == expected, found


def exercise_clashscore2():
  from mmtbx.programs import clashscore2
  fname = "tst_clashscore_dot_merge.pdb"
  with open(fname, "w") as f:
    f.write(pdb_str)
  found = {}
  for cp in [True, False]:
    parser = CCTBXParser(program_class=clashscore2.Program, logger=null_out())
    parser.parse_args([fname, "keep_hydrogens=True", "condensed_probe=%s" % cp])
    task = clashscore2.Program(parser.data_manager, parser.working_phil.extract(),
                               logger=null_out())
    task.validate()
    task.run()
    found[cp] = [c.format_old() for c in task.get_results().results]
  os.remove(fname)
  # condensed_probe=True prints the centre-to-centre overlap for clashscore2 (probe2's
  # condensed JSON has no dot gap), so compare the atom pairs there and the value here.
  pairs = dict((cp, [c.split(":")[0] for c in found[cp]]) for cp in found)
  assert pairs[False] == pairs[True] == [expected[0].split(":")[0]], found
  assert found[False] == expected, found


if (__name__ == "__main__"):
  if (not libtbx.env.has_module(name="probe")):
    print("Skipping tst_clashscore_dot_merge: probe not configured")
    print("OK")
  else:
    t0 = time.time()
    exercise_clashscore()
    exercise_clashscore2()
    print("OK. Time: %8.3f" % (time.time() - t0))
