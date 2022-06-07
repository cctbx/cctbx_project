
from __future__ import absolute_import, division, print_function
from libtbx import easy_run

def exercise():
  pdb_in = """\
CRYST1   25.015   29.415   52.761  89.54  86.10  82.39 P 1
ATOM      1  N   LEU A   4      24.685  46.025  13.868  1.00 21.07           N
ANISOU    1  N   LEU A   4     2950   3310   1747    461   -821   -375       N
ATOM      2  CA  LEU A   4      25.638  44.959  13.585  1.00 20.86           C
ANISOU    2  CA  LEU A   4     2943   3305   1678    371   -612   -550       C
ATOM      3  C   LEU A   4      27.062  45.490  13.642  1.00 20.74           C
ANISOU    3  C   LEU A   4     2921   3288   1670    491   -610   -389       C
ATOM      4  O   LEU A   4      27.332  46.631  13.264  1.00 21.80           O
ANISOU    4  O   LEU A   4     2981   3360   1943    466   -648   -371       O
ATOM      5  CB  LEU A   4      25.385  44.352  12.197  1.00 21.20           C
ANISOU    5  CB  LEU A   4     2920   3379   1758    316   -860   -516       C
ATOM      6  CG  LEU A   4      24.065  43.606  11.957  1.00 22.17           C
ANISOU    6  CG  LEU A   4     2955   3451   2017    289   -654   -556       C
ATOM      7  CD1 LEU A   4      23.889  43.288  10.474  1.00 22.59           C
ANISOU    7  CD1 LEU A   4     3049   3450   2083    277   -581   -559       C
ATOM      8  CD2 LEU A   4      23.972  42.321  12.770  1.00 22.51           C
ANISOU    8  CD2 LEU A   4     3022   3432   2100    219   -528   -600       C
"""
  with open("tmp_rigid_bond_test.pdb", "w") as f:
    f.write(pdb_in)
  result = easy_run.fully_buffered(
    "mmtbx.rigid_bond_test tmp_rigid_bond_test.pdb").raise_if_errors()
  assert (result.return_code == 0)
  #print "\n".join(result.stdout_lines)
  assert ("""  pdb=" N   LEU A   4 " pdb=" CA  LEU A   4 "     24.957""" in
          result.stdout_lines)
  assert ("""    mean = 52.797""" in result.stdout_lines)
  print("OK")

if (__name__ == "__main__"):
  exercise()
