
from __future__ import division
from libtbx import easy_run
import libtbx.load_env
import os

def exercise () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  if (pdb_file is None) :
    print "phenix_regression not available, skipping test"
    return
  result = easy_run.fully_buffered("mmtbx.dssp \"%s\"" % pdb_file
    ).raise_if_errors()
  assert ("\n".join(result.stdout_lines) == """\
HELIX    2   2 ASP A   14  THR A   18  5                                   5
HELIX    6   6 ASP A   37  GLY A   48  1                                  12
HELIX    8   8 SER A   57  GLY A   65  1                                   9
HELIX   10  10 ASN A  119  GLN A  132  1                                  14
HELIX   16  16 GLY A  138  ALA A  152  1                                  15
HELIX   19  19 ASP A  165  VAL A  178  1                                  14
HELIX   22  22 ASP A  181  ARG A  191  1                                  11
HELIX   24  24 SER A  195  GLN A  209  1                                  15
HELIX   29  29 ALA A  216  ALA A  225  1                                  10
HELIX   33  33 SER A  228  GLY A  233  1                                   6
HELIX   35  35 ARG A  235  TYR A  250  1                                  16
HELIX   38  38 SER A  252  ALA A  260  1                                   9
HELIX   41  41 SER A  263  LEU A  275  1                                  13
SHEET    1   1 5 ARG A  13  ASP A  14  0
SHEET    2   1 5 LEU A  27  SER A  30 -1  N  ARG A  29   O  ARG A  13
SHEET    3   1 5 VAL A 156  HIS A 159  1  N  VAL A 156   O  PHE A  28
SHEET    4   1 5 ASP A  51  ASP A  54  1  O  ASP A  51   N  LEU A 157
SHEET    5   1 5 ASP A  74  LEU A  77  1  N  ASP A  74   O  VAL A  52""")
  print "OK"

if (__name__ == "__main__") :
  exercise()
