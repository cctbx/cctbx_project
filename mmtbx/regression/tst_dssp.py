
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
SHEET    4   1 5 ASP A  51  ASP A  54  1  N  ASP A  51   O  LEU A 157
SHEET    5   1 5 ASP A  74  LEU A  77  1  N  ASP A  74   O  VAL A  52""")
  examples_dir = libtbx.env.find_in_repositories(
    relative_path="phenix_examples",
    test=os.path.isdir)
  if (examples_dir is not None) :
    pdb_file_1 = os.path.join(examples_dir, "p9-build", "p9.pdb")
    result = easy_run.fully_buffered("mmtbx.dssp \"%s\"" % pdb_file_1
      ).raise_if_errors()
    pdb_file_2 = os.path.join(examples_dir, "porin-twin", "porin.pdb")
    result = easy_run.fully_buffered("mmtbx.dssp \"%s\"" % pdb_file_2
      ).raise_if_errors()
    assert ("\n".join(result.stdout_lines) == """\
HELIX    5   5 ASP     59  GLY     63  5                                   5
HELIX    7   7 THR     87  VAL     92  1                                   6
HELIX   18  18 ASP    156  VAL    160  5                                   5
HELIX   19  19 ASP    184  ILE    188  5                                   5
SHEET    1   117 ILE     2  TYR    14  0
SHEET    2   117 THR    25  GLU    40 -1  N  VAL    36   O  SER     3
SHEET    3   117 THR    46  ASP    56 -1  N  TRP    55   O  LEU    31
SHEET    4   117 GLN    70  TYR    75 -1  N  SER    74   O  THR    46
SHEET    5   117 VAL    78  GLY    83 -1  N  VAL    82   O  PHE    71
SHEET    6   117 ASN   131  ILE   139 -1  N  THR   136   O  THR    79
SHEET    7   117 VAL   142  ASP   150 -1  N  ASP   150   O  ASN   131
SHEET    8   117 GLU   163  SER   171 -1  N  ASP   169   O  ASN   143
SHEET    9   117 ILE   175  THR   183 -1  N  THR   183   O  PHE   164
SHEET   10   117 ILE   193  LYS   201 -1  N  ALA   199   O  SER   176
SHEET   11   117 GLY   206  ASP   214 -1  N  ASP   214   O  ALA   194
SHEET   12   117 GLN   223  PHE   232 -1  N  ASN   229   O  THR   207
SHEET   13   117 THR   235  ILE   244 -1  N  ASP   243   O  VAL   224
SHEET   14   117 ALA   252  GLN   260 -1  N  ASP   258   O  THR   236
SHEET   15   117 VAL   265  SER   273 -1  N  SER   273   O  TYR   253
SHEET   16   117 THR   279  ASP   288 -1  N  ARG   286   O  LYS   266
SHEET   17   117 TYR     7  VAL    15 -1  N  TYR    14   O  ALA   281""")
  print "OK"

if (__name__ == "__main__") :
  exercise()
