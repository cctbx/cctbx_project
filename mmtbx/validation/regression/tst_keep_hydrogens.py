from __future__ import division
pdb = """
HETATM 1488  C1  GOL A 201      34.477  11.673   5.309  0.74 10.25           C
HETATM 1489  O1  GOL A 201      33.943  11.092   6.476  0.74 11.84           O
HETATM 1490  C2  GOL A 201      34.604  10.654   4.224  0.74  6.14           C
HETATM 1491  O2  GOL A 201      33.332  10.098   3.933  0.74  8.34           O
HETATM 1492  C3  GOL A 201      35.532   9.525   4.602  0.74  6.17           C
HETATM 1493  O3  GOL A 201      35.844   8.778   3.421  0.74  6.93           O
HETATM 1494  HO3 GOL A 201      36.175   7.898   3.672  0.74 12.48           H
HETATM 1495  HO2 GOL A 201      32.969  10.480   3.104  0.74 17.23           H
HETATM 1496  HO1 GOL A 201      34.213  10.157   6.572  0.74  6.66           H
HETATM 1497  H32 GOL A 201      36.444   9.936   5.015  0.74 18.79           H
HETATM 1498  H31 GOL A 201      35.065   8.866   5.328  0.74 10.68           H
HETATM 1499  H2  GOL A 201      34.965  11.134   3.317  0.74  4.28           H
HETATM 1500  H12 GOL A 201      33.837  12.486   4.986  0.74 15.48           H
HETATM 1501  H11 GOL A 201      35.464  12.054   5.504  0.74 29.14           H
HETATM 1544  O  AHOH A 337      32.837  13.298   7.618  0.67  4.89           O
HETATM 1545  O  BHOH A 337      32.593  12.260   7.604  0.33 16.54           O
"""

import sys
from libtbx import easy_run
import StringIO

def run():
  f=file("tst_keep_hydrogens.pdb", "wb")
  f.write(pdb)
  f.close()
  for keep in range(2):
    print keep
    for prog in [
      "phenix.clashscore",
      "phenix.molprobity",
      ]:
      cmd = "%s tst_keep_hydrogens.pdb keep_hydrogens=%s" % (prog, keep)
      print cmd
      er = easy_run.fully_buffered(command=cmd)
      std = StringIO.StringIO()
      er.show_stdout(out=std)
      cs=None
      for line in std.getvalue().splitlines():
        if line.find("clashscore =")>-1:
          cs = float(line.split()[2])
          break
      print 'clashscore',cs
      if cs is None: continue
      if keep: assert cs==71.43, "clashscore is not 71.43 %s" % cs
      else: assert cs==0, "clashscore is not 0 %s" % cs

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
