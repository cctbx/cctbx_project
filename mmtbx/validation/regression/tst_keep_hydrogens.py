from __future__ import absolute_import, division, print_function

# 1yiwH_Tyr-clashOH_Arg54 constructed by Jane
pdb = """
ATOM    415  N   ARG A  54      24.163   6.042   4.912  1.00  8.56           N
ATOM    416  CA  ARG A  54      24.154   7.204   5.783  1.00  7.43           C
ATOM    417  C   ARG A  54      25.587   7.695   6.065  1.00  8.17           C
ATOM    418  O   ARG A  54      26.582   7.006   5.785  1.00  7.99           O
ATOM    419  CB  ARG A  54      23.445   6.888   7.099  1.00  7.83           C
ATOM    420  CG  ARG A  54      21.960   6.614   6.948  1.00  9.29           C
ATOM    421  CD  ARG A  54      21.119   7.876   6.784  1.00 11.64           C
ATOM    422  NE  ARG A  54      19.760   7.560   6.356  1.00 12.47           N
ATOM    423  CZ  ARG A  54      19.322   7.602   5.106  1.00 13.14           C
ATOM    424  NH1 ARG A  54      20.125   7.966   4.106  1.00 14.25           N
ATOM    425  NH2 ARG A  54      18.049   7.296   4.866  1.00 14.31           N
ATOM      0  H   ARG A  54      24.238   5.175   5.406  1.00  8.56           H   new
ATOM      0  HA  ARG A  54      23.605   8.005   5.266  1.00  7.43           H   new
ATOM      0  HB2 ARG A  54      23.924   6.012   7.561  1.00  7.83           H   new
ATOM      0  HB3 ARG A  54      23.583   7.732   7.790  1.00  7.83           H   new
ATOM      0  HG2 ARG A  54      21.803   5.963   6.075  1.00  9.29           H   new
ATOM      0  HG3 ARG A  54      21.606   6.060   7.830  1.00  9.29           H   new
ATOM      0  HD2 ARG A  54      21.087   8.424   7.737  1.00 11.64           H   new
ATOM      0  HD3 ARG A  54      21.592   8.541   6.047  1.00 11.64           H   new
ATOM      0  HE  ARG A  54      19.107   7.290   7.063  1.00 12.47           H   new
ATOM      0 HH11 ARG A  54      21.076   8.213   4.293  1.00 14.25           H   new
ATOM      0 HH12 ARG A  54      19.777   7.992   3.169  1.00 14.25           H   new
ATOM      0 HH21 ARG A  54      17.445   7.041   5.621  1.00 14.31           H   new
ATOM      0 HH22 ARG A  54      17.697   7.321   3.930  1.00 14.31           H   new
ATOM    455  N   TYR A  59      25.004  10.207  12.486  1.00  8.08           N
ATOM    456  CA  TYR A  59      25.033   9.455  13.718  1.00  7.69           C
ATOM    457  C   TYR A  59      25.704  10.209  14.831  1.00  8.15           C
ATOM    458  O   TYR A  59      25.946   9.645  15.894  1.00  8.88           O
ATOM    459  CB  TYR A  59      25.738   8.114  13.502  1.00  6.47           C
ATOM    460  CG  TYR A  59      25.081   7.238  12.459  1.00  7.52           C
ATOM    461  CD1 TYR A  59      25.600   7.151  11.173  1.00  7.62           C
ATOM    462  CD2 TYR A  59      23.944   6.499  12.763  1.00  6.99           C
ATOM    463  CE1 TYR A  59      25.008   6.354  10.215  1.00 10.53           C
ATOM    464  CE2 TYR A  59      23.339   5.696  11.818  1.00  6.39           C
ATOM    465  CZ  TYR A  59      23.871   5.625  10.544  1.00  8.24           C
ATOM    466  OH  TYR A  59      23.274   4.827   9.594  1.00  8.82           O
ATOM      0  H   TYR A  59      25.881  10.214  12.006  1.00  8.08           H   new
ATOM      0  HA  TYR A  59      23.993   9.268  14.023  1.00  7.69           H   new
ATOM      0  HB2 TYR A  59      26.780   8.303  13.204  1.00  6.47           H   new
ATOM      0  HB3 TYR A  59      25.770   7.570  14.458  1.00  6.47           H   new
ATOM      0  HD1 TYR A  59      26.499   7.729  10.914  1.00  7.62           H   new
ATOM      0  HD2 TYR A  59      23.520   6.556  13.776  1.00  6.99           H   new
ATOM      0  HE1 TYR A  59      25.430   6.297   9.200  1.00 10.53           H   new
ATOM      0  HE2 TYR A  59      22.438   5.119  12.075  1.00  6.39           H   new
ATOM      0  HH  TYR A  59      23.201   5.331   8.734  1.00  8.82           H   new
"""

import sys
from libtbx import easy_run
from six.moves import cStringIO as StringIO

def run():
  f=open("tst_keep_hydrogens.pdb", "w")
  f.write(pdb)
  f.close()
  #for keep in range(2):
  for keep in [False,True]:
    print("keep_hydrogens=", str(keep))
    for prog in [
      "phenix.clashscore",
      "phenix.molprobity",
      ]:
      cmd = "%s tst_keep_hydrogens.pdb keep_hydrogens=%s" % (prog, keep)
      print(cmd)
      er = easy_run.fully_buffered(command=cmd)
      std = StringIO()
      er.show_stdout(out=std)
      cs=None
      for line in std.getvalue().splitlines():
        if line.find("clashscore =")>-1:
          cs = float(line.split()[2])
          break
      print('clashscore',cs)
      if cs is None: continue
      if keep: assert cs==44.44, "%s: clashscore is not 44.44 %s" % (prog, cs)
      else: assert cs==0, "%s: clashscore is not 0 %s" % (prog, cs)

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
