from __future__ import division
import sys
from iotbx.cli_parser import run_program
from mmtbx.programs.arginine_geometry import Program

pdb_str = '''CRYST1   38.853   46.881   71.758  90.00  90.00  90.00 P 21 21 21
SCALE1      0.025738  0.000000  0.000000        0.00000
SCALE2      0.000000  0.021331  0.000000        0.00000
SCALE3      0.000000  0.000000  0.013936        0.00000
ATOM     10  N   GLU A 248      -8.310  -2.226  10.254  1.00 17.05           N
ANISOU   10  N   GLU A 248     2878   1993   1606   -179   -365     90       N
ATOM     11  CA  GLU A 248      -8.170  -3.714  10.203  1.00 17.32           C
ANISOU   11  CA  GLU A 248     2932   1996   1654   -158   -409    183       C
ATOM     12  C   GLU A 248      -8.952  -4.335  11.364  1.00 18.33           C
ANISOU   12  C   GLU A 248     3283   2097   1585   -225   -344    222       C
ATOM     13  O   GLU A 248      -9.604  -5.377  11.174  1.00 18.66           O
ANISOU   13  O   GLU A 248     3334   2113   1644   -227   -252    268       O
ATOM     14  CB  GLU A 248      -6.710  -4.119  10.242  1.00 18.66           C
ANISOU   14  CB  GLU A 248     3052   2128   1911   -119   -629    243       C
ATOM     15  CG  GLU A 248      -5.970  -3.729   8.986  1.00 19.08           C
ANISOU   15  CG  GLU A 248     2877   2194   2177    -62   -637    207       C
ATOM     16  CD  GLU A 248      -6.210  -4.566   7.743  1.00 19.74           C
ANISOU   16  CD  GLU A 248     2830   2276   2393    -17   -531    212       C
ATOM     17  OE1 GLU A 248      -5.446  -4.333   6.804  1.00 23.25           O
ANISOU   17  OE1 GLU A 248     3120   2722   2991     18   -536    187       O
ATOM     18  OE2 GLU A 248      -7.190  -5.424   7.676  1.00 17.49           O
ANISOU   18  OE2 GLU A 248     2602   1985   2059    -28   -428    231       O
ATOM     19  N  AARG A 249      -8.924  -3.737  12.553  0.50 19.58           N
ANISOU   19  N  AARG A 249     3635   2254   1550   -290   -374    199       N
ATOM     20  CA AARG A 249      -9.730  -4.282  13.672  0.50 21.22           C
ANISOU   20  CA AARG A 249     4090   2434   1537   -373   -266    232       C
ATOM     21  C  AARG A 249     -11.212  -4.230  13.327  0.50 20.68           C
ANISOU   21  C  AARG A 249     3967   2375   1515   -395     19    165       C
ATOM     22  O  AARG A 249     -11.922  -5.224  13.583  0.50 20.67           O
ANISOU   22  O  AARG A 249     4046   2337   1470   -433    135    221       O
ATOM     23  CB AARG A 249      -9.489  -3.467  14.930  0.50 23.49           C
ANISOU   23  CB AARG A 249     4611   2727   1587   -454   -321    185       C
ATOM     24  CG AARG A 249      -8.072  -3.629  15.421  0.50 25.09           C
ANISOU   24  CG AARG A 249     4885   2911   1738   -445   -644    262       C
ATOM     25  CD AARG A 249      -8.221  -3.512  16.898  0.50 27.81           C
ANISOU   25  CD AARG A 249     5575   3244   1749   -559   -667    266       C
ATOM     26  NE AARG A 249      -8.624  -2.206  17.353  0.50 29.01           N
ANISOU   26  NE AARG A 249     5810   3425   1787   -625   -540    106       N
ATOM     27  CZ AARG A 249      -7.786  -1.191  17.460  0.50 29.90           C
ANISOU   27  CZ AARG A 249     5889   3551   1921   -629   -708     25       C
ATOM     28  NH1AARG A 249      -6.544  -1.279  17.005  0.50 29.11           N
ANISOU   28  NH1AARG A 249     5624   3444   1992   -564   -987     84       N
ATOM     29  NH2AARG A 249      -8.211  -0.073  18.013  0.50 31.92           N
ANISOU   29  NH2AARG A 249     6262   3814   2053   -702   -578   -126       N
ATOM     30  N  BARG A 249      -8.925  -3.732  12.561  0.50 19.84           N
ANISOU   30  N  BARG A 249     3670   2287   1582   -291   -374    198       N
ATOM     31  CA BARG A 249      -9.709  -4.282  13.701  0.50 21.70           C
ANISOU   31  CA BARG A 249     4157   2495   1593   -374   -271    233       C
ATOM     32  C  BARG A 249     -11.206  -4.202  13.383  0.50 20.99           C
ANISOU   32  C  BARG A 249     4017   2415   1545   -399     19    163       C
ATOM     33  O  BARG A 249     -11.925  -5.175  13.683  0.50 21.05           O
ANISOU   33  O  BARG A 249     4113   2385   1500   -439    137    219       O
ATOM     34  CB BARG A 249      -9.372  -3.588  15.023  0.50 24.47           C
ANISOU   34  CB BARG A 249     4757   2847   1695   -456   -351    202       C
ATOM     35  CG BARG A 249      -8.322  -4.340  15.822  0.50 27.11           C
ANISOU   35  CG BARG A 249     5261   3141   1899   -471   -627    332       C
ATOM     36  CD BARG A 249      -8.449  -4.152  17.324  0.50 30.43           C
ANISOU   36  CD BARG A 249     6039   3552   1971   -590   -651    335       C
ATOM     37  NE BARG A 249      -7.980  -2.852  17.724  0.50 31.82           N
ANISOU   37  NE BARG A 249     6257   3760   2073   -627   -736    211       N
ATOM     38  CZ BARG A 249      -7.824  -2.442  18.979  0.50 34.06           C
ANISOU   38  CZ BARG A 249     6849   4043   2050   -736   -813    179       C
ATOM     39  NH1BARG A 249      -8.116  -3.243  19.991  0.50 36.40           N
ANISOU   39  NH1BARG A 249     7464   4312   2056   -822   -809    280       N
ATOM     40  NH2BARG A 249      -7.383  -1.217  19.209  0.50 34.74           N
ANISOU   40  NH2BARG A 249     6937   4149   2112   -768   -887     43       N
ATOM     41  N   SER A 250     -11.652  -3.118  12.750  1.00 19.65           N
ANISOU   41  N   SER A 250     3692   2278   1498   -373    120     55       N
ATOM     42  CA  SER A 250     -13.080  -2.949  12.372  1.00 19.72           C
ANISOU   42  CA  SER A 250     3604   2282   1605   -384    368    -13       C
ATOM     43  C   SER A 250     -13.449  -4.047  11.363  1.00 18.82           C
ANISOU   43  C   SER A 250     3339   2165   1648   -344    381     46       C
ATOM     44  O   SER A 250     -14.476  -4.726  11.556  1.00 18.74           O
ANISOU   44  O   SER A 250     3351   2123   1647   -389    546     49       O
ATOM     45  CB  SER A 250     -13.334  -1.618  11.815  1.00 20.51           C
ANISOU   45  CB  SER A 250     3555   2400   1836   -348    417   -113       C
ATOM     46  OG  SER A 250     -14.727  -1.474  11.599  1.00 23.31           O
ANISOU   46  OG  SER A 250     3816   2734   2306   -359    635   -175       O
END
'''

def main(pdb_file=None):
  if pdb_file is None:
    pdb_file='tst_arginine_geometry.pdb'
    f=open(pdb_file, 'w')
    f.write(pdb_str)
    del f
  results = run_program(program_class=Program, args=[pdb_file])
  assert len(results)==1

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
