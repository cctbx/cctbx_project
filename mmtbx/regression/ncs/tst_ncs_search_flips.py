from __future__ import division
import iotbx.ncs
import iotbx.pdb

test_pdb_1 = """\
CRYST1  272.680  272.680  272.680  90.00  90.00  90.00 F 2 3        96
ATOM     17  N   LEU A   7      -3.960 -35.813 -28.701  1.00 59.07           N
ATOM     18  CA  LEU A   7      -2.543 -35.701 -29.071  1.00 52.98           C
ATOM     19  C   LEU A   7      -1.964 -34.379 -28.552  1.00 07.51           C
ATOM     20  O   LEU A   7      -2.502 -33.305 -28.832  1.00 95.11           O
ATOM     21  CB  LEU A   7      -2.337 -35.825 -30.589  1.00 93.68           C
ATOM     22  CG  LEU A   7      -2.092 -37.246 -31.147  1.00 36.57           C
ATOM     23  CD1 LEU A   7      -3.388 -38.030 -31.314  1.00 45.19           C
ATOM     24  CD2 LEU A   7      -1.342 -37.242 -32.482  1.00 52.49           C
TER
ATOM     17  N   LEU B   7       1.133 -33.422 -31.479  1.00 59.07           N
ATOM     18  CA  LEU B   7       2.550 -33.310 -31.849  1.00 52.98           C
ATOM     19  C   LEU B   7       3.129 -31.988 -31.330  1.00 07.51           C
ATOM     20  O   LEU B   7       2.591 -30.914 -31.610  1.00 95.11           O
ATOM     21  CB  LEU B   7       2.756 -33.434 -33.367  1.00 93.68           C
ATOM     22  CG  LEU B   7       3.001 -34.855 -33.925  1.00 36.57           C
ATOM     24  CD1 LEU B   7       3.751 -34.851 -35.260  1.00 52.49           C
ATOM     23  CD2 LEU B   7       1.705 -35.639 -34.092  1.00 45.19           C
END
"""

test_pdb_2 = """\
CRYST1   53.930  108.200  147.160  90.00  90.00  90.00 P 21 21 21   16
ATOM      1  N   TYR A   1      30.440   8.711   1.306  1.00 64.98           N
ATOM      2  CA  TYR A   1      29.117   9.171   1.763  1.00 67.34           C
ATOM      3  C   TYR A   1      28.007   8.102   1.965  1.00 65.79           C
ATOM      4  O   TYR A   1      26.830   8.428   1.845  1.00 58.29           O
ATOM      5  CB  TYR A   1      29.262  10.057   3.008  1.00 67.25           C
ATOM      6  CG  TYR A   1      29.855  11.415   2.711  1.00 70.44           C
ATOM      8  CD1 TYR A   1      31.185  11.716   3.026  1.00 71.73           C
ATOM      7  CD2 TYR A   1      29.079  12.412   2.112  1.00 68.62           C
ATOM     10  CE1 TYR A   1      31.716  12.971   2.753  1.00 72.09           C
ATOM      9  CE2 TYR A   1      29.603  13.668   1.843  1.00 68.12           C
ATOM     11  CZ  TYR A   1      30.919  13.946   2.160  1.00 69.15           C
ATOM     12  OH  TYR A   1      31.432  15.200   1.886  1.00 61.93           O
TER
ATOM      1  N   TYR B   1      30.440   8.711   6.306  1.00 64.98           N
ATOM      2  CA  TYR B   1      29.117   9.171   6.763  1.00 67.34           C
ATOM      3  C   TYR B   1      28.007   8.102   6.965  1.00 65.79           C
ATOM      4  O   TYR B   1      26.830   8.428   6.845  1.00 58.29           O
ATOM      5  CB  TYR B   1      29.262  10.057   8.008  1.00 67.25           C
ATOM      6  CG  TYR B   1      29.855  11.415   7.711  1.00 70.44           C
ATOM      7  CD1 TYR B   1      29.079  12.412   7.112  1.00 68.62           C
ATOM      8  CD2 TYR B   1      31.185  11.716   8.026  1.00 71.73           C
ATOM      9  CE1 TYR B   1      29.603  13.668   6.843  1.00 68.12           C
ATOM     10  CE2 TYR B   1      31.716  12.971   7.753  1.00 72.09           C
ATOM     11  CZ  TYR B   1      30.919  13.946   7.160  1.00 69.15           C
ATOM     12  OH  TYR B   1      31.432  15.200   6.886  1.00 61.93           O
"""
test_pdb_3 = """\
CRYST1   53.930  108.200  147.160  90.00  90.00  90.00 P 21 21 21   16
ATOM      1  N   TYR A   1      30.440   8.711   1.306  1.00 64.98           N
ATOM      2  CA  TYR A   1      29.117   9.171   1.763  1.00 67.34           C
ATOM      3  C   TYR A   1      28.007   8.102   1.965  1.00 65.79           C
ATOM      4  O   TYR A   1      26.830   8.428   1.845  1.00 58.29           O
ATOM      5  CB  TYR A   1      29.262  10.057   3.008  1.00 67.25           C
ATOM      6  CG  TYR A   1      29.855  11.415   2.711  1.00 70.44           C
ATOM      8  CD1 TYR A   1      31.185  11.716   3.026  1.00 71.73           C
ATOM      7  CD2 TYR A   1      29.079  12.412   2.112  1.00 68.62           C
ATOM      9  CE1 TYR A   1      29.603  13.668   1.843  1.00 68.12           C
ATOM     10  CE2 TYR A   1      31.716  12.971   2.753  1.00 72.09           C
ATOM     11  CZ  TYR A   1      30.919  13.946   2.160  1.00 69.15           C
ATOM     12  OH  TYR A   1      31.432  15.200   1.886  1.00 61.93           O
TER
ATOM      1  N   TYR B   1      30.440   8.711   6.306  1.00 64.98           N
ATOM      2  CA  TYR B   1      29.117   9.171   6.763  1.00 67.34           C
ATOM      3  C   TYR B   1      28.007   8.102   6.965  1.00 65.79           C
ATOM      4  O   TYR B   1      26.830   8.428   6.845  1.00 58.29           O
ATOM      5  CB  TYR B   1      29.262  10.057   8.008  1.00 67.25           C
ATOM      6  CG  TYR B   1      29.855  11.415   7.711  1.00 70.44           C
ATOM      7  CD1 TYR B   1      29.079  12.412   7.112  1.00 68.62           C
ATOM      8  CD2 TYR B   1      31.185  11.716   8.026  1.00 71.73           C
ATOM      9  CE1 TYR B   1      29.603  13.668   6.843  1.00 68.12           C
ATOM     10  CE2 TYR B   1      31.716  12.971   7.753  1.00 72.09           C
ATOM     11  CZ  TYR B   1      30.919  13.946   7.160  1.00 69.15           C
ATOM     12  OH  TYR B   1      31.432  15.200   6.886  1.00 61.93           O
"""

test_pdb_4 = """\
ATOM      1  N   HIS A   3      32.534  13.242   2.403  1.00  0.00           N
ATOM      2  CA  HIS A   3      33.724  12.970   3.202  1.00  0.00           C
ATOM      3  C   HIS A   3      34.993  13.295   2.422  1.00  0.00           C
ATOM      4  O   HIS A   3      35.327  12.618   1.450  1.00  0.00           O
ATOM      5  CB  HIS A   3      33.744  11.503   3.640  1.00  0.00           C
ATOM      6  CG  HIS A   3      32.618  11.130   4.554  1.00  0.00           C
ATOM      7  ND1 HIS A   3      32.586  11.494   5.882  1.00  0.00           N
ATOM      8  CD2 HIS A   3      31.485  10.424   4.330  1.00  0.00           C
ATOM      9  CE1 HIS A   3      31.481  11.029   6.437  1.00  0.00           C
ATOM     10  NE2 HIS A   3      30.795  10.375   5.517  1.00  0.00           N
TER
ATOM      1  N   HIS B   3      34.434   0.922   2.403  1.00  0.00           N
ATOM      2  CA  HIS B   3      35.624   0.650   3.202  1.00  0.00           C
ATOM      3  C   HIS B   3      36.893   0.975   2.422  1.00  0.00           C
ATOM      4  O   HIS B   3      37.227   0.298   1.450  1.00  0.00           O
ATOM      5  CB  HIS B   3      35.644  -0.817   3.640  1.00  0.00           C
ATOM      6  CG  HIS B   3      34.518  -1.190   4.554  1.00  0.00           C
ATOM      7  CD2 HIS B   3      34.486  -0.826   5.882  1.00  0.00           N
ATOM      8  ND1 HIS B   3      33.385  -1.896   4.330  1.00  0.00           C
ATOM      9  NE2 HIS B   3      33.381  -1.291   6.437  1.00  0.00           C
ATOM     10  CE1 HIS B   3      32.695  -1.945   5.517  1.00  0.00           N
"""

test_pdb_5 = """\
ATOM      1  N   ARG A   1      26.061  12.824   1.988  1.00  0.00           N
ATOM      2  CA  ARG A   1      27.253  12.525   2.773  1.00  0.00           C
ATOM      3  C   ARG A   1      28.520  12.882   2.003  1.00  0.00           C
ATOM      4  O   ARG A   1      28.853  12.243   1.005  1.00  0.00           O
ATOM      5  CB  ARG A   1      27.280  11.041   3.156  1.00 10.00           C
ATOM      6  CG  ARG A   1      26.107  10.591   4.022  1.00 10.00           C
ATOM      7  CD  ARG A   1      26.118  11.230   5.409  1.00 10.00           C
ATOM      8  NE  ARG A   1      27.283  10.828   6.201  1.00 10.00           N
ATOM      9  CZ  ARG A   1      27.735  11.441   7.298  1.00 10.00           C
ATOM     10  NH1 ARG A   1      27.146  12.525   7.803  1.00 10.00           N
ATOM     11  NH2 ARG A   1      28.808  10.956   7.908  1.00 10.00           N
ATOM      1  N   ALA A   2      29.223  13.907   2.474  1.00  0.00           N
ATOM      2  CA  ALA A   2      30.455  14.351   1.832  1.00  0.00           C
ATOM      3  C   ALA A   2      31.652  14.171   2.758  1.00  0.00           C
ATOM      4  O   ALA A   2      31.775  14.859   3.772  1.00  0.00           O
ATOM      5  CB  ALA A   2      30.331  15.807   1.408  1.00  0.00           C
ATOM      1  N   HIS A   3      32.534  13.242   2.403  1.00  0.00           N
ATOM      2  CA  HIS A   3      33.724  12.970   3.202  1.00  0.00           C
ATOM      3  C   HIS A   3      34.993  13.295   2.422  1.00  0.00           C
ATOM      4  O   HIS A   3      35.327  12.618   1.450  1.00  0.00           O
ATOM      5  CB  HIS A   3      33.744  11.503   3.640  1.00  0.00           C
ATOM      6  CG  HIS A   3      32.618  11.130   4.554  1.00  0.00           C
ATOM      7  ND1 HIS A   3      32.586  11.494   5.882  1.00  0.00           N
ATOM      8  CD2 HIS A   3      31.485  10.424   4.330  1.00  0.00           C
ATOM      9  CE1 HIS A   3      31.481  11.029   6.437  1.00  0.00           C
ATOM     10  NE2 HIS A   3      30.795  10.375   5.517  1.00  0.00           N
ATOM      1  N   ALA A   4      35.698  14.335   2.856  1.00  0.00           N
ATOM      2  CA  ALA A   4      36.932  14.752   2.201  1.00  0.00           C
ATOM      3  C   ALA A   4      38.127  14.604   3.136  1.00  0.00           C
ATOM      4  O   ALA A   4      38.248  15.329   4.124  1.00  0.00           O
ATOM      5  CB  ALA A   4      36.812  16.192   1.723  1.00  0.00           C
ATOM      1  N   ASP A   5      39.007  13.660   2.818  1.00  0.00           N
ATOM      2  CA  ASP A   5      40.194  13.415   3.630  1.00  0.00           C
ATOM      3  C   ASP A   5      41.467  13.708   2.841  1.00  0.00           C
ATOM      4  O   ASP A   5      41.801  12.995   1.896  1.00  0.00           O
ATOM      5  CB  ASP A   5      40.211  11.966   4.122  1.00  0.00           C
ATOM      6  CG  ASP A   5      41.346  11.691   5.089  1.00  0.00           C
ATOM      7  OD1 ASP A   5      41.256  12.134   6.254  1.00  0.00           O
ATOM      8  OD2 ASP A   5      42.327  11.032   4.685  1.00  0.00           O
ATOM      1  N   ALA A   6      42.172  14.763   3.238  1.00  0.00           N
ATOM      2  CA  ALA A   6      43.409  15.152   2.570  1.00  0.00           C
ATOM      3  C   ALA A   6      44.601  15.036   3.514  1.00  0.00           C
ATOM      4  O   ALA A   6      44.722  15.797   4.474  1.00  0.00           O
ATOM      5  CB  ALA A   6      43.294  16.573   2.039  1.00  0.00           C
ATOM      1  N   GLU A   7      45.480  14.079   3.234  1.00  0.00           N
ATOM      2  CA  GLU A   7      46.665  13.862   4.057  1.00  0.00           C
ATOM      3  C   GLU A   7      47.940  14.122   3.261  1.00  0.00           C
ATOM      4  O   GLU A   7      48.275  13.373   2.344  1.00  0.00           O
ATOM      5  CB  GLU A   7      46.677  12.432   4.604  1.00  0.00           C
ATOM      6  CG  GLU A   7      45.565  12.140   5.599  1.00  0.00           C
ATOM      7  CD  GLU A   7      45.595  10.711   6.103  1.00  0.00           C
ATOM      8  OE1 GLU A   7      46.403   9.912   5.585  1.00  0.00           O
ATOM      9  OE2 GLU A   7      44.809  10.384   7.019  1.00  0.00           O
ATOM      1  N   ALA A   8      48.647  15.189   3.620  1.00  0.00           N
ATOM      2  CA  ALA A   8      49.886  15.550   2.941  1.00  0.00           C
ATOM      3  C   ALA A   8      51.076  15.468   3.892  1.00  0.00           C
ATOM      4  O   ALA A   8      51.196  16.264   4.823  1.00  0.00           O
ATOM      5  CB  ALA A   8      49.776  16.951   2.356  1.00  0.00           C
ATOM      1  N   ALA A  10      55.122  15.615   4.002  1.00  0.00           N
ATOM      2  CA  ALA A  10      56.363  15.948   3.313  1.00  0.00           C
ATOM      3  C   ALA A  10      57.551  15.898   4.269  1.00  0.00           C
ATOM      4  O   ALA A  10      57.671  16.728   5.170  1.00  0.00           O
ATOM      5  CB  ALA A  10      56.258  17.326   2.676  1.00  0.00           C
ATOM      1  N   ASN A  11      58.427  14.919   4.065  1.00  0.00           N
ATOM      2  CA  ASN A  11      59.606  14.759   4.908  1.00  0.00           C
ATOM      3  C   ASN A  11      60.886  14.953   4.102  1.00  0.00           C
ATOM      4  O   ASN A  11      61.222  14.136   3.244  1.00  0.00           O
ATOM      5  CB  ASN A  11      59.609  13.379   5.562  1.00  0.00           C
ATOM      6  CG  ASN A  11      58.532  13.236   6.620  1.00  0.00           C
ATOM      7  OD1 ASN A  11      58.296  14.149   7.410  1.00  0.00           O
ATOM      8  ND2 ASN A  11      57.872  12.083   6.640  1.00  0.00           N
ATOM      1  N   ALA A  12      61.597  16.041   4.383  1.00  0.00           N
ATOM      2  CA  ALA A  12      62.841  16.345   3.686  1.00  0.00           C
ATOM      3  C   ALA A  12      64.025  16.328   4.646  1.00  0.00           C
ATOM      4  O   ALA A  12      64.145  17.191   5.515  1.00  0.00           O
ATOM      5  CB  ALA A  12      62.740  17.698   2.997  1.00  0.00           C
ATOM      1  N   GLN A  13      64.899  15.340   4.481  1.00  0.00           N
ATOM      2  CA  GLN A  13      66.076  15.209   5.332  1.00  0.00           C
ATOM      3  C   GLN A  13      67.359  15.370   4.522  1.00  0.00           C
ATOM      4  O   GLN A  13      67.695  14.521   3.697  1.00  0.00           O
ATOM      5  CB  GLN A  13      66.071  13.849   6.037  1.00  0.00           C
ATOM      6  CG  GLN A  13      67.212  13.651   7.023  1.00  0.00           C
ATOM      7  CD  GLN A  13      67.140  12.317   7.739  1.00  0.00           C
ATOM      8  OE1 GLN A  13      66.251  11.506   7.477  1.00  0.00           O
ATOM      9  NE2 GLN A  13      68.078  12.082   8.650  1.00  0.00           N
ATOM      1  N   ALA A  14      68.071  16.466   4.765  1.00  0.00           N
ATOM      2  CA  ALA A  14      69.318  16.740   4.059  1.00  0.00           C
ATOM      3  C   ALA A  14      70.500  16.757   5.022  1.00  0.00           C
ATOM      4  O   ALA A  14      70.620  17.652   5.859  1.00  0.00           O
ATOM      5  CB  ALA A  14      69.222  18.067   3.320  1.00  0.00           C
ATOM      1  N   LEU A  15      71.372  15.761   4.897  1.00  0.00           N
ATOM      2  CA  LEU A  15      72.547  15.660   5.755  1.00  0.00           C
ATOM      3  C   LEU A  15      73.832  15.788   4.943  1.00  0.00           C
ATOM      4  O   LEU A  15      74.168  14.907   4.151  1.00  0.00           O
ATOM      5  CB  LEU A  15      72.541  14.325   6.508  1.00  0.00           C
ATOM      6  CG  LEU A  15      71.415  14.114   7.526  1.00  0.00           C
ATOM      7  CD1 LEU A  15      71.462  12.699   8.081  1.00  0.00           C
ATOM      8  CD2 LEU A  15      71.487  15.136   8.654  1.00  0.00           C
ATOM      1  N   ALA A  16      74.546  16.890   5.146  1.00  0.00           N
ATOM      2  CA  ALA A  16      75.795  17.135   4.434  1.00  0.00           C
ATOM      3  C   ALA A  16      76.975  17.185   5.399  1.00  0.00           C
ATOM      4  O   ALA A  16      77.095  18.110   6.202  1.00  0.00           O
ATOM      5  CB  ALA A  16      75.704  18.434   3.646  1.00  0.00           C
ATOM      1  N   PHE A  17      77.845  16.184   5.313  1.00  0.00           N
ATOM      2  CA  PHE A  17      79.018  16.111   6.177  1.00  0.00           C
ATOM      3  C   PHE A  17      80.304  16.206   5.364  1.00  0.00           C
ATOM      4  O   PHE A  17      80.640  15.296   4.606  1.00  0.00           O
ATOM      5  CB  PHE A  17      79.005  14.807   6.980  1.00  0.00           C
ATOM      6  CG  PHE A  17      77.889  14.721   7.981  1.00  0.00           C
ATOM      7  CD1 PHE A  17      77.992  15.351   9.210  1.00  0.00           C
ATOM      8  CD2 PHE A  17      76.735  14.009   7.693  1.00  0.00           C
ATOM      9  CE1 PHE A  17      76.966  15.272  10.133  1.00  0.00           C
ATOM     10  CE2 PHE A  17      75.706  13.927   8.612  1.00  0.00           C
ATOM     11  CZ  PHE A  17      75.822  14.560   9.833  1.00  0.00           C
ATOM      1  N   ALA A  18      81.021  17.314   5.528  1.00  0.00           N
ATOM      2  CA  ALA A  18      82.272  17.529   4.810  1.00  0.00           C
ATOM      3  C   ALA A  18      83.450  17.613   5.775  1.00  0.00           C
ATOM      4  O   ALA A  18      83.570  18.567   6.543  1.00  0.00           O
ATOM      5  CB  ALA A  18      82.186  18.797   3.973  1.00  0.00           C
ATOM      1  N   TYR A  19      84.318  16.606   5.729  1.00  0.00           N
ATOM      2  CA  TYR A  19      85.488  16.564   6.598  1.00  0.00           C
ATOM      3  C   TYR A  19      86.777  16.625   5.785  1.00  0.00           C
ATOM      4  O   TYR A  19      87.113  15.687   5.063  1.00  0.00           O
ATOM      5  CB  TYR A  19      85.471  15.292   7.450  1.00  0.00           C
ATOM      6  CG  TYR A  19      84.363  15.258   8.479  1.00  0.00           C
ATOM      7  CD1 TYR A  19      83.662  16.411   8.812  1.00  0.00           C
ATOM      8  CD2 TYR A  19      84.016  14.074   9.116  1.00  0.00           C
ATOM      9  CE1 TYR A  19      82.648  16.386   9.751  1.00  0.00           C
ATOM     10  CE2 TYR A  19      83.002  14.039  10.057  1.00  0.00           C
ATOM     11  CZ  TYR A  19      82.323  15.197  10.369  1.00  0.00           C
ATOM     12  OH  TYR A  19      81.313  15.166  11.305  1.00  0.00           O
ATOM      1  N   ALA A  20      87.496  17.737   5.909  1.00  0.00           N
ATOM      2  CA  ALA A  20      88.749  17.922   5.187  1.00  0.00           C
ATOM      3  C   ALA A  20      89.925  18.039   6.151  1.00  0.00           C
ATOM      4  O   ALA A  20      90.046  19.021   6.883  1.00  0.00           O
ATOM      5  CB  ALA A  20      88.668  19.159   4.303  1.00  0.00           C
ATOM      1  N   VAL A  21      90.791  17.030   6.145  1.00  0.00           N
ATOM      2  CA  VAL A  21      91.959  17.017   7.018  1.00  0.00           C
ATOM      3  C   VAL A  21      93.250  17.045   6.207  1.00  0.00           C
ATOM      4  O   VAL A  21      93.585  16.079   5.521  1.00  0.00           O
ATOM      5  CB  VAL A  21      91.967  15.769   7.925  1.00  0.00           C
ATOM      6  CG1 VAL A  21      93.241  15.722   8.760  1.00  0.00           C
ATOM      7  CG2 VAL A  21      90.735  15.749   8.820  1.00  0.00           C
ATOM      1  N   ALA A  22      93.971  18.159   6.291  1.00  0.00           N
ATOM      2  CA  ALA A  22      95.226  18.315   5.565  1.00  0.00           C
ATOM      3  C   ALA A  22      96.400  18.465   6.527  1.00  0.00           C
ATOM      4  O   ALA A  22      96.521  19.473   7.222  1.00  0.00           O
ATOM      5  CB  ALA A  22      95.150  19.517   4.636  1.00  0.00           C
TER     156      ALA A  22
ATOM      1  N   ARG B   1      27.961   0.504   1.988  1.00  0.00           N
ATOM      2  CA  ARG B   1      29.153   0.205   2.773  1.00  0.00           C
ATOM      3  C   ARG B   1      30.420   0.562   2.003  1.00  0.00           C
ATOM      4  O   ARG B   1      30.753  -0.077   1.005  1.00  0.00           O
ATOM      5  CB  ARG B   1      29.180  -1.279   3.156  1.00 10.00           C
ATOM      6  CG  ARG B   1      28.007  -1.729   4.022  1.00 10.00           C
ATOM      7  CD  ARG B   1      28.018  -1.090   5.409  1.00 10.00           C
ATOM      8  NE  ARG B   1      29.183  -1.492   6.201  1.00 10.00           N
ATOM      9  CZ  ARG B   1      29.635  -0.879   7.298  1.00 10.00           C
ATOM     10  NH2 ARG B   1      29.046   0.205   7.803  1.00 10.00           N
ATOM     11  NH1 ARG B   1      30.708  -1.364   7.908  1.00 10.00           N
ATOM      1  N   ALA B   2      31.123   1.587   2.474  1.00  0.00           N
ATOM      2  CA  ALA B   2      32.355   2.031   1.832  1.00  0.00           C
ATOM      3  C   ALA B   2      33.552   1.851   2.758  1.00  0.00           C
ATOM      4  O   ALA B   2      33.675   2.539   3.772  1.00  0.00           O
ATOM      5  CB  ALA B   2      32.231   3.487   1.408  1.00  0.00           C
ATOM      1  N   HIS B   3      34.434   0.922   2.403  1.00  0.00           N
ATOM      2  CA  HIS B   3      35.624   0.650   3.202  1.00  0.00           C
ATOM      3  C   HIS B   3      36.893   0.975   2.422  1.00  0.00           C
ATOM      4  O   HIS B   3      37.227   0.298   1.450  1.00  0.00           O
ATOM      5  CB  HIS B   3      35.644  -0.817   3.640  1.00  0.00           C
ATOM      6  CG  HIS B   3      34.518  -1.190   4.554  1.00  0.00           C
ATOM      7  CD2 HIS B   3      34.486  -0.826   5.882  1.00  0.00           N
ATOM      8  ND1 HIS B   3      33.385  -1.896   4.330  1.00  0.00           C
ATOM      9  NE2 HIS B   3      33.381  -1.291   6.437  1.00  0.00           C
ATOM     10  CE1 HIS B   3      32.695  -1.945   5.517  1.00  0.00           N
ATOM      1  N   ALA B   4      37.598   2.015   2.856  1.00  0.00           N
ATOM      2  CA  ALA B   4      38.832   2.432   2.201  1.00  0.00           C
ATOM      3  C   ALA B   4      40.027   2.284   3.136  1.00  0.00           C
ATOM      4  O   ALA B   4      40.148   3.009   4.124  1.00  0.00           O
ATOM      5  CB  ALA B   4      38.712   3.872   1.723  1.00  0.00           C
ATOM      1  N   ASP B   5      40.907   1.340   2.818  1.00  0.00           N
ATOM      2  CA  ASP B   5      42.094   1.095   3.630  1.00  0.00           C
ATOM      3  C   ASP B   5      43.367   1.388   2.841  1.00  0.00           C
ATOM      4  O   ASP B   5      43.701   0.675   1.896  1.00  0.00           O
ATOM      5  CB  ASP B   5      42.111  -0.354   4.122  1.00  0.00           C
ATOM      6  CG  ASP B   5      43.246  -0.629   5.089  1.00  0.00           C
ATOM      7  OD2 ASP B   5      43.156  -0.186   6.254  1.00  0.00           O
ATOM      8  OD1 ASP B   5      44.227  -1.288   4.685  1.00  0.00           O
ATOM      1  N   ALA B   6      44.072   2.443   3.238  1.00  0.00           N
ATOM      2  CA  ALA B   6      45.309   2.832   2.570  1.00  0.00           C
ATOM      3  C   ALA B   6      46.501   2.716   3.514  1.00  0.00           C
ATOM      4  O   ALA B   6      46.622   3.477   4.474  1.00  0.00           O
ATOM      5  CB  ALA B   6      45.194   4.253   2.039  1.00  0.00           C
ATOM      1  N   GLU B   7      47.380   1.759   3.234  1.00  0.00           N
ATOM      2  CA  GLU B   7      48.565   1.542   4.057  1.00  0.00           C
ATOM      3  C   GLU B   7      49.840   1.802   3.261  1.00  0.00           C
ATOM      4  O   GLU B   7      50.175   1.053   2.344  1.00  0.00           O
ATOM      5  CB  GLU B   7      48.577   0.112   4.604  1.00  0.00           C
ATOM      6  CG  GLU B   7      47.465  -0.180   5.599  1.00  0.00           C
ATOM      7  CD  GLU B   7      47.495  -1.609   6.103  1.00  0.00           C
ATOM      8  OE2 GLU B   7      48.303  -2.408   5.585  1.00  0.00           O
ATOM      9  OE1 GLU B   7      46.709  -1.936   7.019  1.00  0.00           O
ATOM      1  N   ALA B   8      50.547   2.869   3.620  1.00  0.00           N
ATOM      2  CA  ALA B   8      51.786   3.230   2.941  1.00  0.00           C
ATOM      3  C   ALA B   8      52.976   3.148   3.892  1.00  0.00           C
ATOM      4  O   ALA B   8      53.096   3.944   4.823  1.00  0.00           O
ATOM      5  CB  ALA B   8      51.676   4.631   2.356  1.00  0.00           C
ATOM      1  N   ALA B  10      57.022   3.295   4.002  1.00  0.00           N
ATOM      2  CA  ALA B  10      58.263   3.628   3.313  1.00  0.00           C
ATOM      3  C   ALA B  10      59.451   3.578   4.269  1.00  0.00           C
ATOM      4  O   ALA B  10      59.571   4.408   5.170  1.00  0.00           O
ATOM      5  CB  ALA B  10      58.158   5.006   2.676  1.00  0.00           C
ATOM      1  N   ASN B  11      60.327   2.599   4.065  1.00  0.00           N
ATOM      2  CA  ASN B  11      61.506   2.439   4.908  1.00  0.00           C
ATOM      3  C   ASN B  11      62.786   2.633   4.102  1.00  0.00           C
ATOM      4  O   ASN B  11      63.122   1.816   3.244  1.00  0.00           O
ATOM      5  CB  ASN B  11      61.509   1.059   5.562  1.00  0.00           C
ATOM      6  CG  ASN B  11      60.432   0.916   6.620  1.00  0.00           C
ATOM      7  ND2 ASN B  11      60.196   1.829   7.410  1.00  0.00           O
ATOM      8  OD1 ASN B  11      59.772  -0.237   6.640  1.00  0.00           N
ATOM      1  N   ALA B  12      63.497   3.721   4.383  1.00  0.00           N
ATOM      2  CA  ALA B  12      64.741   4.025   3.686  1.00  0.00           C
ATOM      3  C   ALA B  12      65.925   4.008   4.646  1.00  0.00           C
ATOM      4  O   ALA B  12      66.045   4.871   5.515  1.00  0.00           O
ATOM      5  CB  ALA B  12      64.640   5.378   2.997  1.00  0.00           C
ATOM      1  N   GLN B  13      66.799   3.020   4.481  1.00  0.00           N
ATOM      2  CA  GLN B  13      67.976   2.889   5.332  1.00  0.00           C
ATOM      3  C   GLN B  13      69.259   3.050   4.522  1.00  0.00           C
ATOM      4  O   GLN B  13      69.595   2.201   3.697  1.00  0.00           O
ATOM      5  CB  GLN B  13      67.971   1.529   6.037  1.00  0.00           C
ATOM      6  CG  GLN B  13      69.112   1.331   7.023  1.00  0.00           C
ATOM      7  CD  GLN B  13      69.040  -0.003   7.739  1.00  0.00           C
ATOM      8  NE2 GLN B  13      68.151  -0.814   7.477  1.00  0.00           O
ATOM      9  OE1 GLN B  13      69.978  -0.238   8.650  1.00  0.00           N
ATOM      1  N   ALA B  14      69.971   4.146   4.765  1.00  0.00           N
ATOM      2  CA  ALA B  14      71.218   4.420   4.059  1.00  0.00           C
ATOM      3  C   ALA B  14      72.400   4.437   5.022  1.00  0.00           C
ATOM      4  O   ALA B  14      72.520   5.332   5.859  1.00  0.00           O
ATOM      5  CB  ALA B  14      71.122   5.747   3.320  1.00  0.00           C
ATOM      1  N   LEU B  15      73.272   3.441   4.897  1.00  0.00           N
ATOM      2  CA  LEU B  15      74.447   3.340   5.755  1.00  0.00           C
ATOM      3  C   LEU B  15      75.732   3.468   4.943  1.00  0.00           C
ATOM      4  O   LEU B  15      76.068   2.587   4.151  1.00  0.00           O
ATOM      5  CB  LEU B  15      74.441   2.005   6.508  1.00  0.00           C
ATOM      6  CG  LEU B  15      73.315   1.794   7.526  1.00  0.00           C
ATOM      7  CD2 LEU B  15      73.362   0.379   8.081  1.00  0.00           C
ATOM      8  CD1 LEU B  15      73.387   2.816   8.654  1.00  0.00           C
ATOM      1  N   ALA B  16      76.446   4.570   5.146  1.00  0.00           N
ATOM      2  CA  ALA B  16      77.695   4.815   4.434  1.00  0.00           C
ATOM      3  C   ALA B  16      78.875   4.865   5.399  1.00  0.00           C
ATOM      4  O   ALA B  16      78.995   5.790   6.202  1.00  0.00           O
ATOM      5  CB  ALA B  16      77.604   6.114   3.646  1.00  0.00           C
ATOM      1  N   PHE B  17      79.745   3.864   5.313  1.00  0.00           N
ATOM      2  CA  PHE B  17      80.918   3.791   6.177  1.00  0.00           C
ATOM      3  C   PHE B  17      82.204   3.886   5.364  1.00  0.00           C
ATOM      4  O   PHE B  17      82.540   2.976   4.606  1.00  0.00           O
ATOM      5  CB  PHE B  17      80.905   2.487   6.980  1.00  0.00           C
ATOM      6  CG  PHE B  17      79.789   2.401   7.981  1.00  0.00           C
ATOM      7  CD2 PHE B  17      79.892   3.031   9.210  1.00  0.00           C
ATOM      8  CD1 PHE B  17      78.635   1.689   7.693  1.00  0.00           C
ATOM      9  CE2 PHE B  17      78.866   2.952  10.133  1.00  0.00           C
ATOM     10  CE1 PHE B  17      77.606   1.607   8.612  1.00  0.00           C
ATOM     11  CZ  PHE B  17      77.722   2.240   9.833  1.00  0.00           C
ATOM      1  N   ALA B  18      82.921   4.994   5.528  1.00  0.00           N
ATOM      2  CA  ALA B  18      84.172   5.209   4.810  1.00  0.00           C
ATOM      3  C   ALA B  18      85.350   5.293   5.775  1.00  0.00           C
ATOM      4  O   ALA B  18      85.470   6.247   6.543  1.00  0.00           O
ATOM      5  CB  ALA B  18      84.086   6.477   3.973  1.00  0.00           C
ATOM      1  N   TYR B  19      86.218   4.286   5.729  1.00  0.00           N
ATOM      2  CA  TYR B  19      87.388   4.244   6.598  1.00  0.00           C
ATOM      3  C   TYR B  19      88.677   4.305   5.785  1.00  0.00           C
ATOM      4  O   TYR B  19      89.013   3.367   5.063  1.00  0.00           O
ATOM      5  CB  TYR B  19      87.371   2.972   7.450  1.00  0.00           C
ATOM      6  CG  TYR B  19      86.263   2.938   8.479  1.00  0.00           C
ATOM      7  CD2 TYR B  19      85.562   4.091   8.812  1.00  0.00           C
ATOM      8  CD1 TYR B  19      85.916   1.754   9.116  1.00  0.00           C
ATOM      9  CE2 TYR B  19      84.548   4.066   9.751  1.00  0.00           C
ATOM     10  CE1 TYR B  19      84.902   1.719  10.057  1.00  0.00           C
ATOM     11  CZ  TYR B  19      84.223   2.877  10.369  1.00  0.00           C
ATOM     12  OH  TYR B  19      83.213   2.846  11.305  1.00  0.00           O
ATOM      1  N   ALA B  20      89.396   5.417   5.909  1.00  0.00           N
ATOM      2  CA  ALA B  20      90.649   5.602   5.187  1.00  0.00           C
ATOM      3  C   ALA B  20      91.825   5.719   6.151  1.00  0.00           C
ATOM      4  O   ALA B  20      91.946   6.701   6.883  1.00  0.00           O
ATOM      5  CB  ALA B  20      90.568   6.839   4.303  1.00  0.00           C
ATOM      1  N   VAL B  21      92.691   4.710   6.145  1.00  0.00           N
ATOM      2  CA  VAL B  21      93.859   4.697   7.018  1.00  0.00           C
ATOM      3  C   VAL B  21      95.150   4.725   6.207  1.00  0.00           C
ATOM      4  O   VAL B  21      95.485   3.759   5.521  1.00  0.00           O
ATOM      5  CB  VAL B  21      93.867   3.449   7.925  1.00  0.00           C
ATOM      6  CG2 VAL B  21      95.141   3.402   8.760  1.00  0.00           C
ATOM      7  CG1 VAL B  21      92.635   3.429   8.820  1.00  0.00           C
ATOM      1  N   ALA B  22      95.871   5.839   6.291  1.00  0.00           N
ATOM      2  CA  ALA B  22      97.126   5.995   5.565  1.00  0.00           C
ATOM      3  C   ALA B  22      98.300   6.145   6.527  1.00  0.00           C
ATOM      4  O   ALA B  22      98.421   7.153   7.222  1.00  0.00           O
ATOM      5  CB  ALA B  22      97.050   7.197   4.636  1.00  0.00           C
TER     312      ALA B  22
END
"""


def test_1():
  h = iotbx.pdb.input(source_info=None, lines=test_pdb_1).construct_hierarchy()
  asc = h.atom_selection_cache()
  ncs_inp = iotbx.ncs.input(
      hierarchy=h,
      chain_max_rmsd=0.01)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups) == 1
  # group 1
  assert ncs_groups[0].master_iselection.all_eq(
    asc.selection(string = "chain A").iselection())
  g1_c = ncs_groups[0].copies
  assert len(g1_c)==1
  assert g1_c[0].iselection.all_eq(
    asc.selection(string = "chain B").iselection())

def test_2():
  """ Testing TYR both CD and CE flips"""
  h = iotbx.pdb.input(source_info=None, lines=test_pdb_2).construct_hierarchy()
  asc = h.atom_selection_cache()
  ncs_inp = iotbx.ncs.input(
      hierarchy=h,
      chain_max_rmsd=0.01)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups) == 1
  # group 1
  assert ncs_groups[0].master_iselection.all_eq(
    asc.selection(string = "chain A").iselection())
  g1_c = ncs_groups[0].copies
  assert len(g1_c)==1
  assert g1_c[0].iselection.all_eq(
    asc.selection(string = "chain B").iselection())

def test_3():
  """ Testing TYR only CD and flip. Not sure if such flip is valid
      This test is not working when we use torsion angles because in this
      procedure we don't flip only one atom out of two like in this case"""
  h = iotbx.pdb.input(source_info=None, lines=test_pdb_3).construct_hierarchy()
  asc = h.atom_selection_cache()
  ncs_inp = iotbx.ncs.input(
      hierarchy=h,
      chain_max_rmsd=0.01)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups) == 1
  # group 1
  assert ncs_groups[0].master_iselection.all_eq(
    asc.selection(string = "chain A").iselection())
  g1_c = ncs_groups[0].copies
  assert len(g1_c)==1
  assert g1_c[0].iselection.all_eq(
    asc.selection(string = "chain B").iselection())

def test_4():
  """ HIS is interesting because of different atoms can be flipped """
  h = iotbx.pdb.input(source_info=None, lines=test_pdb_4).construct_hierarchy()
  asc = h.atom_selection_cache()
  ncs_inp = iotbx.ncs.input(
      hierarchy=h,
      chain_max_rmsd=0.01)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups) == 1
  # group 1
  assert ncs_groups[0].master_iselection.all_eq(
    asc.selection(string = "chain A").iselection())
  g1_c = ncs_groups[0].copies
  assert len(g1_c)==1
  assert g1_c[0].iselection.all_eq(
    asc.selection(string = "chain B").iselection())

def test_5():
  """ Testing all possible residues at once. All flipped wherever possible
  in chain B """
  h = iotbx.pdb.input(source_info=None, lines=test_pdb_5).construct_hierarchy()
  asc = h.atom_selection_cache()
  ncs_inp = iotbx.ncs.input(
      hierarchy=h,
      chain_max_rmsd=0.01)
  ncs_groups = ncs_inp.get_ncs_restraints_group_list()
  assert len(ncs_groups) == 1
  # group 1
  assert ncs_groups[0].master_iselection.all_eq(
    asc.selection(string = "chain A").iselection())
  g1_c = ncs_groups[0].copies
  assert len(g1_c)==1
  assert g1_c[0].iselection.all_eq(
    asc.selection(string = "chain B").iselection())


if __name__=='__main__':
  test_1()
  test_2()
  # test_3()
  test_4()
  test_5()
  print "OK"
