"""reduce2: the description's final scores are fresh scores at the final poses,
and no multi-Mover Clique ends below its score after the coarse step.

The search stores each Mover's score when it sets the Mover's state; a Clique
partner that moves afterwards leaves the stored score stale. The description
prints each Mover's score at its final pose against everything around it at the
end, and each Clique's total is the sum of these. A fine or polish move of a
Clique member is accepted only if it raises the Clique's total.

Test:
- 4wsc Lys C 245 NZ with Gln C 319 (and Phe C 219): every printed score equals a
  fresh rescoring, the Clique total their sum, and the Clique's final total is at
  least its total after the coarse step (the "Clique of size N optimized with
  score X" line); the Lys's fine move lowered it before the Clique-total gate.
- 6vw1 Thr A 276 HG1 with Thr A 445 HG1: Thr A 276's stored score is stale
  (15.86, fresh 15.04) after Thr A 445 moves in the polish; every printed score
  equals a fresh rescoring.
"""
from __future__ import absolute_import, division, print_function
from io import StringIO
from iotbx.cli_parser import run_program
from libtbx.test_utils import approx_equal
from mmtbx.programs import reduce2
from mmtbx.probe import Helpers
from mmtbx.reduce import Optimizers

# 4wsc Phe C 219, Lys C 245, Gln C 319
lys_gln_str = """\
CRYST1  135.640  259.390  143.310  90.00 105.12  90.00 P 1 21 1
ATOM   9281  N   PHE C 219      -5.510 -12.352 -31.231  1.00 47.10           N
ATOM   9282  CA  PHE C 219      -4.604 -13.105 -32.094  1.00 49.36           C
ATOM   9283  C   PHE C 219      -3.693 -13.996 -31.243  1.00 37.63           C
ATOM   9284  O   PHE C 219      -3.391 -13.690 -30.090  1.00 31.68           O
ATOM   9285  CB  PHE C 219      -3.782 -12.165 -32.969  1.00 50.42           C
ATOM   9286  CG  PHE C 219      -4.555 -11.593 -34.106  1.00 53.72           C
ATOM   9287  CD1 PHE C 219      -4.933 -12.388 -35.169  1.00 55.35           C
ATOM   9288  CD2 PHE C 219      -4.920 -10.256 -34.105  1.00 64.45           C
ATOM   9289  CE1 PHE C 219      -5.657 -11.861 -36.212  1.00 60.53           C
ATOM   9290  CE2 PHE C 219      -5.641  -9.719 -35.148  1.00 68.56           C
ATOM   9291  CZ  PHE C 219      -6.014 -10.524 -36.203  1.00 66.84           C
ATOM   9479  N   LYS C 245      -6.780  -5.665 -31.988  1.00 62.07           N
ATOM   9480  CA  LYS C 245      -7.180  -6.883 -31.275  1.00 52.07           C
ATOM   9481  C   LYS C 245      -6.128  -7.246 -30.233  1.00 47.04           C
ATOM   9482  O   LYS C 245      -4.956  -6.933 -30.426  1.00 47.22           O
ATOM   9483  CB  LYS C 245      -7.366  -8.039 -32.257  1.00 50.07           C
ATOM   9484  CG  LYS C 245      -8.591  -7.921 -33.151  1.00 58.18           C
ATOM   9485  CD  LYS C 245      -8.826  -9.234 -33.875  1.00 62.94           C
ATOM   9486  CE  LYS C 245     -10.159  -9.269 -34.601  1.00 64.88           C
ATOM   9487  NZ  LYS C 245     -10.397 -10.602 -35.241  1.00 63.79           N
ATOM  10021  N   GLN C 319      -6.118 -17.653 -34.053  1.00 55.81           N
ATOM  10022  CA  GLN C 319      -7.173 -16.807 -33.532  1.00 52.62           C
ATOM  10023  C   GLN C 319      -8.283 -17.609 -32.837  1.00 56.71           C
ATOM  10024  O   GLN C 319      -8.580 -18.742 -33.220  1.00 63.19           O
ATOM  10025  CB  GLN C 319      -7.714 -15.979 -34.690  1.00 51.64           C
ATOM  10026  CG  GLN C 319      -8.318 -14.668 -34.292  1.00 67.66           C
ATOM  10027  CD  GLN C 319      -8.460 -13.739 -35.479  1.00 65.76           C
ATOM  10028  OE1 GLN C 319      -9.202 -12.747 -35.429  1.00 62.32           O
ATOM  10029  NE2 GLN C 319      -7.750 -14.061 -36.564  1.00 54.15           N
END
"""

# 6vw1 Thr A 276, Thr A 445 and everything within 6 A (also in
# tst_reduce2_clique_polish)
pdb_6vw1 = """\
CRYST1   80.435  118.034  112.075  90.00  93.12  90.00 P 1 21 1
ATOM   2029  N   LEU A 267     119.568   8.744 154.595  1.00 78.43           N
ATOM   2030  CA  LEU A 267     118.168   8.844 154.204  1.00 79.08           C
ATOM   2031  C   LEU A 267     117.262   9.428 155.282  1.00 79.01           C
ATOM   2032  O   LEU A 267     116.093   9.705 154.995  1.00 78.08           O
ATOM   2033  CB  LEU A 267     117.643   7.467 153.783  1.00 73.61           C
ATOM   2034  CG  LEU A 267     118.281   6.930 152.501  1.00 74.06           C
ATOM   2035  CD1 LEU A 267     117.674   5.592 152.111  1.00 78.59           C
ATOM   2036  CD2 LEU A 267     118.134   7.944 151.375  1.00 67.22           C
ATOM   2037  N   GLY A 268     117.758   9.627 156.501  1.00 77.56           N
ATOM   2038  CA  GLY A 268     117.000  10.271 157.555  1.00 80.74           C
ATOM   2039  C   GLY A 268     116.670   9.391 158.742  1.00 76.85           C
ATOM   2040  O   GLY A 268     116.356   9.921 159.816  1.00 74.81           O
ATOM   2041  N   ASP A 269     116.721   8.070 158.589  1.00 74.86           N
ATOM   2042  CA  ASP A 269     116.594   7.170 159.728  1.00 78.44           C
ATOM   2043  C   ASP A 269     117.506   5.969 159.497  1.00 78.01           C
ATOM   2044  O   ASP A 269     118.305   5.938 158.556  1.00 81.00           O
ATOM   2045  CB  ASP A 269     115.125   6.783 159.971  1.00 79.45           C
ATOM   2046  CG  ASP A 269     114.490   6.042 158.796  1.00 82.11           C
ATOM   2047  OD1 ASP A 269     115.084   5.072 158.276  1.00 77.35           O
ATOM   2048  OD2 ASP A 269     113.366   6.424 158.406  1.00 80.79           O
ATOM   2075  N   ARG A 273     113.301   4.016 155.752  1.00 79.96           N
ATOM   2076  CA  ARG A 273     111.870   3.753 155.673  1.00 79.38           C
ATOM   2077  C   ARG A 273     111.081   5.038 155.468  1.00 78.54           C
ATOM   2078  O   ARG A 273     110.056   5.033 154.779  1.00 78.91           O
ATOM   2079  CB  ARG A 273     111.400   3.035 156.937  1.00 81.16           C
ATOM   2080  CG  ARG A 273     109.928   2.671 156.935  1.00 86.03           C
ATOM   2081  CD  ARG A 273     109.513   2.057 158.260  1.00 85.54           C
ATOM   2082  NE  ARG A 273     108.071   1.850 158.333  1.00 85.36           N
ATOM   2083  CZ  ARG A 273     107.406   1.619 159.461  1.00 83.49           C
ATOM   2084  NH1 ARG A 273     106.093   1.442 159.433  1.00 85.52           N
ATOM   2085  NH2 ARG A 273     108.054   1.569 160.617  1.00 82.53           N
ATOM   2086  N   PHE A 274     111.540   6.142 156.051  1.00 78.20           N
ATOM   2087  CA  PHE A 274     110.942   7.453 155.841  1.00 78.96           C
ATOM   2088  C   PHE A 274     112.051   8.464 155.597  1.00 78.08           C
ATOM   2089  O   PHE A 274     113.053   8.479 156.318  1.00 80.58           O
ATOM   2090  CB  PHE A 274     110.090   7.885 157.042  1.00 76.81           C
ATOM   2091  CG  PHE A 274     109.010   6.906 157.407  1.00 69.52           C
ATOM   2092  CD1 PHE A 274     107.879   6.776 156.617  1.00 69.53           C
ATOM   2093  CD2 PHE A 274     109.124   6.117 158.540  1.00 63.84           C
ATOM   2094  CE1 PHE A 274     106.883   5.877 156.948  1.00 68.02           C
ATOM   2095  CE2 PHE A 274     108.132   5.216 158.876  1.00 65.03           C
ATOM   2096  CZ  PHE A 274     107.011   5.094 158.078  1.00 66.42           C
ATOM   2097  N   TRP A 275     111.871   9.305 154.577  1.00 82.98           N
ATOM   2098  CA  TRP A 275     112.779  10.413 154.309  1.00 85.43           C
ATOM   2099  C   TRP A 275     112.342  11.698 155.001  1.00 94.04           C
ATOM   2100  O   TRP A 275     112.654  12.797 154.519  1.00102.35           O
ATOM   2101  CB  TRP A 275     112.909  10.635 152.801  1.00 76.67           C
ATOM   2102  CG  TRP A 275     113.405   9.438 152.040  1.00 76.21           C
ATOM   2103  CD1 TRP A 275     113.829   8.247 152.560  1.00 74.63           C
ATOM   2104  CD2 TRP A 275     113.536   9.321 150.619  1.00 80.52           C
ATOM   2105  NE1 TRP A 275     114.208   7.397 151.551  1.00 78.04           N
ATOM   2106  CE2 TRP A 275     114.039   8.033 150.349  1.00 81.74           C
ATOM   2107  CE3 TRP A 275     113.275  10.181 149.548  1.00 80.17           C
ATOM   2108  CZ2 TRP A 275     114.287   7.585 149.053  1.00 80.37           C
ATOM   2109  CZ3 TRP A 275     113.521   9.734 148.263  1.00 81.75           C
ATOM   2110  CH2 TRP A 275     114.021   8.448 148.026  1.00 79.70           C
ATOM   2111  N   THR A 276     111.618  11.581 156.119  1.00 91.54           N
ATOM   2112  CA  THR A 276     111.064  12.751 156.794  1.00 82.49           C
ATOM   2113  C   THR A 276     112.157  13.714 157.244  1.00 76.85           C
ATOM   2114  O   THR A 276     112.013  14.934 157.111  1.00 72.40           O
ATOM   2115  CB  THR A 276     110.219  12.310 157.993  1.00 83.86           C
ATOM   2116  OG1 THR A 276     109.175  11.432 157.551  1.00 83.88           O
ATOM   2117  CG2 THR A 276     109.596  13.513 158.696  1.00 80.89           C
ATOM   2118  N   ASN A 277     113.265  13.186 157.760  1.00 82.04           N
ATOM   2119  CA  ASN A 277     114.321  14.009 158.337  1.00 85.47           C
ATOM   2120  C   ASN A 277     115.249  14.623 157.293  1.00 86.17           C
ATOM   2121  O   ASN A 277     116.277  15.200 157.666  1.00 89.38           O
ATOM   2122  CB  ASN A 277     115.127  13.184 159.342  1.00 90.87           C
ATOM   2123  CG  ASN A 277     114.263  12.629 160.459  1.00 90.01           C
ATOM   2124  OD1 ASN A 277     113.318  13.280 160.907  1.00 88.09           O
ATOM   2125  ND2 ASN A 277     114.576  11.420 160.906  1.00 91.84           N
ATOM   2126  N   LEU A 278     114.921  14.513 156.009  1.00 84.04           N
ATOM   2127  CA  LEU A 278     115.649  15.216 154.963  1.00 81.83           C
ATOM   2128  C   LEU A 278     115.017  16.553 154.607  1.00 86.57           C
ATOM   2129  O   LEU A 278     115.483  17.210 153.673  1.00 88.87           O
ATOM   2130  CB  LEU A 278     115.748  14.354 153.699  1.00 70.37           C
ATOM   2131  CG  LEU A 278     116.411  12.983 153.811  1.00 71.63           C
ATOM   2132  CD1 LEU A 278     116.425  12.298 152.455  1.00 58.68           C
ATOM   2133  CD2 LEU A 278     117.821  13.115 154.363  1.00 67.68           C
ATOM   2134  N   TYR A 279     113.966  16.965 155.321  1.00 82.58           N
ATOM   2135  CA  TYR A 279     113.246  18.175 154.940  1.00 86.94           C
ATOM   2136  C   TYR A 279     114.117  19.417 155.096  1.00 88.94           C
ATOM   2137  O   TYR A 279     114.032  20.345 154.285  1.00 93.06           O
ATOM   2138  CB  TYR A 279     111.968  18.313 155.767  1.00 83.97           C
ATOM   2139  CG  TYR A 279     111.126  19.505 155.374  1.00 86.79           C
ATOM   2140  CD1 TYR A 279     110.529  19.571 154.121  1.00 81.05           C
ATOM   2141  CD2 TYR A 279     110.937  20.568 156.247  1.00 93.91           C
ATOM   2142  CE1 TYR A 279     109.760  20.658 153.751  1.00 85.32           C
ATOM   2143  CE2 TYR A 279     110.169  21.661 155.884  1.00 96.02           C
ATOM   2144  CZ  TYR A 279     109.583  21.701 154.635  1.00 91.92           C
ATOM   2145  OH  TYR A 279     108.818  22.785 154.269  1.00 95.42           O
ATOM   2146  N   SER A 280     114.961  19.455 156.130  1.00 84.50           N
ATOM   2147  CA  SER A 280     115.843  20.603 156.314  1.00 86.28           C
ATOM   2148  C   SER A 280     116.847  20.724 155.177  1.00 91.69           C
ATOM   2149  O   SER A 280     117.305  21.829 154.865  1.00 93.47           O
ATOM   2150  CB  SER A 280     116.572  20.500 157.653  1.00 84.34           C
ATOM   2151  OG  SER A 280     117.303  19.289 157.738  1.00 85.30           O
ATOM   3401  N   LEU A 440     105.805  18.469 149.501  1.00 98.40           N
ATOM   3402  CA  LEU A 440     107.113  17.823 149.512  1.00 91.60           C
ATOM   3403  C   LEU A 440     107.384  17.155 150.853  1.00 95.17           C
ATOM   3404  O   LEU A 440     107.892  16.025 150.901  1.00 94.24           O
ATOM   3405  CB  LEU A 440     108.197  18.853 149.181  1.00 90.48           C
ATOM   3406  CG  LEU A 440     109.654  18.427 148.964  1.00 85.58           C
ATOM   3407  CD1 LEU A 440     110.437  18.352 150.270  1.00 84.95           C
ATOM   3408  CD2 LEU A 440     109.719  17.103 148.225  1.00 81.61           C
ATOM   3409  N   LYS A 441     107.059  17.837 151.956  1.00 97.61           N
ATOM   3410  CA  LYS A 441     107.247  17.213 153.262  1.00 89.44           C
ATOM   3411  C   LYS A 441     106.334  16.003 153.427  1.00 80.46           C
ATOM   3412  O   LYS A 441     106.733  14.996 154.023  1.00 69.90           O
ATOM   3413  CB  LYS A 441     107.021  18.226 154.382  1.00 91.30           C
ATOM   3414  CG  LYS A 441     107.310  17.656 155.763  1.00 95.42           C
ATOM   3415  CD  LYS A 441     107.277  18.726 156.838  1.00107.00           C
ATOM   3416  CE  LYS A 441     107.498  18.119 158.214  1.00108.69           C
ATOM   3417  NZ  LYS A 441     106.385  17.213 158.610  1.00108.68           N
ATOM   3418  N   GLN A 442     105.103  16.083 152.909  1.00 85.05           N
ATOM   3419  CA  GLN A 442     104.217  14.922 152.920  1.00 86.99           C
ATOM   3420  C   GLN A 442     104.795  13.778 152.096  1.00 92.65           C
ATOM   3421  O   GLN A 442     104.707  12.611 152.495  1.00 87.46           O
ATOM   3422  CB  GLN A 442     102.833  15.304 152.388  1.00 91.81           C
ATOM   3423  CG  GLN A 442     102.016  16.225 153.283  1.00104.52           C
ATOM   3424  CD  GLN A 442     101.480  15.524 154.517  1.00107.92           C
ATOM   3425  OE1 GLN A 442     101.640  14.314 154.678  1.00113.98           O
ATOM   3426  NE2 GLN A 442     100.831  16.281 155.391  1.00103.01           N
ATOM   3427  N   ALA A 443     105.388  14.097 150.943  1.00 95.10           N
ATOM   3428  CA  ALA A 443     105.859  13.065 150.028  1.00 97.18           C
ATOM   3429  C   ALA A 443     107.076  12.345 150.582  1.00 94.57           C
ATOM   3430  O   ALA A 443     107.187  11.119 150.450  1.00 91.27           O
ATOM   3431  CB  ALA A 443     106.182  13.677 148.667  1.00102.42           C
ATOM   3432  N   LEU A 444     107.996  13.094 151.203  1.00 91.61           N
ATOM   3433  CA  LEU A 444     109.178  12.489 151.813  1.00 84.59           C
ATOM   3434  C   LEU A 444     108.804  11.318 152.710  1.00 84.85           C
ATOM   3435  O   LEU A 444     109.522  10.314 152.770  1.00 86.90           O
ATOM   3436  CB  LEU A 444     109.954  13.536 152.613  1.00 79.71           C
ATOM   3437  CG  LEU A 444     110.720  14.615 151.848  1.00 83.12           C
ATOM   3438  CD1 LEU A 444     111.252  15.662 152.815  1.00 84.47           C
ATOM   3439  CD2 LEU A 444     111.854  14.000 151.041  1.00 81.49           C
ATOM   3440  N   THR A 445     107.675  11.426 153.409  1.00 84.77           N
ATOM   3441  CA  THR A 445     107.193  10.366 154.283  1.00 78.71           C
ATOM   3442  C   THR A 445     106.345   9.339 153.538  1.00 77.28           C
ATOM   3443  O   THR A 445     106.597   8.135 153.638  1.00 75.57           O
ATOM   3444  CB  THR A 445     106.386  10.967 155.440  1.00 76.37           C
ATOM   3445  OG1 THR A 445     107.211  11.872 156.185  1.00 75.49           O
ATOM   3446  CG2 THR A 445     105.884   9.871 156.365  1.00 72.71           C
ATOM   3447  N   ILE A 446     105.347   9.796 152.779  1.00 80.59           N
ATOM   3448  CA  ILE A 446     104.318   8.889 152.280  1.00 85.02           C
ATOM   3449  C   ILE A 446     104.782   8.163 151.023  1.00 86.46           C
ATOM   3450  O   ILE A 446     104.654   6.938 150.921  1.00 84.51           O
ATOM   3451  CB  ILE A 446     103.000   9.652 152.049  1.00 86.42           C
ATOM   3452  CG1 ILE A 446     102.568  10.351 153.340  1.00 77.33           C
ATOM   3453  CG2 ILE A 446     101.905   8.708 151.582  1.00 81.07           C
ATOM   3454  CD1 ILE A 446     101.263  11.102 153.226  1.00 74.92           C
ATOM   3455  N   VAL A 447     105.318   8.889 150.039  1.00 90.31           N
ATOM   3456  CA  VAL A 447     105.780   8.212 148.830  1.00 88.70           C
ATOM   3457  C   VAL A 447     107.210   7.701 148.974  1.00 88.41           C
ATOM   3458  O   VAL A 447     107.584   6.734 148.294  1.00 96.77           O
ATOM   3459  CB  VAL A 447     105.638   9.116 147.588  1.00 87.55           C
ATOM   3460  CG1 VAL A 447     104.848  10.367 147.918  1.00 94.71           C
ATOM   3461  CG2 VAL A 447     106.989   9.472 146.995  1.00 88.38           C
ATOM   3462  N   GLY A 448     108.010   8.297 149.861  1.00 80.03           N
ATOM   3463  CA  GLY A 448     109.348   7.787 150.102  1.00 75.89           C
ATOM   3464  C   GLY A 448     109.361   6.388 150.686  1.00 78.68           C
ATOM   3465  O   GLY A 448     110.341   5.655 150.527  1.00 84.43           O
ATOM   3466  N   THR A 449     108.285   5.999 151.372  1.00 73.30           N
ATOM   3467  CA  THR A 449     108.218   4.716 152.059  1.00 72.71           C
ATOM   3468  C   THR A 449     107.681   3.584 151.191  1.00 74.77           C
ATOM   3469  O   THR A 449     107.833   2.416 151.565  1.00 72.62           O
ATOM   3470  CB  THR A 449     107.348   4.834 153.320  1.00 69.74           C
ATOM   3471  OG1 THR A 449     107.651   3.757 154.216  1.00 67.41           O
ATOM   3472  CG2 THR A 449     105.870   4.785 152.966  1.00 65.07           C
ATOM   3473  N   LEU A 450     107.065   3.897 150.044  1.00 79.63           N
ATOM   3474  CA  LEU A 450     106.392   2.855 149.268  1.00 82.05           C
ATOM   3475  C   LEU A 450     107.370   1.896 148.598  1.00 82.51           C
ATOM   3476  O   LEU A 450     107.164   0.670 148.700  1.00 88.07           O
ATOM   3477  CB  LEU A 450     105.438   3.498 148.258  1.00 84.21           C
ATOM   3478  CG  LEU A 450     104.311   4.337 148.864  1.00 87.61           C
ATOM   3479  CD1 LEU A 450     103.432   4.935 147.776  1.00 89.02           C
ATOM   3480  CD2 LEU A 450     103.487   3.502 149.833  1.00 86.02           C
TER
END
"""

# run() swaps Optimizers.Optimizer for Recorder while reduce2 runs
_Optimizer = Optimizers.Optimizer

class Recorder(_Optimizer):
  """After the fine step (the end of the search without symmetry copies):
  each Mover's stored score and a fresh score at its pose, with its atoms
  checked at the pose its OptimizerC holds; the names of the Movers of each
  multi-Mover Clique per call, in search order."""
  states = {}
  cliques = []
  def _OptimizeFine(self, optC, cliques, *args):
    _Optimizer._OptimizeFine(self, optC, cliques, *args)
    Recorder.cliques.append([[self._MoverName(self._interactionGraph.vertex_label(v))
      for v in c] for c in cliques if len(c) > 1])
    dot_scorer = Helpers.createDotScorer(self._extraAtomInfo, self._probePhil)
    # _ScoreMover() looks up each Mover's OptimizerC; the base sets it after this
    self._optimizerOf = dict((x, optC) for x in self._movers)
    for m in self._movers:
      coarse, fine = optC.GetCoarseLocation(m), optC.GetFineLocation(m)
      pr, index = m.CoarsePositions(), coarse
      if fine >= 0:
        pr, index = m.FinePositions(coarse), fine
      for a, xyz in zip(pr.atoms, pr.positions[index]):
        assert approx_equal(a.xyz, xyz, eps=1.e-6, out=None), (
          self._MoverName(m), a.name)
      Recorder.states[self._MoverName(m)] = (optC.GetHighScore(m),
        self._ScoreMover(m, dot_scorer))

def run(prefix, pdb_text, extra_args=()):
  with open(prefix + ".pdb", "w") as f:
    f.write(pdb_text)
  Recorder.states = {}
  Recorder.cliques = []
  Optimizers.Optimizer = Recorder
  try:
    run_program(program_class=reduce2.Program, logger=StringIO(), args=[
      prefix + ".pdb", "output.filename=" + prefix + "_out.pdb",
      "output.description_file_name=" + prefix + ".txt", "overwrite=True"]
      + list(extra_args))
  finally:
    Optimizers.Optimizer = _Optimizer
  with open(prefix + ".txt") as f:
    return f.read()

def assert_fresh(desc, fresh):
  """Every printed final score is the fresh score, and every Clique's final
  total the sum of its Movers' fresh scores, to 2 decimals."""
  printed, groups = {}, []
  for l in desc.split("\n"):
    if "Totals:" in l:
      groups.append((l.split("final score ")[1].strip(), []))
    elif "Singleton Movers:" in l:
      groups.append((None, []))
    elif " final score: " in l:
      name = l.split(" Initial score:")[0].strip()
      printed[name] = l.split(" final score: ")[1].split()[0]
      groups[-1][1].append(name)
  assert sorted(printed) == sorted(fresh), (printed, fresh)
  for name in fresh:
    assert printed[name] == "{:.2f}".format(fresh[name]), (
      name, printed[name], fresh[name])
  for total, names in groups:
    if total is not None:
      assert total == "{:.2f}".format(sum(fresh[n] for n in names)), (
        total, names)

def clique_scores(desc, cliques):
  """Per multi-Mover Clique (Mover names) its total after the coarse step and
  its final total, both as printed (2 decimals). cliques: per optimization pass
  (a section of the description) the Cliques' Movers in search order, the
  order of the "Clique of size N optimized" lines."""
  sections = desc.split(" Running Reduce optimization on model index")[1:]
  assert len(sections) == len(cliques), (len(sections), len(cliques))
  ret = []
  for text, names in zip(sections, cliques):
    coarse, totals, group = [], {}, None
    for l in text.split("\n"):
      if l.startswith(" Clique of size ") and "optimized with score" in l:
        coarse.append(float(l.split("optimized with score ")[1].split()[0]))
      elif "Totals:" in l:
        group = [float(l.split("final score ")[1]), []]
      elif "Singleton Movers:" in l:
        group = None
      elif group is not None and " final score: " in l:
        group[1].append(l.split(" Initial score:")[0].strip())
        totals[frozenset(group[1])] = group[0]
    assert len(coarse) == len(names), (coarse, names)
    for c, n in zip(coarse, names):
      ret.append((n, c, totals[frozenset(n)]))
  return ret

def assert_not_below_coarse(desc, cliques):
  for names, coarse, final in clique_scores(desc, cliques):
    assert final >= coarse, (names, coarse, final)

def test_stale_partner():
  desc = run("tst_reduce2_fresh_scores_4wsc", lys_gln_str,
    ["add_flip_movers=True"])
  states = Recorder.states
  assert_fresh(desc, dict((n, s[1]) for n, s in states.items()))
  scores = clique_scores(desc, Recorder.cliques)
  assert len(scores) == 1, scores
  assert_not_below_coarse(desc, Recorder.cliques)

def test_stale_stored_score():
  desc = run("tst_reduce2_fresh_scores_6vw1", pdb_6vw1, ["add_flip_movers=True"])
  states = Recorder.states
  # The case must stay: the stored score differs from the fresh one.
  stored, fresh = states["SingleHydrogenRotator at chain A THR 276 HG1"]
  assert abs(stored - fresh) > 0.5, (stored, fresh)
  assert_fresh(desc, dict((n, s[1]) for n, s in states.items()))

if __name__ == "__main__":
  test_stale_partner()
  test_stale_stored_score()
  print("OK")
