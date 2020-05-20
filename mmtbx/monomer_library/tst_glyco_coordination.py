from __future__ import absolute_import, division, print_function
from libtbx import easy_run

pdbs = {
  'linking_BDP.pdb' : '''
CRYST1   30.794   82.051   32.278  90.00 117.89  90.00 P 1 21 1
SCALE1      0.032474  0.000000  0.017187        0.00000
SCALE2      0.000000  0.012188  0.000000        0.00000
SCALE3      0.000000  0.000000  0.035052        0.00000
ATOM      1  NZ  LYS A  76     -10.288 -11.233  -9.001  1.00 11.47           N
ATOM      2  HD3 LYS A  76      -8.299 -12.103  -7.539  1.00 13.18           H
ATOM      3  HZ1 LYS A  76     -10.482 -12.102  -8.722  0.00 11.31           H
ATOM      4  HZ3 LYS A  76      -9.860 -11.203  -9.840  0.00 11.31           H
ATOM      5  OH  TYR A  83       2.764 -11.498  13.475  1.00 14.89           O
ATOM      6  HH  TYR A  83       3.087 -11.482  14.250  0.00 12.06           H
ATOM      7  N   ALA A 102       3.459 -14.630  10.397  1.00  9.86           N
ATOM      8  H   ALA A 102       3.898 -13.982  10.750  1.00  9.75           H
ATOM      9  H   ALA A 103       3.966 -16.501  11.919  1.00 11.43           H
ATOM     10  CD  GLU A 131      -8.631  -9.699 -12.199  1.00 11.23           C
ATOM     11  OE2 GLU A 131      -8.579 -10.186 -11.066  1.00 15.02           O
ATOM     12  HG3 GLU A 131      -7.187  -9.899 -13.559  1.00  8.96           H
TER
HETATM   13  C3  NAG A1175       3.340 -15.138  17.816  1.00 18.88           C
HETATM   14  C7  NAG A1175       5.232 -15.904  20.454  1.00 23.06           C
HETATM   15  O3  NAG A1175       4.706 -15.282  17.480  1.00 19.56           O
HETATM   16  O4  NAG A1175       2.467 -14.929  15.523  1.00 18.37           O
HETATM   17  O7  NAG A1175       5.969 -15.238  21.198  1.00 27.13           O
HETATM   18  H3  NAG A1175       3.124 -14.101  18.073  1.00 19.18           H
HETATM   19  H82 NAG A1175       6.014 -16.957  18.773  1.00 21.76           H
HETATM   20  HO4 NAG A1175       3.094 -14.171  15.694  0.00 20.13           H
HETATM   21  C1  BDP A1176       5.315 -14.207  16.731  1.00 18.20           C
HETATM   22  C2  BDP A1176       6.818 -14.233  17.016  1.00 18.52           C
HETATM   23  C3  BDP A1176       7.628 -13.300  16.107  1.00 17.62           C
HETATM   24  C4  BDP A1176       7.296 -13.523  14.622  1.00 15.67           C
HETATM   25  C5  BDP A1176       5.772 -13.482  14.483  1.00 15.86           C
HETATM   26  C6  BDP A1176       5.352 -13.910  13.096  1.00 14.73           C
HETATM   27  O2  BDP A1176       7.014 -13.855  18.368  1.00 21.01           O
HETATM   28  O3  BDP A1176       9.012 -13.490  16.376  1.00 18.59           O
HETATM   29  O4  BDP A1176       7.901 -12.567  13.765  1.00 13.68           O
HETATM   30  O5  BDP A1176       5.100 -14.402  15.349  1.00 17.16           O
HETATM   31  O6A BDP A1176       5.577 -15.096  12.735  1.00 15.20           O
HETATM   32  O6B BDP A1176       4.762 -13.054  12.384  1.00 15.71           O
HETATM   33  H1  BDP A1176       4.902 -13.247  17.038  1.00 18.05           H
HETATM   34  H2  BDP A1176       7.167 -15.256  16.878  1.00 18.76           H
HETATM   35  H3  BDP A1176       7.358 -12.285  16.370  1.00 17.48           H
HETATM   36  HB  BDP A1176       6.144 -13.679  18.812  0.00 18.97           H
HETATM   37  H4  BDP A1176       7.635 -14.523  14.351  1.00 15.66           H
HETATM   38  H5  BDP A1176       5.414 -12.471  14.677  1.00 15.52           H
HETATM   39  HC  BDP A1176       9.118 -14.232  17.068  0.00 17.73           H
HETATM   40  C1  NAG A1177       8.992 -12.956  12.918  1.00 13.58           C
HETATM   41  C2  NAG A1177       9.213 -11.887  11.849  1.00 12.17           C
HETATM   42  N2  NAG A1177       7.954 -11.763  11.092  1.00 12.56           N
HETATM   43  O5  NAG A1177      10.174 -13.136  13.662  1.00 14.08           O
HETATM   44  H1  NAG A1177       8.740 -13.896  12.428  1.00 13.21           H
HETATM   45  H2  NAG A1177       9.424 -10.949  12.357  1.00 12.00           H
HETATM   46  H82 NAG A1177       5.243 -11.105  10.751  1.00 10.52           H
HETATM   47  HN2 NAG A1177       7.399 -12.601  11.034  1.00 11.60           H
''',
  'linking_FUC.pdb' : '''
CRYST1   72.333   98.984  152.019  90.00  90.00  90.00 P 21 21 21
SCALE1      0.013825  0.000000  0.000000        0.00000
SCALE2      0.000000  0.010103  0.000000        0.00000
SCALE3      0.000000  0.000000  0.006578        0.00000
HETATM    1  C2  GAL E 207     -50.952  -3.199  42.155  0.90 39.70           C
ANISOU    1  C2  GAL E 207     5050   5702   4333   -619   -389    745       C
HETATM    2  O2  GAL E 207     -51.145  -3.682  40.855  0.90 37.04           O
ANISOU    2  O2  GAL E 207     4034   3948   6093  -1230     63   -597       O
HETATM    3  C1  FUC E 210     -50.167  -4.605  40.258  0.90 40.07           C
ANISOU    3  C1  FUC E 210     5310   4344   5571   -330    -46  -1039       C
HETATM    4  C2  FUC E 210     -50.885  -5.345  39.138  0.90 39.95           C
ANISOU    4  C2  FUC E 210     4770   4247   6162   -993    207  -1570       C
HETATM    5  C3  FUC E 210     -51.328  -4.387  38.002  0.90 37.97           C
ANISOU    5  C3  FUC E 210     3982   4753   5691  -1411   -306  -2427       C
HETATM    6  C4  FUC E 210     -50.128  -3.526  37.490  0.90 42.09           C
ANISOU    6  C4  FUC E 210     4378   5641   5972  -1302   -138  -1850       C
HETATM    7  C5  FUC E 210     -49.452  -2.849  38.740  0.90 38.93           C
ANISOU    7  C5  FUC E 210     4139   4684   5967  -1869   -447  -1235       C
HETATM    8  C6  FUC E 210     -48.234  -1.996  38.345  0.90 34.08           C
ANISOU    8  C6  FUC E 210     2501   5771   4677   -766    279   -903       C
HETATM    9  O2  FUC E 210     -52.027  -5.944  39.685  0.90 45.56           O
ANISOU    9  O2  FUC E 210     4712   5115   7484   -531   1136  -1747       O
HETATM   10  O3  FUC E 210     -51.952  -5.101  36.892  0.90 48.86           O
ANISOU   10  O3  FUC E 210     4944   7622   6000  -2188   -157  -3300       O
HETATM   11  O4  FUC E 210     -49.187  -4.308  36.634  0.90 47.16           O
ANISOU   11  O4  FUC E 210     3538   8928   5451  -2181      3  -3117       O
HETATM   12  O5  FUC E 210     -49.036  -3.825  39.759  0.90 43.21           O
ANISOU   12  O5  FUC E 210     5470   4074   6872   -230   -266  -1185       O
''',
  'linking_FUL.pdb' : '''
CRYST1   65.168   79.104   68.912  90.00  99.54  90.00 P 1 21 1
SCALE1      0.015345  0.000000  0.002579        0.00000
SCALE2      0.000000  0.012642  0.000000        0.00000
SCALE3      0.000000  0.000000  0.014715        0.00000
HETATM   22  C1  FUL A   1       3.616   6.912  40.821  1.00 17.72           C
HETATM   23  C2  FUL A   1       3.218   7.956  39.796  1.00 16.30           C
HETATM   24  C3  FUL A   1       3.320   9.339  40.423  1.00 18.19           C
HETATM   25  C4  FUL A   1       4.718   9.554  41.013  1.00 18.04           C
HETATM   26  C5  FUL A   1       5.092   8.396  41.931  1.00 19.32           C
HETATM   27  C6  FUL A   1       6.539   8.468  42.405  1.00 17.45           C
HETATM   28  O2  FUL A   1       1.868   7.723  39.365  1.00 21.28           O
HETATM   29  O3  FUL A   1       3.038  10.352  39.442  1.00 16.14           O
HETATM   30  O4  FUL A   1       5.707   9.719  39.980  1.00 17.14           O
HETATM   31  O5  FUL A   1       4.942   7.163  41.219  1.00 17.82           O
HETATM   32  C1  GAL A   2       2.069   3.824  40.258  1.00 23.45           C
HETATM   33  C2  GAL A   2       3.046   4.660  41.077  1.00 21.87           C
HETATM   34  C3  GAL A   2       4.160   3.797  41.666  1.00 19.09           C
HETATM   35  C4  GAL A   2       3.600   2.543  42.324  1.00 24.16           C
HETATM   36  C5  GAL A   2       2.790   1.808  41.283  1.00 23.39           C
HETATM   37  C6  GAL A   2       2.199   0.520  41.847  1.00 30.62           C
HETATM   38  O2  GAL A   2       3.591   5.625  40.203  1.00 21.90           O
HETATM   39  O3  GAL A   2       4.915   4.525  42.617  1.00 21.17           O
HETATM   40  O4  GAL A   2       2.717   2.879  43.372  1.00 26.99           O
HETATM   41  O5  GAL A   2       1.700   2.637  40.920  1.00 25.10           O
HETATM   42  O6  GAL A   2       3.225  -0.291  42.368  1.00 41.60           O
HETATM   43  C2  NAG A   3      -2.260   3.717  38.368  1.00 37.89           C
HETATM   44  C3  NAG A   3      -1.264   3.867  39.503  1.00 33.20           C
HETATM   45  C4  NAG A   3       0.142   4.145  38.966  1.00 25.84           C
HETATM   46  C5  NAG A   3       0.131   5.190  37.840  1.00 22.60           C
HETATM   47  C6  NAG A   3       1.505   5.377  37.199  1.00 23.56           C
HETATM   48  O3  NAG A   3      -1.258   2.696  40.294  1.00 31.55           O
HETATM   49  O4  NAG A   3       0.912   4.638  40.032  1.00 25.82           O
HETATM   50  O5  NAG A   3      -0.838   4.866  36.859  1.00 26.45           O
HETATM   51  O6  NAG A   3       2.039   4.148  36.747  1.00 24.74           O
''',
'linking_BMA.pdb' : '''
CRYST1   54.153   75.989   80.970 105.69 107.32 106.93 P 1
SCALE1      0.018466  0.005621  0.008645        0.00000
SCALE2      0.000000  0.013756  0.005844        0.00000
SCALE3      0.000000  0.000000  0.014056        0.00000
HETATM  127  C1  NAG A 601      18.073   8.873  24.744  1.00 14.91           C
HETATM  128  C2  NAG A 601      17.516   7.798  23.848  1.00 12.61           C
HETATM  129  C3  NAG A 601      18.605   7.334  22.886  1.00 13.69           C
HETATM  130  C4  NAG A 601      19.200   8.525  22.137  1.00 16.16           C
HETATM  131  C5  NAG A 601      19.611   9.620  23.119  1.00 15.35           C
HETATM  132  C6  NAG A 601      20.081  10.890  22.445  1.00 17.36           C
HETATM  133  C7  NAG A 601      15.747   6.206  24.449  1.00 15.82           C
HETATM  134  C8  NAG A 601      15.353   5.083  25.347  1.00 15.25           C
HETATM  135  N2  NAG A 601      16.983   6.685  24.621  1.00 13.76           N
HETATM  136  O3  NAG A 601      18.035   6.409  21.970  1.00 15.17           O
HETATM  137  O4  NAG A 601      20.387   8.143  21.441  1.00 15.20           O
HETATM  138  O5  NAG A 601      18.496   9.974  23.952  1.00 14.87           O
HETATM  139  O6  NAG A 601      19.017  11.504  21.736  1.00 27.11           O
HETATM  140  O7  NAG A 601      14.977   6.665  23.612  1.00 15.56           O
HETATM  141  H1  NAG A 601      18.837   8.518  25.237  1.00 17.89           H
HETATM  142  H2  NAG A 601      16.792   8.185  23.320  1.00 15.13           H
HETATM  143  H3  NAG A 601      19.311   6.889  23.392  1.00 16.43           H
HETATM  144  H4  NAG A 601      18.546   8.878  21.505  1.00 19.39           H
HETATM  145  H5  NAG A 601      20.330   9.280  23.684  1.00 18.42           H
HETATM  146  H61 NAG A 601      20.801  10.675  21.822  1.00 20.83           H
HETATM  147  H62 NAG A 601      20.415  11.508  23.122  1.00 20.83           H
HETATM  148  H81 NAG A 601      15.443   5.363  26.278  1.00 18.30           H
HETATM  149  H82 NAG A 601      15.932   4.315  25.181  1.00 18.30           H
HETATM  150  H83 NAG A 601      14.425   4.835  25.171  1.00 18.30           H
HETATM  151  HN2 NAG A 601      17.508   6.313  25.267  1.00 16.51           H
HETATM  152  HO3 NAG A 601      18.620   5.763  21.799  1.00 18.20           H
HETATM  153  HO6 NAG A 601      19.250  12.332  21.513  1.00 32.53           H
HETATM  154  C1  NAG A 602      20.291   7.982  20.029  1.00 15.66           C
HETATM  155  C2  NAG A 602      21.694   8.153  19.475  1.00 17.28           C
HETATM  156  C3  NAG A 602      21.716   7.782  18.005  1.00 15.28           C
HETATM  157  C4  NAG A 602      21.137   6.392  17.799  1.00 15.84           C
HETATM  158  C5  NAG A 602      19.761   6.316  18.435  1.00 15.39           C
HETATM  159  C6  NAG A 602      19.131   4.952  18.330  1.00 16.01           C
HETATM  160  C7  NAG A 602      23.130   9.806  20.573  1.00 19.18           C
HETATM  161  C8  NAG A 602      23.534  11.246  20.642  1.00 25.12           C
HETATM  162  N2  NAG A 602      22.191   9.504  19.668  1.00 18.24           N
HETATM  163  O3  NAG A 602      23.055   7.851  17.521  1.00 19.36           O
HETATM  164  O4  NAG A 602      21.007   6.158  16.402  1.00 17.40           O
HETATM  165  O5  NAG A 602      19.866   6.637  19.832  1.00 14.98           O
HETATM  166  O6  NAG A 602      19.984   3.975  18.908  1.00 18.37           O
HETATM  167  O7  NAG A 602      23.626   8.958  21.303  1.00 22.72           O
HETATM  168  H1  NAG A 602      19.683   8.584  19.558  1.00 18.79           H
HETATM  169  H2  NAG A 602      22.284   7.537  19.949  1.00 20.74           H
HETATM  170  H3  NAG A 602      21.172   8.425  17.512  1.00 18.34           H
HETATM  171  H4  NAG A 602      21.726   5.723  18.196  1.00 19.01           H
HETATM  172  H5  NAG A 602      19.177   6.968  18.004  1.00 18.47           H
HETATM  173  H61 NAG A 602      18.276   4.954  18.800  1.00 19.22           H
HETATM  174  H62 NAG A 602      18.983   4.736  17.390  1.00 19.22           H
HETATM  175  H81 NAG A 602      22.749  11.796  20.825  1.00 30.15           H
HETATM  176  H82 NAG A 602      23.928  11.515  19.790  1.00 30.15           H
HETATM  177  H83 NAG A 602      24.189  11.366  21.356  1.00 30.15           H
HETATM  178  HN2 NAG A 602      21.854  10.173  19.147  1.00 21.88           H
HETATM  179  HO3 NAG A 602      23.059   8.216  16.711  1.00 23.23           H
HETATM  180  HO6 NAG A 602      19.539   3.212  19.002  1.00 22.05           H
HETATM  181  C1  BMA A 603      21.666   4.994  15.926  1.00 16.37           C
HETATM  182  C2  BMA A 603      20.957   4.639  14.627  1.00 15.23           C
HETATM  183  C3  BMA A 603      21.802   3.580  13.908  1.00 15.63           C
HETATM  184  C4  BMA A 603      23.218   4.011  13.748  1.00 16.11           C
HETATM  185  C5  BMA A 603      23.792   4.383  15.139  1.00 17.26           C
HETATM  186  C6  BMA A 603      25.211   4.922  15.067  1.00 17.68           C
HETATM  187  O2  BMA A 603      20.835   5.765  13.752  1.00 18.00           O
HETATM  188  O3  BMA A 603      21.247   3.245  12.617  1.00 15.82           O
HETATM  189  O4  BMA A 603      23.960   2.929  13.176  1.00 18.33           O
HETATM  190  O5  BMA A 603      22.968   5.428  15.696  1.00 17.11           O
HETATM  191  O6  BMA A 603      25.237   6.022  14.189  1.00 22.16           O
HETATM  192  H1  BMA A 603      21.685   4.111  16.593  1.00 19.64           H
HETATM  193  H2  BMA A 603      19.965   4.228  14.868  1.00 18.27           H
HETATM  194  H3  BMA A 603      21.797   2.650  14.491  1.00 18.76           H
HETATM  195  H4  BMA A 603      23.244   4.908  13.110  1.00 19.33           H
HETATM  196  H5  BMA A 603      23.785   3.498  15.797  1.00 20.72           H
HETATM  197  H61 BMA A 603      25.873   4.112  14.724  1.00 21.22           H
HETATM  198  H62 BMA A 603      25.511   5.202  16.086  1.00 21.22           H
HETATM  199  HO2 BMA A 603      19.884   5.942  13.673  1.00 21.60           H
HETATM  200  HO4 BMA A 603      24.121   3.188  12.255  1.00 22.00           H
HETATM  201  C1  MAN A 604      26.507   6.653  14.048  1.00 44.19           C
HETATM  202  C2  MAN A 604      26.028   7.844  13.187  1.00 52.01           C
HETATM  203  C3  MAN A 604      25.606   7.365  11.770  1.00 44.32           C
HETATM  204  C4  MAN A 604      26.691   6.484  11.139  1.00 37.58           C
HETATM  205  C5  MAN A 604      26.986   5.304  12.075  1.00 43.40           C
HETATM  206  C6  MAN A 604      28.053   4.355  11.521  1.00 42.65           C
HETATM  207  O2  MAN A 604      27.085   8.780  12.982  1.00 49.10           O
HETATM  208  O3  MAN A 604      25.279   8.451  10.896  1.00 60.52           O
HETATM  209  O4  MAN A 604      26.256   5.987   9.869  1.00 52.68           O
HETATM  210  O5  MAN A 604      27.432   5.840  13.374  1.00 32.14           O
HETATM  211  O6  MAN A 604      29.263   5.091  11.348  1.00 54.50           O
HETATM  212  H1  MAN A 604      27.067   7.026  14.919  1.00 53.03           H
HETATM  213  H2  MAN A 604      25.175   8.320  13.687  1.00 62.41           H
HETATM  214  H3  MAN A 604      24.696   6.763  11.862  1.00 53.19           H
HETATM  215  H4  MAN A 604      27.612   7.083  11.036  1.00 45.09           H
HETATM  216  H5  MAN A 604      26.062   4.728  12.213  1.00 52.08           H
HETATM  217  H61 MAN A 604      28.184   3.523  12.230  1.00 51.18           H
HETATM  218  H62 MAN A 604      27.681   3.947  10.572  1.00 51.18           H
HETATM  219  HO2 MAN A 604      27.279   9.223  13.817  1.00 58.92           H
HETATM  220  HO3 MAN A 604      25.135   8.045  10.025  1.00 72.62           H
HETATM  221  HO4 MAN A 604      26.982   6.182   9.256  1.00 63.22           H
HETATM  222  C1  MAN A 605      20.779   1.936  12.559  1.00 15.17           C
HETATM  223  C2  MAN A 605      20.527   1.587  11.091  1.00 14.93           C
HETATM  224  C3  MAN A 605      19.321   2.377  10.553  1.00 15.06           C
HETATM  225  C4  MAN A 605      18.102   2.184  11.429  1.00 16.79           C
HETATM  226  C5  MAN A 605      18.481   2.544  12.880  1.00 16.62           C
HETATM  227  C6  MAN A 605      17.368   2.242  13.879  1.00 31.84           C
HETATM  228  O2  MAN A 605      20.150   0.231  10.926  1.00 21.91           O
HETATM  229  O3  MAN A 605      19.038   1.972   9.209  1.00 16.78           O
HETATM  230  O4  MAN A 605      17.034   3.060  11.020  1.00 17.30           O
HETATM  231  O5  MAN A 605      19.610   1.767  13.304  1.00 18.47           O
HETATM  232  O6  MAN A 605      17.620   3.039  15.033  1.00 36.31           O
HETATM  233  H1  MAN A 605      21.481   1.210  12.998  1.00 18.20           H
HETATM  234  H2  MAN A 605      21.426   1.824  10.505  1.00 17.91           H
HETATM  235  H3  MAN A 605      19.567   3.444  10.523  1.00 18.08           H
HETATM  236  H4  MAN A 605      17.796   1.124  11.393  1.00 20.15           H
HETATM  237  H5  MAN A 605      18.718   3.615  12.935  1.00 19.94           H
HETATM  238  H61 MAN A 605      16.400   2.486  13.416  1.00 38.21           H
HETATM  239  H62 MAN A 605      17.393   1.167  14.102  1.00 38.21           H
HETATM  240  HO3 MAN A 605      19.304   2.723   8.654  1.00 20.14           H
HETATM  241  HO4 MAN A 605      16.415   2.497  10.530  1.00 20.76           H
''',
  }

links = {
  'linking_BDP.pdb' : '''
  Links applied
    BETA1-3
      " NAG A1175 " - " BDP A1176 "
    BETA1-4
      " BDP A1176 " - " NAG A1177 "''',
  'linking_FUC.pdb' : '''
  Links applied
    BETA1-2
      " GAL E 207 " - " FUC E 210 "
       ~> Even though FUC is an alpha isomer, a beta linkage is required...''',
  'linking_FUL.pdb' : '''
  Links applied
    ALPHA1-2
      " FUL A   1 " - " GAL A   2 "
       ~> Even though GAL is a beta isomer, an alpha linkage is required...
    BETA1-4
      " GAL A   2 " - " NAG A   3 "''',
  'linking_BMA.pdb' : '''
  Links applied
    ALPHA1-3
      " BMA A 603 " - " MAN A 605 "
    ALPHA1-6
      " BMA A 603 " - " MAN A 604 "
    BETA1-4
      " NAG A 601 " - " NAG A 602 "
      " NAG A 602 " - " BMA A 603 "''',
  }

def run(only_i=None):
  try: only_i=int(only_i)
  except ValueError: only_i=None
  except TypeError: only_i=None
  cifs = ""
  for pdb in pdbs:
    f=open(pdb, "w")
    f.write(pdbs[pdb])
    f.close()
    if pdb.endswith(".cif"): cifs += " %s" % pdb
    links[pdb] = links[pdb].splitlines()
  j=0
  for pdb in sorted(pdbs):
    #break
    if pdb.endswith(".cif"): continue
    if pdb.endswith(".params"): continue
    print('pdb',pdb)
    j+=1
    if only_i is not None and only_i!=j: continue
    # log_filename = "%s.log" % (pdb)
    cmd = "phenix.pdb_interpretation %s write_geo_file=False" % pdb
    result = easy_run.fully_buffered(cmd).raise_if_errors()
    assert (result.return_code == 0)
    for line in result.stdout_lines:
      if (line in links[pdb]):
        links[pdb].remove(line)
    if links[pdb]:
      raise RuntimeError("Missing expected log output %s" % pdb)
    print("OK")

if __name__=="__main__":
  import sys
  run(*tuple(sys.argv[1:]))
