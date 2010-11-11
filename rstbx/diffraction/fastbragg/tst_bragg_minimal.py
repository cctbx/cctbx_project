from scitbx.array_family import flex
from scitbx.matrix import sqr
from rstbx.diffraction.fastbragg import foo
from rstbx.diffraction.fastbragg import detector,camera,crystal
from rstbx.diffraction.fastbragg import fast_bragg_simulation
import libtbx.load_env

pdb_lines = """HEADER    ELECTRON TRANSPORT                      17-JAN-08   3BZ1
CRYST1  127.692  225.403  306.106  90.00  90.00  90.00 P 21 21 21    4
ATOM    100  CH2 TRP A  20      56.430  40.246  32.394  1.00 81.85           C
ATOM    200  N   PHE A  33      53.226  47.627  40.654  1.00 70.99           N
ATOM    300  CD1 ILE A  46      45.639  29.364  40.352  1.00 64.65           C
ATOM    400  CB  ILE A  60      31.998  29.624  61.196  1.00 48.31           C
ATOM    500  CZ  TYR A  73      33.394  17.375  51.997  1.00 59.10           C
ATOM    600  N   ILE A  89      40.139  30.719  66.768  1.00 57.41           N
ATOM    700  CB  ALA A 100      48.660  18.955  63.036  1.00 50.19           C
ATOM    800  OH  TYR A 112      35.813  37.284  55.881  1.00 49.58           O
ATOM    900  O   CYS A 125      47.637  48.097  42.161  1.00 69.98           O
ATOM   1000  CB  ARG A 136      53.053  54.518  31.934  1.00 62.41           C
ATOM   1100  O   SER A 148      38.270  53.004  48.759  1.00 59.64           O
ATOM   1200  OH  TYR A 161      29.848  43.600  64.133  1.00 51.65           O
ATOM   1300  CG1 ILE A 176      25.474  34.051  46.196  1.00 49.30           C
ATOM   1400  C   GLU A 189      24.280  45.586  66.763  1.00 54.83           C
ATOM   1500  CB  LEU A 200      25.906  59.092  53.794  1.00 57.65           C
ATOM   1600  CD2 HIS A 215      26.386  60.725  30.820  1.00 66.14           C
ATOM   1700  O   GLU A 229      15.446  53.147  11.386  1.00 94.97           O
ATOM   1800  OE1 GLN A 241      31.691  60.361  15.701  1.00 85.33           O
ATOM   1900  CG  TYR A 254      14.970  67.624  27.886  1.00 80.00           C
ATOM   2000  CE1 PHE A 265      26.289  68.506  37.590  1.00 88.41           C
ATOM   2100  N   TRP A 278      29.660  62.229  43.659  1.00 63.59           N
ATOM   2200  CG1 ILE A 290      29.934  48.592  56.874  1.00 56.42           C
ATOM   2300  C   ASN A 303      13.257  53.636  73.305  1.00 53.67           C
ATOM   2400  CG2 THR A 316       7.855  57.090  68.127  1.00 58.80           C
ATOM   2500  CB  GLU A 329      19.451  40.098  69.217  1.00 49.38           C
ATOM   2600  CA  LEU A 341      25.291  40.222  77.420  1.00 58.62           C
ATOM   2700  CD2 HIS B   9      10.231  39.622  13.536  1.00 61.64           C
ATOM   2800  O   HIS B  23      -1.741  36.161  13.565  1.00 73.51           O
ATOM   2900  C   MET B  37      -3.628  21.765  27.157  1.00 60.16           C
ATOM   3000  N   PRO B  50     -14.550  15.131  42.592  1.00 59.82           N
ATOM   3100  CE2 PHE B  61      -3.211  29.588  31.858  1.00 36.72           C
ATOM   3200  CG  TRP B  75     -15.698  12.830  23.209  1.00 67.94           C
ATOM   3300  CD  PRO B  88     -20.869  13.773  26.882  1.00 88.90           C
ATOM   3400  CG1 ILE B 101      -4.955  20.401  18.664  1.00 55.70           C
ATOM   3500  CD2 HIS B 114      -3.854  35.717   4.227  1.00 76.12           C
ATOM   3600  CD  ARG B 124       1.132  49.272  -3.409  1.00 75.98           C
ATOM   3700  CD  LYS B 137     -13.446  52.508   4.977  1.00 81.73           C
ATOM   3800  CB  CYS B 150     -19.609  42.068  22.286  1.00 71.10           C
ATOM   3900  N   PRO B 164     -25.221  27.178  25.947  1.00 78.19           N
ATOM   4000  CB  SER B 177     -24.093  23.363  37.543  1.00 69.60           C
ATOM   4100  CZ  PHE B 190     -23.890  49.302  34.563  1.00 82.38           C
ATOM   4200  C   ALA B 204     -24.194  46.989  25.481  1.00 68.99           C
ATOM   4300  CD2 LEU B 218     -21.126  54.934   7.329  1.00 70.45           C
ATOM   4400  CG  ARG B 230       2.163  61.215  14.361  1.00 77.19           C
ATOM   4500  CB  ALA B 244      -6.668  48.852  25.158  1.00 60.23           C
ATOM   4600  CZ2 TRP B 257      -3.347  36.591  46.798  1.00 52.77           C
ATOM   4700  C   THR B 271      -8.867  33.510  49.162  1.00 53.31           C
ATOM   4800  CD  GLN B 281     -13.383  29.920  62.401  1.00 67.13           C
ATOM   4900  CA  SER B 294     -33.026  21.366  66.399  1.00 69.25           C
ATOM   5000  OE2 GLU B 307     -14.800  13.540  54.306  1.00 68.45           O
ATOM   5100  CD  PRO B 319      -5.922  31.146  48.057  1.00 47.82           C
ATOM   5200  N   ASP B 334       0.902  19.258  50.198  1.00 44.55           N
ATOM   5300  CA  ARG B 347     -10.219  20.652  76.188  1.00 61.50           C
ATOM   5400  CG  ARG B 358      -7.285  29.120  62.704  1.00 50.92           C
ATOM   5500  CB  LEU B 370      -4.691  26.163  73.320  1.00 50.59           C
ATOM   5600  CE2 PHE B 383      13.890  27.916  70.647  1.00 33.21           C
ATOM   5700  N   GLN B 395      -1.283  20.760  81.492  1.00 56.45           N
ATOM   5800  N   GLN B 409     -14.069   9.629  68.475  1.00 59.62           N
ATOM   5900  O   ALA B 421      -2.576  22.480  66.140  1.00 57.75           O
ATOM   6000  C   ASP B 433      -0.960  11.201  49.784  1.00 55.09           C
ATOM   6100  N   SER B 446      -2.763  27.862  44.195  1.00 50.25           N
ATOM   6200  C   PHE B 458       2.089  42.286  31.803  1.00 53.50           C
ATOM   6300  CA  HIS B 469       5.917  52.473  21.275  1.00 65.93           C
ATOM   6400  C   GLY B 481      18.369  63.066  15.911  1.00 86.80           C
ATOM   6500  C   GLU C  29      42.829  84.391  42.151  1.00 96.99           C
ATOM   6600  N   LEU C  42      49.063  79.537  47.972  1.00 80.65           N
ATOM   6700  C   HIS C  56      47.160  76.079  63.881  1.00 75.29           C
ATOM   6800  CD1 LEU C  69      40.241  76.971  80.264  1.00 73.46           C
ATOM   6900  CB  MET C  81      37.685  59.675  83.702  1.00 66.57           C
ATOM   7000  N   THR C  94      47.149  59.964  83.079  1.00 72.30           N
ATOM   7100  CG2 THR C 108      39.735  68.134  90.385  1.00 84.58           C
ATOM   7200  CA  SER C 121      49.614  80.118  74.233  1.00 71.01           C
ATOM   7300  CB  ARG C 135      50.883  92.723  59.014  1.00 88.70           C
ATOM   7400  CA  GLY C 148      62.945  81.923  53.170  1.00 86.41           C
ATOM   7500  OG1 THR C 159      69.683  71.849  55.898  1.00 76.62           O
ATOM   7600  CG  LEU C 173      64.970  62.992  74.939  1.00 75.86           C
ATOM   7700  CE2 TYR C 186      46.358  53.614  86.921  1.00 54.85           C
ATOM   7800  O   THR C 200      65.438  52.726  87.800  1.00 81.03           O
ATOM   7900  CZ  TYR C 212      57.265  45.691  74.323  1.00 64.32           C
ATOM   8000  CG2 VAL C 225      48.991  48.489  68.817  1.00 54.23           C
ATOM   8100  CG  TRP C 239      66.638  53.123  72.325  1.00 79.64           C
ATOM   8200  CB  ILE C 252      71.570  69.941  60.611  1.00 88.29           C
ATOM   8300  CD1 PHE C 264      60.724  60.118  48.299  1.00 79.04           C
ATOM   8400  CG  LEU C 276      48.684  64.549  53.259  1.00 64.89           C
ATOM   8500  CG2 VAL C 290      41.832  54.897  72.981  1.00 57.48           C
ATOM   8600  CZ  PHE C 301      44.591  57.265  86.195  1.00 70.05           C
ATOM   8700  C   THR C 316      30.583  36.329  87.074  1.00 59.01           C
ATOM   8800  CG2 VAL C 328      35.267  26.607  88.411  1.00 63.86           C
ATOM   8900  CB  ARG C 343      41.831  23.983  76.148  1.00 51.02           C
ATOM   9000  CA  ARG C 357      38.290  34.957  71.104  1.00 61.60           C
ATOM   9100  C   GLU C 367      46.540  37.293  87.483  1.00 64.97           C
ATOM   9200  O   ILE C 380      37.697  30.576  92.067  1.00 67.63           O
ATOM   9300  O   ARG C 391      35.039  45.120  86.966  1.00 60.78           O
ATOM   9400  CD1 LEU C 404      27.268  58.766  66.001  1.00 57.48           C
ATOM   9500  CB  PHE C 419      27.493  55.491  79.223  1.00 67.79           C
ATOM   9600  CD1 PHE C 431      40.201  52.507  60.945  1.00 57.40           C
ATOM   9700  CG  TRP C 443      43.756  64.474  45.968  1.00 75.60           C
ATOM   9800  CA  LYS C 457      46.056  69.932  33.682  1.00 75.97           C
ATOM   9900  CG  PRO C 470      53.283  49.954  20.292  1.00 81.53           C
ATOM  10000  CG  TRP D  21      -0.801  81.734  29.348  1.00 77.02           C
ATOM  10100  CB  TRP D  32      -0.735  70.313  28.373  1.00 71.40           C
ATOM  10200  C   LEU D  45       3.971  67.431  50.438  1.00 56.33           C
ATOM  10300  NE1 TRP D  58      -3.948  56.069  71.824  1.00 63.14           N
ATOM  10400  O   CYS D  71       9.769  64.233  65.991  1.00 66.26           O
ATOM  10500  CG  MET D  85      -9.271  51.951  57.279  1.00 61.82           C
ATOM  10600  CD  GLN D  98     -20.334  58.870  58.371  1.00 79.48           C
ATOM  10700  O   TRP D 111      -1.668  57.883  50.989  1.00 55.61           O
ATOM  10800  N   GLY D 124       1.991  62.768  35.633  1.00 65.92           N
ATOM  10900  CG  LEU D 135       7.056  71.817  23.809  1.00 70.97           C
ATOM  11000  CA  PRO D 149      10.850  54.606  39.091  1.00 55.93           C
ATOM  11100  CB  PRO D 161      -0.735  47.345  49.987  1.00 48.29           C
ATOM  11200  CZ  PHE D 173       3.153  56.984  49.943  1.00 50.71           C
ATOM  11300  CZ  PHE D 185       7.181  43.691  50.494  1.00 37.45           C
ATOM  11400  C   HIS D 197      17.164  41.521  43.801  1.00 50.10           C
ATOM  11500  CA  ILE D 213      33.634  55.063  37.092  1.00 69.26           C
ATOM  11600  N   GLY D 226      43.193  75.028  28.981  1.00 80.84           N
ATOM  11700  CG  GLN D 239      28.634  70.306  20.768  1.00 96.48           C
ATOM  11800  C   PHE D 252      42.989  51.474  30.264  1.00 69.10           C
ATOM  11900  CA  LYS D 264      30.303  56.876  23.582  1.00 65.51           C
ATOM  12000  N   VAL D 274      21.376  49.633  32.352  1.00 59.67           N
ATOM  12100  CA  GLY D 288       5.063  38.350  43.091  1.00 51.26           C
ATOM  12200  CB  SER D 300      17.107  20.052  47.341  1.00 49.18           C
ATOM  12300  OE1 GLU D 312      25.457  24.792  61.967  1.00 56.19           O
ATOM  12400  O   GLY D 324      11.625  37.142  65.784  1.00 50.64           O
ATOM  12500  CE1 HIS D 336      -2.635  43.943  61.172  1.00 62.69           C
ATOM  12600  CG  ARG D 348      13.912  33.212  72.962  1.00 50.25           C
ATOM  12700  CA  ASP E  12      18.164  91.898  42.777  1.00 89.88           C
ATOM  12800  CG  HIS E  23       9.724  89.213  52.411  1.00 85.44           C
ATOM  12900  CB  LEU E  36      -0.117  81.291  67.900  1.00 70.26           C
ATOM  13000  CG2 THR E  49      -1.179  64.296  67.834  1.00 61.18           C
ATOM  13100  CB  ARG E  61      -2.323  49.122  79.882  1.00 88.20           C
ATOM  13200  CD  LYS E  73     -17.743  60.195  61.429  1.00 74.64           C
ATOM  13300  N   VAL F  11      15.993  96.516  32.978  1.00 95.73           N
ATOM  13400  O   ALA F  22       5.710  78.633  45.555  1.00 73.32           O
ATOM  13500  C   GLY F  35       6.641  77.242  64.377  1.00 67.22           C
ATOM  13600  CA  ARG H   4      -6.132  42.027  -7.866  1.00 90.52           C
ATOM  13700  ND2 ASN H  15      -6.361  49.254   1.776  1.00 81.21           N
ATOM  13800  CA  LEU H  30     -18.790  65.525  21.049  1.00 73.95           C
ATOM  13900  CG  LEU H  42     -14.264  53.453  33.427  1.00 79.57           C
ATOM  14000  CD1 ILE H  54     -20.487  50.335  50.523  1.00 68.77           C
ATOM  14100  CA  GLU I   2      58.092  23.554  57.876  1.00 71.08           C
ATOM  14200  CG  PHE I  14      61.361  36.617  43.358  1.00 66.86           C
ATOM  14300  C   GLY I  26      62.908  55.772  38.054  1.00 78.38           C
ATOM  14400  O   PRO J   9      24.648  91.273  51.531  1.00 78.21           O
ATOM  14500  CB  VAL J  23      19.540  80.157  65.363  1.00 69.47           C
ATOM  14600  CD1 LEU J  36       5.079  66.010  73.040  1.00 74.63           C
ATOM  14700  CE2 PHE K  18      34.702  85.604  85.570  1.00 81.28           C
ATOM  14800  CA  PHE K  32      37.105  79.841  61.828  1.00 68.56           C
ATOM  14900  C   PHE K  45      35.499  84.515  44.358  1.00 84.69           C
ATOM  15000  CB  VAL L  10      19.568  40.277  10.480  1.00 81.89           C
ATOM  15100  N   LEU L  23      25.370  37.709  29.295  1.00 48.48           N
ATOM  15200  CG  PHE L  35      15.700  26.715  36.587  1.00 57.45           C
ATOM  15300  CB  ALA M  10      14.779  24.043  33.437  1.00 50.93           C
ATOM  15400  CD1 ILE M  23      27.117  29.969  19.939  1.00 68.61           C
ATOM  15500  OG1 THR O  30      59.852  28.238  88.273  1.00 86.55           O
ATOM  15600  CB  LYS O  44      52.660  18.801  77.543  1.00 63.87           C
ATOM  15700  CG2 ILE O  58      49.231   3.465 101.214  1.00 92.18           C
ATOM  15800  CA  LEU O  71      45.849  12.711  81.045  1.00 67.23           C
ATOM  15900  CE  LYS O  83      29.412  -3.685  52.037  1.00 91.64           C
ATOM  16000  CA  LEU O  96      43.140  15.192  64.970  1.00 64.05           C
ATOM  16100  CD  GLN O 108      40.291  18.746  93.673  1.00 84.07           C
ATOM  16200  C   VAL O 122      38.086  10.051  93.162  1.00 68.13           C
ATOM  16300  CA  GLN O 135      38.800  12.567  60.896  1.00 54.37           C
ATOM  16400  C   VAL O 148      36.165  10.783  83.952  1.00 59.64           C
ATOM  16500  CA  ILE O 162      53.767   2.333  94.433  1.00 84.28           C
ATOM  16600  CG2 VAL O 174      32.148   7.788  78.609  1.00 58.15           C
ATOM  16700  N   GLY O 187      16.598  24.287  64.445  1.00 52.49           N
ATOM  16800  CD  PRO O 201       5.762  22.790  55.974  1.00 49.53           C
ATOM  16900  CG  LYS O 214      16.607  16.684  70.039  1.00 46.98           C
ATOM  17000  CA  VAL O 227      52.726   3.511  85.940  1.00 79.62           C
ATOM  17100  CA  GLU O 242      44.541   3.458  77.560  1.00 70.49           C
ATOM  17200  O   GLU O 255      26.874   6.548  65.762  1.00 58.69           O
ATOM  17300  CA  ALA O 267      49.743  12.180  89.104  1.00 75.32           C
ATOM  17400  CA  PHE T   8      30.107  24.156  33.778  1.00 58.16           C
ATOM  17500  CE2 PHE T  19      37.073  30.806  23.807  1.00 66.58           C
ATOM  17600  O   LYS T  31      38.930  53.026   9.904  1.00104.23           O
ATOM  17700  O   GLY U  48       2.470  14.492  84.844  1.00 53.19           O
ATOM  17800  OD1 ASN U  61      17.430  27.269  88.022  1.00 64.67           O
ATOM  17900  CG  PRO U  73      20.509  34.023  98.785  1.00 68.36           C
ATOM  18000  CG  GLU U  86      13.983   7.114  98.623  1.00 75.48           C
ATOM  18100  N   ARG U 100      16.728  26.887 112.613  1.00 79.81           N
ATOM  18200  CA  HIS U 111       5.985  23.384  97.418  1.00 66.26           C
ATOM  18300  OE1 GLU U 123      25.495  13.972  83.507  1.00 74.96           O
ATOM  18400  O   GLU V  28      -1.415  59.687  86.032  1.00 69.71           O
ATOM  18500  N   GLY V  42       1.728  48.438 106.910  1.00 77.83           N
ATOM  18600  O   GLU V  54      16.395  61.529  91.352  1.00 67.27           O
ATOM  18700  C   HIS V  67      18.618  48.373  86.041  1.00 51.86           C
ATOM  18800  CA  ARG V  81       9.654  43.737  90.640  1.00 59.19           C
ATOM  18900  CA  ASN V  94       8.577  49.578 100.936  1.00 64.63           C
ATOM  19000  OG1 THR V 106      19.762  41.352 102.868  1.00 66.61           O
ATOM  19100  CA  PRO V 119      27.043  49.726  96.368  1.00 53.71           C
ATOM  19200  NE  ARG V 131      34.037  57.790  98.066  1.00 77.69           N
ATOM  19300  N   ILE V 145      11.653  52.943  92.598  1.00 57.53           N
ATOM  19400  CA  GLY V 157      10.349  43.420  85.728  1.00 62.22           C
ATOM  19500  C   ILE y  25      27.680  83.976  70.831  1.00 93.60           C
ATOM  19600  C   LEU y  39      36.215  96.397  55.160  1.00 88.39           C
ATOM  19700  C   LYS X  17     -23.484  62.051  48.307  1.00 73.05           C
ATOM  19800  C   GLY X  31     -12.949  75.783  36.869  1.00 78.68           C
ATOM  19900  N   LYS X  45       0.108  87.147  24.605  1.00 94.37           N
ATOM  20000  CB  UNK Y  15      -2.576  94.042  54.728  1.00100.98           C
ATOM  20100  C   PHE Z   5      40.041  89.364  89.872  1.00 91.94           C
ATOM  20200  O   VAL Z  18      39.978  94.277  69.809  1.00 76.44           O
ATOM  20300  CB  ASP Z  32      47.282 103.853  60.501  1.00100.60           C
ATOM  20400  N   SER Z  44      48.888  94.580  73.401  1.00 86.76           N
ATOM  20500  CD2 LEU Z  57      44.301  89.773  91.041  1.00 78.59           C
HETATM20600  C2  CLA A 362      31.764  47.668  47.431  1.00 60.44           C
HETATM20700  NB  CLA A 363      29.710  41.445  46.930  1.00 45.58           N
HETATM20800  C7  CLA D 364      16.133  64.229  46.053  1.00 76.65           C
HETATM20900  CBB PHO D 355      14.716  56.516  46.747  1.00 56.43           C
HETATM21000  C14 CLA A 366      66.544  40.085  54.813  1.00 90.68           C
HETATM21100  C4C CLA B 511     -26.404  52.016  32.721  1.00 96.16           C
HETATM21200  C19 CLA B 512      -3.098  50.381  41.085  1.00 83.96           C
HETATM21300  C1D CLA B 514      -5.234  38.699  29.202  1.00 71.34           C
HETATM21400  CHC CLA B 516     -20.209  32.957  21.111  1.00 89.96           C
HETATM21500  CAD CLA B 517       2.554  21.339  34.550  1.00 55.30           C
HETATM21600  C3A CLA B 519     -10.472  56.508  18.414  1.00 84.57           C
HETATM21700  O2D CLA B 520     -10.650  46.894  16.518  1.00 72.71           O
HETATM21800  CGA CLA B 522      -3.932  39.345  15.732  1.00 83.81           C
HETATM21900  C4  CLA B 523       8.852  43.002  22.934  1.00 70.48           C
HETATM22000  C2B CLA B 525     -14.006  41.766   8.133  1.00 97.68           C
HETATM22100  C9  CLA B 526     -11.824  25.329  -0.631  1.00106.01           C
HETATM22200  CBB CLA C 475      39.928  67.338  64.645  1.00 59.03           C
HETATM22300  C14 CLA C 476      50.357  62.165  64.457  1.00 87.49           C
HETATM22400  C4C CLA C 478      53.106  55.767  51.571  1.00 72.93           C
HETATM22500  C19 CLA C 479      59.949  48.251  61.645  1.00100.03           C
HETATM22600  C1D CLA C 481      42.030  67.462  48.801  1.00 85.05           C
HETATM22700  CHC CLA C 483      44.314  69.269  63.901  1.00 50.61           C
HETATM22800  CAD CLA C 484      43.342  83.454  49.029  1.00 82.21           C
HETATM22900  C3A CLA C 486      61.049  84.971  60.287  1.00 97.06           C
HETATM23000 FE   HEM V 164      21.702  46.175  93.090  1.00 67.31          FE
HETATM23100  C4  PL9 A 367      23.534  65.392  32.049  1.00 87.59           C
HETATM23200  C14 BCR A 369      49.813  31.115  42.316  1.00 69.53           C
HETATM23300  C38 BCR B 528       4.937  19.657  26.930  1.00 79.68           C
HETATM23400  C14 BCR H 107     -21.140  55.848  30.373  1.00 89.39           C
HETATM23500  C38 BCR K 112      25.295  79.179  53.876  1.00 74.98           C
HETATM23600  C14 BCR Z 116      52.931  81.451  77.107  1.00 77.87           C
HETATM23700  C3D DGD A 370      53.637  32.456  61.804  1.00 85.34           C
HETATM23800  C28 LHG A 371      35.937  62.504  44.119  1.00 65.86           C
HETATM23900  C3G DGD C 490      27.108  67.075  73.164  1.00 77.16           C
HETATM24000  C6  LMG D 359       9.071  67.359  68.690  1.00 71.69           C
HETATM24100  C4  LMG B 531      19.088  52.863  20.315  1.00 79.55           C
HETATM24200  C4  LMG D 360      31.370  49.703  21.577  1.00 74.79           C
HETATM24300  O47 SQD T 213      39.292  32.149  16.094  1.00 89.97           O
HETATM24400  C5  LMG B 534      -5.220   3.042  32.449  1.00 88.99           C
HETATM24500  C16 LMG M 217      31.439  29.591  17.567  1.00 80.15           C
HETATM24600  C32 LMG C 492      29.124  71.911  68.040  1.00 81.16           C
HETATM24700  C1  LHG A 374      26.048  76.110  36.215  1.00 94.79           C
HETATM24800  C13 SQD F 224      -4.639  79.397  39.277  1.00 94.00           C
HETATM24900  C5B LMT M 226      18.726  10.738  43.575  1.00 89.50           C
HETATM25000  C12 LMT B 535     -16.566  32.829  12.754  1.00 71.46           C
HETATM25100  C11 LMT B 536     -11.459  65.368  27.690  1.00 83.60           C
"""

def amplitudes_from_pdb(resolution,algorithm=None,anomalous=False):
  from iotbx import pdb
  pdb_inp = pdb.input(source_info=None,lines = pdb_lines)
  xray_structure = pdb_inp.xray_structure_simple()
  primitive_xray_structure = xray_structure.primitive_setting()
  P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
  fcalc = P1_primitive_xray_structure.structure_factors(
    d_min=resolution, anomalous_flag=anomalous, algorithm=algorithm).f_calc()
  return fcalc.amplitudes()

def tst_all():
  F = foo()
  assert F == (1,2,3,4)

  size=1516
  D = detector(slowpixels=1516,fastpixels=1516,pixel_size=0.00011)
  D.set_region_of_interest(0,int(0.6*size),int(0.4*size),int(0.6*size))
  D.set_oversampling(1)

  C = camera()
  C.distance = 0.18166
  C.Ybeam = 0.08338
  C.Zbeam = 0.08338
  C.lambda0 = 6.2E-10
  C.dispersion = 0.002
  C.dispsteps = 4
  C.hdivrange = 0
  C.vdivrange = 0
  C.hdivstep = 1
  C.vdivstep = 1
  C.source_distance = 10.
  C.fluence = 1.E24

  Amat = sqr((127.6895065259495, 0.6512077339887,-0.4403031342553,
               -1.1449112128916,225.3922539826207,1.8393136632579,
               1.0680694468752,-2.4923062985132,306.0953037195841))

  PSII = amplitudes_from_pdb(8.,"fft",True)

  from cctbx import crystal_orientation
  X = crystal()
  X.orientation = crystal_orientation.crystal_orientation(
                  Amat,crystal_orientation.basis_type.direct)
  X.miller = PSII.indices()
  X.amplitudes = PSII.data()
  X.Na = 6; X.Nb = 6; X.Nc = 6

  SIM = fast_bragg_simulation()
  SIM.set_detector(D)
  SIM.set_camera(C)
  SIM.set_crystal(X)
  SIM.sweep_over_detector()
  data = D.raw
  scale_factor = 55000./flex.max(data)
  #print "scale_factor",scale_factor
  fileout="intimage_001.img"
  SIM.to_smv_format(fileout=fileout,
                    intfile_scale = scale_factor, saturation=40000)
  import os
  assert os.path.isfile(fileout)
  os.remove(fileout)

  #simulation is complete, now we'll autoindex the image fragment and verify
  # that the indexed cell is similar to the input cell.

  if (not libtbx.env.has_module("annlib")):
    print "Skipping some tests: annlib not available."
    return
  # 1. Analysis of the image to identify the Bragg peak centers.
  #    This step uses an inefficient algorithm and implementation and
  #    is most time consuming; but the code is only for testing, not production
  from rstbx.diffraction.fastbragg.tst_utils_clustering import specific_libann_cluster
  M=specific_libann_cluster(scale_factor*data,intensity_cutoff = 25,distance_cutoff=17)
  # M is a dictionary of peak intensities indexed by pixel coordinates

  # 2. Now autoindex the pattern
  from rstbx.diffraction.fastbragg.tst_utils_clustering import index_wrapper
  SIM.C = C
  SIM.D = D
  ai,ref_uc = index_wrapper(M.keys(), SIM, PSII)
  tst_uc = ai.getOrientation().unit_cell()
  #print ref_uc  # (127.692, 225.403, 306.106, 90, 90, 90)
  #print tst_uc  # (106.432, 223.983, 303.102, 90.3185, 91.5998, 90.5231)

  # 3. Final assertion.  In the given orientation,
  #  the unit cell A vector is into the beam and is not well sampled,
  #  so tolerances have to be fairly relaxed, 5%.
  #  Labelit does better with the target_cell restraint, but this improved
  #  algorithm is not used here for the test
  assert ref_uc.is_similar_to(tst_uc, relative_length_tolerance=0.20,
                                     absolute_angle_tolerance= 2.0)


if __name__=="__main__":
  tst_all()
  print "OK"
