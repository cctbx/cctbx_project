from __future__ import absolute_import, division, print_function
import time
from mmtbx.validation import validate_ligands
from mmtbx.validation.validate_ligands import master_params_str
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal

pdb_str_1 = """
CRYST1   26.971   23.398   30.626  90.00  90.00  90.00 P 1
SCALE1      0.037077  0.000000  0.000000        0.00000
SCALE2      0.000000  0.042739  0.000000        0.00000
SCALE3      0.000000  0.000000  0.032652        0.00000
ATOM    230  CB  PHE A  37     -11.978  68.427 -40.130  1.00 30.99           C
ATOM    231  CG  PHE A  37     -13.053  67.383 -40.232  1.00 32.04           C
ATOM    232  CD1 PHE A  37     -13.355  66.792 -41.449  1.00 35.89           C
ATOM    233  CD2 PHE A  37     -13.792  67.026 -39.114  1.00 30.41           C
ATOM    234  CE1 PHE A  37     -14.383  65.855 -41.558  1.00 39.57           C
ATOM    235  CE2 PHE A  37     -14.825  66.094 -39.199  1.00 30.06           C
ATOM    236  CZ  PHE A  37     -15.122  65.505 -40.426  1.00 35.49           C
ATOM    307  CB  VAL A  46     -18.628  66.498 -36.433  1.00 36.27           C
ATOM    308  CG1 VAL A  46     -17.310  66.316 -35.689  1.00 41.89           C
ATOM    309  CG2 VAL A  46     -19.693  65.584 -35.848  1.00 38.02           C
ATOM    342  CB  PHE A  50      -7.537  67.137 -39.092  1.00 17.63           C
ATOM    343  CG  PHE A  50      -7.384  67.839 -40.418  1.00 19.74           C
ATOM    344  CD1 PHE A  50      -6.411  67.439 -41.328  1.00 21.42           C
ATOM    345  CD2 PHE A  50      -8.188  68.924 -40.738  1.00 19.01           C
ATOM    346  CE1 PHE A  50      -6.241  68.116 -42.543  1.00 23.86           C
ATOM    347  CE2 PHE A  50      -8.030  69.609 -41.944  1.00 21.42           C
ATOM    348  CZ  PHE A  50      -7.049  69.200 -42.850  1.00 20.02           C
ATOM    364  CB  TRP A  52      -3.350  65.049 -42.868  1.00 19.63           C
ATOM    365  CG  TRP A  52      -3.229  64.112 -44.010  1.00 22.71           C
ATOM    366  CD1 TRP A  52      -2.722  64.391 -45.251  1.00 23.99           C
ATOM    367  CD2 TRP A  52      -3.643  62.743 -44.045  1.00 25.36           C
ATOM    368  NE1 TRP A  52      -2.800  63.278 -46.056  1.00 22.55           N
ATOM    369  CE2 TRP A  52      -3.364  62.253 -45.342  1.00 23.11           C
ATOM    370  CE3 TRP A  52      -4.230  61.882 -43.109  1.00 30.70           C
ATOM    371  CZ2 TRP A  52      -3.643  60.938 -45.721  1.00 33.42           C
ATOM    372  CZ3 TRP A  52      -4.508  60.574 -43.489  1.00 36.99           C
ATOM    373  CH2 TRP A  52      -4.218  60.118 -44.784  1.00 36.75           C
ATOM    458  CB  LEU A  62      -5.758  60.896 -39.665  1.00 23.50           C
ATOM    459  CG  LEU A  62      -6.916  61.887 -39.482  1.00 27.86           C
ATOM    460  CD1 LEU A  62      -6.520  63.234 -40.050  1.00 29.46           C
ATOM    461  CD2 LEU A  62      -8.175  61.370 -40.176  1.00 25.44           C
ATOM    474  CB  LEU A  64      -9.315  64.303 -35.082  1.00 26.03           C
ATOM    475  CG  LEU A  64      -9.424  63.795 -36.518  1.00 25.34           C
ATOM    476  CD1 LEU A  64     -10.622  62.864 -36.629  1.00 33.60           C
ATOM    477  CD2 LEU A  64      -9.566  64.972 -37.469  1.00 25.46           C
ATOM    532  OE1 GLU A  72     -12.523  65.769 -32.644  1.00 48.55           O
ATOM    539  CG  LEU A  73     -10.375  59.956 -34.221  1.00 37.03           C
ATOM    540  CD1 LEU A  73      -9.368  59.507 -35.269  1.00 31.81           C
ATOM    541  CD2 LEU A  73     -11.727  59.325 -34.515  1.00 30.45           C
ATOM    755  CB  ILE A 101     -12.754  56.211 -40.650  1.00 32.70           C
ATOM    756  CG1 ILE A 101     -11.239  56.390 -40.490  1.00 39.17           C
ATOM    757  CG2 ILE A 101     -13.443  57.571 -40.510  1.00 34.03           C
ATOM    758  CD1 ILE A 101     -10.836  57.124 -39.223  1.00 40.02           C
ATOM    781  CA  THR A 105     -17.017  57.403 -35.722  1.00 34.56           C
ATOM    782  C   THR A 105     -18.439  57.668 -36.183  1.00 31.79           C
ATOM    783  O   THR A 105     -19.019  58.696 -35.843  1.00 33.82           O
ATOM    784  CB  THR A 105     -16.093  58.393 -36.466  1.00 35.10           C
ATOM    785  OG1 THR A 105     -16.364  58.347 -37.870  1.00 25.12           O
ATOM    786  CG2 THR A 105     -14.635  58.044 -36.223  1.00 35.53           C
ATOM    795  CB  MET A 107     -18.589  60.322 -40.140  1.00 40.28           C
ATOM    796  CG  MET A 107     -19.438  61.415 -39.493  1.00 49.13           C
ATOM    797  SD  MET A 107     -18.470  62.795 -38.823  1.00 53.68           S
ATOM    798  CE  MET A 107     -18.197  62.131 -37.271  1.00 48.80           C
ATOM    807  CA  ILE A 109     -15.736  59.538 -44.760  1.00 34.00           C
ATOM    810  CB  ILE A 109     -14.662  59.947 -43.731  1.00 32.38           C
ATOM    811  CG1 ILE A 109     -15.322  60.670 -42.546  1.00 39.13           C
ATOM    812  CG2 ILE A 109     -13.619  60.833 -44.387  1.00 30.68           C
ATOM    813  CD1 ILE A 109     -16.068  61.936 -42.922  1.00 38.08           C
ATOM    824  O   LEU A 111      -9.007  58.605 -48.267  1.00 27.74           O
ATOM    826  CG  LEU A 111      -9.184  57.521 -45.004  1.00 29.09           C
ATOM    827  CD1 LEU A 111      -9.009  58.934 -44.480  1.00 25.47           C
ATOM    844  CA  LEU A 114      -9.615  61.831 -50.012  1.00 25.43           C
ATOM    845  C   LEU A 114      -8.086  61.843 -49.971  1.00 24.51           C
ATOM    846  O   LEU A 114      -7.461  62.897 -50.102  1.00 30.03           O
ATOM    847  CB  LEU A 114     -10.188  61.983 -48.602  1.00 29.73           C
ATOM    848  CG  LEU A 114     -10.792  63.312 -48.118  1.00 38.28           C
ATOM    849  CD1 LEU A 114     -11.278  64.169 -49.282  1.00 30.02           C
ATOM    850  CD2 LEU A 114     -11.948  62.984 -47.172  1.00 23.21           C
ATOM    851  N   ARG A 115      -7.475  60.678 -49.804  1.00 19.10           N
ATOM    852  CA  ARG A 115      -6.021  60.601 -49.771  1.00 26.17           C
ATOM    853  C   ARG A 115      -5.413  61.078 -51.083  1.00 24.58           C
ATOM    854  O   ARG A 115      -4.250  61.474 -51.129  1.00 28.35           O
ATOM    855  CB  ARG A 115      -5.555  59.172 -49.473  1.00 25.43           C
ATOM    856  CG  ARG A 115      -5.841  58.743 -48.034  1.00 37.74           C
ATOM    857  CD  ARG A 115      -4.840  57.723 -47.542  1.00 33.28           C
ATOM    885  N   ILE A 118      -7.064  65.279 -51.419  1.00 14.85           N
ATOM    886  CA  ILE A 118      -6.470  66.275 -50.537  1.00 15.64           C
ATOM    889  CB  ILE A 118      -6.228  65.720 -49.115  1.00 18.62           C
ATOM    890  CG1 ILE A 118      -6.064  66.891 -48.135  1.00 18.82           C
ATOM    891  CG2 ILE A 118      -4.955  64.853 -49.084  1.00 15.18           C
ATOM    892  CD1 ILE A 118      -6.051  66.491 -46.671  1.00 22.26           C
ATOM   1362  CG  LEU A 178     -11.177  68.527 -45.827  1.00 25.53           C
ATOM   1363  CD1 LEU A 178     -10.397  67.589 -46.722  1.00 27.97           C
ATOM   1364  CD2 LEU A 178     -10.911  68.189 -44.355  1.00 24.30           C
HETATM 1447  C1  PG5 A 201      -6.715  62.650 -46.134  1.00 65.81           C
HETATM 1448  C2  PG5 A 201      -7.998  63.566 -44.320  1.00 67.25           C
HETATM 1449  C3  PG5 A 201      -9.335  64.007 -44.092  1.00 72.55           C
HETATM 1450  C4  PG5 A 201     -10.886  63.195 -42.372  1.00 72.63           C
HETATM 1451  C5  PG5 A 201     -11.575  63.550 -41.117  1.00 72.20           C
HETATM 1452  C6  PG5 A 201     -12.731  62.045 -39.328  1.00 70.18           C
HETATM 1453  C7  PG5 A 201     -13.987  62.396 -38.666  1.00 68.83           C
HETATM 1454  C8  PG5 A 201     -15.334  62.240 -36.581  1.00 68.39           C
HETATM 1455  O1  PG5 A 201      -7.873  63.394 -45.750  1.00 67.20           O
HETATM 1456  O2  PG5 A 201      -9.815  64.102 -42.646  1.00 71.78           O
HETATM 1457  O3  PG5 A 201     -12.779  62.573 -40.682  1.00 69.06           O
HETATM 1458  O4  PG5 A 201     -14.165  61.906 -37.203  1.00 71.25           O
HETATM 1469  O   HOH A 312      -2.851  61.885 -48.431  1.00 33.49           O
HETATM 1497  O   HOH A 342     -14.478  61.372 -30.793  1.00 39.31           O
"""

pdb_str_2 = '''
CRYST1   44.666   33.974   46.691  90.00  90.00  90.00 P 1
SCALE1      0.022388  0.000000  0.000000        0.00000
SCALE2      0.000000  0.029434  0.000000        0.00000
SCALE3      0.000000  0.000000  0.021417        0.00000
ATOM      1  N   ALA A   1      20.064 -26.291   0.046  1.00  6.09           N
ANISOU    1  N   ALA A   1     1005    711    597     62    -25     96       N
ATOM      2  CA  ALA A   1      21.152 -26.951   0.767  1.00  6.82           C
ANISOU    2  CA  ALA A   1     1214    700    679    158   -101    161       C
ATOM      3  C   ALA A   1      21.867 -26.032   1.762  1.00  6.66           C
ANISOU    3  C   ALA A   1     1031    833    668     97    -32    147       C
ATOM      4  O   ALA A   1      22.541 -26.551   2.653  1.00  8.52           O
ANISOU    4  O   ALA A   1     1495    989    753    114   -325    229       O
ATOM      5  CB  ALA A   1      22.184 -27.500  -0.229  1.00  8.10           C
ANISOU    5  CB  ALA A   1     1318   1031    728    400    -97    138       C
ATOM      6  N  ASER A   2      21.784 -24.715   1.626  0.64  6.65           N
ANISOU    6  N  ASER A   2      908    873    746     67   -338    145       N
ATOM      7  CA ASER A   2      22.466 -23.794   2.534  0.64  6.81           C
ANISOU    7  CA ASER A   2      781    869    936    -95   -250     93       C
ATOM      8  C  ASER A   2      21.522 -22.817   3.228  0.64  7.07           C
ANISOU    8  C  ASER A   2      832   1080    773   -170    -40     -8       C
ATOM      9  O  ASER A   2      21.992 -21.784   3.761  0.64  8.65           O
ANISOU    9  O  ASER A   2      797   1255   1236   -163    -33   -274       O
ATOM     10  CB ASER A   2      23.532 -22.969   1.784  0.64  9.30           C
ANISOU   10  CB ASER A   2      718   1161   1654    -41     68     -8       C
ATOM     11  OG ASER A   2      24.587 -23.797   1.419  0.64  9.73           O
ANISOU   11  OG ASER A   2      843   1391   1463    123    -63   -198       O
ATOM     12  N  BSER A   2      21.645 -24.736   1.517  0.36  7.69           N
ANISOU   12  N  BSER A   2     1224    684   1014     97    291    115       N
ATOM     13  CA BSER A   2      22.383 -23.726   2.275  0.36  7.53           C
ANISOU   13  CA BSER A   2      909    934   1017    244   -208    141       C
ATOM     14  C  BSER A   2      21.537 -22.892   3.226  0.36  6.35           C
ANISOU   14  C  BSER A   2      760    834    817    -36   -325     57       C
ATOM     15  O  BSER A   2      22.068 -22.056   3.980  0.36  7.61           O
ANISOU   15  O  BSER A   2      940    939   1014   -276   -257     -4       O
ATOM     16  CB BSER A   2      23.111 -22.858   1.207  0.36  8.82           C
ANISOU   16  CB BSER A   2      798    957   1598     -3    112     16       C
ATOM     17  OG BSER A   2      23.874 -23.715   0.380  0.36 12.89           O
ANISOU   17  OG BSER A   2     1061   2019   1817     72    270   -413       O
ATOM     18  N   ASN A   3      20.220 -23.070   3.223  1.00  6.74           N
ANISOU   18  N   ASN A   3      777   1004    779   -106   -141    -48       N
ATOM     19  CA  ASN A   3      19.257 -22.226   3.915  1.00  6.97           C
ANISOU   19  CA  ASN A   3      846   1000    801    -35   -127    -59       C
ATOM     20  C   ASN A   3      18.270 -23.004   4.725  1.00  8.43           C
ANISOU   20  C   ASN A   3      911   1411    882    -80    -48     90       C
ATOM     21  O   ASN A   3      17.296 -22.364   5.219  1.00 10.52           O
ANISOU   21  O   ASN A   3     1121   1810   1068    223     63    131       O
ATOM     22  CB  ASN A   3      18.498 -21.312   2.907  1.00  6.84           C
ANISOU   22  CB  ASN A   3      745   1002    853    -46    -62   -115       C
ATOM     23  CG  ASN A   3      19.364 -20.233   2.409  1.00  6.98           C
ANISOU   23  CG  ASN A   3      890    852    911    -55   -124   -186       C
ATOM     24  OD1 ASN A   3      19.488 -19.190   3.069  1.00 10.15           O
ANISOU   24  OD1 ASN A   3     1320   1066   1471   -229    130   -580       O
ATOM     25  ND2 ASN A   3      20.003 -20.475   1.265  1.00  6.72           N
ANISOU   25  ND2 ASN A   3      995    684    873   -198    -61    -65       N
ATOM     26  OXT ASN A   3      18.421 -24.209   4.852  1.00 11.11           O
ANISOU   26  OXT ASN A   3     1492   1291   1439   -209    424     45       O
TER
HETATM   27 CA    CA A   1     -10.015 -24.840   3.178  1.00  4.60          Ca
ANISOU   27 CA    CA A   1      577    628    544     34     81     11      Ca
HETATM   28  C  ABEN A   2      -2.188 -15.645 -18.318  0.56  4.74           C
ANISOU   28  C  ABEN A   2      725    387    690     11     36     96       C
HETATM   29  C1 ABEN A   2      -1.864 -14.876 -17.107  0.56  5.08           C
ANISOU   29  C1 ABEN A   2      634    542    755     95    101     64       C
HETATM   30  C2 ABEN A   2      -1.559 -15.529 -15.913  0.56  5.32           C
ANISOU   30  C2 ABEN A   2      737    539    746     54    -14    -71       C
HETATM   31  C3 ABEN A   2      -1.287 -14.833 -14.751  0.56  6.27           C
ANISOU   31  C3 ABEN A   2      913    606    863    109    -70   -181       C
HETATM   32  C4 ABEN A   2      -1.265 -13.455 -14.805  0.56  7.69           C
ANISOU   32  C4 ABEN A   2     1297    640    986    247   -169   -287       C
HETATM   33  C5 ABEN A   2      -1.473 -12.777 -15.977  0.56  8.11           C
ANISOU   33  C5 ABEN A   2     1388    670   1023    205   -175   -134       C
HETATM   34  C6 ABEN A   2      -1.797 -13.467 -17.152  0.56  6.61           C
ANISOU   34  C6 ABEN A   2      961    537   1013     88   -106   -103       C
HETATM   35  N1 ABEN A   2      -1.908 -16.952 -18.415  0.56  4.87           N
ANISOU   35  N1 ABEN A   2      568    432    850     73     61    120       N
HETATM   36  N2 ABEN A   2      -2.750 -15.113 -19.412  0.56  5.75           N
ANISOU   36  N2 ABEN A   2      851    507    825    188     21    163       N
HETATM   37  C  BBEN A   2      -2.149 -15.627 -18.289  0.44  5.73           C
ANISOU   37  C  BBEN A   2      682    690    807     71    -33    -98       C
HETATM   38  C1 BBEN A   2      -1.870 -14.923 -17.020  0.44  5.50           C
ANISOU   38  C1 BBEN A   2      777    488    824     12    -79   -124       C
HETATM   39  C2 BBEN A   2      -0.938 -15.394 -16.117  0.44  5.87           C
ANISOU   39  C2 BBEN A   2      685    625    919     10     -8     21       C
HETATM   40  C3 BBEN A   2      -0.697 -14.731 -14.934  0.44  6.98           C
ANISOU   40  C3 BBEN A   2      995    682    974    -81   -172     74       C
HETATM   41  C4 BBEN A   2      -1.449 -13.622 -14.603  0.44  8.17           C
ANISOU   41  C4 BBEN A   2     1035   1144    927    126     14   -161       C
HETATM   42  C5 BBEN A   2      -2.368 -13.128 -15.507  0.44  7.70           C
ANISOU   42  C5 BBEN A   2     1169    678   1078    219    -95   -348       C
HETATM   43  C6 BBEN A   2      -2.577 -13.756 -16.763  0.44  6.82           C
ANISOU   43  C6 BBEN A   2      940    595   1057    135   -140   -150       C
HETATM   44  N1 BBEN A   2      -1.607 -16.841 -18.438  0.44  5.07           N
ANISOU   44  N1 BBEN A   2      572    654    702    -27     32   -187       N
HETATM   45  N2 BBEN A   2      -2.983 -15.074 -19.170  0.44  5.46           N
ANISOU   45  N2 BBEN A   2      770    521    785    -11   -172   -163       N
HETATM   46  O1  SO4 A   3      -1.816 -12.816 -11.138  0.65  8.54           O
HETATM   47  O2  SO4 A   3      -0.023 -11.767 -10.034  0.65 11.23           O
HETATM   48  O3  SO4 A   3      -2.003 -10.425 -10.555  0.65 10.57           O
HETATM   49  O4  SO4 A   3      -0.441 -11.129 -12.259  0.65 13.14           O
HETATM   50  S   SO4 A   3      -1.151 -11.482 -11.029  0.65  7.42           S
HETATM   51  O1  SO4 A   4      13.061 -11.888 -31.472  0.48 14.55           O
ANISOU   51  O1  SO4 A   4     1714   1932   1884    366    830    561       O
HETATM   52  O2  SO4 A   4      11.726 -13.037 -29.761  0.48 10.45           O
ANISOU   52  O2  SO4 A   4     1483   1272   1216    302    392    327       O
HETATM   53  O3  SO4 A   4      10.828 -11.094 -30.953  0.48 10.27           O
ANISOU   53  O3  SO4 A   4     1547   1072   1285    320    402    154       O
HETATM   54  O4  SO4 A   4      12.579 -10.820 -29.293  0.48 14.33           O
ANISOU   54  O4  SO4 A   4     2083   1487   1874     52     33    388       O
HETATM   55  S   SO4 A   4      12.076 -11.776 -30.446  0.48 12.03           S
ANISOU   55  S   SO4 A   4     1400   1230   1939    310    477    589       S
HETATM   56  C1  GOL A   5      22.489  -3.691  -9.702  0.67 58.70           C
HETATM   57  C2  GOL A   5      21.482  -4.248 -10.663  0.67110.35           C
HETATM   58  C3  GOL A   5      20.064  -4.076 -10.144  0.67111.50           C
HETATM   59  O1  GOL A   5      23.893  -3.768 -10.129  0.67 99.94           O
HETATM   60  O2  GOL A   5      21.641  -3.526 -11.889  0.67114.91           O
HETATM   61  O3  GOL A   5      19.180  -4.727 -11.061  0.67 86.14           O
HETATM   62  O   HOH A   6      -9.233 -22.998   4.376  0.98  5.96           O
ANISOU   62  O   HOH A   6      791    803    669    -21     61   -177       O
HETATM   63  O   HOH A   7     -10.079 -26.946   2.058  0.97  5.27           O
ANISOU   63  O   HOH A   7      743    599    662    -80    112      1       O
HETATM   64  O   HOH A   8      -1.938 -21.274  -8.865  1.00  4.16           O
ANISOU   64  O   HOH A   8      477    534    569    -13     18   -121       O
HETATM   65  O   HOH A   9      11.300 -16.164 -16.305  0.98  5.35           O
ANISOU   65  O   HOH A   9      793    561    680     46    -95    -87       O
END
'''

filenames = ['test_1.pdb', 'test_2.pdb']
pdb_strings = [pdb_str_1, pdb_str_2]


def write_model_files(filenames, pdb_strings):
  for filename, pdb_str in zip(filenames, pdb_strings):
    f = open("%s" % filename, "w")
    f.write(pdb_str)
    f.close()


def tst_get_adps(vl_manager):
  '''
  Test getting ADPs of ligands and surrounding atoms
  '''
  n_iso_answer = (0,0,0,0,0)
  n_aniso_answer = (9,9,5,5,6)
  #print(vl_manager)
  for id_tuple, ligand_dict in vl_manager.items():
    #print(ligand_dict)
    for altloc, lr in ligand_dict.items():
      adps = lr.get_adps()
      id_str = lr.id_str
      if (id_str.strip() == 'A BEN    2 A'):
        assert(adps.n_iso == 0)
        assert(adps.n_aniso == 9)
        assert(adps.n_above_100 == 0)
        assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
          [4.7, 8.1, 6.0], eps=0.1)
        assert(approx_equal(
          [adps.b_min_within, adps.b_max_within, adps.b_mean_within],
          [6.1,12.9,7.8], eps=0.1))
      if (id_str.strip() == 'A BEN    2 B'):
        assert(adps.n_iso == 0)
        assert(adps.n_aniso == 9)
        assert(adps.n_above_100 == 0)
        assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
          [5.1, 8.2, 6.4], eps=0.1)
        assert(approx_equal(
          [adps.b_min_within, adps.b_max_within, adps.b_mean_within],
          [6.1,12.9,7.8], eps=0.1))
      if (id_str.strip() == 'A SO4    3'):
        assert(adps.n_iso == 5)
        assert(adps.n_aniso == 0)
        assert(adps.n_above_100 == 0)
        assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
          [7.4,13.1,10.2], eps=0.1)
        assert(approx_equal(
          [adps.b_min_within, adps.b_max_within, adps.b_mean_within],
          [6.3, 11.1, 7.8], eps=0.1))
      if (id_str.strip() == 'A SO4    4'):
        assert(adps.n_iso == 0)
        assert(adps.n_aniso == 5)
        assert(adps.b_min_within is None)
        assert(adps.n_above_100 == 0)
        assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
          [10.3,14.6,12.3], eps=0.1)
      if (id_str.strip() == 'A GOL    5'):
        assert(adps.n_iso == 6)
        assert(adps.n_aniso == 0)
        assert(adps.b_min_within is None)
        assert(adps.n_above_100 == 3)
        assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
          [58.7,114.9,96.9], eps=0.1)


def run_test1():
  '''
  Test if iselection for ligand PG5 (chain A resseq 201) is correct.
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_1.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input = pdb_inp)

  params = iotbx.phil.parse(
    input_string=master_params_str, process_includes=True).extract()
  # do not place H atoms for this test
  #params.validate_ligands.place_hydrogens = False

  fn = filenames[0]

  vl_manager = validate_ligands.manager(
    model = model,
    model_fn = fn,
    fmodel = None,
    params = params.validate_ligands,
    log   = null_out)
  vl_manager.run()

  tst_get_ligands(vl_manager = vl_manager)
  tst_get_overlaps(vl_manager = vl_manager)


def tst_get_ligands(vl_manager):
  '''
  Test finding ligand
  '''
  assert (len(vl_manager) == 1)
  # test iselection
  for id_tuple, ligand_dict in vl_manager.items():
    assert (id_tuple == ('', 'A', ' 201'))
    lr = ligand_dict['']
    assert (list(lr.isel) == [84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95])


def tst_get_overlaps(vl_manager):
  '''
  Test nonbonded overlaps
  '''
  for id_tuple, ligand_dict in vl_manager.items():
    for altloc, lr in ligand_dict.items():
      clashes_result = lr.get_overlaps()
      print(clashes_result.clashscore)
      assert(clashes_result.n_clashes == 5)
      assert approx_equal(clashes_result.clashscore, 31.6, eps=1.0)
# anaconda
#(['pdb=" HE3 MET A 107 "', 'pdb=" H81 PG5 A 201 "'], 17, 54, 2.0370952358689647, 2.44, '', None),
#(['pdb=" CE  MET A 107 "', 'pdb=" C8  PG5 A 201 "'], 15, 34, 2.946989989803154, 3.4, '', None),
#(['pdb=" CE  MET A 107 "', 'pdb=" H83 PG5 A 201 "'], 15, 56, 2.4839921497460486, 2.92, '', None)
#
# MAC
#(['pdb=" CE  MET A 107 "', 'pdb=" C8  PG5 A 201 "'], 16, 35, 2.946989989803154, 3.4, '', None),
#(['pdb=" HE3 MET A 107 "', 'pdb=" H83 PG5 A 201 "'], 18, 57, 2.026073542594147, 2.44, '', None),
#(['pdb=" CE  MET A 107 "', 'pdb=" H81 PG5 A 201 "'], 16, 55, 2.4973179613337146, 2.92, '', None)
#
# Readyset gives different names to H atoms.


def tst_get_occupancies(vl_manager):
  '''
  Test occupancy determination
  '''
  assert (len(vl_manager) == 4)
  id_tuple_answer = [('', 'A', '   2'), ('', 'A', '   3'), ('', 'A', '   4'), ('', 'A', '   5')]
  ligand_dict_length_answer = [2, 1, 1, 1]
  occupancy_answer = []
  for id_tuple, id_tuple_answer, length_answer in zip(vl_manager.keys(), id_tuple_answer, ligand_dict_length_answer):
    ligand_dict = vl_manager[id_tuple]
    assert (id_tuple == id_tuple_answer)
    assert (len(ligand_dict) == length_answer)
    for altloc, lr in ligand_dict.items():
      occs = lr.get_occupancies()
      id_str = lr.id_str
      if (id_str.strip() == 'A BEN    2 A'):
        assert(occs.occ_mean == 0.56)
      if (id_str.strip() == 'A BEN    2 B'):
        assert(occs.occ_mean == 0.44)
      if (id_str.strip() == 'A SO4    3'):
        assert(occs.occ_mean == 0.65)
      if (id_str.strip() == 'A SO4    4'):
        assert(occs.occ_mean == 0.48)
      if (id_str.strip() == 'A GOL    5'):
        assert(occs.occ_mean == 0.67)


def run_test2():
  '''
  Test
  - occupancy determination for ligands
  - adp determination for ligands and neighbors
  Tests are combined to decrease computing time (restraints manager is slow).
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str_2.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input = pdb_inp)
  params = iotbx.phil.parse(
    input_string=master_params_str, process_includes=True).extract()
  fn = filenames[1]
  vl_manager = validate_ligands.manager(
    model = model,
    model_fn = fn,
    fmodel = None,
    params = params.validate_ligands,
    log   = null_out)
  vl_manager.run()

  tst_get_occupancies(vl_manager = vl_manager)
  tst_get_adps(vl_manager = vl_manager)


def run():
  write_model_files(filenames = filenames, pdb_strings = pdb_strings)
  run_test1()
  run_test2()


if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f" % (time.time()-t0))
