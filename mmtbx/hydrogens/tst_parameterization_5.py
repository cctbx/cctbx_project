from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from six.moves import zip

#-----------------------------------------------------------------------------
# This test checks the parameterization of H atoms for multiple conformations
# These examples are from pdb and initially failed
#-----------------------------------------------------------------------------

def exercise(pdb_str, type_list_known):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(make_restraints=True)
  pdb_hierarchy = model.get_hierarchy()
  sites_cart = model.get_sites_cart()
  atoms = pdb_hierarchy.atoms()

  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()

  h_para = riding_h_manager.h_parameterization

  diagnostics = riding_h_manager.diagnostics(
    sites_cart = sites_cart,
    threshold  = 0.05)
  h_distances   = diagnostics.h_distances
  type_list     = diagnostics.type_list

  number_h = model.get_hd_selection().count(True)
  number_h_para = len(h_para) - h_para.count(None)

  assert (number_h_para == number_h), 'Not all H atoms are parameterized'

  for ih in h_distances:
    labels = atoms[ih].fetch_labels()
    assert (h_distances[ih] < 0.2), \
      'distance too large: %s  atom: %s (%s) residue: %s ' \
      % (h_para[ih].htype, atoms[ih].name, ih, labels.resseq.strip())

  for type1, type2 in zip(type_list, type_list_known):
    assert (type1 == type2)
    #print "'%s'," % type1,

# Several fragments from pdb with double conformations which caused crashes
# There are 68 H atoms in the pdb_string, check if all of them are recognized
# Note: one H atom (VAL 7 HA) is bound to two CA atoms at once
pdb_str1 = """\
CRYST1   27.832   30.931   25.475  90.00  90.00  90.00 P 1
SCALE1      0.035930  0.000000  0.000000        0.00000
SCALE2      0.000000  0.032330  0.000000        0.00000
SCALE3      0.000000  0.000000  0.039254        0.00000
ATOM      1  CB  SER A   1      10.278  16.843  12.582  1.00  6.41           C
ANISOU    1  CB  SER A   1      638    591   1207   -162    148    106       C
ATOM      2  OG  SER A   1       9.084  17.462  12.244  1.00  7.60           O
ANISOU    2  OG  SER A   1      680    955   1251   -244   -164    339       O
ATOM      3  HB2 SER A   1      10.836  16.758  11.793  1.00  6.41           H
ATOM      4  HB3 SER A   1      10.100  15.944  12.901  1.00  6.41           H
ATOM      5  HG  SER A   1       8.783  17.124  11.536  1.00  7.60           H
ATOM      6  N  ASER A   1      11.457  18.906  13.130  0.40  5.10           N
ANISOU    6  N  ASER A   1      500    632    808   -107     58    104       N
ATOM      7  CA ASER A   1      11.030  17.638  13.665  0.40  5.95           C
ANISOU    7  CA ASER A   1      582    671   1007   -168    135    158       C
ATOM      8  C  ASER A   1      12.200  16.876  14.296  0.40  5.79           C
ANISOU    8  C  ASER A   1      661    588    952   -116    206    135       C
ATOM      9  O  ASER A   1      12.094  16.315  15.382  0.40  6.83           O
ANISOU    9  O  ASER A   1      874    705   1015    -39    314    258       O
ATOM     10  HA ASER A   1      10.350  17.812  14.335  0.40  5.95           H
ATOM     11  N  BSER A   1      11.557  18.906  13.130  0.60  5.10           N
ANISOU   11  N  BSER A   1      500    632    808   -107     58    104       N
ATOM     12  CA BSER A   1      11.130  17.638  13.665  0.60  5.95           C
ANISOU   12  CA BSER A   1      582    671   1007   -168    135    158       C
ATOM     13  C  BSER A   1      12.300  16.876  14.296  0.60  5.79           C
ANISOU   13  C  BSER A   1      661    588    952   -116    206    135       C
ATOM     14  O  BSER A   1      12.194  16.315  15.382  0.60  6.83           O
ANISOU   14  O  BSER A   1      874    705   1015    -39    314    258       O
ATOM     15  HA BSER A   1      10.450  17.812  14.335  0.60  5.95           H
ATOM     16  CA  SER A   2       7.255  21.911  13.410  1.00  5.95           C
ANISOU   16  CA  SER A   2      582    671   1007   -168    135    158       C
ATOM     17  OG  SER A   2       5.325  21.752  11.965  1.00  7.60           O
ANISOU   17  OG  SER A   2      680    955   1251   -244   -164    339       O
ATOM     18  HA  SER A   2       6.575  22.085  14.080  1.00  5.95           H
ATOM     19  HG  SER A   2       5.027  21.411  11.258  1.00  7.60           H
ATOM     20  N  ASER A   2       7.582  23.179  12.875  0.40  5.10           N
ANISOU   20  N  ASER A   2      500    632    808   -107     58    104       N
ATOM     21  C  ASER A   2       8.325  21.149  14.041  0.40  5.79           C
ANISOU   21  C  ASER A   2      661    588    952   -116    206    135       C
ATOM     22  O  ASER A   2       8.219  20.588  15.127  0.40  6.83           O
ANISOU   22  O  ASER A   2      874    705   1015    -39    314    258       O
ATOM     23  CB ASER A   2       6.403  21.116  12.327  0.40  6.41           C
ANISOU   23  CB ASER A   2      638    591   1207   -162    148    106       C
ATOM     24  HB2ASER A   2       6.961  21.031  11.538  0.40  6.41           H
ATOM     25  HB3ASER A   2       6.225  20.217  12.646  0.40  6.41           H
ATOM     26  N  BSER A   2       7.682  23.179  12.875  0.60  5.10           N
ANISOU   26  N  BSER A   2      500    632    808   -107     58    104       N
ATOM     27  C  BSER A   2       8.425  21.149  14.041  0.60  5.79           C
ANISOU   27  C  BSER A   2      661    588    952   -116    206    135       C
ATOM     28  O  BSER A   2       8.319  20.588  15.127  0.60  6.83           O
ANISOU   28  O  BSER A   2      874    705   1015    -39    314    258       O
ATOM     29  CB BSER A   2       6.503  21.116  12.327  0.60  6.41           C
ANISOU   29  CB BSER A   2      638    591   1207   -162    148    106       C
ATOM     30  HB2BSER A   2       7.070  21.011  11.547  0.60  6.41           H
ATOM     31  HB3BSER A   2       6.305  20.224  12.655  0.60  6.41           H
ATOM     32  OG  SER A   3       8.348  24.487  15.352  1.00  7.60           O
ANISOU   32  OG  SER A   3      680    955   1251   -244   -164    339       O
ATOM     33  HG  SER A   3       8.047  24.149  14.644  1.00  7.60           H
ATOM     34  N  ASER A   3      10.721  25.931  16.238  0.40  5.10           N
ANISOU   34  N  ASER A   3      500    632    808   -107     58    104       N
ATOM     35  CA ASER A   3      10.294  24.663  16.773  0.40  5.95           C
ANISOU   35  CA ASER A   3      582    671   1007   -168    135    158       C
ATOM     36  C  ASER A   3      11.464  23.901  17.404  0.40  5.79           C
ANISOU   36  C  ASER A   3      661    588    952   -116    206    135       C
ATOM     37  O  ASER A   3      11.358  23.340  18.490  0.40  6.83           O
ANISOU   37  O  ASER A   3      874    705   1015    -39    314    258       O
ATOM     38  CB ASER A   3       9.542  23.868  15.690  0.40  6.41           C
ANISOU   38  CB ASER A   3      638    591   1207   -162    148    106       C
ATOM     39  HA ASER A   3       9.614  24.837  17.443  0.40  5.95           H
ATOM     40  HB2ASER A   3      10.100  23.783  14.901  0.40  6.41           H
ATOM     41  HB3ASER A   3       9.364  22.969  16.009  0.40  6.41           H
ATOM     42  N  BSER A   3      10.821  25.931  16.238  0.60  5.10           N
ANISOU   42  N  BSER A   3      500    632    808   -107     58    104       N
ATOM     43  CA BSER A   3      10.394  24.663  16.773  0.60  5.95           C
ANISOU   43  CA BSER A   3      582    671   1007   -168    135    158       C
ATOM     44  C  BSER A   3      11.564  23.901  17.404  0.60  5.79           C
ANISOU   44  C  BSER A   3      661    588    952   -116    206    135       C
ATOM     45  O  BSER A   3      11.458  23.340  18.490  0.60  6.83           O
ANISOU   45  O  BSER A   3      874    705   1015    -39    314    258       O
ATOM     46  CB BSER A   3       9.642  23.868  15.690  0.60  6.41           C
ANISOU   46  CB BSER A   3      638    591   1207   -162    148    106       C
ATOM     47  HA BSER A   3       9.714  24.837  17.443  0.60  5.95           H
ATOM     48  HB2BSER A   3      10.209  23.763  14.910  0.60  6.41           H
ATOM     49  HB3BSER A   3       9.444  22.976  16.018  0.60  6.41           H
ATOM     50  CA  TYR A   4      22.290  10.873   7.693  1.00 13.00           C
ANISOU   50  CA  TYR A   4     2514   1244   1180    -28    116    139       C
ATOM     51  C   TYR A   4      22.074   9.526   7.008  1.00 11.81           C
ANISOU   51  C   TYR A   4     2332   1026   1130   -105     88    190       C
ATOM     52  O   TYR A   4      21.538   8.599   7.605  1.00 12.59           O
ANISOU   52  O   TYR A   4     2489   1154   1140    -72    111    219       O
ATOM     53  N   VAL A   5      22.517   9.437   5.761  1.00 10.80           N
ANISOU   53  N   VAL A   5     2053    884   1167    -29     80    167       N
ATOM     54  H   VAL A   5      22.832  10.103   5.318  1.00 12.96           H
ATOM     55  CA AVAL A   5      22.516   8.185   5.022  0.40 10.37           C
ANISOU   55  CA AVAL A   5     1783    902   1255     67   -148    -14       C
ATOM     56  CA BVAL A   5      22.616   8.185   5.022  0.60 10.37           C
ANISOU   56  CA BVAL A   5     1783    902   1255     67   -148    -14       C
ATOM     57  N   SER A   6      17.758  20.026  12.546  1.00  4.20           N
ANISOU   57  N   SER A   6      578    469    547    -59     10   -119       N
ATOM     58  CA  SER A   6      17.692  20.410  11.150  1.00  4.34           C
ANISOU   58  CA  SER A   6      576    472    600     -2      0      0       C
ATOM     59  C   SER A   6      18.920  19.904  10.378  1.00  3.90           C
ANISOU   59  C   SER A   6      554    387    543    -50    -15    -34       C
ATOM     60  O   SER A   6      18.798  19.501   9.222  1.00  4.24           O
ANISOU   60  O   SER A   6      634    443    534    -77    -67    -59       O
ATOM     62  HA  SER A   6      16.905  19.937  10.839  1.00  4.34           H
ATOM     63  CB ASER A   6      17.651  21.951  11.144  0.60  5.87           C
ANISOU   63  CB ASER A   6      843    568    821    229    -39    -32       C
ATOM     64  OG ASER A   6      17.494  22.382   9.839  0.60  6.95           O
ANISOU   64  OG ASER A   6     1107    702    833    257    200    132       O
ATOM     65  HB2ASER A   6      16.919  22.271  11.695  0.60  4.80           H
ATOM     66  HB3ASER A   6      18.469  22.311  11.522  0.60  4.80           H
ATOM     67  HG ASER A   6      16.967  23.036   9.821  0.60  4.95           H
ATOM     68  CB BSER A   6      17.234  21.814  10.776  0.20  4.89           C
ANISOU   68  CB BSER A   6      429    509    922     98    134    102       C
ATOM     69  OG BSER A   6      18.212  22.672  11.327  0.20  5.81           O
ANISOU   69  OG BSER A   6      757    477    971     72     15   -142       O
ATOM     70  HB2BSER A   6      17.175  21.920   9.814  0.20  4.80           H
ATOM     71  HB3BSER A   6      16.354  22.005  11.137  0.20  4.80           H
ATOM     72  HG BSER A   6      18.020  23.469  11.146  0.20  4.95           H
ATOM     73  CB CSER A   6      17.677  21.930  10.974  0.20  4.80           C
ANISOU   73  CB CSER A   6      432    322   1070   -132    173      4       C
ATOM     74  OG CSER A   6      16.492  22.377  11.598  0.20  4.95           O
ANISOU   74  OG CSER A   6      511    563    806    114    -22    -75       O
ATOM     75  HB2CSER A   6      18.459  22.335  11.381  0.20  4.80           H
ATOM     76  HB3CSER A   6      17.688  22.172  10.035  0.20  4.80           H
ATOM     77  HG CSER A   6      16.620  22.446  12.425  0.20  4.95           H
ATOM     78  N   VAL A   7      17.941  13.419  15.481  1.00 10.80           N
ANISOU   78  N   VAL A   7     2053    884   1167    -29     80    167       N
ATOM     79  C   VAL A   7      19.247  12.139  13.940  1.00  9.56           C
ANISOU   79  C   VAL A   7     1874    685   1073      3   -185     85       C
ATOM     80  O   VAL A   7      19.671  13.155  13.422  1.00 10.25           O
ANISOU   80  O   VAL A   7     2013    634   1247     25      7     95       O
ATOM     81  HA  VAL A   7      17.920  11.415  15.343  0.44 12.11           H
ATOM     82  CA AVAL A   7      17.940  12.167  14.742  0.44 10.37           C
ANISOU   82  CA AVAL A   7     1783    902   1255     67   -148    -14       C
ATOM     83  CB AVAL A   7      16.693  12.086  13.853  0.44 11.26           C
ANISOU   83  CB AVAL A   7     1571    963   1743     96   -176     62       C
ATOM     84  CG1AVAL A   7      16.690  13.235  12.923  0.44 11.06           C
ANISOU   84  CG1AVAL A   7     1704   1072   1425    125   -105     37       C
ATOM     85  CG2AVAL A   7      16.608  10.757  13.097  0.44 10.02           C
ANISOU   85  CG2AVAL A   7     1510    782   1516    -92     19    232       C
ATOM     86  HB AVAL A   7      15.905  12.157  14.414  0.44 13.51           H
ATOM     87 HG11AVAL A   7      16.677  14.057  13.438  0.44 13.27           H
ATOM     88 HG12AVAL A   7      15.902  13.184  12.361  0.44 13.27           H
ATOM     89 HG13AVAL A   7      17.490  13.200  12.376  0.44 13.27           H
ATOM     90 HG21AVAL A   7      15.789  10.303  13.349  0.44 12.03           H
ATOM     91 HG22AVAL A   7      17.375  10.211  13.331  0.44 12.03           H
ATOM     92 HG23AVAL A   7      16.608  10.936  12.143  0.44 12.03           H
ATOM     93  CA BVAL A   7      17.935  12.170  14.724  0.56 10.10           C
ANISOU   93  CA BVAL A   7     1914    779   1143    -16      5    185       C
ATOM     94  CB BVAL A   7      16.746  12.075  13.715  0.56 11.08           C
ANISOU   94  CB BVAL A   7     1662   1363   1184     84    -25     58       C
ATOM     95  CG1BVAL A   7      16.699  10.688  13.068  0.56 12.75           C
ANISOU   95  CG1BVAL A   7     1479   1600   1766    255   -607   -238       C
ATOM     96  CG2BVAL A   7      15.398  12.392  14.382  0.56 11.66           C
ANISOU   96  CG2BVAL A   7     1508   1829   1094    214   -319     28       C
ATOM     97  HB BVAL A   7      16.885  12.725  13.009  0.56 13.29           H
ATOM     98 HG11BVAL A   7      17.532  10.533  12.595  0.56 15.30           H
ATOM     99 HG12BVAL A   7      15.954  10.654  12.448  0.56 15.30           H
ATOM    100 HG13BVAL A   7      16.583  10.020  13.762  0.56 15.30           H
ATOM    101 HG21BVAL A   7      14.999  13.157  13.939  0.56 13.99           H
ATOM    102 HG22BVAL A   7      15.550  12.594  15.318  0.56 13.99           H
ATOM    103 HG23BVAL A   7      14.817  11.620  14.299  0.56 13.99           H
ATOM    104  N   SER A   8      14.526  19.060  18.223  1.00  5.10           N
ANISOU  104  N   SER A   8      500    632    808   -107     58    104       N
ATOM    105  CA  SER A   8      14.099  17.792  18.758  1.00  5.95           C
ANISOU  105  CA  SER A   8      582    671   1007   -168    135    158       C
ATOM    106  C   SER A   8      15.269  17.030  19.389  1.00  5.79           C
ANISOU  106  C   SER A   8      661    588    952   -116    206    135       C
ATOM    107  O   SER A   8      15.163  16.469  20.475  1.00  6.83           O
ANISOU  107  O   SER A   8      874    705   1015    -39    314    258       O
ATOM    108  HA  SER A   8      13.419  17.966  19.428  1.00  5.95           H
ATOM    109  CB ASER A   8      13.570  16.905  17.583  0.52 10.56           C
ANISOU  109  CB ASER A   8     1137    993   1883   -635   -478    170       C
ATOM    110  OG ASER A   8      13.140  15.675  18.052  0.52 13.20           O
ANISOU  110  OG ASER A   8     1612   1164   2238   -898    145    -29       O
ATOM    111  HB2ASER A   8      12.840  17.358  17.133  0.52  6.41           H
ATOM    112  HB3ASER A   8      14.272  16.778  16.926  0.52  6.41           H
ATOM    113  HG ASER A   8      12.347  15.739  18.322  0.52  7.60           H
ATOM    114  CB BSER A   8      13.347  16.997  17.675  0.47  6.41           C
ANISOU  114  CB BSER A   8      638    591   1207   -162    148    106       C
ATOM    115  OG BSER A   8      12.153  17.616  17.337  0.47  7.60           O
ANISOU  115  OG BSER A   8      680    955   1251   -244   -164    339       O
ATOM    116  HB2BSER A   8      13.905  16.912  16.886  0.47  6.41           H
ATOM    117  HB3BSER A   8      13.169  16.098  17.994  0.47  6.41           H
ATOM    118  HG BSER A   8      11.852  17.278  16.629  0.47  7.60           H
ATOM    119  C  AASN A   9      19.204   6.973  12.924  0.46 10.01           C
ANISOU  119  C  AASN A   9     1311   1084   1408    -15   -324    -40       C
ATOM    120  C  BASN A   9      19.123   6.961  12.995  0.54  9.67           C
ANISOU  120  C  BASN A   9     1289   1143   1241   -122   -106    160       C
ATOM    121  N   SER A  10      19.921   6.127  13.648  1.00  9.28           N
ANISOU  121  N   SER A  10     1274    902   1350    -21    -87    -74       N
ATOM    122  H   SER A  10      20.776   6.182  13.716  1.00 11.14           H
ATOM    123  CA ASER A  10      19.327   5.020  14.363  0.63  9.94           C
ANISOU  123  CA ASER A  10     1406    774   1598     87   -335   -196       C
ATOM    124  CA BSER A  10      19.344   5.000  14.369  0.37  9.47           C
ANISOU  124  CA BSER A  10     1420    764   1413     89    320   -137       C
ATOM    125  N   TYR A  11       9.640   8.356   6.598  1.00 15.00           N
ATOM    126  C   TYR A  11      11.549   9.847   6.197  1.00 15.00           C
ATOM    127  O   TYR A  11      11.823   9.939   5.000  1.00 15.00           O
ATOM    128  HA  TYR A  11       9.553  10.281   6.071  1.00 15.00           H
ATOM    129  CA ATYR A  11      10.100   9.738   6.660  0.50 15.00           C
ATOM    130  CB ATYR A  11       9.960  10.286   8.081  0.50 15.00           C
ATOM    131  CG ATYR A  11       8.637   9.956   8.735  0.50 15.00           C
ATOM    132  CD1ATYR A  11       7.436  10.254   8.105  0.50 15.00           C
ATOM    133  CD2ATYR A  11       8.589   9.348   9.982  0.50 15.00           C
ATOM    134  CE1ATYR A  11       6.224   9.955   8.699  0.50 15.00           C
ATOM    135  CE2ATYR A  11       7.383   9.044  10.584  0.50 15.00           C
ATOM    136  CZ ATYR A  11       6.203   9.350   9.938  0.50 15.00           C
ATOM    137  OH ATYR A  11       5.000   9.050  10.533  0.50 15.00           O
ATOM    138  HB2ATYR A  11      10.664   9.911   8.632  0.50 15.00           H
ATOM    139  HB3ATYR A  11      10.044  11.252   8.053  0.50 15.00           H
ATOM    140  HD1ATYR A  11       7.446  10.662   7.270  0.50 15.00           H
ATOM    141  HD2ATYR A  11       9.383   9.140  10.420  0.50 15.00           H
ATOM    142  HE1ATYR A  11       5.427  10.160   8.265  0.50 15.00           H
ATOM    143  HE2ATYR A  11       7.366   8.637  11.419  0.50 15.00           H
ATOM    144  HH ATYR A  11       5.131   8.687  11.279  0.50 15.00           H
ATOM    145  CA BTYR A  11      10.100   9.738   6.660  0.50 15.00           C
ATOM    146  CB BTYR A  11       9.960  10.286   8.081  0.50 15.00           C
ATOM    147  CG BTYR A  11      10.474  11.698   8.244  0.50 15.00           C
ATOM    148  CD1BTYR A  11       9.764  12.780   7.741  0.50 15.00           C
ATOM    149  CD2BTYR A  11      11.671  11.950   8.902  0.50 15.00           C
ATOM    150  CE1BTYR A  11      10.230  14.072   7.888  0.50 15.00           C
ATOM    151  CE2BTYR A  11      12.146  13.239   9.053  0.50 15.00           C
ATOM    152  CZ BTYR A  11      11.422  14.296   8.545  0.50 15.00           C
ATOM    153  OH BTYR A  11      11.890  15.582   8.693  0.50 15.00           O
ATOM    154  HB2BTYR A  11       9.021  10.283   8.326  0.50 15.00           H
ATOM    155  HB3BTYR A  11      10.459   9.717   8.687  0.50 15.00           H
ATOM    156  HD1BTYR A  11       8.960  12.632   7.297  0.50 15.00           H
ATOM    157  HD2BTYR A  11      12.162  11.239   9.246  0.50 15.00           H
ATOM    158  HE1BTYR A  11       9.744  14.787   7.545  0.50 15.00           H
ATOM    159  HE2BTYR A  11      12.949  13.393   9.496  0.50 15.00           H
ATOM    160  HH BTYR A  11      12.620  15.575   9.109  0.50 15.00           H
TER
END
"""

type_list_known1 = ['2tetra', '2tetra', 'alg1b', '3neigbs', '3neigbs',
  '3neigbs', 'alg1b', '2tetra', '2tetra', '2tetra', '2tetra', 'alg1b',
  '3neigbs', '2tetra', '2tetra', '3neigbs', '2tetra', '2tetra', 'flat_2neigbs',
  '3neigbs', '2tetra', '2tetra', 'alg1b', '2tetra', '2tetra', 'alg1b',
  '2tetra', '2tetra', 'alg1b', '3neigbs', '3neigbs', 'prop', 'prop', 'prop',
  'prop', 'prop', 'prop', '3neigbs', 'prop', 'prop', 'prop', 'prop', 'prop',
  'prop', '3neigbs', '2tetra', '2tetra', 'alg1b', '2tetra', '2tetra', 'alg1b',
  'flat_2neigbs', '3neigbs', '2tetra', '2tetra', 'flat_2neigbs', 'flat_2neigbs',
  'flat_2neigbs', 'flat_2neigbs', 'alg1b', '2tetra', '2tetra', 'flat_2neigbs',
  'flat_2neigbs', 'flat_2neigbs', 'flat_2neigbs', 'alg1b']

pdb_str2 = """
CRYST1   13.142   13.841   12.524  90.00  90.00  90.00 P 1
SCALE1      0.076092  0.000000  0.000000        0.00000
SCALE2      0.000000  0.072249  0.000000        0.00000
SCALE3      0.000000  0.000000  0.079847        0.00000
ATOM    285  CB  ALA A  25      21.105  -3.928  25.422  1.00 11.78           C
ATOM    286  N  AALA A  25      18.694  -4.221  25.276  0.50  9.87           N
ATOM    287  CA AALA A  25      19.812  -3.247  25.214  0.50 10.86           C
ATOM    288  C  AALA A  25      19.902  -2.496  23.898  0.50  9.98           C
ATOM    289  O  AALA A  25      20.720  -1.565  23.768  0.50  7.99           O
ATOM    291  HA AALA A  25      19.625  -2.606  25.918  0.50 11.03           H
ATOM    292  HB1AALA A  25      21.824  -3.278  25.378  0.50 11.78           H
ATOM    293  HB2AALA A  25      21.109  -4.356  26.292  0.50 11.78           H
ATOM    294  HB3AALA A  25      21.233  -4.598  24.732  0.50 11.78           H
ATOM    295  N  BALA A  25      18.682  -4.224  25.328  0.50 10.05           N
ATOM    296  CA BALA A  25      19.756  -3.242  25.276  0.50 11.02           C
ATOM    297  C  BALA A  25      19.678  -2.389  23.992  0.50  9.97           C
ATOM    298  O  BALA A  25      20.024  -1.203  24.029  0.50 11.29           O
ATOM    300  HA BALA A  25      19.650  -2.631  26.022  0.50 11.03           H
ATOM    301  HB1BALA A  25      21.811  -3.264  25.386  0.50 11.78           H
ATOM    302  HB2BALA A  25      21.142  -4.393  26.273  0.50 11.78           H
ATOM    303  HB3BALA A  25      21.224  -4.566  24.701  0.50 11.78           H
"""

type_list_known2 = ['3neigbs', 'prop', 'prop', 'prop',
  '3neigbs', 'prop', 'prop', 'prop']


# This fragment is from PDB model 1qjh
# A residue (Arg 71 in 1qjh) is close to its symmetry mate.
# The geo file lists two asu bond restraints for atom HH22 (there should be
# only one, IMHO!).
# The fact that HH22 has no simple bond proxy, but angle proxies, caused
# an error (asu covalent bond restraints are currently ignored).

pdb_str3 = """
CRYST1   49.701   49.701   73.255  90.00  90.00  90.00 P 42 21 2     8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.020120  0.000000  0.000000        0.00000
SCALE2      0.000000  0.020120  0.000000        0.00000
SCALE3      0.000000  0.000000  0.013651        0.00000
ATOM    575  N   ARG A  71       4.405  19.140  21.765  1.00 56.18           N
ATOM    576  CA  ARG A  71       4.384  20.597  21.643  1.00 50.81           C
ATOM    577  C   ARG A  71       5.782  21.217  21.647  1.00 45.46           C
ATOM    578  O   ARG A  71       5.933  22.440  21.582  1.00 46.55           O
ATOM    579  CB  ARG A  71       3.558  21.214  22.773  1.00 51.01           C
ATOM    580  CG  ARG A  71       2.166  20.635  22.956  1.00 52.13           C
ATOM    581  CD  ARG A  71       1.358  20.625  21.661  1.00 60.92           C
ATOM    582  NE  ARG A  71       1.134  21.965  21.122  1.00 63.41           N
ATOM    583  CZ  ARG A  71       0.470  22.927  21.755  1.00 62.11           C
ATOM    584  NH1 ARG A  71      -0.038  22.699  22.955  1.00 60.26           N
ATOM    585  NH2 ARG A  71       0.310  24.119  21.188  1.00 62.26           N
ATOM      0  HA  ARG A  71       3.979  20.793  20.783  1.00 50.81           H
ATOM      0  HB2 ARG A  71       4.047  21.109  23.604  1.00 51.01           H
ATOM      0  HB3 ARG A  71       3.476  22.167  22.610  1.00 51.01           H
ATOM      0  HG2 ARG A  71       2.238  19.729  23.295  1.00 52.13           H
ATOM      0  HG3 ARG A  71       1.690  21.151  23.626  1.00 52.13           H
ATOM      0  HD2 ARG A  71       1.822  20.088  20.999  1.00 60.92           H
ATOM      0  HD3 ARG A  71       0.502  20.199  21.822  1.00 60.92           H
ATOM      0  HE  ARG A  71       1.453  22.143  20.343  1.00 63.41           H
ATOM      0 HH11 ARG A  71       0.062  21.929  23.325  1.00 60.26           H
ATOM      0 HH12 ARG A  71      -0.468  23.321  23.365  1.00 60.26           H
ATOM      0 HH21 ARG A  71       0.637  24.271  20.407  1.00 62.26           H
ATOM      0 HH22 ARG A  71      -0.120  24.738  21.602  1.00 62.26           H
"""

type_list_known3 = ['3neigbs', '2tetra', '2tetra', '2tetra', '2tetra',
 '2tetra', '2tetra', 'flat_2neigbs', 'alg1a', 'alg1a']

def exercise3(pdb_str, type_list_known):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  params = mmtbx.model.manager.get_default_pdb_interpretation_scope().extract()
  params.pdb_interpretation.allow_polymer_cross_special_position=True
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(pdb_interpretation_params=params, make_restraints=True)
  pdb_hierarchy = model.get_hierarchy()
  sites_cart = model.get_sites_cart()
  atoms = pdb_hierarchy.atoms()

  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()

  h_para = riding_h_manager.h_parameterization

  diagnostics = riding_h_manager.diagnostics(
    sites_cart = sites_cart,
    threshold  = 0.05)
  h_distances   = diagnostics.h_distances
  type_list     = diagnostics.type_list

  number_h = model.get_hd_selection().count(True)
  number_h_para = len(h_para) - h_para.count(None)

  assert (number_h_para == number_h-2), 'Not all H atoms are parameterized'

  for ih in h_distances:
    labels = atoms[ih].fetch_labels()
    assert (h_distances[ih] < 0.2), \
      'distance too large: %s  atom: %s (%s) residue: %s ' \
      % (h_para[ih].htype, atoms[ih].name, ih, labels.resseq.strip())

  for type1, type2 in zip(type_list, type_list_known):
    assert (type1 == type2)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise(pdb_str = pdb_str1, type_list_known = type_list_known1)
  exercise(pdb_str = pdb_str2, type_list_known = type_list_known2)
  exercise3(pdb_str = pdb_str3, type_list_known = type_list_known3)

  print("OK. Time: %8.3f"%(time.time()-t0))
