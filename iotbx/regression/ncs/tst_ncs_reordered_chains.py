from __future__ import absolute_import, division, print_function
import time
from tst_ncs_user_selections import get_ncs_groups

#
# Example of shortcoming of NCS verification procedure. The error is due to
# exotic order of chains inside NCS groups.
#
# The failure shows that not something is broken, but is an example of the corner
# case with which current code cannot deal with.
#


pdb_str="""
CRYST1  106.713   90.134   63.871  90.00  90.00  90.00 P 1
SCALE1      0.009371  0.000000  0.000000        0.00000
SCALE2      0.000000  0.011095  0.000000        0.00000
SCALE3      0.000000  0.000000  0.015657        0.00000
ATOM      1  N   SER A   3      13.654  19.632   7.555  1.00 34.27           N
ATOM      2  CA  SER A   3      12.388  18.913   7.485  1.00 34.27           C
ATOM      3  C   SER A   3      12.618  17.449   7.145  1.00 34.27           C
ATOM      4  O   SER A   3      13.438  17.137   6.282  1.00 34.27           O
ATOM      5  CB  SER A   3      11.461  19.548   6.451  1.00 34.27           C
ATOM      6  OG  SER A   3      11.917  19.290   5.135  1.00 34.27           O
ATOM      7  N   PRO A   4      11.898  16.555   7.831  1.00 33.04           N
ATOM      8  CA  PRO A   4      12.068  15.118   7.555  1.00 33.04           C
ATOM      9  C   PRO A   4      11.762  14.729   6.120  1.00 33.04           C
ATOM     10  O   PRO A   4      12.428  13.841   5.571  1.00 33.04           O
ATOM     11  CB  PRO A   4      11.101  14.461   8.549  1.00 33.04           C
ATOM     12  CG  PRO A   4      10.135  15.537   8.917  1.00 33.04           C
ATOM     13  CD  PRO A   4      10.919  16.800   8.900  1.00 33.04           C
TER
ATOM     14  N   SER B   3      49.210  85.134  17.642  1.00 20.33           N
ATOM     15  CA  SER B   3      48.236  84.429  18.467  1.00 20.33           C
ATOM     16  C   SER B   3      48.154  82.969  18.058  1.00 20.33           C
ATOM     17  O   SER B   3      48.185  82.668  16.869  1.00 20.33           O
ATOM     18  CB  SER B   3      46.860  85.086  18.357  1.00 20.33           C
ATOM     19  OG  SER B   3      45.877  84.345  19.060  1.00 20.33           O
ATOM     20  N   PRO B   4      48.033  82.073  19.053  1.00 17.85           N
ATOM     21  CA  PRO B   4      48.012  80.628  18.778  1.00 17.85           C
ATOM     22  C   PRO B   4      47.212  80.215  17.552  1.00 17.85           C
ATOM     23  O   PRO B   4      47.712  79.456  16.715  1.00 17.85           O
ATOM     24  CB  PRO B   4      47.390  80.040  20.052  1.00 17.85           C
ATOM     25  CG  PRO B   4      47.304  81.176  21.047  1.00 17.85           C
ATOM     26  CD  PRO B   4      48.041  82.344  20.496  1.00 17.85           C
TER
ATOM     27  N   LEU C  27       5.385  30.101  36.788  1.00 41.94           N
ATOM     28  CA  LEU C  27       5.992  30.923  37.828  1.00 41.94           C
ATOM     29  C   LEU C  27       6.816  32.057  37.227  1.00 41.94           C
ATOM     30  O   LEU C  27       6.274  33.081  36.813  1.00 41.94           O
ATOM     31  CB  LEU C  27       6.870  30.069  38.742  1.00 41.94           C
ATOM     32  CG  LEU C  27       6.182  28.907  39.459  1.00 41.94           C
ATOM     33  CD1 LEU C  27       7.173  28.162  40.339  1.00 41.94           C
ATOM     34  CD2 LEU C  27       5.000  29.405  40.276  1.00 41.94           C
ATOM     35  N   LEU C  28       8.132  31.860  37.178  1.00 43.81           N
ATOM     36  CA  LEU C  28       9.055  32.874  36.690  1.00 43.81           C
ATOM     37  C   LEU C  28      10.305  32.191  36.157  1.00 43.81           C
ATOM     38  O   LEU C  28      10.553  31.011  36.416  1.00 43.81           O
ATOM     39  CB  LEU C  28       9.439  33.868  37.792  1.00 43.81           C
ATOM     40  CG  LEU C  28       8.481  35.000  38.154  1.00 43.81           C
ATOM     41  CD1 LEU C  28       8.934  35.675  39.438  1.00 43.81           C
ATOM     42  CD2 LEU C  28       8.405  36.006  37.021  1.00 43.81           C
TER
ATOM     43  O5'  DA K   1      24.122  38.760  12.604  1.00 49.90           O
ATOM     44  C5'  DA K   1      24.158  37.777  11.573  1.00 49.90           C
ATOM     45  C4'  DA K   1      24.551  38.399  10.243  1.00 49.90           C
ATOM     46  O4'  DA K   1      23.742  39.580  10.004  1.00 49.90           O
ATOM     47  C3'  DA K   1      26.001  38.862  10.146  1.00 49.90           C
ATOM     48  O3'  DA K   1      26.480  38.693   8.827  1.00 49.90           O
ATOM     49  C2'  DA K   1      25.908  40.332  10.513  1.00 49.90           C
ATOM     50  C1'  DA K   1      24.571  40.717   9.897  1.00 49.90           C
ATOM     51  N9   DA K   1      23.923  41.829  10.585  1.00 49.90           N
ATOM     52  C8   DA K   1      23.464  41.847  11.872  1.00 49.90           C
ATOM     53  N7   DA K   1      22.930  42.991  12.226  1.00 49.90           N
ATOM     54  C5   DA K   1      23.050  43.778  11.094  1.00 49.90           C
ATOM     55  C6   DA K   1      22.670  45.105  10.816  1.00 49.90           C
ATOM     56  N6   DA K   1      22.070  45.899  11.707  1.00 49.90           N
ATOM     57  N1   DA K   1      22.933  45.584   9.582  1.00 49.90           N
ATOM     58  C2   DA K   1      23.533  44.785   8.693  1.00 49.90           C
ATOM     59  N3   DA K   1      23.936  43.525   8.839  1.00 49.90           N
ATOM     60  C4   DA K   1      23.656  43.076  10.072  1.00 49.90           C
TER
ATOM     61  O5'  DG L   2      76.950  13.593  57.985  1.00 91.25           O
ATOM     62  C5'  DG L   2      76.703  14.675  58.871  1.00 91.25           C
ATOM     63  C4'  DG L   2      76.734  16.000  58.130  1.00 91.25           C
ATOM     64  O4'  DG L   2      77.717  15.944  57.064  1.00 91.25           O
ATOM     65  C3'  DG L   2      75.435  16.402  57.445  1.00 91.25           C
ATOM     66  O3'  DG L   2      75.390  17.820  57.381  1.00 91.25           O
ATOM     67  C2'  DG L   2      75.629  15.800  56.050  1.00 91.25           C
ATOM     68  C1'  DG L   2      77.071  16.215  55.835  1.00 91.25           C
ATOM     69  N9   DG L   2      77.814  15.539  54.772  1.00 91.25           N
ATOM     70  C8   DG L   2      77.703  14.240  54.328  1.00 91.25           C
ATOM     71  N7   DG L   2      78.562  13.946  53.377  1.00 91.25           N
ATOM     72  C5   DG L   2      79.288  15.124  53.202  1.00 91.25           C
ATOM     73  C6   DG L   2      80.359  15.443  52.317  1.00 91.25           C
ATOM     74  O6   DG L   2      80.916  14.723  51.477  1.00 91.25           O
ATOM     75  N1   DG L   2      80.789  16.763  52.490  1.00 91.25           N
ATOM     76  C2   DG L   2      80.253  17.654  53.392  1.00 91.25           C
ATOM     77  N2   DG L   2      80.791  18.886  53.421  1.00 91.25           N
ATOM     78  N3   DG L   2      79.264  17.365  54.210  1.00 91.25           N
ATOM     79  C4   DG L   2      78.835  16.097  54.062  1.00 91.25           C
TER
ATOM     80  N   SER D   3      93.227  70.323   7.415  1.00 34.27           N
ATOM     81  CA  SER D   3      94.496  71.037   7.347  1.00 34.27           C
ATOM     82  C   SER D   3      94.273  72.500   6.999  1.00 34.27           C
ATOM     83  O   SER D   3      93.457  72.811   6.131  1.00 34.27           O
ATOM     84  CB  SER D   3      95.425  70.394   6.319  1.00 34.27           C
ATOM     85  OG  SER D   3      94.975  70.647   5.000  1.00 34.27           O
ATOM     86  N   PRO D   4      94.993  73.395   7.684  1.00 33.04           N
ATOM     87  CA  PRO D   4      94.830  74.831   7.400  1.00 33.04           C
ATOM     88  C   PRO D   4      95.143  75.212   5.965  1.00 33.04           C
ATOM     89  O   PRO D   4      94.482  76.100   5.409  1.00 33.04           O
ATOM     90  CB  PRO D   4      95.795  75.489   8.395  1.00 33.04           C
ATOM     91  CG  PRO D   4      96.756  74.411   8.772  1.00 33.04           C
ATOM     92  CD  PRO D   4      95.967  73.151   8.758  1.00 33.04           C
TER
ATOM     93  N   SER E   3      57.388   5.000  17.660  1.00 20.33           N
ATOM     94  CA  SER E   3      58.361   5.705  18.486  1.00 20.33           C
ATOM     95  C   SER E   3      58.450   7.163  18.071  1.00 20.33           C
ATOM     96  O   SER E   3      58.425   7.459  16.880  1.00 20.33           O
ATOM     97  CB  SER E   3      59.735   5.043  18.385  1.00 20.33           C
ATOM     98  OG  SER E   3      60.718   5.783  19.088  1.00 20.33           O
ATOM     99  N   PRO E   4      58.571   8.063  19.062  1.00 17.85           N
ATOM    100  CA  PRO E   4      58.598   9.507  18.780  1.00 17.85           C
ATOM    101  C   PRO E   4      59.405   9.911  17.556  1.00 17.85           C
ATOM    102  O   PRO E   4      58.911  10.668  16.713  1.00 17.85           O
ATOM    103  CB  PRO E   4      59.217  10.098  20.054  1.00 17.85           C
ATOM    104  CG  PRO E   4      59.295   8.967  21.055  1.00 17.85           C
ATOM    105  CD  PRO E   4      58.556   7.799  20.506  1.00 17.85           C
TER
ATOM    106  N   LEU F  27     101.339  59.959  36.729  1.00 41.94           N
ATOM    107  CA  LEU F  27     100.725  59.144  37.770  1.00 41.94           C
ATOM    108  C   LEU F  27      99.899  58.010  37.171  1.00 41.94           C
ATOM    109  O   LEU F  27     100.439  56.983  36.764  1.00 41.94           O
ATOM    110  CB  LEU F  27      99.847  60.006  38.677  1.00 41.94           C
ATOM    111  CG  LEU F  27     100.536  61.168  39.391  1.00 41.94           C
ATOM    112  CD1 LEU F  27      99.544  61.921  40.264  1.00 41.94           C
ATOM    113  CD2 LEU F  27     101.713  60.670  40.215  1.00 41.94           C
ATOM    114  N   LEU F  28      98.584  58.212  37.116  1.00 43.81           N
ATOM    115  CA  LEU F  28      97.659  57.199  36.629  1.00 43.81           C
ATOM    116  C   LEU F  28      96.414  57.884  36.088  1.00 43.81           C
ATOM    117  O   LEU F  28      96.169  59.066  36.340  1.00 43.81           O
ATOM    118  CB  LEU F  28      97.267  56.212  37.734  1.00 43.81           C
ATOM    119  CG  LEU F  28      98.220  55.078  38.105  1.00 43.81           C
ATOM    120  CD1 LEU F  28      97.759  54.411  39.390  1.00 43.81           C
ATOM    121  CD2 LEU F  28      98.296  54.066  36.977  1.00 43.81           C
TER
ATOM    122  O5'  DA I   1      82.668  51.257  12.510  1.00 49.90           O
ATOM    123  C5'  DA I   1      82.640  52.236  11.474  1.00 49.90           C
ATOM    124  C4'  DA I   1      82.250  51.609  10.145  1.00 49.90           C
ATOM    125  O4'  DA I   1      83.055  50.424   9.915  1.00 49.90           O
ATOM    126  C3'  DA I   1      80.798  51.151  10.045  1.00 49.90           C
ATOM    127  O3'  DA I   1      80.325  51.316   8.723  1.00 49.90           O
ATOM    128  C2'  DA I   1      80.885  49.682  10.419  1.00 49.90           C
ATOM    129  C1'  DA I   1      82.223  49.289   9.810  1.00 49.90           C
ATOM    130  N9   DA I   1      82.864  48.178  10.506  1.00 49.90           N
ATOM    131  C8   DA I   1      83.317  48.165  11.795  1.00 49.90           C
ATOM    132  N7   DA I   1      83.846  47.020  12.156  1.00 49.90           N
ATOM    133  C5   DA I   1      83.727  46.228  11.027  1.00 49.90           C
ATOM    134  C6   DA I   1      84.104  44.899  10.757  1.00 49.90           C
ATOM    135  N6   DA I   1      84.697  44.107  11.654  1.00 49.90           N
ATOM    136  N1   DA I   1      83.844  44.415   9.524  1.00 49.90           N
ATOM    137  C2   DA I   1      83.250  45.212   8.629  1.00 49.90           C
ATOM    138  N3   DA I   1      82.851  46.474   8.768  1.00 49.90           N
ATOM    139  C4   DA I   1      83.128  46.928  10.000  1.00 49.90           C
TER
ATOM    140  O5'  DG J   2      29.751  76.829  57.560  1.00 91.25           O
ATOM    141  C5'  DG J   2      29.991  75.750  58.452  1.00 91.25           C
ATOM    142  C4'  DG J   2      29.958  74.422  57.717  1.00 91.25           C
ATOM    143  O4'  DG J   2      28.979  74.477  56.647  1.00 91.25           O
ATOM    144  C3'  DG J   2      31.258  74.012  57.039  1.00 91.25           C
ATOM    145  O3'  DG J   2      31.298  72.594  56.982  1.00 91.25           O
ATOM    146  C2'  DG J   2      31.072  74.608  55.640  1.00 91.25           C
ATOM    147  C1'  DG J   2      29.629  74.198  55.422  1.00 91.25           C
ATOM    148  N9   DG J   2      28.893  74.871  54.352  1.00 91.25           N
ATOM    149  C8   DG J   2      29.010  76.168  53.903  1.00 91.25           C
ATOM    150  N7   DG J   2      28.156  76.461  52.947  1.00 91.25           N
ATOM    151  C5   DG J   2      27.427  75.285  52.775  1.00 91.25           C
ATOM    152  C6   DG J   2      26.358  74.965  51.887  1.00 91.25           C
ATOM    153  O6   DG J   2      25.807  75.684  51.041  1.00 91.25           O
ATOM    154  N1   DG J   2      25.923  73.648  52.064  1.00 91.25           N
ATOM    155  C2   DG J   2      26.452  72.759  52.972  1.00 91.25           C
ATOM    156  N2   DG J   2      25.909  71.529  53.005  1.00 91.25           N
ATOM    157  N3   DG J   2      27.438  73.048  53.793  1.00 91.25           N
ATOM    158  C4   DG J   2      27.873  74.314  53.641  1.00 91.25           C
TER
END
"""

def exercise_1():
  """
  The desired outcome of the ncs_search procedure is that user's selections are fully preserved.
  The current procedure results in dropping chains I,J and K,L from selections.
  Chains match as following:
  A <--> D
  B <--> E
  C <--> F
  I <--> K
  J <--> L
  The reason for this is that in input pdb the chains order require the procedure to do the match:
    m m m           m m  (master)
    A-B-C-K-L-D-E-F-I-J
          c c c c c      (copy)
  where in master chains I,J are following protein ABC and in copy chains K,L are preceding
  protein DEF.
  """

  phil_str="""
  ncs_group {
    reference = chain A or chain B or chain C or chain I or chain J
    selection = chain F or chain D or chain E or chain K or chain L
  }
"""
  ncs_groups = get_ncs_groups(phil_str, pdb_str)
  # import iotbx.pdb
  # h = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
  # ncs_groups._show(hierarchy=h, brief=False)
  assert len(ncs_groups) == 1
  assert ncs_groups[0].master_str_selection == "chain A or chain B or chain C or chain I or chain J"
  assert ncs_groups[0].copies[0].str_selection == "chain F or chain D or chain E or chain K or chain L"

if (__name__ == "__main__"):
  t0=time.time()

  exercise_1()

  print("Time: %6.4f"%(time.time()-t0))
  print("OK")
