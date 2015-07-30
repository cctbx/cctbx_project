from __future__ import division

from libtbx.utils import null_out
from cStringIO import StringIO
from mmtbx.secondary_structure.find_ss_from_ca import find_secondary_structure
from libtbx import test_utils

def remove_blank(text):
  return text.replace(" ","").replace("\n","")

std_text="""
ATOM      2  CA  THRAa   3     186.743 125.884 251.259  1.00100.00           C
ATOM      5  CA  ASNAa   4     189.629 123.742 252.763  1.00100.00           C
ATOM      8  CA  SERAa   5     191.072 126.112 255.320  1.00100.00           C
ATOM     11  CA  ASPAa   6     192.080 124.928 258.848  1.00100.00           C
ATOM     14  CA  PHEAa   7     189.384 124.585 261.530  1.00100.00           C
ATOM     17  CA  VALAa   8     189.248 124.466 265.315  1.00100.00           C
ATOM     20  CA  VALAa   9     187.059 122.294 267.547  1.00100.00           C
ATOM     23  CA  ILEAa  10     185.534 123.893 270.679  1.00100.00           C
ATOM     26  CA  LYSAa  11     183.570 122.134 273.450  1.00100.00           C
ATOM     29  CA  ALAAa  12     181.897 124.298 276.085  1.00100.00           C
ATOM     32  CA  LEUAa  13     182.733 123.145 279.601  1.00100.00           C
ATOM     35  CA  GLUAa  14     180.241 125.609 281.156  1.00100.00           C
ATOM     38  CA  ASPAa  15     177.155 127.540 279.985  1.00100.00           C
ATOM     41  CA  GLYAa  16     177.637 130.843 278.162  1.00100.00           C
ATOM     44  CA  VALAa  17     180.958 130.212 276.395  1.00100.00           C
ATOM     47  CA  ASNAa  18     181.477 132.715 273.547  1.00100.00           C
ATOM     50  CA  VALAa  19     183.320 131.753 270.320  1.00100.00           C
ATOM     53  CA  ILEAa  20     184.043 135.156 268.674  1.00100.00           C
ATOM     56  CA  GLYAa  21     185.054 135.558 264.994  1.00100.00           C
ATOM     59  CA  LEUAa  22     187.345 138.529 264.419  1.00100.00           C
ATOM     62  CA  THRAa  23     187.310 140.218 261.033  1.00100.00           C
ATOM     65  CA  ARGAa  24     189.831 139.523 258.335  1.00100.00           C
ATOM     68  CA  GLYAa  25     191.359 142.673 256.805  1.00100.00           C
ATOM     71  CA  ALAAa  26     192.794 146.041 257.837  1.00100.00           C
ATOM     74  CA  ASPAa  27     190.126 146.289 260.564  1.00100.00           C
ATOM     77  CA  THRAa  28     189.912 143.928 263.570  1.00100.00           C
ATOM     80  CA  ARGAa  29     186.413 143.856 265.033  1.00100.00           C
ATOM     83  CA  PHEAa  30     183.873 141.240 266.091  1.00100.00           C
ATOM     86  CA  HISAa  31     181.625 140.079 263.343  1.00100.00           C
ATOM     89  CA  HISAa  32     179.931 137.209 265.203  1.00100.00           C
ATOM     92  CA  SERAa  33     179.805 135.702 268.677  1.00100.00           C
ATOM     95  CA  GLUAa  34     178.501 132.109 268.857  1.00100.00           C
ATOM     98  CA  CYSAa  35     177.222 131.284 272.342  1.00100.00           C
ATOM    101  CA  LEUAa  36     177.646 127.700 273.502  1.00100.00           C
ATOM    104  CA  ASPAa  37     175.969 125.990 276.438  1.00100.00           C
ATOM    107  CA  LYSAa  38     177.682 123.298 278.488  1.00100.00           C
ATOM    110  CA  GLYAa  39     178.623 120.300 276.385  1.00100.00           C
ATOM    113  CA  GLUAa  40     177.892 121.761 272.941  1.00100.00           C
ATOM    116  CA  VALAa  41     180.597 121.439 270.276  1.00100.00           C
ATOM    119  CA  LEUAa  42     181.492 123.998 267.594  1.00100.00           C
ATOM    122  CA  ILEAa  43     183.793 123.155 264.645  1.00100.00           C
ATOM    125  CA  ALAAa  44     184.701 126.388 262.889  1.00100.00           C
ATOM    128  CA  GLNAa  45     186.987 127.209 259.959  1.00100.00           C
ATOM    131  CA  PHEAa  46     189.115 130.161 259.157  1.00100.00           C
ATOM    134  CA  THRAa  47     187.356 131.901 256.203  1.00100.00           C
ATOM    137  CA  GLUAa  48     187.180 134.953 253.965  1.00100.00           C
ATOM    140  CA  HISAa  49     185.578 136.805 256.905  1.00100.00           C
ATOM    143  CA  THRAa  50     187.343 135.292 259.938  1.00100.00           C
ATOM    146  CA  SERAa  51     191.129 135.327 260.339  1.00100.00           C
ATOM    149  CA  ALAAa  52     191.231 135.094 264.170  1.00100.00           C
ATOM    152  CA  ILEAa  53     188.989 133.390 266.744  1.00100.00           C
ATOM    155  CA  LYSAa  54     188.770 134.368 270.428  1.00100.00           C
ATOM    158  CA  VALAa  55     187.303 131.970 273.016  1.00100.00           C
ATOM    161  CA  ARGAa  56     185.817 133.382 276.214  1.00100.00           C
ATOM    164  CA  GLYAa  57     184.672 131.065 278.997  1.00100.00           C
ATOM    167  CA  LYSAa  58     185.698 127.553 280.004  1.00100.00           C
ATOM    170  CA  ALAAa  59     186.172 125.294 276.966  1.00100.00           C
ATOM    173  CA  TYRAa  60     188.258 122.444 275.620  1.00100.00           C
ATOM    176  CA  ILEAa  61     189.863 123.277 272.265  1.00100.00           C
ATOM    179  CA  GLNAa  62     191.492 121.098 269.577  1.00100.00           C
ATOM    182  CA  THRAa  63     193.550 122.431 266.653  1.00100.00           C
ATOM    185  CA  ARGAa  64     196.271 121.116 264.358  1.00100.00           C
ATOM    188  CA  HISAa  65     198.826 122.305 266.995  1.00100.00           C
ATOM    191  CA  GLYAa  66     197.443 120.330 269.914  1.00100.00           C
ATOM    194  CA  VALAa  67     194.865 120.679 272.646  1.00100.00           C
ATOM    197  CA  ILEAa  68     194.232 123.486 275.120  1.00100.00           C
ATOM    200  CA  GLUAa  69     191.576 124.693 277.564  1.00100.00           C
ATOM    203  CA  SERAa  70     190.301 128.219 277.907  1.00100.00           C
ATOM    206  CA  GLUAa  71     189.167 129.249 281.377  1.00100.00           C
ATOM    209  CA  GLYAa  72     186.003 131.073 282.428  1.00100.00           C
"""

helices_text="""
ATOM      2  CA  ALA A   1      11.323  32.055  11.635  1.00 40.00           C
ATOM      7  CA  ALA A   2       8.288  29.768  10.916  1.00 40.00           C
ATOM     12  CA  ALA A   3      10.313  27.854   8.231  1.00 40.00           C
ATOM     17  CA  ALA A   4      13.089  27.116  10.822  1.00 40.00           C
ATOM     22  CA  ALA A   5      10.573  25.488  13.298  1.00 40.00           C
ATOM     27  CA  ALA A   6       9.258  23.514  10.260  1.00 40.00           C
ATOM     32  CA  ALA A   7      12.788  22.543   8.962  1.00 40.00           C
ATOM     37  CA  ALA A   8      13.846  21.459  12.515  1.00 40.00           C
ATOM     42  CA  ALA A   9      10.716  19.261  12.994  1.00 40.00           C
ATOM     47  CA  ALA A  10      11.063  17.985   9.357  1.00 40.00           C
ATOM     52  CA  ALA A  11      14.754  17.018   9.967  1.00 40.00           C
ATOM     57  CA  ALA A  12      13.721  15.483  13.371  1.00 40.00           C
ATOM     62  CA  ALA A  13      10.821  13.516  11.708  1.00 40.00           C
ATOM     67  CA  ALA A  14      13.246  12.367   8.939  1.00 40.00           C
ATOM     72  CA  ALA A  15      15.847  11.407  11.629  1.00 40.00           C
ATOM     77  CA  ALA A  16      13.099   9.317  13.370  1.00 40.00           C
ATOM      2  CA  ALA B   2       1.733  -3.620  -2.296  1.00  1.00
ATOM      7  CA  ALA B   3      -1.902  -4.065  -1.341  1.00  1.00
ATOM     12  CA  ALA B   4      -2.941  -0.441  -1.685  1.00  1.00
ATOM     17  CA  ALA B   5      -0.320   0.578  -4.218  1.00  1.00
ATOM     22  CA  ALA B   6       0.221  -2.836  -5.759  1.00  1.00
ATOM     27  CA  ALA B   7      -3.192  -4.271  -4.973  1.00  1.00
ATOM     32  CA  ALA B   8      -5.081  -0.993  -4.849  1.00  1.00
ATOM     37  CA  ALA B   9      -2.802   0.969  -7.148  1.00  1.00
ATOM     42  CA  ALA B  10      -1.460  -1.967  -9.123  1.00  1.00
ATOM     47  CA  ALA B  11      -4.418  -4.277  -8.632  1.00  1.00
ATOM     52  CA  ALA B  12      -7.044  -1.601  -8.116  1.00  1.00
ATOM     57  CA  ALA B  13      -5.323   1.151 -10.064  1.00  1.00
ATOM     62  CA  ALA B  14      -3.322  -1.073 -12.383  1.00  1.00
ATOM     67  CA  ALA B  15      -5.629  -4.072 -12.291  1.00  1.00
ATOM     72  CA  ALA B  16      -8.822  -2.205 -11.488  1.00  1.00
ATOM     77  CA  ALA B  17      -7.833   1.122 -12.996  1.00  1.00
ATOM     82  CA  ALA B  18      -5.368  -0.211 -15.540  1.00  1.00
ATOM     87  CA  ALA B  19      -6.878  -3.661 -15.920  1.00  1.00
ATOM     92  CA  ALA B  20     -10.423  -2.748 -14.958  1.00  1.00
ATOM     97  CA  ALA B  21     -10.280   0.896 -15.972  1.00  1.00
ATOM    102  CA  ALA B  22      -7.582   0.562 -18.606  1.00  1.00
ATOM      2  CA  ALA C   2       1.202  -3.661  -1.646  1.00  1.00           C
ATOM      7  CA  ALA C   3      -1.466  -2.408  -4.020  1.00  1.00           C
ATOM     12  CA  ALA C   4       1.288  -2.503  -6.614  1.00  1.00           C
ATOM     17  CA  ALA C   5       0.312  -6.139  -7.010  1.00  1.00           C
ATOM     22  CA  ALA C   6      -2.284  -4.816  -9.426  1.00  1.00           C
ATOM     27  CA  ALA C   7       0.502  -5.008 -11.981  1.00  1.00           C
ATOM     32  CA  ALA C   8      -0.579  -8.614 -12.375  1.00  1.00           C
ATOM     37  CA  ALA C   9      -3.100  -7.225 -14.833  1.00  1.00           C
ATOM     42  CA  ALA C  10      -0.285  -7.514 -17.347  1.00  1.00           C
ATOM     47  CA  ALA C  11      -1.470 -11.087 -17.740  1.00  1.00           C
ATOM     52  CA  ALA C  12      -3.913  -9.634 -20.239  1.00  1.00           C
ATOM     57  CA  ALA C  13      -1.074 -10.021 -22.713  1.00  1.00           C
ATOM     62  CA  ALA C  14      -2.362 -13.558 -23.106  1.00  1.00           C
ATOM     67  CA  ALA C  15      -4.725 -12.045 -25.646  1.00  1.00           C
ATOM     72  CA  ALA C  16      -1.865 -12.529 -28.077  1.00  1.00           C
ATOM     77  CA  ALA C  17      -3.254 -16.028 -28.473  1.00  1.00           C
ATOM     82  CA  ALA C  18      -5.534 -14.456 -31.052  1.00  1.00           C
ATOM     87  CA  ALA C  19      -2.657 -15.038 -33.442  1.00  1.00           C
ATOM     92  CA  ALA C  20      -4.146 -18.495 -33.840  1.00  1.00           C
ATOM     97  CA  ALA C  21      -6.342 -16.867 -36.458  1.00  1.00           C
ATOM    102  CA  ALA C  22      -3.451 -17.549 -38.805  1.00  1.00           C
"""

one_full_helix_text="""
ATOM      1  N   ALA A  15      29.207  -2.952  12.868  1.00 16.39           N
ATOM      2  CA  ALA A  15      27.822  -3.418  12.724  1.00 17.10           C
ATOM      3  C   ALA A  15      27.023  -3.016  13.951  1.00 16.98           C
ATOM      4  O   ALA A  15      25.872  -2.551  13.769  1.00 16.78           O
ATOM      5  N   LEU A  16      27.570  -3.117  15.127  1.00 15.97           N
ATOM      6  CA  LEU A  16      26.958  -2.649  16.351  1.00 18.20           C
ATOM      7  C   LEU A  16      26.614  -1.169  16.344  1.00 20.28           C
ATOM      8  O   LEU A  16      25.599  -0.734  16.933  1.00 18.32           O
ATOM      9  N   ILE A  17      27.514  -0.365  15.791  1.00 20.97           N
ATOM     10  CA  ILE A  17      27.343   1.056  15.618  1.00 20.41           C
ATOM     11  C   ILE A  17      26.081   1.392  14.758  1.00 18.17           C
ATOM     12  O   ILE A  17      25.380   2.240  15.282  1.00 16.46           O
ATOM     13  N   SER A  18      25.930   0.759  13.657  1.00 16.97           N
ATOM     14  CA  SER A  18      24.825   0.827  12.744  1.00 19.98           C
ATOM     15  C   SER A  18      23.499   0.405  13.438  1.00 18.89           C
ATOM     16  O   SER A  18      22.557   1.165  13.352  1.00 18.37           O
ATOM     17  N   TRP A  19      23.512  -0.661  14.161  1.00 17.71           N
ATOM     18  CA  TRP A  19      22.492  -1.085  15.081  1.00 15.72           C
ATOM     19  C   TRP A  19      22.083   0.004  16.012  1.00 18.02           C
ATOM     20  O   TRP A  19      20.820   0.244  16.160  1.00 16.93           O
ATOM     21  N   ILE A  20      22.930   0.594  16.823  1.00 14.82           N
ATOM     22  CA  ILE A  20      22.628   1.633  17.766  1.00 15.67           C
ATOM     23  C   ILE A  20      21.917   2.819  17.080  1.00 17.51           C
ATOM     24  O   ILE A  20      20.942   3.365  17.655  1.00 17.70           O
ATOM     25  N   LYS A  21      22.464   3.177  15.957  1.00 16.61           N
ATOM     26  CA  LYS A  21      21.888   4.236  15.157  1.00 19.84           C
ATOM     27  C   LYS A  21      20.436   3.910  14.752  1.00 21.02           C
ATOM     28  O   LYS A  21      19.685   4.899  14.971  1.00 22.80           O
"""
one_helix_beginning_text="""
ATOM      2  CA  ALA A   1      11.323  32.055  11.635  1.00 40.00           C
ATOM      7  CA  ALA A   2       8.288  29.768  10.916  1.00 40.00           C
ATOM     12  CA  ALA A   3      10.313  27.854   8.231  1.00 40.00           C
ATOM     17  CA  ALA A   4      13.089  27.116  10.822  1.00 40.00           C
ATOM     22  CA  ALA A   5      10.573  25.488  13.298  1.00 40.00           C
ATOM     27  CA  ALA A   6       9.258  23.514  10.260  1.00 40.00           C
ATOM     32  CA  ALA A   7      12.788  22.543   8.962  1.00 40.00           C
ATOM     37  CA  ALA A   8      13.846  21.459  12.515  1.00 40.00           C
"""
one_helix_end_text="""
ATOM     42  CA  ALA A   9      10.716  19.261  12.994  1.00 40.00           C
ATOM     47  CA  ALA A  10      11.063  17.985   9.357  1.00 40.00           C
ATOM     52  CA  ALA A  11      14.754  17.018   9.967  1.00 40.00           C
ATOM     57  CA  ALA A  12      13.721  15.483  13.371  1.00 40.00           C
ATOM     62  CA  ALA A  13      10.821  13.516  11.708  1.00 40.00           C
ATOM     67  CA  ALA A  14      13.246  12.367   8.939  1.00 40.00           C
ATOM     72  CA  ALA A  15      15.847  11.407  11.629  1.00 40.00           C
ATOM     77  CA  ALA A  16      13.099   9.317  13.370  1.00 40.00           C
"""
one_helix_middle_text="""
ATOM     22  CA  ALA A   5      10.573  25.488  13.298  1.00 40.00           C
ATOM     27  CA  ALA A   6       9.258  23.514  10.260  1.00 40.00           C
ATOM     32  CA  ALA A   7      12.788  22.543   8.962  1.00 40.00           C
ATOM     37  CA  ALA A   8      13.846  21.459  12.515  1.00 40.00           C
ATOM     42  CA  ALA A   9      10.716  19.261  12.994  1.00 40.00           C
ATOM     47  CA  ALA A  10      11.063  17.985   9.357  1.00 40.00           C
ATOM     52  CA  ALA A  11      14.754  17.018   9.967  1.00 40.00           C
ATOM     57  CA  ALA A  12      13.721  15.483  13.371  1.00 40.00           C
ATOM     62  CA  ALA A  13      10.821  13.516  11.708  1.00 40.00           C
"""
one_helix_text="""
ATOM      2  CA  ALA A   1      11.323  32.055  11.635  1.00 40.00           C
ATOM      7  CA  ALA A   2       8.288  29.768  10.916  1.00 40.00           C
ATOM     12  CA  ALA A   3      10.313  27.854   8.231  1.00 40.00           C
ATOM     17  CA  ALA A   4      13.089  27.116  10.822  1.00 40.00           C
ATOM     22  CA  ALA A   5      10.573  25.488  13.298  1.00 40.00           C
ATOM     27  CA  ALA A   6       9.258  23.514  10.260  1.00 40.00           C
ATOM     32  CA  ALA A   7      12.788  22.543   8.962  1.00 40.00           C
ATOM     37  CA  ALA A   8      13.846  21.459  12.515  1.00 40.00           C
ATOM     42  CA  ALA A   9      10.716  19.261  12.994  1.00 40.00           C
ATOM     47  CA  ALA A  10      11.063  17.985   9.357  1.00 40.00           C
ATOM     52  CA  ALA A  11      14.754  17.018   9.967  1.00 40.00           C
ATOM     57  CA  ALA A  12      13.721  15.483  13.371  1.00 40.00           C
ATOM     62  CA  ALA A  13      10.821  13.516  11.708  1.00 40.00           C
ATOM     67  CA  ALA A  14      13.246  12.367   8.939  1.00 40.00           C
ATOM     72  CA  ALA A  15      15.847  11.407  11.629  1.00 40.00           C
ATOM     77  CA  ALA A  16      13.099   9.317  13.370  1.00 40.00           C
"""

two_helix_text="""
ATOM      2  CA  GLY A   1      43.603 -11.488  24.325  1.00 35.57
ATOM      6  CA  ILE A   2      44.200  -8.183  22.475  1.00 27.55
ATOM     14  CA  GLY A   3      43.999 -10.264  19.329  1.00 21.05
ATOM     18  CA  ALA A   4      40.378 -11.260  20.106  1.00 21.80
ATOM     23  CA  VAL A   5      39.355  -7.658  21.083  1.00 19.34
ATOM     30  CA  LEU A   6      41.062  -6.432  17.957  1.00 17.59
ATOM     38  CA  LYS A   7      39.079  -8.646  15.636  1.00 22.55
ATOM     47  CA  VAL A   8      35.792  -7.369  17.211  1.00 20.52
ATOM     69  CA  THR A  11      34.132  -6.405  12.343  1.00 24.14
ATOM     76  CA  GLY A  12      31.584  -6.595  15.140  1.00 24.17
ATOM     80  CA  LEU A  13      31.923  -2.919  16.364  1.00 23.24
ATOM     88  CA  PRO A  14      31.026  -1.278  13.030  1.00 17.52
ATOM     95  CA  ALA A  15      27.822  -3.418  12.724  1.00 17.10
ATOM    100  CA  LEU A  16      26.958  -2.649  16.351  1.00 18.20
ATOM    108  CA  ILE A  17      27.343   1.056  15.618  1.00 20.41
ATOM    116  CA  SER A  18      24.825   0.827  12.744  1.00 19.98
ATOM    122  CA  TRP A  19      22.492  -1.085  15.081  1.00 15.72
ATOM    136  CA  ILE A  20      22.628   1.633  17.766  1.00 15.67
ATOM    144  CA  LYS A  21      21.888   4.236  15.157  1.00 19.84
ATOM    153  CA  ARG A  22      18.740   2.273  14.020  1.00 20.38
ATOM    164  CA  LYS A  23      17.500   1.928  17.550  1.00 22.62
ATOM    173  CA  ARG A  24      18.059   5.674  18.276  1.00 27.11
ATOM    184  CA  GLN A  25      15.836   6.730  15.339  1.00 37.50
ATOM    193  CA  GLN A  26      13.132   4.360  16.583  1.00 46.66
"""

two_chain_text="""
ATOM    375  N   TYR A  50       6.211 -13.569   8.292  1.00 10.98           N
ATOM    376  CA  TYR A  50       7.318 -12.627   8.487  1.00 10.61           C
ATOM    377  C   TYR A  50       6.766 -11.320   9.020  1.00  9.44           C
ATOM    378  O   TYR A  50       5.608 -10.952   8.800  1.00 10.44           O
ATOM    379  CB  TYR A  50       8.112 -12.427   7.166  1.00 11.77           C
ATOM    380  CG  TYR A  50       8.994 -13.649   6.910  1.00 14.04           C
ATOM    381  CD1 TYR A  50       8.495 -14.805   6.394  1.00 14.79           C
ATOM    382  CD2 TYR A  50      10.334 -13.645   7.280  1.00 15.78           C
ATOM    383  CE1 TYR A  50       9.309 -15.912   6.201  1.00 16.35           C
ATOM    384  CE2 TYR A  50      11.185 -14.715   7.104  1.00 17.74           C
ATOM    385  CZ  TYR A  50      10.649 -15.862   6.545  1.00 17.61           C
ATOM    386  OH  TYR A  50      11.444 -16.974   6.380  1.00 22.90           O
ATOM    387  N   ILE A  51       7.626 -10.601   9.716  1.00  9.43           N
ATOM    388  CA  ILE A  51       7.269  -9.357  10.407  1.00  8.71           C
ATOM    389  C   ILE A  51       7.916  -8.211   9.657  1.00  8.94           C
ATOM    390  O   ILE A  51       9.112  -7.923   9.789  1.00  9.49           O
ATOM    391  CB  ILE A  51       7.636  -9.391  11.898  1.00  9.57           C
ATOM    392  CG1 ILE A  51       6.877 -10.489  12.634  1.00 11.51           C
ATOM    393  CG2 ILE A  51       7.340  -8.025  12.509  1.00  9.99           C
ATOM    394  CD1 ILE A  51       7.189 -11.913  12.447  1.00 13.04           C
ATOM    395  N   TYR A  52       7.103  -7.549   8.818  1.00  8.87           N
ATOM    396  CA  TYR A  52       7.523  -6.401   8.024  1.00  8.95           C
ATOM    397  C   TYR A  52       7.566  -5.196   8.933  1.00  8.79           C
ATOM    398  O   TYR A  52       6.546  -4.870   9.560  1.00  9.43           O
ATOM    399  CB  TYR A  52       6.599  -6.227   6.809  1.00  9.82           C
ATOM    400  CG  TYR A  52       6.769  -7.352   5.826  1.00  9.77           C
ATOM    401  CD1 TYR A  52       6.171  -8.583   6.037  1.00 10.04           C
ATOM    402  CD2 TYR A  52       7.581  -7.227   4.699  1.00 11.49           C
ATOM    403  CE1 TYR A  52       6.330  -9.627   5.177  1.00 11.90           C
ATOM    404  CE2 TYR A  52       7.751  -8.274   3.800  1.00 12.84           C
ATOM    405  CZ  TYR A  52       7.116  -9.468   4.045  1.00 12.86           C
ATOM    406  OH  TYR A  52       7.270 -10.516   3.193  1.00 14.57           O
ATOM    407  N   THR A  53       8.737  -4.586   9.055  1.00  8.30           N
ATOM    408  CA  THR A  53       8.996  -3.643  10.133  1.00  8.29           C
ATOM    409  C   THR A  53       9.384  -2.299   9.550  1.00  8.26           C
ATOM    410  O   THR A  53      10.386  -2.196   8.794  1.00 10.05           O
ATOM    411  CB  THR A  53      10.098  -4.225  11.054  1.00  9.01           C
ATOM    412  OG1 THR A  53       9.695  -5.497  11.548  1.00  9.73           O
ATOM    413  CG2 THR A  53      10.340  -3.320  12.250  1.00 10.75           C
ATOM    414  N   TYR A  54       8.595  -1.292   9.883  1.00  8.50           N
ATOM    415  CA  TYR A  54       8.688   0.046   9.368  1.00  8.54           C
ATOM    416  C   TYR A  54       9.130   0.977  10.480  1.00  9.00           C
ATOM    417  O   TYR A  54       8.469   1.106  11.537  1.00  9.89           O
ATOM    418  CB  TYR A  54       7.350   0.517   8.786  1.00  9.33           C
ATOM    419  CG  TYR A  54       6.866  -0.378   7.660  1.00  9.22           C
ATOM    420  CD1 TYR A  54       6.113  -1.523   7.899  1.00  9.26           C
ATOM    421  CD2 TYR A  54       7.207  -0.083   6.368  1.00  9.45           C
ATOM    422  CE1 TYR A  54       5.683  -2.351   6.899  1.00  9.56           C
ATOM    423  CE2 TYR A  54       6.783  -0.929   5.368  1.00 10.19           C
ATOM    424  CZ  TYR A  54       6.039  -2.054   5.609  1.00  9.31           C
ATOM    425  OH  TYR A  54       5.593  -2.884   4.606  1.00 10.62           O
ATOM    426  N   ARG A  55      10.273   1.635  10.317  1.00  8.99           N
ATOM    427  CA  ARG A  55      10.779   2.540  11.331  1.00  9.30           C
ATOM    428  C   ARG A  55      10.006   3.855  11.284  1.00  9.04           C
ATOM    429  O   ARG A  55       9.809   4.438  10.215  1.00 10.80           O
ATOM    430  CB  ARG A  55      12.255   2.778  11.090  1.00  9.60           C
ATOM    431  CG  ARG A  55      13.089   1.503  11.251  1.00 11.09           C
ATOM    432  CD  ARG A  55      14.531   1.724  10.963  1.00 12.60           C
ATOM    433  NE  ARG A  55      15.393   0.536  10.991  1.00 12.79           N
ATOM    434  CZ  ARG A  55      15.853  -0.067   9.907  1.00 11.24           C
ATOM    435  NH1 ARG A  55      15.601   0.267   8.635  1.00 13.70           N
ATOM    436  NH2 ARG A  55      16.683  -1.110  10.084  1.00 12.83           N
ATOM    437  N   VAL A  56       9.574   4.315  12.463  1.00  9.17           N
ATOM    438  CA  VAL A  56       8.706   5.468  12.645  1.00  9.47           C
ATOM    439  C   VAL A  56       9.399   6.439  13.577  1.00  9.55           C
ATOM    440  O   VAL A  56       9.819   6.065  14.671  1.00 11.35           O
ATOM    441  CB  VAL A  56       7.334   5.020  13.148  1.00 10.23           C
ATOM    442  CG1 VAL A  56       6.464   6.185  13.566  1.00 13.31           C
ATOM    443  CG2 VAL A  56       6.623   4.148  12.100  1.00 11.48           C
ATOM    444  N   SER A  57       9.486   7.688  13.122  1.00 10.07           N
ATOM    445  CA  SER A  57      10.220   8.739  13.797  1.00 11.56           C
ATOM    446  C   SER A  57       9.413  10.020  13.948  1.00 10.72           C
ATOM    447  O   SER A  57       8.594  10.296  13.076  1.00 10.93           O
ATOM    448  CB  SER A  57      11.541   9.078  13.005  1.00 13.02           C
ATOM    449  OG  SER A  57      12.275   7.871  12.757  1.00 16.49           O
ATOM    877  N   LEU B 278      13.003 -13.579  11.217  1.00 14.16           N
ATOM    878  CA  LEU B 278      11.645 -13.391  10.746  1.00 13.77           C
ATOM    879  C   LEU B 278      11.240 -11.932  10.544  1.00 12.73           C
ATOM    880  O   LEU B 278      10.118 -11.639  10.140  1.00 19.80           O
ATOM    881  CB  LEU B 278      10.621 -14.010  11.721  1.00 17.07           C
ATOM    882  CG  LEU B 278      10.887 -15.612  11.784  1.00 20.72           C
ATOM    883  CD1 LEU B 278       9.743 -16.058  12.768  1.00 26.78           C
ATOM    884  CD2 LEU B 278      10.513 -16.229  10.192  1.00 22.51           C
ATOM    885  N   THR B 279      12.124 -11.019  10.844  1.00 10.85           N
ATOM    886  CA  THR B 279      11.888  -9.580  10.658  1.00 10.91           C
ATOM    887  C   THR B 279      12.501  -9.102   9.356  1.00 11.12           C
ATOM    888  O   THR B 279      13.671  -9.407   9.078  1.00 12.89           O
ATOM    889  CB  THR B 279      12.517  -8.801  11.813  1.00 10.60           C
ATOM    890  OG1 THR B 279      11.721  -8.935  12.994  1.00 11.57           O
ATOM    891  CG2 THR B 279      12.676  -7.301  11.552  1.00 12.09           C
ATOM    892  N   ILE B 280      11.699  -8.402   8.572  1.00 10.77           N
ATOM    893  CA  ILE B 280      12.101  -7.808   7.306  1.00 10.70           C
ATOM    894  C   ILE B 280      11.929  -6.297   7.452  1.00  9.56           C
ATOM    895  O   ILE B 280      10.802  -5.815   7.480  1.00 11.05           O
ATOM    896  CB  ILE B 280      11.283  -8.353   6.135  1.00 11.04           C
ATOM    897  CG1 ILE B 280      11.443  -9.854   5.977  1.00 12.99           C
ATOM    898  CG2 ILE B 280      11.674  -7.630   4.848  1.00 13.91           C
ATOM    899  CD1 ILE B 280      10.627 -10.530   4.933  1.00 15.29           C
ATOM    900  N   TYR B 281      13.017  -5.563   7.612  1.00  9.95           N
ATOM    901  CA  TYR B 281      12.918  -4.104   7.622  1.00 10.00           C
ATOM    902  C   TYR B 281      12.476  -3.631   6.240  1.00 10.83           C
ATOM    903  O   TYR B 281      13.016  -4.065   5.216  1.00 13.82           O
ATOM    904  CB  TYR B 281      14.242  -3.422   8.111  1.00 11.29           C
ATOM    905  CG  TYR B 281      14.359  -3.643   9.634  1.00 11.38           C
ATOM    906  CD1 TYR B 281      15.128  -4.632  10.167  1.00 12.02           C
ATOM    907  CD2 TYR B 281      13.622  -2.828  10.483  1.00 12.75           C
ATOM    908  CE1 TYR B 281      15.185  -4.790  11.579  1.00 12.84           C
ATOM    909  CE2 TYR B 281      13.646  -2.973  11.857  1.00 13.88           C
ATOM    910  CZ  TYR B 281      14.427  -3.954  12.373  1.00 14.15           C
ATOM    911  OH  TYR B 281      14.404  -4.084  13.775  1.00 16.24           O
ATOM    912  N   ALA B 282      11.500  -2.767   6.208  1.00 10.01           N
ATOM    913  CA  ALA B 282      10.815  -2.397   4.984  1.00 11.13           C
ATOM    914  C   ALA B 282      10.577  -0.889   4.965  1.00 10.30           C
ATOM    915  O   ALA B 282      10.663  -0.156   5.937  1.00 10.83           O
ATOM    916  CB  ALA B 282       9.496  -3.180   4.899  1.00 14.45           C
"""
def tst_00():
  print "Finding sheets, splitting and merging...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(std_text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())
  records=fss.helix_strand_segments.get_all_pdb_records()
  import iotbx.pdb.secondary_structure as ioss
  annotation=ioss.annotation.from_records(records=flex.split_lines(records))
  f=StringIO()
  print >>f, "New records: \n",annotation.as_pdb_str()
  spl=annotation.split_sheets()
  print >>f, "After split_sheets: \n",spl.as_pdb_str()
  merged=spl.merge_sheets()
  print >>f, "After merge_sheets: \n",merged.as_pdb_str()
  print >>f, "\nSpl:\n",spl.as_pdb_str()
  assert merged.is_same_as(annotation)
  print >>f, "\nComparing merged and spl:"
  print >>f, "\nMerged:\n",merged.as_pdb_str()
  print >>f, "\nSpl:\n",spl.as_pdb_str()
  print >>f, "\nFINAL PDB selections:\n",merged.as_atom_selections()
  assert merged.is_same_as(spl)
  found_text=f.getvalue()

  expected_text="""
New records:
SHEET    1   1 3 HISAa  32  LEUAa  36  0
SHEET    2   1 3 VALAa  17  LEUAa  22 -1  N  GLYAa  21   O  HISAa  32
SHEET    3   1 3 ALAAa  52  VALAa  55 -1  N  LYSAa  54   O  ILEAa  20
SHEET    1   2 4 GLUAa  40  GLNAa  45  0
SHEET    2   2 4 PHEAa   7  ALAAa  12 -1  N  ALAAa  12   O  GLUAa  40
SHEET    3   2 4 LYSAa  58  THRAa  63 -1  N  GLNAa  62   O  VALAa   9
SHEET    4   2 4 GLYAa  66  GLUAa  71 -1  N  SERAa  70   O  ALAAa  59
After split_sheets:
SHEET    1   1 2 HISAa  32  LEUAa  36  0
SHEET    2   1 2 VALAa  17  LEUAa  22 -1  N  GLYAa  21   O  HISAa  32
SHEET    1   2 2 VALAa  17  LEUAa  22  0
SHEET    2   2 2 ALAAa  52  VALAa  55 -1  N  LYSAa  54   O  ILEAa  20
SHEET    1   3 2 GLUAa  40  GLNAa  45  0
SHEET    2   3 2 PHEAa   7  ALAAa  12 -1  N  ALAAa  12   O  GLUAa  40
SHEET    1   4 2 PHEAa   7  ALAAa  12  0
SHEET    2   4 2 LYSAa  58  THRAa  63 -1  N  GLNAa  62   O  VALAa   9
SHEET    1   5 2 LYSAa  58  THRAa  63  0
SHEET    2   5 2 GLYAa  66  GLUAa  71 -1  N  SERAa  70   O  ALAAa  59
After merge_sheets:
SHEET    1   1 2 HISAa  32  LEUAa  36  0
SHEET    2   1 2 VALAa  17  LEUAa  22 -1  N  GLYAa  21   O  HISAa  32
SHEET    2   2 2 ALAAa  52  VALAa  55 -1  N  LYSAa  54   O  ILEAa  20
SHEET    1   3 2 GLUAa  40  GLNAa  45  0
SHEET    2   3 2 PHEAa   7  ALAAa  12 -1  N  ALAAa  12   O  GLUAa  40
SHEET    2   4 2 LYSAa  58  THRAa  63 -1  N  GLNAa  62   O  VALAa   9
SHEET    2   5 2 GLYAa  66  GLUAa  71 -1  N  SERAa  70   O  ALAAa  59

Spl:
SHEET    1   1 2 HISAa  32  LEUAa  36  0
SHEET    2   1 2 VALAa  17  LEUAa  22 -1  N  GLYAa  21   O  HISAa  32
SHEET    1   2 2 VALAa  17  LEUAa  22  0
SHEET    2   2 2 ALAAa  52  VALAa  55 -1  N  LYSAa  54   O  ILEAa  20
SHEET    1   3 2 GLUAa  40  GLNAa  45  0
SHEET    2   3 2 PHEAa   7  ALAAa  12 -1  N  ALAAa  12   O  GLUAa  40
SHEET    1   4 2 PHEAa   7  ALAAa  12  0
SHEET    2   4 2 LYSAa  58  THRAa  63 -1  N  GLNAa  62   O  VALAa   9
SHEET    1   5 2 LYSAa  58  THRAa  63  0
SHEET    2   5 2 GLYAa  66  GLUAa  71 -1  N  SERAa  70   O  ALAAa  59

Comparing merged and spl:

Merged:
SHEET    1   1 2 HISAa  32  LEUAa  36  0
SHEET    2   1 2 VALAa  17  LEUAa  22 -1  N  GLYAa  21   O  HISAa  32
SHEET    2   2 2 ALAAa  52  VALAa  55 -1  N  LYSAa  54   O  ILEAa  20
SHEET    1   3 2 GLUAa  40  GLNAa  45  0
SHEET    2   3 2 PHEAa   7  ALAAa  12 -1  N  ALAAa  12   O  GLUAa  40
SHEET    2   4 2 LYSAa  58  THRAa  63 -1  N  GLNAa  62   O  VALAa   9
SHEET    2   5 2 GLYAa  66  GLUAa  71 -1  N  SERAa  70   O  ALAAa  59

Spl:
SHEET    1   1 2 HISAa  32  LEUAa  36  0
SHEET    2   1 2 VALAa  17  LEUAa  22 -1  N  GLYAa  21   O  HISAa  32
SHEET    1   2 2 VALAa  17  LEUAa  22  0
SHEET    2   2 2 ALAAa  52  VALAa  55 -1  N  LYSAa  54   O  ILEAa  20
SHEET    1   3 2 GLUAa  40  GLNAa  45  0
SHEET    2   3 2 PHEAa   7  ALAAa  12 -1  N  ALAAa  12   O  GLUAa  40
SHEET    1   4 2 PHEAa   7  ALAAa  12  0
SHEET    2   4 2 LYSAa  58  THRAa  63 -1  N  GLNAa  62   O  VALAa   9
SHEET    1   5 2 LYSAa  58  THRAa  63  0
SHEET    2   5 2 GLYAa  66  GLUAa  71 -1  N  SERAa  70   O  ALAAa  59


FINAL PDB selections:
["chain 'Aa' and resid 32  through 36 ", "chain 'Aa' and resid 17  through 22 ", "chain 'Aa' and resid 52  through 55 ", "chain 'Aa' and resid 40  through 45 ", "chain 'Aa' and resid 7  through 12 ", "chain 'Aa' and resid 58  through 63 ", "chain 'Aa' and resid 66  through 71 "]

  """
  if remove_blank(found_text)!=remove_blank(expected_text):
    print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
    raise AssertionError, "FAILED"


  print "OK"


def tst_01():
  print "Finding helices...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(two_helix_text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())

  expected_text="""
Model 1  N: 8  Start: 1 End: 8
Class:  Alpha helix  N: 8 Start: 1 End: 8  Rise: 1.56 A Dot: 0.98

Model 2  N: 16  Start: 11 End: 26
Class:  Alpha helix  N: 16 Start: 11 End: 26  Rise: 1.58 A Dot: 0.98

FINAL PDB RECORDS:
HELIX    1   1 GLY A    1  VAL A    8  1                                   8
HELIX    2   2 THR A   11  GLN A   26  1                                  16

FINAL PDB selections:
" ( chain 'A' and resid 1 through 8 )  or  ( chain 'A' and resid 11 through 26 ) "

"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  #assert not test_utils.show_diff(found_text, expected_text)
  if remove_blank(found_text)!=remove_blank(expected_text):
    print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
    raise AssertionError, "FAILED"
  print "OK"

def tst_02():
  text="""
ATOM      2  CA  GLY A   1      43.603 -11.488  24.325  1.00 35.57
ATOM      6  CA  ILE A   2      44.200  -8.183  22.475  1.00 27.55
ATOM     14  CA  GLY A   3      43.999 -10.264  19.329  1.00 21.05
ATOM     18  CA  ALA A   4      40.378 -11.260  20.106  1.00 21.80
ATOM     23  CA  VAL A   5      39.355  -7.658  21.083  1.00 19.34
ATOM     30  CA  LEU A   6      41.062  -6.432  17.957  1.00 17.59
ATOM     38  CA  LYS A   7      39.079  -8.646  15.636  1.00 22.55
ATOM     47  CA  VAL A   8      35.792  -7.369  17.211  1.00 20.52
ATOM     69  CA  THR A  11      34.132  -6.405  12.343  1.00 24.14
ATOM     76  CA  GLY A  12      31.584  -6.595  15.140  1.00 24.17
ATOM     80  CA  LEU A  13      31.923  -2.919  16.364  1.00 23.24
ATOM     88  CA  PRO A  14      31.026  -1.278  13.030  1.00 17.52
ATOM     95  CA  ALA A  15      27.822  -3.418  12.724  1.00 17.10
ATOM    100  CA  LEU A  16      26.958  -2.649  16.351  1.00 18.20
ATOM    108  CA  ILE A  17      27.343   1.056  15.618  1.00 20.41
ATOM    116  CA  SER A  18      24.825   0.827  12.744  1.00 19.98
ATOM    122  CA  TRP A  19      22.492  -1.085  15.081  1.00 15.72
ATOM    136  CA  ILE A  20      22.628   1.633  17.766  1.00 15.67
ATOM    144  CA  LYS A  21      21.888   4.236  15.157  1.00 19.84
ATOM    153  CA  ARG A  22      18.740   2.273  14.020  1.00 20.38
ATOM    164  CA  LYS A  23      17.500   1.928  17.550  1.00 22.62
ATOM    173  CA  ARG A  24      18.059   5.674  18.276  1.00 27.11
ATOM    184  CA  GLN A  25      15.836   6.730  15.339  1.00 37.50
ATOM    193  CA  GLN A  26      13.132   4.360  16.583  1.00 46.66
"""
  print "Finding helices...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())

  expected_text="""
Model 1  N: 8  Start: 1 End: 8
Class:  Alpha helix  N: 8 Start: 1 End: 8  Rise: 1.56 A Dot: 0.98

Model 2  N: 16  Start: 11 End: 26
Class:  Alpha helix  N: 16 Start: 11 End: 26  Rise: 1.58 A Dot: 0.98

FINAL PDB RECORDS:
HELIX    1   1 GLY A    1  VAL A    8  1                                   8
HELIX    2   2 THR A   11  GLN A   26  1                                  16



FINAL PDB selections:
" ( chain 'A' and resid 1  through 8 )  or  ( chain 'A' and resid 11  through 26 ) "
"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  #assert not test_utils.show_diff(found_text, expected_text)
  if remove_blank(found_text)!=remove_blank(expected_text):
    print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
    raise AssertionError, "FAILED"
  print "OK"

def tst_03():
  print "Finding alpha,3-10 and pi helices...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(helices_text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())

  expected_text="""
Model 1  N: 16  Start: 1 End: 16
Class:  Alpha helix  N: 16 Start: 1 End: 16  Rise: 1.51 A Dot: 0.98

Model 2  N: 21  Start: 2 End: 22
Class:     Pi helix  N: 21 Start: 2 End: 22  Rise: 0.96 A Dot: 0.98

Model 3  N: 21  Start: 2 End: 22
Class:   3-10 helix  N: 20 Start: 2 End: 21  Rise: 1.99 A Dot: 1.00

FINAL PDB RECORDS:
HELIX    1   1 ALA A    1  ALA A   16  1                                  16
HELIX    1   1 ALA C    2  ALA C   21  5                                  20
HELIX    1   1 ALA B    2  ALA B   22  3                                  21



FINAL PDB selections:
" ( chain 'A' and resid 1  through 16 )  or  ( chain 'C' and resid 2  through 21 )  or  ( chain 'B' and resid 2  through 22 ) "
"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  #assert not test_utils.show_diff(found_text, expected_text)
  if remove_blank(found_text)!=remove_blank(expected_text):
    print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
    raise AssertionError, "FAILED"
  print "OK"

def tst_04():
  text="""
ATOM      2  CA  THRAa   3     186.743 125.884 251.259  1.00100.00           C
ATOM      5  CA  ASNAa   4     189.629 123.742 252.763  1.00100.00           C
ATOM      8  CA  SERAa   5     191.072 126.112 255.320  1.00100.00           C
ATOM     11  CA  ASPAa   6     192.080 124.928 258.848  1.00100.00           C
ATOM     14  CA  PHEAa   7     189.384 124.585 261.530  1.00100.00           C
ATOM     17  CA  VALAa   8     189.248 124.466 265.315  1.00100.00           C
ATOM     20  CA  VALAa   9     187.059 122.294 267.547  1.00100.00           C
ATOM     23  CA  ILEAa  10     185.534 123.893 270.679  1.00100.00           C
ATOM     26  CA  LYSAa  11     183.570 122.134 273.450  1.00100.00           C
ATOM     29  CA  ALAAa  12     181.897 124.298 276.085  1.00100.00           C
ATOM     32  CA  LEUAa  13     182.733 123.145 279.601  1.00100.00           C
ATOM     35  CA  GLUAa  14     180.241 125.609 281.156  1.00100.00           C
ATOM     38  CA  ASPAa  15     177.155 127.540 279.985  1.00100.00           C
ATOM     41  CA  GLYAa  16     177.637 130.843 278.162  1.00100.00           C
ATOM     44  CA  VALAa  17     180.958 130.212 276.395  1.00100.00           C
ATOM     47  CA  ASNAa  18     181.477 132.715 273.547  1.00100.00           C
ATOM     50  CA  VALAa  19     183.320 131.753 270.320  1.00100.00           C
ATOM     53  CA  ILEAa  20     184.043 135.156 268.674  1.00100.00           C
ATOM     56  CA  GLYAa  21     185.054 135.558 264.994  1.00100.00           C
ATOM     59  CA  LEUAa  22     187.345 138.529 264.419  1.00100.00           C
ATOM     62  CA  THRAa  23     187.310 140.218 261.033  1.00100.00           C
ATOM     65  CA  ARGAa  24     189.831 139.523 258.335  1.00100.00           C
ATOM     68  CA  GLYAa  25     191.359 142.673 256.805  1.00100.00           C
ATOM     71  CA  ALAAa  26     192.794 146.041 257.837  1.00100.00           C
ATOM     74  CA  ASPAa  27     190.126 146.289 260.564  1.00100.00           C
ATOM     77  CA  THRAa  28     189.912 143.928 263.570  1.00100.00           C
ATOM     80  CA  ARGAa  29     186.413 143.856 265.033  1.00100.00           C
ATOM     83  CA  PHEAa  30     183.873 141.240 266.091  1.00100.00           C
ATOM     86  CA  HISAa  31     181.625 140.079 263.343  1.00100.00           C
ATOM     89  CA  HISAa  32     179.931 137.209 265.203  1.00100.00           C
ATOM     92  CA  SERAa  33     179.805 135.702 268.677  1.00100.00           C
ATOM     95  CA  GLUAa  34     178.501 132.109 268.857  1.00100.00           C
ATOM     98  CA  CYSAa  35     177.222 131.284 272.342  1.00100.00           C
ATOM    101  CA  LEUAa  36     177.646 127.700 273.502  1.00100.00           C
ATOM    104  CA  ASPAa  37     175.969 125.990 276.438  1.00100.00           C
ATOM    107  CA  LYSAa  38     177.682 123.298 278.488  1.00100.00           C
ATOM    110  CA  GLYAa  39     178.623 120.300 276.385  1.00100.00           C
ATOM    113  CA  GLUAa  40     177.892 121.761 272.941  1.00100.00           C
ATOM    116  CA  VALAa  41     180.597 121.439 270.276  1.00100.00           C
ATOM    119  CA  LEUAa  42     181.492 123.998 267.594  1.00100.00           C
ATOM    122  CA  ILEAa  43     183.793 123.155 264.645  1.00100.00           C
ATOM    125  CA  ALAAa  44     184.701 126.388 262.889  1.00100.00           C
ATOM    128  CA  GLNAa  45     186.987 127.209 259.959  1.00100.00           C
ATOM    131  CA  PHEAa  46     189.115 130.161 259.157  1.00100.00           C
ATOM    134  CA  THRAa  47     187.356 131.901 256.203  1.00100.00           C
ATOM    137  CA  GLUAa  48     187.180 134.953 253.965  1.00100.00           C
ATOM    140  CA  HISAa  49     185.578 136.805 256.905  1.00100.00           C
ATOM    143  CA  THRAa  50     187.343 135.292 259.938  1.00100.00           C
ATOM    146  CA  SERAa  51     191.129 135.327 260.339  1.00100.00           C
ATOM    149  CA  ALAAa  52     191.231 135.094 264.170  1.00100.00           C
ATOM    152  CA  ILEAa  53     188.989 133.390 266.744  1.00100.00           C
ATOM    155  CA  LYSAa  54     188.770 134.368 270.428  1.00100.00           C
ATOM    158  CA  VALAa  55     187.303 131.970 273.016  1.00100.00           C
ATOM    161  CA  ARGAa  56     185.817 133.382 276.214  1.00100.00           C
ATOM    164  CA  GLYAa  57     184.672 131.065 278.997  1.00100.00           C
ATOM    167  CA  LYSAa  58     185.698 127.553 280.004  1.00100.00           C
ATOM    170  CA  ALAAa  59     186.172 125.294 276.966  1.00100.00           C
ATOM    173  CA  TYRAa  60     188.258 122.444 275.620  1.00100.00           C
ATOM    176  CA  ILEAa  61     189.863 123.277 272.265  1.00100.00           C
ATOM    179  CA  GLNAa  62     191.492 121.098 269.577  1.00100.00           C
ATOM    182  CA  THRAa  63     193.550 122.431 266.653  1.00100.00           C
ATOM    185  CA  ARGAa  64     196.271 121.116 264.358  1.00100.00           C
ATOM    188  CA  HISAa  65     198.826 122.305 266.995  1.00100.00           C
ATOM    191  CA  GLYAa  66     197.443 120.330 269.914  1.00100.00           C
ATOM    194  CA  VALAa  67     194.865 120.679 272.646  1.00100.00           C
ATOM    197  CA  ILEAa  68     194.232 123.486 275.120  1.00100.00           C
ATOM    200  CA  GLUAa  69     191.576 124.693 277.564  1.00100.00           C
ATOM    203  CA  SERAa  70     190.301 128.219 277.907  1.00100.00           C
ATOM    206  CA  GLUAa  71     189.167 129.249 281.377  1.00100.00           C
ATOM    209  CA  GLYAa  72     186.003 131.073 282.428  1.00100.00           C
"""
  print "Finding sheets...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())

  expected_text="""
Model 1  N: 70  Start: 3 End: 72
Class:  Beta strand  N: 10 Start: 3 End: 12  Rise: 3.32 A Dot: 0.88
Class:  Beta strand  N: 9 Start: 16 End: 24  Rise: 3.24 A Dot: 0.97
Class:  Beta strand  N: 4 Start: 27 End: 30  Rise: 3.34 A Dot: 0.95
Class:  Beta strand  N: 6 Start: 31 End: 36  Rise: 3.29 A Dot: 0.99
Class:  Beta strand  N: 8 Start: 40 End: 47  Rise: 3.30 A Dot: 0.96
Class:  Beta strand  N: 5 Start: 51 End: 55  Rise: 3.41 A Dot: 1.00
Class:  Beta strand  N: 6 Start: 58 End: 63  Rise: 3.41 A Dot: 0.96
Class:  Beta strand  N: 7 Start: 66 End: 72  Rise: 3.41 A Dot: 0.98

FINAL PDB RECORDS:
SHEET    1   1 3 HISAa  32  LEUAa  36  0
SHEET    2   1 3 VALAa  17  LEUAa  22 -1  N  GLYAa  21   O  HISAa  32
SHEET    3   1 3 ALAAa  52  VALAa  55 -1  N  LYSAa  54   O  ILEAa  20
SHEET    1   2 4 GLUAa  40  GLNAa  45  0
SHEET    2   2 4 PHEAa   7  ALAAa  12 -1  N  ALAAa  12   O  GLUAa  40
SHEET    3   2 4 LYSAa  58  THRAa  63 -1  N  GLNAa  62   O  VALAa   9
SHEET    4   2 4 GLYAa  66  GLUAa  71 -1  N  SERAa  70   O  ALAAa  59



FINAL PDB selections:
" ( chain 'Aa' and resid 32  through 36 )  or  ( chain 'Aa' and resid 17  through 22 )  or  ( chain 'Aa' and resid 52  through 55 )  or  ( chain 'Aa' and resid 40  through 45 )  or  ( chain 'Aa' and resid 7  through 12 )  or  ( chain 'Aa' and resid 58  through 63 )  or  ( chain 'Aa' and resid 66  through 71 ) "
"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  #assert not test_utils.show_diff(found_text, expected_text)
  if remove_blank(found_text)!=remove_blank(expected_text):
    print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
    raise AssertionError, "FAILED"
  print "OK"

def tst_05():
  print "Finding sheets with separate chains...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(two_chain_text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())

  expected_text="""
Model 1  N: 8  Start: 50 End: 57
Class:  Beta strand  N: 8 Start: 50 End: 57  Rise: 3.21 A Dot: 0.97

Model 2  N: 5  Start: 278 End: 282
Class:  Beta strand  N: 5 Start: 278 End: 282  Rise: 3.16 A Dot: 0.98

FINAL PDB RECORDS:
SHEET    1   1 2 TYR A  50  TYR A  54  0
SHEET    2   1 2 LEU B 278  ALA B 282  1  N  ILE B 280   O  ILE A  51



FINAL PDB selections:
" ( chain 'A' and resid 50 through 54 )  or  ( chain 'B' and resid 278 through 282 ) "

"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  if remove_blank(found_text)!=remove_blank(expected_text):
    print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
    raise AssertionError, "FAILED"
  print "OK"

def tst_06():
  text="""
ATOM      8  CA  GLY A   2      24.485  19.185   6.248  1.00 11.14           C
HETATM   15  CA  23F A   3      26.939  16.455   5.194  1.00  9.61           C
ATOM     33  CA  ALA A   4      29.149  18.888   3.424  1.00  9.96           C
HETATM   43  CA  23F A   5      30.573  19.304   6.910  1.00  6.42           C
HETATM   61  CA  23F A   6      32.558  16.167   6.280  1.00  6.41           C
ATOM     79  CA  ALA A   7      35.089  18.339   4.563  1.00  6.26           C
HETATM   89  CA  23F A   8      36.195  19.092   8.094  1.00  6.38           C
HETATM  107  CA  23F A   9      38.283  15.914   7.621  1.00  7.78           C
ATOM    125  CA  ALA A  10      40.789  18.180   5.892  1.00  8.66           C
ATOM    135  CA  GLY A  11      41.608  19.716   9.325  1.00 10.78           C
ATOM    142  CA  GLY A  12      44.498  17.479   9.975  1.00 17.00           C
ATOM    149  CA  GLY A  13      43.927  17.193  13.603  1.00 13.58           C
ATOM    156  CA  GLY A  14      41.242  17.379  16.363  1.00 11.14           C
HETATM  163  CA  23F A  15      39.608  20.319  14.616  1.00  7.70           C
ATOM    181  CA  ALA A  16      38.402  17.853  12.023  1.00  7.08           C
ATOM    191  CA  LEU A  17      35.810  16.973  14.649  1.00  6.22           C
HETATM  210  CA  23F A  18      34.098  20.219  13.633  1.00  6.81           C
ATOM    228  CA  ALA A  19      32.642  18.019  10.889  1.00  6.28           C
ATOM    238  CA  LEU A  20      30.139  16.927  13.574  1.00  6.81           C
HETATM  257  CA  23F A  21      28.460  20.242  12.654  1.00  8.80           C
ATOM    275  CA  ALA A  22      27.017  18.382   9.700  1.00  7.89           C
"""
  print "Finding sheets with unusual residues...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())

  expected_text="""
Model 1  N: 21  Start: 2 End: 22
Class:  Alpha helix  N: 10 Start: 3 End: 12  Rise: 2.00 A Dot: 0.98
Class:  Alpha helix  N: 10 Start: 13 End: 22  Rise: 1.96 A Dot: 0.97
FINAL PDB RECORDS:
HELIX    1   1 23F A    3  GLY A   12  1                                  10
HELIX    2   2 GLY A   13  ALA A   22  1                                  10
FINAL PDB selections:
" ( chain 'A' and resid 3 through 12 )  or  ( chain 'A' and resid 13 through 22 ) "
"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  if remove_blank(found_text)!=remove_blank(expected_text):
    print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
    raise AssertionError, "FAILED"
  print "OK"

def tst_07():
  text="""
ATOM    651  CA BPRO E   1      14.350   6.490 -29.205  0.50 16.99           C
ATOM    658  CA BPRO E   2      12.612   6.495 -25.864  0.50 14.37           C
ATOM    666  CA AGLY E   3      12.816   7.962 -32.315  0.55 19.29           C
ATOM    667  CA BGLY E   3      13.074   9.621 -23.839  0.45 12.65           C
ATOM    674  CA APRO E   4      14.350   6.490 -29.205  0.50 16.99           C
ATOM    675  CA BPRO E   4      15.262  10.063 -20.808  0.50 15.34           C
ATOM    688  CA APRO E   5      12.612   6.495 -25.864  0.50 14.37           C
ATOM    689  CA BPRO E   5      14.316   8.840 -17.372  0.50 10.66           C
ATOM    702  CA AGLY E   6      13.074   9.621 -23.839  0.50 12.65           C
ATOM    703  CA BGLY E   6      11.932  10.884 -15.276  0.50  7.99           C
ATOM    710  CA APRO E   7      15.262  10.063 -20.808  0.50 15.34           C
ATOM    711  CA BPRO E   7      13.150  12.796 -12.241  0.50  9.73           C  """
  print "Finding sheets with alt confs where there is no A for first res...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())

  expected_text="""
Model 1  N: 5  Start: 3 End: 7
Class:  Beta strand  N: 4 Start: 4 End: 7  Rise: 3.27 A Dot: 0.91"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  if remove_blank(found_text)!=remove_blank(expected_text):
    print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
    raise AssertionError, "FAILED"
  print "OK"

def tst_08():
  print "Checking similar annotations"

  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(one_helix_text)).construct_hierarchy()

  ann_one_helix_beg=get_annotation(one_helix_beginning_text)
  ann_one_helix_middle=get_annotation(one_helix_middle_text)
  ann_one_helix=get_annotation(one_helix_text)
  for h1 in ann_one_helix_beg.helices:
    for h2 in ann_one_helix_beg.helices:
       print "Should be same:",h1.is_similar_to(other=h2,hierarchy=hierarchy)
       assert h1.is_similar_to(other=h2,hierarchy=hierarchy)
  for maximum_length_difference in [4,8]:
    for minimum_overlap in [6,10]:
      for h1 in ann_one_helix_beg.helices:
        for h2 in ann_one_helix.helices:
           value=h1.is_similar_to(other=h2,hierarchy=hierarchy,
             maximum_length_difference=maximum_length_difference,
             minimum_overlap=minimum_overlap)
           print "Comparison:",value
           assert (value and maximum_length_difference==8 and
              minimum_overlap==6) or not value


  # Now strands and sheets

  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(std_text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())
  records=fss.helix_strand_segments.get_all_pdb_records()
  import iotbx.pdb.secondary_structure as ioss


  s1_records="""
SHEET    1   1 3 HISAa  32  LEUAa  36  0
SHEET    2   1 3 VALAa  17  LEUAa  22 -1  N  GLYAa  21   O  HISAa  32
SHEET    3   1 3 ALAAa  52  VALAa  55 -1  N  LYSAa  54   O  ILEAa  20
"""

  s1_similar_records="""
SHEET    1   1 3 HISAa  32  LEUAa  36  0
SHEET    2   1 3 VALAa  18  LEUAa  22 -1  N  VALAa  19   O  GLUAa  34
SHEET    3   1 3 ALAAa  50  VALAa  57 -1  N  LYSAa  54   O  ILEAa  20
"""

  s1_diff_records="""
SHEET    1   2 4 GLUAa  40  GLNAa  45  0
SHEET    2   2 4 PHEAa   7  ALAAa  12 -1  N  ALAAa  12   O  GLUAa  40
SHEET    3   2 4 LYSAa  58  THRAa  63 -1  N  GLNAa  62   O  VALAa   9
"""
  s1_reverse_diff_records="""
SHEET    1   2 4 GLUAa  40  GLNAa  45  0
SHEET    2   2 4 PHEAa   7  ALAAa  12  1  N  ALAAa  12   O  GLUAa  40
SHEET    3   2 4 LYSAa  58  THRAa  63 -1  N  GLNAa  62   O  VALAa   9
"""

  s1_similar_reg_records="""
SHEET    1   1 3 HISAa  32  LEUAa  36  0
SHEET    2   1 3 VALAa  17  LEUAa  22 -1  N  ILEAa  20   O  SERAa  33
SHEET    3   1 3 ALAAa  52  VALAa  55 -1  N  LYSAa  54   O  ILEAa  20
"""

  s1_similar_backwards="""
SHEET    2   1 3 VALAa  18  LEUAa  22  0
SHEET    3   1 3 ALAAa  50  VALAa  57 -1  O  ILEAa  20   N  LYSAa  54
SHEET    1   1 3 HISAa  32  LEUAa  36 -1  O  SERAa  33   N  ILEAa  20
"""

  s1_similar_backwards_2="""
SHEET    2   1 3 VALAa  18  LEUAa  22  0
SHEET    3   1 3 ALAAa  50  VALAa  57 -1  O  ILEAa  20   N  LYSAa  54
SHEET    1   1 3 HISAa  32  LEUAa  36 -1  O  SERAa  33   N  ILEAa  20
"""

  s1_full=ioss.annotation.from_records(records=flex.split_lines(s1_records))
  s1_similar=ioss.annotation.from_records(records=flex.split_lines(s1_similar_records))
  s1_diff=ioss.annotation.from_records(records=flex.split_lines(s1_diff_records))
  s1_reverse_diff=ioss.annotation.from_records(records=flex.split_lines(s1_reverse_diff_records))

  print "\nChecking similar strands:",
  f=StringIO()
  for s1 in s1_full.sheets:
    for str1 in s1.strands:
      for s2 in s1_similar.sheets:
        for str2 in s2.strands:
          value=str1.is_similar_to(other=str2,hierarchy=hierarchy,
             maximum_length_difference=4)
          print str1.as_atom_selections(),str2.as_atom_selections(),value
          print >>f,str1.as_atom_selections(),str2.as_atom_selections(),value
  assert f.getvalue()=="""chain 'Aa' and resid 32  through 36  chain 'Aa' and resid 32  through 36  True
chain 'Aa' and resid 32  through 36  chain 'Aa' and resid 18  through 22  False
chain 'Aa' and resid 32  through 36  chain 'Aa' and resid 50  through 57  False
chain 'Aa' and resid 17  through 22  chain 'Aa' and resid 32  through 36  False
chain 'Aa' and resid 17  through 22  chain 'Aa' and resid 18  through 22  True
chain 'Aa' and resid 17  through 22  chain 'Aa' and resid 50  through 57  False
chain 'Aa' and resid 52  through 55  chain 'Aa' and resid 32  through 36  False
chain 'Aa' and resid 52  through 55  chain 'Aa' and resid 18  through 22  False
chain 'Aa' and resid 52  through 55  chain 'Aa' and resid 50  through 57  True
"""

  print "\nChecking different strands:",
  for s1 in s1_full.sheets:
    for str1 in s1.strands:
      for s2 in s1_diff.sheets:
        for str2 in s2.strands:
          print str1.as_atom_selections(),str2.as_atom_selections(),\
            str1.is_similar_to(other=str2,hierarchy=hierarchy,
             maximum_length_difference=4)

  print "\nChecking similar sheets:",
  for s1 in s1_full.sheets:
    for s2 in s1_similar.sheets:
      value=s1.is_similar_to(other=s2,hierarchy=hierarchy,
            maximum_length_difference=4)
      print value
      assert value

  print "\nChecking non-similar sheets:",
  for s1 in s1_full.sheets:
    for s2 in s1_diff.sheets:
      value=s1.is_similar_to(other=s2,hierarchy=hierarchy,
            maximum_length_difference=4)
      print value
      assert not value

  print "\nChecking similar overall annotations:",
  value=s1_full.is_similar_to(other=s1_similar,hierarchy=hierarchy,
            maximum_length_difference=4)
  print value
  assert value

  print "\nChecking different overall annotations:",
  value=s1_full.is_similar_to(other=s1_diff,hierarchy=hierarchy,
            maximum_length_difference=4)
  print value
  assert not value

  print "\nChecking different overall directions:",
  value=s1_full.is_similar_to(other=s1_reverse_diff,hierarchy=hierarchy,
            maximum_length_difference=4)
  print value
  assert not value


  # parallel strands
  print "\n\nChecking parallel strands..."

  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(two_chain_text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())
  records=fss.helix_strand_segments.get_all_pdb_records()
  import iotbx.pdb.secondary_structure as ioss


  s2_records="""
FINAL PDB RECORDS:
SHEET    1   1 2 TYR A  50  TYR A  54  0
SHEET    2   1 2 LEU B 278  ALA B 282  1  N  ILE B 280   O  ILE A  51
"""

  s2_similar_records="""
FINAL PDB RECORDS:
SHEET    1   1 2 TYR A  50  TYR A  54  0
SHEET    2   1 2 LEU B 278  ALA B 282  1  N  ALA B 282   O  THR A  53
"""

  s2_similar_records_2="""
FINAL PDB RECORDS:
SHEET    1   1 2 TYR A  50  TYR A  54  0
SHEET    2   1 2 LEU B 278  ALA B 282  1  O  LEU B 278   N  ILE A  51
"""

  s2_different_records="""
FINAL PDB RECORDS:
SHEET    1   1 2 TYR A  50  TYR A  54  0
SHEET    2   1 2 LEU B 278  ALA B 282  1  O  ILE B 280   N  ILE A  51
"""

  s2_full=ioss.annotation.from_records(records=flex.split_lines(s2_records))
  s2_similar=ioss.annotation.from_records(records=flex.split_lines(s2_similar_records))
  s2_similar_2=ioss.annotation.from_records(records=flex.split_lines(s2_similar_records_2))
  s2_different=ioss.annotation.from_records(records=flex.split_lines(s2_different_records))

  print "\nChecking similar overall annotations (offset by 2):",
  value=s2_full.is_similar_to(other=s2_similar,hierarchy=hierarchy,
            maximum_length_difference=4)
  print value
  assert value

  print "\nChecking similar overall annotations (switch N/O):",
  value=s2_full.is_similar_to(other=s2_similar_2,hierarchy=hierarchy,
            maximum_length_difference=4)
  print value
  assert value

  print "\nChecking different overall annotations (offset by 1):",
  value=s2_full.is_similar_to(other=s2_different,hierarchy=hierarchy,
            maximum_length_difference=4)
  print value
  assert not value

  print "\nOK"



def tst_09():
  print "Comparing sheets and helices...",

  helix_1="""
HELIX    1   1 GLY A    1  VAL A    8  1                                   8
"""
  helix_2="""
HELIX    2   2 THR A   11  GLN A   26  1                                  16
"""
  sheet_1="""
SHEET    1   1 3 HISAa  32  LEUAa  36  0
SHEET    2   1 3 VALAa  17  LEUAa  22 -1  N  GLYAa  21   O  HISAa  32
SHEET    3   1 3 ALAAa  52  VALAa  55  1  N  LYSAa  54   O  ILEAa  20
"""
  sheet_2="""
SHEET    1   2 4 GLUAa  40  GLNAa  45  0
SHEET    2   2 4 PHEAa   7  ALAAa  12 -1  N  ALAAa  12   O  GLUAa  40
SHEET    3   2 4 LYSAa  58  THRAa  63  1  N  GLNAa  62   O  VALAa   9
SHEET    4   2 4 GLYAa  66  GLUAa  71 -1  N  SERAa  70   O  ALAAa  59
"""
  import iotbx.pdb.secondary_structure as ioss
  from cctbx.array_family import flex

  h1=ioss.annotation.from_records(records=flex.split_lines(helix_1))
  h2=ioss.annotation.from_records(records=flex.split_lines(helix_2))
  s1=ioss.annotation.from_records(records=flex.split_lines(sheet_1))
  s2=ioss.annotation.from_records(records=flex.split_lines(sheet_2))
  assert h1.is_same_as(h1)
  assert h2.is_same_as(h2)
  assert not h1.is_same_as(h2)
  assert not h1.is_same_as(s1)
  assert not s1.is_same_as(s2)
  assert s1.is_same_as(s1)
  for a in s1.sheets:
    for b in s1.sheets:
      assert (a==b and a.is_same_as(b)) or (not a.is_same_as(b))
      for sa in a.strands:
        for sb in b.strands:
          assert (sa==sb and sa.is_same_as(sb)) or (not sa.is_same_as(sb))
  print "OK"


def get_annotation(text):
  import iotbx.pdb
  from cctbx.array_family import flex
  import iotbx.pdb.secondary_structure as ioss
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())
  records=fss.helix_strand_segments.get_all_pdb_records()
  return ioss.annotation.from_records(records=flex.split_lines(records))

def tst_10():

  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(one_full_helix_text+two_chain_text)).construct_hierarchy()

  text_helix_1="""
HELIX    1   1 ALA A   15  LYS A   21  1                                   7
"""
  text_helix_2="""
HELIX    1   1 LEU A   16  LYS A   21  1                                   6
"""
  text_sheet_1="""
SHEET    1   1 2 ILE A  51  TYR A  54  0
SHEET    2   1 2 THR B 279  ALA B 282  1  N  ILE B 280   O  ILE A  51
"""
  text_sheet_2="""
SHEET    1   1 2 TYR A  50  TYR A  54  0
SHEET    2   1 2 LEU B 278  ALA B 282  1  N  ILE B 280   O  ILE A  51
"""
  import iotbx.pdb.secondary_structure as ioss
  from cctbx.array_family import flex
  h1=ioss.annotation.from_records(records=flex.split_lines(text_helix_1))
  h2=ioss.annotation.from_records(records=flex.split_lines(text_helix_2))
  s1=ioss.annotation.from_records(records=flex.split_lines(text_sheet_1))
  s2=ioss.annotation.from_records(records=flex.split_lines(text_sheet_2))

  hs1=ioss.annotation.from_records(records=flex.split_lines(text_helix_1+text_sheet_1))
  hs2=ioss.annotation.from_records(records=flex.split_lines(text_helix_2+text_sheet_2))

  print "\nCombining annotations"

  ann_all=get_annotation(one_full_helix_text+two_chain_text)
  print "\nFull annotation\n",ann_all.as_pdb_str()

  print "\nCombining annotations from two parts"
  ann_one_full_helix=get_annotation(one_full_helix_text)
  print "\nAnnotation for helix:\n",ann_one_full_helix.as_pdb_str()
  ann_two_chain=get_annotation(two_chain_text)
  print "\nAnnotation for two chains:\n",ann_two_chain.as_pdb_str()

  ann_combined=ann_one_full_helix.combine_annotations(other=ann_two_chain,
    hierarchy=hierarchy)
  if ann_combined:
    print "Combined: \n",ann_combined.as_pdb_str()

  print "\nCombining annotations from overlapping helix annotations"
  print "\nHelix 1 and 2: \n",h1.as_pdb_str(),"\n",h2.as_pdb_str()
  ann_combined=h1.combine_annotations(other=h2,hierarchy=hierarchy)
  print "\nCombined: \n",ann_combined.as_pdb_str()

  print "\nCombining annotations from overlapping strand annotations"
  print "\nStrand 1 and 2: \n",s1.as_pdb_str(),"\n",s2.as_pdb_str()
  ann_combined=s1.combine_annotations(other=s2,hierarchy=hierarchy,
     minimum_overlap=3,out=sys.stdout)
  print "\nCombined: \n",ann_combined.as_pdb_str()

  print "\nCombining annotations from overlapping strand and helix annotations"
  print "\nStrand/helix 1:\n",hs1.as_pdb_str(),"Strand/helix 2:\n",hs2.as_pdb_str()
  ann_combined=hs1.combine_annotations(other=hs2,hierarchy=hierarchy,
     minimum_overlap=3,out=sys.stdout)
  print "\nCombined: \n",ann_combined.as_pdb_str()


  print "OK"

def tst_11():

  print "\nCounting H-bonds"

  import iotbx.pdb
  import iotbx.pdb.secondary_structure as ioss
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(two_chain_text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())
  ann=fss.get_annotation()
  print ann.as_pdb_str()

  print "Good H-bonds: %d  Poor H-Bonds: %d" %(
         fss.number_of_good_h_bonds,
         fss.number_of_poor_h_bonds,)
  assert fss.number_of_good_h_bonds==4 and fss.number_of_poor_h_bonds==0

  print "\nCounting H-bonds using ioss.annotation:"

  number_of_good_h_bonds,number_of_poor_h_bonds=ann.count_h_bonds(
    hierarchy=hierarchy)
  print "Good H-bonds: %d  Poor H-Bonds: %d" %(
         fss.number_of_good_h_bonds,
         fss.number_of_poor_h_bonds,)
  assert fss.number_of_good_h_bonds==4 and fss.number_of_poor_h_bonds==0

  print "\nCounting residues in secondary structure:",
  print ann.count_residues(hierarchy=hierarchy)
  assert ann.count_residues(hierarchy=hierarchy)==10

  print "\nCounting H-bonds in helix:"

  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(one_full_helix_text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())
  ann=fss.get_annotation()
  print ann.as_pdb_str()

  print "\nH-bonds with cutoff=3.5 (default):\n"
  number_of_good_h_bonds,number_of_poor_h_bonds=ann.count_h_bonds(
    hierarchy=hierarchy)
  print "Good H-bonds: %d  Poor H-Bonds: %d" %(
         number_of_good_h_bonds,
         number_of_poor_h_bonds,)
  assert number_of_good_h_bonds==3 and number_of_poor_h_bonds==0

  print "\nH-bonds with cutoff=3.0:\n"
  number_of_good_h_bonds,number_of_poor_h_bonds=ann.count_h_bonds(
    hierarchy=hierarchy,max_h_bond_length=3.0)
  print "Good H-bonds: %d  Poor H-Bonds: %d" %(
         number_of_good_h_bonds,
         number_of_poor_h_bonds,)
  assert number_of_good_h_bonds==1 and number_of_poor_h_bonds==2

  print "\nCount number of residues in secondary structure:",
  print ann.count_residues(hierarchy=hierarchy)
  assert ann.count_residues(hierarchy=hierarchy) ==7

  print "\nH-bonds in mixed helix/strand"

  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(two_chain_text+one_full_helix_text)
         ).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())
  ann=fss.get_annotation()
  print ann.as_pdb_str()

  print "\nH-bonds with cutoff=3.0 :\n"
  number_of_good_h_bonds,number_of_poor_h_bonds=ann.count_h_bonds(
    hierarchy=hierarchy,max_h_bond_length=3.0)
  print "Good H-bonds: %d  Poor H-Bonds: %d" %(
         number_of_good_h_bonds,
         number_of_poor_h_bonds,)
  assert number_of_good_h_bonds==5 and number_of_poor_h_bonds==2

  print "\nCount number of residues in secondary structure:",
  print ann.count_residues(hierarchy=hierarchy)
  assert ann.count_residues(hierarchy=hierarchy) ==17

def tst_12():

  text="""
ATOM      8  CA  GLY A   2      24.485  19.185   6.248  1.00 11.14           C
HETATM   15  CA  23F A   3      26.939  16.455   5.194  1.00  9.61           C
ATOM     33  CA  ALA A   4      29.149  18.888   3.424  1.00  9.96           C
HETATM   43  CA  23F A   5      30.573  19.304   6.910  1.00  6.42           C
HETATM   61  CA  23F A   6      32.558  16.167   6.280  1.00  6.41           C
ATOM     79  CA  ALA A   7      35.089  18.339   4.563  1.00  6.26           C
HETATM   89  CA  23F A   8      36.195  19.092   8.094  1.00  6.38           C
HETATM  107  CA  23F A   9      38.283  15.914   7.621  1.00  7.78           C
ATOM    125  CA  ALA A  10      40.789  18.180   5.892  1.00  8.66           C
ATOM    135  CA  GLY A  11      41.608  19.716   9.325  1.00 10.78           C
ATOM    142  CA  GLY A  12      44.498  17.479   9.975  1.00 17.00           C
ATOM    149  CA  GLY A  13      43.927  17.193  13.603  1.00 13.58           C
ATOM    156  CA  GLY A  14      41.242  17.379  16.363  1.00 11.14           C
HETATM  163  CA  23F A  15      39.608  20.319  14.616  1.00  7.70           C
ATOM    181  CA  ALA A  16      38.402  17.853  12.023  1.00  7.08           C
ATOM    191  CA  LEU A  17      35.810  16.973  14.649  1.00  6.22           C
HETATM  210  CA  23F A  18      34.098  20.219  13.633  1.00  6.81           C
ATOM    228  CA  ALA A  19      32.642  18.019  10.889  1.00  6.28           C
ATOM    238  CA  LEU A  20      30.139  16.927  13.574  1.00  6.81           C
HETATM  257  CA  23F A  21      28.460  20.242  12.654  1.00  8.80           C
ATOM    275  CA  ALA A  22      27.017  18.382   9.700  1.00  7.89           C
"""
if __name__=="__main__":
  import sys
  tst_00()
  tst_01()
  tst_02()
  tst_03()
  tst_04()
  tst_05()
  tst_06()
  tst_07()
  tst_08()
  tst_09()
  tst_10()
  tst_11()
  tst_12()
  print "OK"
