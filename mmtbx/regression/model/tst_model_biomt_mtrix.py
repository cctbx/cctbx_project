from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
import time

"""
Test multiplication of hierarchy and SS annotations in different combinations
of MTRIX and BIOMT records presence.
"""

single_mtrix_txt = """
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.479787 -0.038259 -0.876550        0.00000
MTRIX2   2 -0.530698  0.782918 -0.324654        0.00000
MTRIX3   2  0.698688  0.620947  0.355330        0.00000
"""


mtrix_txt = """
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.479787 -0.038259 -0.876550        0.00000
MTRIX2   2 -0.530698  0.782918 -0.324654        0.00000
MTRIX3   2  0.698688  0.620947  0.355330        0.00000
MTRIX1   3 -0.361936 -0.592602 -0.719600        0.00000
MTRIX2   3 -0.896947  0.431671  0.095646        0.00000
MTRIX3   3  0.253950  0.680060 -0.687769        0.00000
"""

biomt_txt = """
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2  0.500000 -0.809017  0.309017        0.00000
REMARK 350   BIOMT2   2  0.809017  0.309017 -0.500000        0.00000
REMARK 350   BIOMT3   2  0.309017  0.500000  0.809017        0.00000
REMARK 350   BIOMT1   3 -0.309017 -0.500000  0.809017        0.00000
REMARK 350   BIOMT2   3  0.500000 -0.809017 -0.309017        0.00000
REMARK 350   BIOMT3   3  0.809017  0.309017  0.500000        0.00000
"""

ss_txt = """
HELIX    6   6 ARG A  316  LEU A  318  5                                   3
HELIX    7   7 SER A  335  ASN A  341  1                                   7
SHEET    1   E 2 TYR A 305  SER A 308  0
SHEET    2   E 2 GLN A 311  GLU A 314 -1  O  ARG A 313   N  PHE A 306

"""

# 300 atoms
atoms_txt = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
ATOM   2065  N   GLY A 304       3.950 -35.449 102.015  1.00 21.30           N
ATOM   2066  CA  GLY A 304       4.631 -35.764 103.257  1.00 19.87           C
ATOM   2067  C   GLY A 304       6.074 -36.196 103.097  1.00 19.37           C
ATOM   2068  O   GLY A 304       6.642 -36.823 103.994  1.00 18.69           O
ATOM   2069  N   TYR A 305       6.673 -35.863 101.957  1.00 19.21           N
ATOM   2070  CA  TYR A 305       8.065 -36.210 101.688  1.00 18.15           C
ATOM   2071  C   TYR A 305       8.975 -35.002 101.819  1.00 18.33           C
ATOM   2072  O   TYR A 305       8.673 -33.920 101.314  1.00 19.14           O
ATOM   2073  CB  TYR A 305       8.202 -36.803 100.289  1.00 17.45           C
ATOM   2074  CG  TYR A 305       7.826 -38.258 100.228  1.00 18.23           C
ATOM   2075  CD1 TYR A 305       8.761 -39.250 100.515  1.00 18.89           C
ATOM   2076  CD2 TYR A 305       6.526 -38.647  99.917  1.00 17.60           C
ATOM   2077  CE1 TYR A 305       8.411 -40.595 100.493  1.00 18.65           C
ATOM   2078  CE2 TYR A 305       6.166 -39.982  99.895  1.00 18.16           C
ATOM   2079  CZ  TYR A 305       7.112 -40.954 100.182  1.00 18.93           C
ATOM   2080  OH  TYR A 305       6.756 -42.284 100.149  1.00 20.67           O
ATOM   2081  N   PHE A 306      10.097 -35.200 102.498  1.00 18.08           N
ATOM   2082  CA  PHE A 306      11.061 -34.133 102.707  1.00 18.50           C
ATOM   2083  C   PHE A 306      12.458 -34.585 102.325  1.00 19.27           C
ATOM   2084  O   PHE A 306      12.737 -35.780 102.239  1.00 19.20           O
ATOM   2085  CB  PHE A 306      11.071 -33.720 104.181  1.00 18.35           C
ATOM   2086  CG  PHE A 306       9.784 -33.115 104.650  1.00 18.59           C
ATOM   2087  CD1 PHE A 306       9.582 -31.744 104.580  1.00 18.76           C
ATOM   2088  CD2 PHE A 306       8.765 -33.919 105.141  1.00 18.71           C
ATOM   2089  CE1 PHE A 306       8.379 -31.180 104.994  1.00 19.56           C
ATOM   2090  CE2 PHE A 306       7.558 -33.366 105.557  1.00 20.09           C
ATOM   2091  CZ  PHE A 306       7.365 -31.993 105.483  1.00 20.22           C
ATOM   2092  N   MET A 307      13.330 -33.615 102.088  1.00 20.68           N
ATOM   2093  CA  MET A 307      14.717 -33.899 101.773  1.00 21.16           C
ATOM   2094  C   MET A 307      15.423 -33.667 103.103  1.00 22.68           C
ATOM   2095  O   MET A 307      15.811 -32.544 103.417  1.00 23.48           O
ATOM   2096  CB  MET A 307      15.243 -32.920 100.729  1.00 20.61           C
ATOM   2097  CG  MET A 307      16.666 -33.197 100.265  1.00 21.45           C
ATOM   2098  SD  MET A 307      16.835 -34.726  99.306  1.00 21.53           S
ATOM   2099  CE  MET A 307      18.018 -35.567 100.270  1.00 21.98           C
ATOM   2100  N   SER A 308      15.551 -34.726 103.900  1.00 23.90           N
ATOM   2101  CA  SER A 308      16.196 -34.629 105.204  1.00 25.50           C
ATOM   2102  C   SER A 308      17.711 -34.736 105.044  1.00 26.51           C
ATOM   2103  O   SER A 308      18.287 -35.817 105.177  1.00 26.25           O
ATOM   2104  CB  SER A 308      15.681 -35.736 106.127  1.00 25.33           C
ATOM   2105  OG  SER A 308      16.066 -35.503 107.471  1.00 28.15           O
ATOM   2106  N   ASN A 309      18.344 -33.602 104.754  1.00 27.44           N
ATOM   2107  CA  ASN A 309      19.787 -33.533 104.556  1.00 28.19           C
ATOM   2108  C   ASN A 309      20.235 -34.192 103.260  1.00 28.98           C
ATOM   2109  O   ASN A 309      20.106 -33.612 102.183  1.00 30.58           O
ATOM   2110  CB  ASN A 309      20.522 -34.170 105.737  1.00 28.18           C
ATOM   2111  CG  ASN A 309      20.631 -33.238 106.923  1.00 28.56           C
ATOM   2112  OD1 ASN A 309      21.308 -32.212 106.855  1.00 28.58           O
ATOM   2113  ND2 ASN A 309      19.963 -33.585 108.017  1.00 27.94           N
ATOM   2114  N   ASP A 310      20.753 -35.410 103.369  1.00 28.99           N
ATOM   2115  CA  ASP A 310      21.249 -36.144 102.212  1.00 28.59           C
ATOM   2116  C   ASP A 310      20.246 -37.125 101.612  1.00 27.76           C
ATOM   2117  O   ASP A 310      20.386 -37.526 100.457  1.00 28.39           O
ATOM   2118  CB  ASP A 310      22.514 -36.904 102.598  1.00 30.83           C
ATOM   2119  CG  ASP A 310      22.271 -37.900 103.717  1.00 33.72           C
ATOM   2120  OD1 ASP A 310      21.917 -37.469 104.838  1.00 32.95           O
ATOM   2121  OD2 ASP A 310      22.426 -39.118 103.473  1.00 34.80           O
ATOM   2122  N   GLN A 311      19.238 -37.512 102.388  1.00 26.40           N
ATOM   2123  CA  GLN A 311      18.239 -38.461 101.910  1.00 24.95           C
ATOM   2124  C   GLN A 311      16.812 -37.956 101.941  1.00 23.86           C
ATOM   2125  O   GLN A 311      16.499 -36.941 102.561  1.00 23.90           O
ATOM   2126  CB  GLN A 311      18.268 -39.739 102.741  1.00 27.17           C
ATOM   2127  CG  GLN A 311      19.513 -40.558 102.642  1.00 31.00           C
ATOM   2128  CD  GLN A 311      19.343 -41.880 103.340  1.00 32.73           C
ATOM   2129  OE1 GLN A 311      18.986 -41.928 104.519  1.00 32.73           O
ATOM   2130  NE2 GLN A 311      19.591 -42.969 102.618  1.00 35.40           N
ATOM   2131  N   ILE A 312      15.948 -38.706 101.272  1.00 22.64           N
ATOM   2132  CA  ILE A 312      14.529 -38.409 101.221  1.00 20.51           C
ATOM   2133  C   ILE A 312      13.903 -39.134 102.404  1.00 21.02           C
ATOM   2134  O   ILE A 312      14.209 -40.300 102.653  1.00 21.21           O
ATOM   2135  CB  ILE A 312      13.886 -38.956  99.933  1.00 19.13           C
ATOM   2136  CG1 ILE A 312      14.452 -38.229  98.715  1.00 18.32           C
ATOM   2137  CG2 ILE A 312      12.375 -38.824 100.008  1.00 17.77           C
ATOM   2138  CD1 ILE A 312      13.918 -38.745  97.395  1.00 17.16           C
ATOM   2139  N   ARG A 313      13.042 -38.442 103.140  1.00 21.29           N
ATOM   2140  CA  ARG A 313      12.363 -39.049 104.275  1.00 21.69           C
ATOM   2141  C   ARG A 313      10.872 -38.764 104.223  1.00 22.24           C
ATOM   2142  O   ARG A 313      10.445 -37.705 103.766  1.00 22.50           O
ATOM   2143  CB  ARG A 313      12.932 -38.539 105.601  1.00 21.74           C
ATOM   2144  CG  ARG A 313      14.216 -39.223 106.028  1.00 23.27           C
ATOM   2145  CD  ARG A 313      14.488 -39.007 107.511  1.00 24.27           C
ATOM   2146  NE  ARG A 313      15.647 -39.768 107.970  1.00 26.04           N
ATOM   2147  CZ  ARG A 313      16.906 -39.483 107.652  1.00 26.93           C
ATOM   2148  NH1 ARG A 313      17.177 -38.443 106.873  1.00 26.76           N
ATOM   2149  NH2 ARG A 313      17.895 -40.244 108.103  1.00 27.18           N
ATOM   2150  N   GLU A 314      10.085 -39.729 104.686  1.00 23.21           N
ATOM   2151  CA  GLU A 314       8.637 -39.591 104.709  1.00 23.33           C
ATOM   2152  C   GLU A 314       8.256 -39.096 106.107  1.00 23.56           C
ATOM   2153  O   GLU A 314       8.950 -39.370 107.084  1.00 23.95           O
ATOM   2154  CB  GLU A 314       7.990 -40.946 104.405  1.00 23.60           C
ATOM   2155  CG  GLU A 314       6.517 -40.906 104.006  1.00 25.54           C
ATOM   2156  CD  GLU A 314       5.571 -40.837 105.196  1.00 27.64           C
ATOM   2157  OE1 GLU A 314       5.803 -41.568 106.184  1.00 27.37           O
ATOM   2158  OE2 GLU A 314       4.586 -40.068 105.137  1.00 27.35           O
ATOM   2159  N   ARG A 315       7.162 -38.349 106.187  1.00 23.65           N
ATOM   2160  CA  ARG A 315       6.672 -37.790 107.443  1.00 22.98           C
ATOM   2161  C   ARG A 315       6.921 -38.651 108.687  1.00 22.60           C
ATOM   2162  O   ARG A 315       7.532 -38.196 109.653  1.00 22.85           O
ATOM   2163  CB  ARG A 315       5.173 -37.516 107.323  1.00 23.23           C
ATOM   2164  CG  ARG A 315       4.735 -36.167 107.842  1.00 24.91           C
ATOM   2165  CD  ARG A 315       3.231 -36.143 108.070  1.00 28.73           C
ATOM   2166  NE  ARG A 315       2.444 -36.339 106.853  1.00 31.17           N
ATOM   2167  CZ  ARG A 315       2.358 -35.454 105.863  1.00 34.22           C
ATOM   2168  NH1 ARG A 315       3.018 -34.301 105.937  1.00 36.07           N
ATOM   2169  NH2 ARG A 315       1.596 -35.715 104.805  1.00 32.99           N
ATOM   2170  N   ARG A 316       6.447 -39.893 108.657  1.00 22.00           N
ATOM   2171  CA  ARG A 316       6.578 -40.806 109.792  1.00 21.99           C
ATOM   2172  C   ARG A 316       7.984 -41.091 110.302  1.00 23.09           C
ATOM   2173  O   ARG A 316       8.149 -41.529 111.439  1.00 23.65           O
ATOM   2174  CB  ARG A 316       5.886 -42.136 109.475  1.00 21.22           C
ATOM   2175  CG  ARG A 316       4.373 -42.042 109.402  1.00 23.38           C
ATOM   2176  CD  ARG A 316       3.836 -42.668 108.123  1.00 26.73           C
ATOM   2177  NE  ARG A 316       3.995 -44.121 108.085  1.00 28.42           N
ATOM   2178  CZ  ARG A 316       4.479 -44.789 107.040  1.00 28.61           C
ATOM   2179  NH1 ARG A 316       4.857 -44.132 105.954  1.00 29.34           N
ATOM   2180  NH2 ARG A 316       4.571 -46.113 107.074  1.00 27.85           N
ATOM   2181  N   ASP A 317       9.000 -40.846 109.485  1.00 25.17           N
ATOM   2182  CA  ASP A 317      10.363 -41.127 109.913  1.00 26.92           C
ATOM   2183  C   ASP A 317      11.076 -39.946 110.566  1.00 26.95           C
ATOM   2184  O   ASP A 317      12.094 -40.123 111.234  1.00 27.91           O
ATOM   2185  CB  ASP A 317      11.184 -41.652 108.733  1.00 31.20           C
ATOM   2186  CG  ASP A 317      12.583 -42.075 109.139  1.00 35.93           C
ATOM   2187  OD1 ASP A 317      13.435 -41.187 109.360  1.00 40.17           O
ATOM   2188  OD2 ASP A 317      12.831 -43.295 109.250  1.00 37.90           O
ATOM   2189  N   LEU A 318      10.550 -38.741 110.378  1.00 26.49           N
ATOM   2190  CA  LEU A 318      11.155 -37.560 110.985  1.00 26.24           C
ATOM   2191  C   LEU A 318      10.717 -37.537 112.448  1.00 26.95           C
ATOM   2192  O   LEU A 318       9.800 -36.810 112.828  1.00 27.99           O
ATOM   2193  CB  LEU A 318      10.693 -36.301 110.252  1.00 24.34           C
ATOM   2194  CG  LEU A 318      11.071 -36.288 108.768  1.00 23.72           C
ATOM   2195  CD1 LEU A 318      10.394 -35.132 108.065  1.00 24.49           C
ATOM   2196  CD2 LEU A 318      12.580 -36.193 108.629  1.00 24.02           C
ATOM   2197  N   THR A 319      11.392 -38.343 113.259  1.00 27.73           N
ATOM   2198  CA  THR A 319      11.076 -38.486 114.676  1.00 27.93           C
ATOM   2199  C   THR A 319      11.764 -37.518 115.637  1.00 28.67           C
ATOM   2200  O   THR A 319      11.364 -37.407 116.795  1.00 29.01           O
ATOM   2201  CB  THR A 319      11.379 -39.919 115.137  1.00 27.93           C
ATOM   2202  OG1 THR A 319      12.746 -40.235 114.841  1.00 27.70           O
ATOM   2203  CG2 THR A 319      10.472 -40.911 114.417  1.00 26.51           C
ATOM   2204  N   THR A 320      12.797 -36.826 115.169  1.00 29.46           N
ATOM   2205  CA  THR A 320      13.513 -35.872 116.011  1.00 29.47           C
ATOM   2206  C   THR A 320      12.739 -34.564 116.150  1.00 30.76           C
ATOM   2207  O   THR A 320      12.574 -34.041 117.251  1.00 31.05           O
ATOM   2208  CB  THR A 320      14.900 -35.571 115.435  1.00 29.84           C
ATOM   2209  OG1 THR A 320      15.722 -36.738 115.556  1.00 28.30           O
ATOM   2210  CG2 THR A 320      15.548 -34.398 116.168  1.00 30.39           C
ATOM   2211  N   SER A 321      12.274 -34.035 115.023  1.00 31.46           N
ATOM   2212  CA  SER A 321      11.505 -32.799 115.014  1.00 31.46           C
ATOM   2213  C   SER A 321      10.184 -33.067 114.316  1.00 30.56           C
ATOM   2214  O   SER A 321      10.104 -33.934 113.448  1.00 32.00           O
ATOM   2215  CB  SER A 321      12.270 -31.700 114.277  1.00 33.39           C
ATOM   2216  OG  SER A 321      13.478 -31.389 114.950  1.00 37.51           O
ATOM   2217  N   VAL A 322       9.149 -32.325 114.694  1.00 28.32           N
ATOM   2218  CA  VAL A 322       7.829 -32.503 114.103  1.00 25.27           C
ATOM   2219  C   VAL A 322       7.760 -31.974 112.671  1.00 24.50           C
ATOM   2220  O   VAL A 322       8.001 -30.794 112.424  1.00 24.31           O
ATOM   2221  CB  VAL A 322       6.756 -31.795 114.945  1.00 24.16           C
ATOM   2222  CG1 VAL A 322       5.380 -32.058 114.362  1.00 24.15           C
ATOM   2223  CG2 VAL A 322       6.830 -32.275 116.380  1.00 23.55           C
ATOM   2224  N   PRO A 323       7.429 -32.850 111.708  1.00 23.60           N
ATOM   2225  CA  PRO A 323       7.325 -32.477 110.294  1.00 22.62           C
ATOM   2226  C   PRO A 323       6.166 -31.512 110.044  1.00 22.04           C
ATOM   2227  O   PRO A 323       5.076 -31.683 110.594  1.00 21.04           O
ATOM   2228  CB  PRO A 323       7.097 -33.816 109.596  1.00 22.18           C
ATOM   2229  CG  PRO A 323       7.725 -34.803 110.523  1.00 23.36           C
ATOM   2230  CD  PRO A 323       7.282 -34.306 111.868  1.00 23.35           C
ATOM   2231  N   PRO A 324       6.389 -30.483 109.214  1.00 21.26           N
ATOM   2232  CA  PRO A 324       5.345 -29.504 108.905  1.00 20.28           C
ATOM   2233  C   PRO A 324       4.265 -30.149 108.052  1.00 20.15           C
ATOM   2234  O   PRO A 324       4.394 -31.297 107.620  1.00 21.43           O
ATOM   2235  CB  PRO A 324       6.093 -28.424 108.128  1.00 20.19           C
ATOM   2236  CG  PRO A 324       7.499 -28.565 108.593  1.00 22.05           C
ATOM   2237  CD  PRO A 324       7.676 -30.054 108.644  1.00 21.65           C
ATOM   2238  N   VAL A 325       3.207 -29.394 107.797  1.00 18.93           N
ATOM   2239  CA  VAL A 325       2.106 -29.878 106.986  1.00 17.74           C
ATOM   2240  C   VAL A 325       1.533 -28.648 106.265  1.00 18.27           C
ATOM   2241  O   VAL A 325       1.462 -27.564 106.845  1.00 19.21           O
ATOM   2242  CB  VAL A 325       1.065 -30.589 107.891  1.00 15.23           C
ATOM   2243  CG1 VAL A 325       0.310 -29.581 108.721  1.00 16.64           C
ATOM   2244  CG2 VAL A 325       0.144 -31.438 107.066  1.00 16.39           C
ATOM   2245  N   ALA A 326       1.160 -28.802 104.998  1.00 18.66           N
ATOM   2246  CA  ALA A 326       0.636 -27.680 104.211  1.00 19.12           C
ATOM   2247  C   ALA A 326      -0.783 -27.267 104.588  1.00 19.46           C
ATOM   2248  O   ALA A 326      -1.755 -27.780 104.035  1.00 22.10           O
ATOM   2249  CB  ALA A 326       0.699 -28.016 102.723  1.00 17.55           C
ATOM   2250  N   LEU A 327      -0.898 -26.324 105.516  1.00 17.63           N
ATOM   2251  CA  LEU A 327      -2.205 -25.855 105.961  1.00 16.12           C
ATOM   2252  C   LEU A 327      -2.471 -24.440 105.467  1.00 16.88           C
ATOM   2253  O   LEU A 327      -1.590 -23.795 104.901  1.00 18.47           O
ATOM   2254  CB  LEU A 327      -2.279 -25.897 107.487  1.00 13.74           C
ATOM   2255  CG  LEU A 327      -1.897 -27.240 108.115  1.00 11.17           C
ATOM   2256  CD1 LEU A 327      -1.930 -27.126 109.624  1.00  9.68           C
ATOM   2257  CD2 LEU A 327      -2.840 -28.330 107.633  1.00  9.69           C
ATOM   2258  N   THR A 328      -3.692 -23.962 105.683  1.00 16.51           N
ATOM   2259  CA  THR A 328      -4.076 -22.623 105.254  1.00 15.88           C
ATOM   2260  C   THR A 328      -3.155 -21.559 105.881  1.00 15.66           C
ATOM   2261  O   THR A 328      -3.011 -21.475 107.101  1.00 14.78           O
ATOM   2262  CB  THR A 328      -5.577 -22.379 105.580  1.00 14.98           C
ATOM   2263  OG1 THR A 328      -5.835 -20.975 105.690  1.00 15.31           O
ATOM   2264  CG2 THR A 328      -5.968 -23.098 106.862  1.00 16.50           C
ATOM   2265  N   ALA A 329      -2.535 -20.758 105.015  1.00 15.68           N
ATOM   2266  CA  ALA A 329      -1.570 -19.718 105.391  1.00 15.70           C
ATOM   2267  C   ALA A 329      -2.002 -18.573 106.305  1.00 16.13           C
ATOM   2268  O   ALA A 329      -3.165 -18.169 106.317  1.00 15.82           O
ATOM   2269  CB  ALA A 329      -0.949 -19.136 104.123  1.00 15.02           C
ATOM   2270  N   THR A 330      -1.028 -18.050 107.054  1.00 17.52           N
ATOM   2271  CA  THR A 330      -1.220 -16.930 107.984  1.00 18.76           C
ATOM   2272  C   THR A 330       0.039 -16.068 108.050  1.00 18.94           C
ATOM   2273  O   THR A 330       1.134 -16.520 107.713  1.00 18.98           O
ATOM   2274  CB  THR A 330      -1.486 -17.400 109.427  1.00 18.55           C
ATOM   2275  OG1 THR A 330      -2.582 -18.314 109.441  1.00 25.08           O
ATOM   2276  CG2 THR A 330      -1.826 -16.215 110.316  1.00 17.33           C
ATOM   2277  N   LYS A 331      -0.130 -14.826 108.494  1.00 18.76           N
ATOM   2278  CA  LYS A 331       0.981 -13.893 108.648  1.00 18.12           C
ATOM   2279  C   LYS A 331       1.276 -13.801 110.140  1.00 18.19           C
ATOM   2280  O   LYS A 331       2.396 -13.505 110.551  1.00 19.51           O
ATOM   2281  CB  LYS A 331       0.595 -12.500 108.149  1.00 17.57           C
ATOM   2282  CG  LYS A 331       0.218 -12.405 106.686  1.00 18.42           C
ATOM   2283  CD  LYS A 331      -0.306 -11.007 106.356  1.00 18.16           C
ATOM   2284  CE  LYS A 331       0.701  -9.919 106.725  1.00 16.46           C
ATOM   2285  NZ  LYS A 331       0.187  -8.552 106.431  1.00 15.24           N
ATOM   2286  N   LEU A 332       0.247 -14.064 110.937  1.00 18.16           N
ATOM   2287  CA  LEU A 332       0.315 -14.001 112.393  1.00 18.94           C
ATOM   2288  C   LEU A 332       1.112 -15.152 113.006  1.00 20.72           C
ATOM   2289  O   LEU A 332       1.034 -15.408 114.207  1.00 21.12           O
ATOM   2290  CB  LEU A 332      -1.111 -13.996 112.944  1.00 17.77           C
ATOM   2291  CG  LEU A 332      -2.055 -13.055 112.181  1.00 15.45           C
ATOM   2292  CD1 LEU A 332      -3.501 -13.375 112.509  1.00 14.04           C
ATOM   2293  CD2 LEU A 332      -1.723 -11.617 112.522  1.00 12.75           C
ATOM   2294  N   ASN A 333       1.880 -15.838 112.166  1.00 23.44           N
ATOM   2295  CA  ASN A 333       2.703 -16.971 112.578  1.00 23.93           C
ATOM   2296  C   ASN A 333       4.150 -16.545 112.728  1.00 24.70           C
ATOM   2297  O   ASN A 333       4.853 -16.992 113.632  1.00 24.54           O
ATOM   2298  CB  ASN A 333       2.646 -18.062 111.515  1.00 26.61           C
ATOM   2299  CG  ASN A 333       1.858 -19.250 111.956  1.00 29.55           C
ATOM   2300  OD1 ASN A 333       2.146 -19.843 112.991  1.00 33.22           O
ATOM   2301  ND2 ASN A 333       0.853 -19.617 111.173  1.00 34.06           N
ATOM   2302  N   GLN A 334       4.586 -15.696 111.803  1.00 24.83           N
ATOM   2303  CA  GLN A 334       5.948 -15.190 111.768  1.00 24.33           C
ATOM   2304  C   GLN A 334       6.154 -14.072 112.782  1.00 23.29           C
ATOM   2305  O   GLN A 334       5.292 -13.209 112.949  1.00 23.17           O
ATOM   2306  CB  GLN A 334       6.257 -14.671 110.362  1.00 26.56           C
ATOM   2307  CG  GLN A 334       7.356 -15.429 109.623  1.00 29.28           C
ATOM   2308  CD  GLN A 334       7.077 -16.915 109.515  1.00 28.54           C
ATOM   2309  OE1 GLN A 334       6.021 -17.330 109.038  1.00 28.36           O
ATOM   2310  NE2 GLN A 334       8.031 -17.725 109.954  1.00 29.34           N
ATOM   2311  N   SER A 335       7.298 -14.094 113.459  1.00 21.22           N
ATOM   2312  CA  SER A 335       7.616 -13.065 114.439  1.00 19.84           C
ATOM   2313  C   SER A 335       8.243 -11.891 113.700  1.00 19.90           C
ATOM   2314  O   SER A 335       8.584 -12.005 112.522  1.00 19.68           O
ATOM   2315  CB  SER A 335       8.600 -13.597 115.479  1.00 19.73           C
ATOM   2316  OG  SER A 335       9.870 -13.836 114.897  1.00 20.23           O
ATOM   2317  N   ALA A 336       8.401 -10.767 114.393  1.00 19.82           N
ATOM   2318  CA  ALA A 336       8.988  -9.579 113.787  1.00 18.90           C
ATOM   2319  C   ALA A 336      10.408  -9.856 113.309  1.00 19.44           C
ATOM   2320  O   ALA A 336      10.807  -9.413 112.232  1.00 18.92           O
ATOM   2321  CB  ALA A 336       8.987  -8.436 114.780  1.00 18.18           C
ATOM   2322  N   SER A 337      11.169 -10.592 114.112  1.00 20.65           N
ATOM   2323  CA  SER A 337      12.543 -10.918 113.756  1.00 21.77           C
ATOM   2324  C   SER A 337      12.601 -11.896 112.588  1.00 21.99           C
ATOM   2325  O   SER A 337      13.536 -11.852 111.790  1.00 22.49           O
ATOM   2326  CB  SER A 337      13.280 -11.495 114.963  1.00 21.52           C
ATOM   2327  OG  SER A 337      12.584 -12.609 115.487  1.00 27.59           O
ATOM   2328  N   ASN A 338      11.610 -12.779 112.484  1.00 22.03           N
ATOM   2329  CA  ASN A 338      11.584 -13.729 111.378  1.00 22.50           C
ATOM   2330  C   ASN A 338      11.427 -12.965 110.068  1.00 22.06           C
ATOM   2331  O   ASN A 338      12.155 -13.213 109.108  1.00 21.25           O
ATOM   2332  CB  ASN A 338      10.430 -14.727 111.517  1.00 25.28           C
ATOM   2333  CG  ASN A 338      10.682 -15.778 112.586  1.00 27.64           C
ATOM   2334  OD1 ASN A 338      11.812 -16.237 112.773  1.00 28.04           O
ATOM   2335  ND2 ASN A 338       9.620 -16.181 113.280  1.00 28.11           N
ATOM   2336  N   ASN A 339      10.477 -12.032 110.038  1.00 21.49           N
ATOM   2337  CA  ASN A 339      10.232 -11.229 108.844  1.00 20.40           C
ATOM   2338  C   ASN A 339      11.511 -10.562 108.357  1.00 19.57           C
ATOM   2339  O   ASN A 339      11.764 -10.496 107.156  1.00 20.36           O
ATOM   2340  CB  ASN A 339       9.186 -10.147 109.115  1.00 20.32           C
ATOM   2341  CG  ASN A 339       7.869 -10.714 109.585  1.00 20.71           C
ATOM   2342  OD1 ASN A 339       7.471 -11.806 109.178  1.00 22.61           O
ATOM   2343  ND2 ASN A 339       7.174  -9.969 110.435  1.00 19.61           N
ATOM   2344  N   LEU A 340      12.312 -10.065 109.294  1.00 17.69           N
ATOM   2345  CA  LEU A 340      13.561  -9.405 108.942  1.00 16.22           C
ATOM   2346  C   LEU A 340      14.507 -10.326 108.185  1.00 16.58           C
ATOM   2347  O   LEU A 340      15.286  -9.866 107.349  1.00 18.23           O
ATOM   2348  CB  LEU A 340      14.251  -8.872 110.198  1.00 13.20           C
ATOM   2349  CG  LEU A 340      13.589  -7.648 110.832  1.00 10.55           C
ATOM   2350  CD1 LEU A 340      14.276  -7.322 112.136  1.00  9.92           C
ATOM   2351  CD2 LEU A 340      13.663  -6.465 109.877  1.00  9.35           C
ATOM   2352  N   ASN A 341      14.435 -11.624 108.468  1.00 16.36           N
ATOM   2353  CA  ASN A 341      15.295 -12.597 107.798  1.00 15.90           C
ATOM   2354  C   ASN A 341      14.594 -13.269 106.626  1.00 15.59           C
ATOM   2355  O   ASN A 341      14.917 -14.403 106.269  1.00 15.34           O
ATOM   2356  CB  ASN A 341      15.754 -13.670 108.782  1.00 15.76           C
ATOM   2357  CG  ASN A 341      16.460 -13.089 109.977  1.00 16.78           C
ATOM   2358  OD1 ASN A 341      17.302 -12.200 109.841  1.00 19.51           O
ATOM   2359  ND2 ASN A 341      16.132 -13.593 111.161  1.00 16.09           N
ATOM   2360  N   ALA A 342      13.640 -12.566 106.027  1.00 14.50           N
ATOM   2361  CA  ALA A 342      12.894 -13.105 104.900  1.00 13.77           C
ATOM   2362  C   ALA A 342      13.788 -13.386 103.695  1.00 14.02           C
ATOM   2363  O   ALA A 342      14.703 -12.618 103.392  1.00 14.49           O
ATOM   2364  CB  ALA A 342      11.782 -12.143 104.511  1.00 12.29           C
"""

def exercise_single_mtrix():
  inp = iotbx.pdb.input(lines=single_mtrix_txt+ss_txt+atoms_txt, source_info=None)
  model = mmtbx.model.manager(
    model_input = inp)
  # print (model.model_as_pdb())
  assert model.get_number_of_atoms() == 600, model.get_number_of_atoms()
  assert model.get_hierarchy().atoms_size() == 600
  assert model.get_xray_structure().scatterers().size() == 600
  ss_ann = model.get_ss_annotation()
  # print ss_ann.as_pdb_str()
  assert ss_ann.get_n_helices() == 4
  assert ss_ann.get_n_sheets() == 2

def exercise_mtrix():
  inp = iotbx.pdb.input(lines=mtrix_txt+ss_txt+atoms_txt, source_info=None)
  model = mmtbx.model.manager(
    model_input = inp)
  assert model.get_number_of_atoms() == 900, model.get_number_of_atoms()
  assert model.get_hierarchy().atoms_size() == 900
  assert model.get_xray_structure().scatterers().size() == 900
  ss_ann = model.get_ss_annotation()
  # print ss_ann.as_pdb_str()
  assert ss_ann.get_n_helices() == 6
  assert ss_ann.get_n_sheets() == 3


def exercise_biomt():
  inp = iotbx.pdb.input(lines=biomt_txt+ss_txt+atoms_txt, source_info=None)
  model = mmtbx.model.manager(
    model_input = inp)
  assert model.get_number_of_atoms() == 300, model.get_number_of_atoms()
  assert model.get_hierarchy().atoms_size() == 300
  assert model.get_xray_structure().scatterers().size() == 300
  ss_ann = model.get_ss_annotation()
  assert ss_ann.get_n_helices() == 2
  assert ss_ann.get_n_sheets() == 1
  model.expand_with_BIOMT_records()
  assert model.get_number_of_atoms() == 900, model.get_number_of_atoms()
  assert model.get_hierarchy().atoms_size() == 900
  assert model.get_xray_structure().scatterers().size() == 900, model.get_xray_structure().scatterers().size()
  ss_ann = model.get_ss_annotation()
  assert ss_ann.get_n_helices() == 6
  assert ss_ann.get_n_sheets() == 3

def exercise_both():
  inp = iotbx.pdb.input(lines=mtrix_txt+biomt_txt+ss_txt+atoms_txt, source_info=None)
  model = mmtbx.model.manager(
    model_input = inp)
  assert model.get_number_of_atoms() == 900, model.get_number_of_atoms()
  assert model.get_hierarchy().atoms_size() == 900
  assert model.get_xray_structure().scatterers().size() == 900
  ss_ann = model.get_ss_annotation()
  # print ss_ann.as_pdb_str()
  # print "="*30
  assert ss_ann.get_n_helices() == 6
  assert ss_ann.get_n_sheets() == 3
  model.expand_with_BIOMT_records()
  assert model.get_number_of_atoms() == 2700, model.get_number_of_atoms()
  assert model.get_hierarchy().atoms_size() == 2700
  assert model.get_xray_structure().scatterers().size() == 2700, model.get_xray_structure().scatterers().size()
  ss_ann = model.get_ss_annotation()
  # print ss_ann.as_pdb_str()
  assert ss_ann.get_n_helices() == 18
  assert ss_ann.get_n_sheets() == 9
  return


if (__name__ == "__main__"):
  t0 = time.time()
  exercise_single_mtrix()
  exercise_mtrix()
  exercise_biomt()
  exercise_both()
  print("Total time: %8.3f"%(time.time() - t0))
  print("OK.")
