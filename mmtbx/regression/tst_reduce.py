##################################################################################
# This is a test program to validate that the Python wrapping of Reduce worked.
#

#                Copyright 2021-2023  Richardson Lab at Duke University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, nested_scopes, generators, division
from __future__ import absolute_import
from mmtbx.reduce import Movers
from mmtbx.reduce import InteractionGraph
from mmtbx.reduce import Optimizers
from mmtbx.programs import reduce2
from iotbx.cli_parser import run_program
from six.moves import cStringIO as StringIO
import os.path as op
import os
import math

def RunReduceTests():

  #========================================================================
  # Call the test functions for all of our files.

  print('Testing Movers objects')
  ret = Movers.Test()
  if len(ret) != 0:
    print("Failure:",ret)
  assert (len(ret) == 0)

  print('Testing InteractionGraph objects')
  ret = InteractionGraph.Test()
  if len(ret) != 0:
    print("Failure:",ret)
  assert (len(ret) == 0)

  print('Testing Optimizers')
  ret = Optimizers.Test()
  if len(ret) != 0:
    print("Failure:",ret)
  assert (len(ret) == 0), 'len(ret)=%s' % len(ret)
  return ret

# Each test case has a name, a raw PDB file, a chain ID, a residue ID, an atom name,
# a list of possible locations for the atom (for sets of 3 hydrogens, any 60-degree
# rotation is equivalent in score), a maximum distance threshold, and a list of extra
# arguments to make to Reduce2 for this case.
testCases = [

  ["7c31_single_hydrogen_rotator",
   """\
CRYST1   27.854   27.854   99.605  90.00  90.00  90.00 P 43          8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.035901  0.000000  0.000000        0.00000
SCALE2      0.000000  0.035901  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010040        0.00000
ATOM     68  N   SER A   5     -31.155  49.742   0.887  1.00 10.02           N
ATOM     69  CA  SER A   5     -32.274  48.937   0.401  1.00  9.76           C
ATOM     70  C   SER A   5     -33.140  49.851  -0.454  1.00  9.47           C
ATOM     71  O   SER A   5     -33.502  50.939  -0.012  1.00 10.76           O
ATOM     72  CB  SER A   5     -33.086  48.441   1.599  1.00 12.34           C
ATOM     73  OG  SER A   5     -34.118  47.569   1.179  1.00 19.50           O
ATOM    758  N   VAL B   2     -34.430  42.959   3.043  1.00 19.95           N
ATOM    759  CA  VAL B   2     -32.994  42.754   2.904  0.50 18.09           C
ATOM    760  C   VAL B   2     -32.311  44.083   2.629  1.00 17.65           C
ATOM    761  O   VAL B   2     -32.869  44.976   1.975  1.00 18.97           O
ATOM    762  CB  VAL B   2     -32.703  41.738   1.778  0.50 19.70           C
ATOM    763  CG1 VAL B   2     -33.098  42.307   0.419  0.50 19.61           C
ATOM    764  CG2 VAL B   2     -31.224  41.334   1.755  0.50 16.49           C
TER    1447      CYS B  47
END
""",
   "A",
   5,
   "HG",
   [ (-33.906, 46.773, 1.342)
   ],
   0.1,
   []
  ],

  ["1xso_amide_unflipped",
   """\
CRYST1   73.450   68.940   58.760  90.00  90.00  90.00 P 21 21 21    8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.013615  0.000000  0.000000        0.00000
SCALE2      0.000000  0.014505  0.000000        0.00000
SCALE3      0.000000  0.000000  0.017018        0.00000
MTRIX1   1  0.987830 -0.010290 -0.155220        1.43777    1
MTRIX2   1 -0.008360 -0.999880  0.013070       33.15571    1
MTRIX3   1 -0.155330 -0.011610 -0.987790       17.07902    1
ATOM    147  N   GLN A  22      32.268  17.697  -8.014  1.00 15.32           N
ATOM    148  CA  GLN A  22      33.617  18.261  -8.189  1.00 14.46           C
ATOM    149  C   GLN A  22      34.546  17.130  -8.669  1.00 18.34           C
ATOM    150  O   GLN A  22      34.525  16.056  -8.055  1.00 20.98           O
ATOM    151  CB  GLN A  22      34.120  18.672  -6.771  1.00 15.81           C
ATOM    152  CG  GLN A  22      35.452  19.443  -6.847  1.00 17.90           C
ATOM    153  CD  GLN A  22      35.878  19.981  -5.501  1.00 14.30           C
ATOM    154  OE1 GLN A  22      35.529  19.359  -4.472  1.00 18.17           O
ATOM    155  NE2 GLN A  22      36.611  21.049  -5.457  1.00 20.19           N
ATOM    173  N   GLU A  25      40.443  19.021  -7.568  1.00 47.81           N
ATOM    174  CA  GLU A  25      40.160  20.153  -6.631  1.00 42.14           C
ATOM    175  C   GLU A  25      39.675  21.334  -7.469  1.00 36.95           C
ATOM    176  O   GLU A  25      39.955  22.496  -7.223  1.00 45.86           O
ATOM    177  CB  GLU A  25      41.460  20.522  -5.912  0.00 41.53           C
ATOM    178  CG  GLU A  25      42.518  19.432  -5.947  0.00 41.54           C
ATOM    179  CD  GLU A  25      43.919  19.983  -6.125  0.00 41.41           C
ATOM    180  OE1 GLU A  25      44.311  20.257  -7.279  0.00 41.37           O
ATOM    181  OE2 GLU A  25      44.629  20.145  -5.110  0.00 41.34           O
ATOM    182  N   GLY A  25A     38.940  21.006  -8.547  1.00 35.16           N
ATOM    183  CA  GLY A  25A     38.466  22.062  -9.438  1.00 30.84           C
ATOM    184  C   GLY A  25A     37.135  22.649  -9.008  1.00 26.52           C
ATOM    185  O   GLY A  25A     36.678  22.555  -7.874  1.00 24.92           O
ATOM    191  N   VAL A  27      33.241  23.494  -8.492  1.00 13.48           N
ATOM    192  CA  VAL A  27      32.081  22.686  -8.152  1.00 14.56           C
ATOM    193  C   VAL A  27      30.869  23.108  -9.008  1.00 14.24           C
ATOM    194  O   VAL A  27      30.552  24.277  -9.114  1.00 15.24           O
ATOM    195  CB  VAL A  27      31.716  22.923  -6.662  1.00 13.92           C
ATOM    196  CG1 VAL A  27      30.493  22.048  -6.301  1.00 19.33           C
ATOM    197  CG2 VAL A  27      32.934  22.549  -5.797  1.00 15.63           C
ATOM    749  N   SER A 103      34.060  24.070  -2.055  1.00 11.00           N
ATOM    750  CA  SER A 103      34.714  22.821  -2.367  1.00 11.78           C
ATOM    751  C   SER A 103      33.932  21.688  -1.675  1.00 11.18           C
ATOM    752  O   SER A 103      33.182  21.910  -0.755  1.00 10.91           O
ATOM    753  CB  SER A 103      36.150  22.787  -1.786  1.00 14.12           C
ATOM    754  OG  SER A 103      36.858  21.681  -2.392  1.00 15.49           O
ATOM    755  N   LEU A 104      34.225  20.450  -2.085  1.00 12.60           N
ATOM    756  CA  LEU A 104      33.727  19.278  -1.440  1.00 13.49           C
ATOM    757  C   LEU A 104      34.844  18.562  -0.681  1.00 14.51           C
ATOM    758  O   LEU A 104      34.606  17.475  -0.159  1.00 16.82           O
ATOM    759  CB  LEU A 104      33.019  18.367  -2.402  1.00 13.27           C
ATOM    760  CG  LEU A 104      31.817  18.935  -3.135  1.00 13.24           C
ATOM    761  CD1 LEU A 104      31.200  17.904  -4.040  1.00 16.51           C
ATOM    762  CD2 LEU A 104      30.821  19.574  -2.182  1.00 22.63           C
ATOM    763  N   LYS A 105      36.033  19.165  -0.594  1.00 14.92           N
ATOM    764  CA  LYS A 105      37.139  18.617   0.145  1.00 16.30           C
ATOM    765  C   LYS A 105      38.020  19.765   0.678  1.00 17.78           C
ATOM    766  O   LYS A 105      37.916  20.873   0.200  1.00 20.07           O
ATOM    767  CB  LYS A 105      37.943  17.627  -0.649  1.00 26.20           C
ATOM    768  CG  LYS A 105      38.587  18.132  -1.914  1.00 38.16           C
ATOM    769  CD  LYS A 105      38.295  17.264  -3.120  0.00 45.15           C
ATOM    770  CE  LYS A 105      38.612  15.811  -2.887  1.00 51.58           C
ATOM    771  NZ  LYS A 105      39.344  15.183  -4.045  1.00 74.91           N
TER    1093      PRO A 151
END
""",
   "A",
   22,
   "OE1",
   [ (35.529, 19.359, -4.472)
   ],
   0.1,
   []
  ],

  ["1xso_amide_flipped",
   """\
CRYST1   73.450   68.940   58.760  90.00  90.00  90.00 P 21 21 21    8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.013615  0.000000  0.000000        0.00000
SCALE2      0.000000  0.014505  0.000000        0.00000
SCALE3      0.000000  0.000000  0.017018        0.00000
MTRIX1   1  0.987830 -0.010290 -0.155220        1.43777    1
MTRIX2   1 -0.008360 -0.999880  0.013070       33.15571    1
MTRIX3   1 -0.155330 -0.011610 -0.987790       17.07902    1
ATOM    924  N   GLY A 127      23.876  47.136  -0.439  1.00 18.38           N
ATOM    925  CA  GLY A 127      22.692  47.748   0.105  1.00 17.97           C
ATOM    926  C   GLY A 127      22.873  49.006   0.862  1.00 18.02           C
ATOM    927  O   GLY A 127      21.908  49.714   1.225  1.00 24.20           O
ATOM    932  N   ASN A 129      24.271  49.439   3.964  1.00 14.98           N
ATOM    933  CA  ASN A 129      24.135  49.362   5.428  1.00 15.72           C
ATOM    934  C   ASN A 129      24.971  48.198   5.914  1.00 13.04           C
ATOM    935  O   ASN A 129      25.459  47.398   5.129  1.00 14.62           O
ATOM    936  CB  ASN A 129      22.691  49.260   5.878  1.00 15.14           C
ATOM    937  CG  ASN A 129      21.911  48.176   5.160  1.00 19.11           C
ATOM    938  OD1 ASN A 129      22.371  47.040   5.114  1.00 15.87           O
ATOM    939  ND2 ASN A 129      20.811  48.560   4.517  1.00 30.20           N
ATOM    948  N   GLU A 131      24.275  45.568   7.628  1.00 14.32           N
ATOM    949  CA  GLU A 131      23.665  44.242   7.404  1.00 15.35           C
ATOM    950  C   GLU A 131      23.963  43.701   6.022  1.00 13.39           C
ATOM    951  O   GLU A 131      24.150  42.499   5.830  1.00 13.14           O
ATOM    952  CB  GLU A 131      22.201  44.278   7.753  1.00 18.42           C
ATOM    953  CG  GLU A 131      21.474  42.964   7.698  1.00 23.30           C
ATOM    954  CD  GLU A 131      21.905  42.014   8.811  1.00 27.07           C
ATOM    955  OE1 GLU A 131      22.005  42.457   9.969  1.00 37.65           O
ATOM    956  OE2 GLU A 131      22.106  40.834   8.513  1.00 32.51           O
ATOM    957  N   SER A 132      24.016  44.584   5.024  1.00 13.86           N
ATOM    958  CA  SER A 132      24.353  44.193   3.664  1.00 12.64           C
ATOM    959  C   SER A 132      25.696  43.468   3.641  1.00 12.67           C
ATOM    960  O   SER A 132      25.800  42.483   2.891  1.00 13.76           O
ATOM    961  CB  SER A 132      24.320  45.391   2.737  1.00 14.43           C
ATOM    962  OG  SER A 132      24.685  45.044   1.405  1.00 13.58           O
ATOM    991  N   ASN A 137      20.989  41.281   3.329  1.00 12.56           N
ATOM    992  CA  ASN A 137      19.607  41.651   3.597  1.00 13.15           C
ATOM    993  C   ASN A 137      18.639  40.676   2.991  1.00 10.75           C
ATOM    994  O   ASN A 137      17.511  41.040   2.611  1.00 15.89           O
ATOM    995  CB  ASN A 137      19.280  43.097   3.159  1.00 13.87           C
ATOM    996  CG  ASN A 137      19.949  44.125   4.040  1.00 16.63           C
ATOM    997  OD1 ASN A 137      19.570  44.245   5.254  1.00 19.77           O
ATOM    998  ND2 ASN A 137      20.916  44.887   3.583  1.00 16.13           N
TER    1093      PRO A 151
END
""",
   "A",
   129,
   "OD1",
   [ (22.371, 47.040, 5.114)
   ],
   0.1,
   []
  ],

  ["1xso_histidine_unflipped",
   """\
CRYST1   73.450   68.940   58.760  90.00  90.00  90.00 P 21 21 21    8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.013615  0.000000  0.000000        0.00000
SCALE2      0.000000  0.014505  0.000000        0.00000
SCALE3      0.000000  0.000000  0.017018        0.00000
MTRIX1   1  0.987830 -0.010290 -0.155220        1.43777    1
MTRIX2   1 -0.008360 -0.999880  0.013070       33.15571    1
MTRIX3   1 -0.155330 -0.011610 -0.987790       17.07902    1
ATOM    103  N   VAL A  17      15.777  20.067  -3.917  1.00 11.00           N
ATOM    104  CA  VAL A  17      16.936  19.211  -3.933  1.00 12.41           C
ATOM    105  C   VAL A  17      18.055  19.934  -4.702  1.00  9.37           C
ATOM    106  O   VAL A  17      17.836  20.513  -5.742  1.00 13.60           O
ATOM    107  CB  VAL A  17      16.635  17.858  -4.617  1.00 16.16           C
ATOM    108  CG1 VAL A  17      17.933  16.983  -4.557  1.00 21.33           C
ATOM    109  CG2 VAL A  17      15.541  17.118  -3.839  1.00 19.65           C
ATOM    117  N   HIS A  19      22.094  19.284  -6.044  1.00 11.05           N
ATOM    118  CA  HIS A  19      23.225  18.356  -6.269  1.00 10.79           C
ATOM    119  C   HIS A  19      24.527  19.153  -6.473  1.00 11.77           C
ATOM    120  O   HIS A  19      24.547  20.205  -7.093  1.00 13.93           O
ATOM    121  CB  HIS A  19      22.906  17.530  -7.527  1.00 18.16           C
ATOM    122  CG  HIS A  19      21.858  16.494  -7.322  1.00 26.17           C
ATOM    123  ND1 HIS A  19      22.149  15.200  -6.888  1.00 28.19           N
ATOM    124  CD2 HIS A  19      20.525  16.560  -7.452  1.00 32.51           C
ATOM    125  CE1 HIS A  19      21.056  14.541  -6.758  1.00 32.39           C
ATOM    126  NE2 HIS A  19      20.028  15.310  -7.102  1.00 34.60           N
HETATM 2433  O   HOH B 215      21.490  11.727  -5.002  1.00 34.66           O
HETATM 2436  O   HOH B 218      16.653  14.508   9.188  1.00 35.04           O
HETATM 2498  O   HOH B 280       8.420  -1.071  13.530  1.00 50.43           O
HETATM 2525  O   HOH B 307      37.167  -8.451   4.794  1.00 59.87           O
HETATM 2538  O   HOH B 320      14.429  18.062  16.293  1.00 67.01           O
TER    1093      PRO A 151
END
""",
   "A",
   19,
   "ND1",
   [ (22.149, 15.200, -6.888)
   ],
   0.1,
   []
  ],

  ["1xso_histidine_flipped",
   """\
CRYST1   73.450   68.940   58.760  90.00  90.00  90.00 P 21 21 21
SCALE1      0.013615  0.000000  0.000000        0.00000
SCALE2      0.000000  0.014505  0.000000        0.00000
SCALE3      0.000000  0.000000  0.017018        0.00000
ATOM      1  N   VAL A  17      15.777  20.067  -3.917  1.00 11.00           N
ATOM      2  CA  VAL A  17      16.936  19.211  -3.933  1.00 12.41           C
ATOM      3  C   VAL A  17      18.055  19.934  -4.702  1.00  9.37           C
ATOM      4  O   VAL A  17      17.836  20.513  -5.742  1.00 13.60           O
ATOM      5  CB  VAL A  17      16.635  17.858  -4.617  1.00 16.16           C
ATOM      6  CG1 VAL A  17      17.933  16.983  -4.557  1.00 21.33           C
ATOM      7  CG2 VAL A  17      15.541  17.118  -3.839  1.00 19.65           C
ATOM      8  H   VAL A  17      15.824  20.729  -4.464  1.00 11.00           H
ATOM      9  HA  VAL A  17      17.211  19.025  -3.022  1.00 12.41           H
ATOM     10  HB  VAL A  17      16.348  18.006  -5.532  1.00 16.16           H
ATOM     11 HG11 VAL A  17      18.178  16.845  -3.629  1.00 21.33           H
ATOM     12 HG12 VAL A  17      18.647  17.445  -5.023  1.00 21.33           H
ATOM     13 HG13 VAL A  17      17.758  16.129  -4.983  1.00 21.33           H
ATOM     14 HG21 VAL A  17      14.738  17.662  -3.828  1.00 19.65           H
ATOM     15 HG22 VAL A  17      15.848  16.963  -2.932  1.00 19.65           H
ATOM     16 HG23 VAL A  17      15.360  16.271  -4.276  1.00 19.65           H
ATOM     17  N   HIS A  19      22.094  19.284  -6.044  1.00 11.05           N
ATOM     18  CA  HIS A  19      23.225  18.356  -6.269  1.00 10.79           C
ATOM     19  C   HIS A  19      24.527  19.153  -6.473  1.00 11.77           C
ATOM     20  O   HIS A  19      24.547  20.205  -7.093  1.00 13.93           O
ATOM     21  CB  HIS A  19      22.920  17.515  -7.520  1.00 18.16           C
ATOM     22  CG  HIS A  19      22.000  16.372  -7.273  1.00 26.17           C
ATOM     23  ND1 HIS A  19      20.613  16.482  -7.381  1.00 28.19           N
ATOM     24  CD2 HIS A  19      22.252  15.112  -6.890  1.00 32.51           C
ATOM     25  CE1 HIS A  19      20.075  15.354  -7.090  1.00 32.39           C
ATOM     26  NE2 HIS A  19      21.021  14.479  -6.766  1.00 34.60           N
ATOM     27  H   HIS A  19      22.262  20.099  -6.260  1.00 11.05           H
ATOM     28  HA  HIS A  19      23.356  17.763  -5.513  1.00 10.79           H
ATOM     29  HB2 HIS A  19      23.752  17.154  -7.863  1.00 18.16           H
ATOM     30  HB3 HIS A  19      22.507  18.087  -8.186  1.00 18.16           H
ATOM     31  HD2 HIS A  19      23.088  14.735  -6.737  1.00 32.51           H
ATOM     32  HE1 HIS A  19      19.162  15.177  -7.104  1.00 32.39           H
ATOM     33  HE2 HIS A  19      20.894  13.664  -6.522  1.00 34.60           H
TER
HETATM   34  O   HOH B 215      21.490  11.727  -5.002  1.00 34.66           O
HETATM   35  O   HOH B 218      16.653  14.508   9.188  1.00 35.04           O
HETATM   36  O   HOH B 280       8.420  -1.071  13.530  1.00 50.43           O
HETATM   37  O   HOH B 307      37.167  -8.451   4.794  1.00 59.87           O
HETATM   38  O   HOH B 320      14.429  18.062  16.293  1.00 67.01           O
END
""",
   "A",
   19,
   "ND1",
   [ (22.149, 15.200, -6.888)
   ],
   0.1,
   []
  ],

  ["1xso_histidine_uncertain",
   """\
CRYST1   73.450   68.940   58.760  90.00  90.00  90.00 P 21 21 21    8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.013615  0.000000  0.000000        0.00000
SCALE2      0.000000  0.014505  0.000000        0.00000
SCALE3      0.000000  0.000000  0.017018        0.00000
MTRIX1   1  0.987830 -0.010290 -0.155220        1.43777    1
MTRIX2   1 -0.008360 -0.999880  0.013070       33.15571    1
MTRIX3   1 -0.155330 -0.011610 -0.987790       17.07902    1
ATOM   1115  N   VAL B   5      26.079  16.885  14.997  1.00  8.83           N
ATOM   1116  CA  VAL B   5      24.686  16.588  15.374  1.00  8.78           C
ATOM   1117  C   VAL B   5      23.969  16.005  14.172  1.00  9.25           C
ATOM   1118  O   VAL B   5      24.335  16.257  13.044  1.00 10.86           O
ATOM   1119  CB  VAL B   5      23.982  17.879  15.867  1.00 10.17           C
ATOM   1120  CG1 VAL B   5      23.763  18.899  14.740  1.00 14.21           C
ATOM   1121  CG2 VAL B   5      22.685  17.620  16.633  1.00 14.81           C
ATOM   1196  N   VAL B  17      17.298  13.044  18.182  1.00 10.96           N
ATOM   1197  CA  VAL B  17      18.457  13.879  18.007  1.00  9.04           C
ATOM   1198  C   VAL B  17      19.658  13.118  18.600  1.00 11.44           C
ATOM   1199  O   VAL B  17      19.557  12.560  19.692  1.00 14.98           O
ATOM   1200  CB  VAL B  17      18.288  15.251  18.678  1.00 13.50           C
ATOM   1201  CG1 VAL B  17      19.502  16.123  18.558  1.00 17.18           C
ATOM   1202  CG2 VAL B  17      17.090  16.002  18.108  1.00 20.25           C
ATOM   1210  N   HIS B  19      23.916  13.702  19.282  1.00  9.93           N
ATOM   1211  CA  HIS B  19      25.062  14.561  19.425  1.00  9.29           C
ATOM   1212  C   HIS B  19      26.366  13.749  19.441  1.00  8.50           C
ATOM   1213  O   HIS B  19      26.409  12.668  20.031  1.00 12.76           O
ATOM   1214  CB  HIS B  19      24.996  15.427  20.688  1.00 12.81           C
ATOM   1215  CG  HIS B  19      23.810  16.312  20.701  1.00 19.55           C
ATOM   1216  ND1 HIS B  19      23.897  17.622  20.251  1.00 24.11           N
ATOM   1217  CD2 HIS B  19      22.556  16.137  21.208  1.00 25.44           C
ATOM   1218  CE1 HIS B  19      22.723  18.181  20.377  1.00 27.18           C
ATOM   1219  NE2 HIS B  19      21.897  17.308  20.959  1.00 29.89           N
HETATM 2413  O   HOH B 195      10.371  -1.926  15.288  1.00 30.69           O
HETATM 2486  O   HOH B 268      10.982   5.736  29.100  1.00 46.97           O
HETATM 2528  O   HOH B 310      19.641  17.849  22.382  1.00 61.47           O
TER    1093      PRO A 151
END
""",
   "B",
   19,
   "ND1",
   [ (23.897, 17.622, 20.251)
   ],
   0.1,
   []
  ],

  ["1xso_histidine_large_unflip_preference",
   """\
CRYST1   73.450   68.940   58.760  90.00  90.00  90.00 P 21 21 21    8
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.013615  0.000000  0.000000        0.00000
SCALE2      0.000000  0.014505  0.000000        0.00000
SCALE3      0.000000  0.000000  0.017018        0.00000
MTRIX1   1  0.987830 -0.010290 -0.155220        1.43777    1
MTRIX2   1 -0.008360 -0.999880  0.013070       33.15571    1
MTRIX3   1 -0.155330 -0.011610 -0.987790       17.07902    1
ATOM   1115  N   VAL B   5      26.079  16.885  14.997  1.00  8.83           N
ATOM   1116  CA  VAL B   5      24.686  16.588  15.374  1.00  8.78           C
ATOM   1117  C   VAL B   5      23.969  16.005  14.172  1.00  9.25           C
ATOM   1118  O   VAL B   5      24.335  16.257  13.044  1.00 10.86           O
ATOM   1119  CB  VAL B   5      23.982  17.879  15.867  1.00 10.17           C
ATOM   1120  CG1 VAL B   5      23.763  18.899  14.740  1.00 14.21           C
ATOM   1121  CG2 VAL B   5      22.685  17.620  16.633  1.00 14.81           C
ATOM   1196  N   VAL B  17      17.298  13.044  18.182  1.00 10.96           N
ATOM   1197  CA  VAL B  17      18.457  13.879  18.007  1.00  9.04           C
ATOM   1198  C   VAL B  17      19.658  13.118  18.600  1.00 11.44           C
ATOM   1199  O   VAL B  17      19.557  12.560  19.692  1.00 14.98           O
ATOM   1200  CB  VAL B  17      18.288  15.251  18.678  1.00 13.50           C
ATOM   1201  CG1 VAL B  17      19.502  16.123  18.558  1.00 17.18           C
ATOM   1202  CG2 VAL B  17      17.090  16.002  18.108  1.00 20.25           C
ATOM   1210  N   HIS B  19      23.916  13.702  19.282  1.00  9.93           N
ATOM   1211  CA  HIS B  19      25.062  14.561  19.425  1.00  9.29           C
ATOM   1212  C   HIS B  19      26.366  13.749  19.441  1.00  8.50           C
ATOM   1213  O   HIS B  19      26.409  12.668  20.031  1.00 12.76           O
ATOM   1214  CB  HIS B  19      24.996  15.427  20.688  1.00 12.81           C
ATOM   1215  CG  HIS B  19      23.810  16.312  20.701  1.00 19.55           C
ATOM   1216  ND1 HIS B  19      23.897  17.622  20.251  1.00 24.11           N
ATOM   1217  CD2 HIS B  19      22.556  16.137  21.208  1.00 25.44           C
ATOM   1218  CE1 HIS B  19      22.723  18.181  20.377  1.00 27.18           C
ATOM   1219  NE2 HIS B  19      21.897  17.308  20.959  1.00 29.89           N
HETATM 2413  O   HOH B 195      10.371  -1.926  15.288  1.00 30.69           O
HETATM 2486  O   HOH B 268      10.982   5.736  29.100  1.00 46.97           O
HETATM 2528  O   HOH B 310      19.641  17.849  22.382  1.00 61.47           O
TER    1093      PRO A 151
END
""",
   "B",
   19,
   "ND1",
   [ (23.897, 17.622, 20.251)
   ],
   0.1,
   ["non_flip_preference=1000"]
  ],

["1ehz_aromatic_methyl_rotator",
   """\
CRYST1   54.981   33.389   61.921  90.00  90.20  90.00 P 1 21 1      2
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.018188  0.000000  0.000063        0.00000
SCALE2      0.000000  0.029950  0.000000        0.00000
SCALE3      0.000000  0.000000  0.016150        0.00000
ATOM    133  P     U A   7      67.463  47.074  35.969  1.00 44.37           P
ATOM    134  OP1   U A   7      68.318  46.756  34.822  1.00 48.09           O
ATOM    135  OP2   U A   7      67.945  47.948  37.077  1.00 45.68           O
ATOM    136  O5'   U A   7      66.104  47.724  35.455  1.00 40.88           O
ATOM    137  C5'   U A   7      65.285  47.024  34.459  1.00 37.89           C
ATOM    138  C4'   U A   7      64.055  47.852  34.101  1.00 35.74           C
ATOM    139  O4'   U A   7      63.297  48.107  35.326  1.00 38.13           O
ATOM    140  C3'   U A   7      64.317  49.197  33.459  1.00 36.87           C
ATOM    141  O3'   U A   7      63.402  49.394  32.378  1.00 37.45           O
ATOM    142  C2'   U A   7      64.097  50.171  34.624  1.00 36.55           C
ATOM    143  O2'   U A   7      63.595  51.417  34.246  1.00 35.54           O
ATOM    144  C1'   U A   7      63.015  49.475  35.442  1.00 37.23           C
ATOM    145  N1    U A   7      63.056  49.858  36.864  1.00 36.91           N
ATOM    146  C2    U A   7      62.011  50.628  37.343  1.00 34.52           C
ATOM    147  O2    U A   7      61.087  50.966  36.653  1.00 34.65           O
ATOM    148  N3    U A   7      62.112  50.993  38.659  1.00 37.03           N
ATOM    149  C4    U A   7      63.131  50.684  39.541  1.00 40.15           C
ATOM    150  O4    U A   7      63.105  51.143  40.699  1.00 36.62           O
ATOM    151  C5    U A   7      64.179  49.865  38.971  1.00 36.52           C
ATOM    152  C6    U A   7      64.106  49.490  37.691  1.00 36.25           C
HETATM 1061  P   5MC A  49      65.638  54.464  33.461  1.00 31.61           P
HETATM 1062  OP1 5MC A  49      66.252  55.251  34.624  1.00 26.54           O
HETATM 1063  OP2 5MC A  49      65.643  52.989  33.540  1.00 33.43           O
HETATM 1064  O5' 5MC A  49      64.126  54.972  33.441  1.00 31.88           O
HETATM 1065  C5' 5MC A  49      63.204  54.605  32.392  1.00 32.68           C
HETATM 1066  C4' 5MC A  49      61.796  55.006  32.810  1.00 30.80           C
HETATM 1067  O4' 5MC A  49      61.292  54.110  33.848  1.00 32.40           O
HETATM 1068  C3' 5MC A  49      61.745  56.381  33.468  1.00 31.44           C
HETATM 1069  O3' 5MC A  49      61.666  57.398  32.483  1.00 31.25           O
HETATM 1070  C2' 5MC A  49      60.485  56.287  34.312  1.00 34.62           C
HETATM 1071  O2' 5MC A  49      59.414  56.469  33.417  1.00 35.41           O
HETATM 1072  C1' 5MC A  49      60.568  54.849  34.825  1.00 33.10           C
HETATM 1073  N1  5MC A  49      61.331  54.787  36.102  1.00 32.34           N
HETATM 1074  C2  5MC A  49      60.737  55.340  37.271  1.00 32.42           C
HETATM 1075  O2  5MC A  49      59.568  55.832  37.197  1.00 28.24           O
HETATM 1076  N3  5MC A  49      61.420  55.348  38.428  1.00 29.50           N
HETATM 1077  C4  5MC A  49      62.626  54.785  38.503  1.00 32.97           C
HETATM 1078  N4  5MC A  49      63.188  54.744  39.702  1.00 28.75           N
HETATM 1079  C5  5MC A  49      63.270  54.215  37.358  1.00 34.81           C
HETATM 1080  C6  5MC A  49      62.583  54.243  36.161  1.00 33.69           C
HETATM 1081  CM5 5MC A  49      64.655  53.581  37.465  1.00 34.52           C
HETATM 1692  O   HOH A 131      65.908  53.779  40.692  1.00 47.68           O
HETATM 1783  O   HOH A 675      68.361  55.414  38.613  1.00 58.70           O
TER    1653        A A  76
END
""",
   "A",
   49,
   "HM51",
   [ (64.602, 52.635, 37.674), (65.188, 53.997, 38.160), (65.153, 53.663, 36.637)
   ],
   0.1,
   []
  ]

]

def RunRegressionTests():
  for tc in testCases:
    name = tc[0]
    pdb_raw = tc[1]
    chain = tc[2]
    resID = tc[3]
    atomName = tc[4]
    positions = tc[5]
    maxDist = tc[6]
    extraArgs = tc[7]

    print('Testing regression on', name)

    pdb_file = "./deleteme.pdb"
    with open(pdb_file, "w") as f:
      f.write(pdb_raw)
    if (op.exists("./deletemeFH.pdb")):
      os.remove("./deletemeFH.pdb")
    if (op.exists("./deleteme_description.txt")):
      os.remove("./deleteme_description.txt")
    out = StringIO()
    try:
      # Run the program
      args = [pdb_file, "add_flip_movers=True", "output.description_file=./deleteme_description.txt"]
      args.extend(extraArgs)
      results = run_program(program_class=reduce2.Program, logger=out, args=args)
      # Check the position of the atom to see if it is near enough to one of the expected locations.
      found = False
      for c in results.model.chains():
        if c.id == chain:
          for rg in c.residue_groups():
            if rg.resseq_as_int() == resID:
              for atom in rg.atoms():
                if atom.name.strip().upper() == atomName:
                  found = True
                  loc = atom.xyz
                  closeEnough = False
                  for pos in positions:
                    dist = math.sqrt( (loc[0]-pos[0])*(loc[0]-pos[0]) +
                      (loc[1]-pos[1])*(loc[1]-pos[1]) +
                      (loc[2]-pos[2])*(loc[2]-pos[2])
                    )
                    if dist <= maxDist:
                      closeEnough = True
                  if not closeEnough:
                    return ("Atom "+chain+" "+str(resID)+" "+atomName+" in "+name+" too far from expected locations: " +
                            "Found at " + str(loc) + " but expected one of " + str(positions)
                           )
      if not found:
        return "Did not find atom "+atomName+" in chain "+chain+" residue "+str(resID)+" of "+name
    except Exception as e:
      return "Exception when running reduce2: "+str(e)
    if (op.exists(pdb_file)):
      os.remove(pdb_file)
    if (op.exists("./deletemeFH.pdb")):
      os.remove("./deletemeFH.pdb")
    if (op.exists("./deleteme_description.txt")):
      os.remove("./deleteme_description.txt")

  return ""

if __name__ == '__main__':

  ret = RunReduceTests()
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)

  ret = RunRegressionTests()
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)
