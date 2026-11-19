"""
probe2's dot output must match fingerprints recorded from the Python dot loop that ContactDB
replaced, and clashscore2's per-atom-pair contact records must match the per-dot JSON parse.
"""
from __future__ import absolute_import, division, print_function
import hashlib, re, time
import iotbx.cli_parser
import iotbx.phil
import iotbx.pdb
import libtbx.load_env
import mmtbx.model
from iotbx.data_manager import DataManager
from libtbx.utils import null_out
from mmtbx.probe import Helpers
from mmtbx.programs import probe2
from mmtbx.validation import clashscore2

# Alternates on a backbone and a side chain, partial-occupancy waters (from tst_clashscore_2).
# HOH 248 overlaps Val 27 CG2: a clash that does not depend on reduce2's H placement.
pdb_str = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASN A   1      12.388  -5.203  29.298  1.00 13.23           N
ATOM      2  CA  ASN A   1      12.398  -3.931  28.500  1.00  6.85           C
ATOM      3  CB  ASN A   1      13.493  -3.033  29.043  1.00  7.32           C
ATOM      4  CG  ASN A   1      13.414  -1.559  28.658  1.00  7.93           C
ATOM      5  OD1 ASN A   1      14.388  -0.782  28.842  1.00 10.99           O
ATOM      6  ND2 ASN A   1      12.315  -1.201  28.134  1.00 11.12           N
ATOM      7  C  AASN A   1      12.384  -4.463  27.091  0.45  6.81           C
ATOM      8  O  AASN A   1      12.645  -5.635  26.750  0.45  8.28           O
ATOM      9  C  BASN A   1      12.891  -4.205  27.064  0.55  9.54           C
ATOM     10  O  BASN A   1      14.004  -4.638  26.750  0.55 12.04           O
ATOM    180  N   GLU A  25      17.119   0.317  19.886  1.00  4.46           N
ATOM    181  CA  GLU A  25      17.951   0.073  21.066  1.00  4.05           C
ATOM    182  C   GLU A  25      17.083  -0.029  22.303  1.00  3.92           C
ATOM    183  O   GLU A  25      16.053   0.642  22.379  1.00  4.03           O
ATOM    184  CB  GLU A  25      18.993   1.197  21.334  1.00  5.58           C
ATOM    185  CG  GLU A  25      20.050   1.356  20.259  1.00  6.70           C
ATOM    186  CD AGLU A  25      21.014   2.440  20.770  0.47  7.87           C
ATOM    187  OE1AGLU A  25      22.134   2.188  21.301  0.47  9.47           O
ATOM    188  OE2AGLU A  25      20.526   3.611  20.664  0.47  7.81           O
ATOM    189  CD BGLU A  25      20.877   2.608  20.405  0.53  8.00           C
ATOM    190  OE1BGLU A  25      21.950   2.605  19.769  0.53  8.34           O
ATOM    191  OE2BGLU A  25      20.372   3.528  21.137  0.53  7.48           O
ATOM    192  N   VAL A  26      17.559  -0.774  23.291  1.00  4.51           N
ATOM    193  CA  VAL A  26      16.931  -0.812  24.603  1.00  4.35           C
ATOM    194  C   VAL A  26      17.394   0.427  25.384  1.00  4.53           C
ATOM    195  O   VAL A  26      18.565   0.603  25.717  1.00  6.63           O
ATOM    196  CB  VAL A  26      17.285  -2.120  25.322  1.00  6.28           C
ATOM    197  CG1 VAL A  26      16.710  -2.090  26.737  1.00  8.38           C
ATOM    198  CG2 VAL A  26      16.793  -3.328  24.541  1.00  7.99           C
ATOM    199  N   VAL A  27      16.430   1.306  25.641  1.00  4.09           N
ATOM    200  CA  VAL A  27      16.630   2.593  26.282  1.00  4.04           C
ATOM    201  C   VAL A  27      15.546   2.821  27.325  1.00  3.81           C
ATOM    202  O   VAL A  27      14.515   2.145  27.331  1.00  4.70           O
ATOM    203  CB  VAL A  27      16.645   3.734  25.238  1.00  4.32           C
ATOM    204  CG1 VAL A  27      17.867   3.681  24.325  1.00  5.54           C
ATOM    205  CG2 VAL A  27      15.361   3.732  24.400  1.00  4.82           C
ATOM    206  N   THR A  28      15.797   3.825  28.171  1.00  3.67           N
ATOM    207  CA  THR A  28      14.793   4.268  29.126  1.00  3.78           C
ATOM    208  C   THR A  28      14.774   5.796  29.120  1.00  3.63           C
ATOM    209  O   THR A  28      15.846   6.424  29.117  1.00  4.10           O
ATOM    210  CB  THR A  28      15.041   3.745  30.540  1.00  4.19           C
ATOM    211  OG1 THR A  28      16.342   4.127  31.000  1.00  4.27           O
ATOM    212  CG2 THR A  28      14.895   2.227  30.627  1.00  5.21           C
ATOM    213  N   PRO A  29      13.578   6.399  29.122  1.00  3.53           N
ATOM    214  CA  PRO A  29      12.262   5.767  29.064  1.00  3.65           C
ATOM    215  C   PRO A  29      11.992   5.133  27.695  1.00  3.56           C
ATOM    216  O   PRO A  29      12.719   5.286  26.720  1.00  5.32           O
ATOM    217  CB  PRO A  29      11.293   6.907  29.372  1.00  4.74           C
ATOM    218  CG  PRO A  29      12.017   8.133  28.814  1.00  4.83           C
ATOM    219  CD  PRO A  29      13.477   7.862  29.201  1.00  4.02           C
ATOM    220  N   MET A  30      10.875   4.415  27.628  1.00  4.82           N
ATOM    221  CA  MET A  30      10.472   3.706  26.429  1.00  3.92           C
ATOM    222  C   MET A  30       9.962   4.670  25.386  1.00  3.81           C
ATOM    223  O   MET A  30       9.205   5.595  25.690  1.00  4.84           O
ATOM    224  CB  MET A  30       9.331   2.734  26.807  1.00  4.57           C
ATOM    225  CG  MET A  30       8.873   1.861  25.654  1.00  4.88           C
ATOM    226  SD  MET A  30      10.130   0.613  25.308  1.00  5.54           S
ATOM    227  CE AMET A  30       9.590   0.158  23.659  0.80  5.92           C
ATOM    228  CE BMET A  30       9.842  -0.500  26.625  0.20 14.35           C
ATOM    229  N   GLY A  31      10.375   4.449  24.142  1.00  3.74           N
ATOM    230  CA  GLY A  31       9.846   5.184  23.024  1.00  4.22           C
ATOM    231  C   GLY A  31       8.887   4.357  22.184  1.00  3.46           C
ATOM    232  O   GLY A  31       8.350   3.328  22.642  1.00  3.89           O
ATOM    233  N   ILE A  32       8.655   4.794  20.950  1.00  3.51           N
ATOM    234  CA  ILE A  32       7.877   3.976  20.038  1.00  3.38           C
ATOM    235  C   ILE A  32       8.559   2.602  19.944  1.00  3.44           C
ATOM    236  O   ILE A  32       9.782   2.531  19.754  1.00  3.63           O
ATOM    237  CB  ILE A  32       7.807   4.651  18.658  1.00  3.57           C
ATOM    238  CG1 ILE A  32       7.052   5.977  18.750  1.00  4.18           C
ATOM    239  CG2 ILE A  32       7.139   3.697  17.669  1.00  4.59           C
ATOM    240  CD1 ILE A  32       7.076   6.794  17.464  1.00  4.77           C
ATOM    241  N   PRO A  33       7.817   1.492  20.034  1.00  3.75           N
ATOM    242  CA  PRO A  33       8.459   0.197  19.927  1.00  4.43           C
ATOM    243  C   PRO A  33       9.162   0.012  18.587  1.00  3.94           C
ATOM    244  O   PRO A  33       8.668   0.406  17.535  1.00  4.50           O
ATOM    245  CB  PRO A  33       7.289  -0.810  20.093  1.00  6.13           C
ATOM    246  CG  PRO A  33       6.309  -0.027  20.992  1.00  6.43           C
ATOM    247  CD  PRO A  33       6.394   1.383  20.414  1.00  4.84           C
ATOM    248  N   ALA A  34      10.310  -0.693  18.648  1.00  4.36           N
ATOM    249  CA  ALA A  34      11.123  -0.882  17.448  1.00  4.68           C
ATOM    250  C   ALA A  34      10.374  -1.543  16.299  1.00  4.72           C
ATOM    251  O   ALA A  34      10.729  -1.285  15.141  1.00  5.19           O
ATOM    252  CB  ALA A  34      12.370  -1.700  17.814  1.00  7.22           C
TER
HETATM  156  O   HOH A 173       9.659  -2.091  27.019  0.29  6.50           O
HETATM  191  O   HOH A 208      24.121   2.918  18.439  0.20  7.03           O
HETATM  217  O   HOH A 234      10.755   0.622  29.147  0.23  8.84           O
HETATM  230  O   HOH A 247      23.962   1.544  18.100  0.26  4.79           O
HETATM  231  O   HOH A 248      13.860   5.233  22.899  1.00  5.00           O
END
"""

# exercise_dots checks probe2's dots against the Python loop, not H placement, so it runs on
# pdb_str with H already placed (reduce_hydrogen.place_hydrogens). HOH 249 takes an H-bond
# from Val 26 N-H.
pdb_str_h = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASN A   1      12.388  -5.203  29.298  1.00 13.23           N
ATOM      2  CA  ASN A   1      12.398  -3.931  28.500  1.00  6.85           C
ATOM      3  CB  ASN A   1      13.493  -3.033  29.043  1.00  7.32           C
ATOM      4  CG  ASN A   1      13.414  -1.559  28.658  1.00  7.93           C
ATOM      5  OD1 ASN A   1      14.388  -0.782  28.842  1.00 10.99           O
ATOM      6  ND2 ASN A   1      12.315  -1.201  28.134  1.00 11.12           N
ATOM      7  H1  ASN A   1      12.231  -5.015  30.154  1.00 13.23           H
ATOM      8  H2  ASN A   1      11.749  -5.741  28.990  1.00 13.23           H
ATOM      9  H3  ASN A   1      13.178  -5.607  29.223  1.00 13.23           H
ATOM     10  C  AASN A   1      12.384  -4.463  27.091  0.45  6.81           C
ATOM     11  O  AASN A   1      12.645  -5.635  26.750  0.45  8.28           O
ATOM     12  HA AASN A   1      11.675  -3.286  28.541  0.45  6.85           H
ATOM     13  HB3AASN A   1      14.345  -3.366  28.720  0.45  7.32           H
ATOM     14  HB2AASN A   1      13.466  -3.076  30.012  0.45  7.32           H
ATOM     15 HD21AASN A   1      12.200  -0.384  27.892  0.45 11.12           H
ATOM     16 HD22AASN A   1      11.685  -1.776  28.021  0.45 11.12           H
ATOM     17  C  BASN A   1      12.891  -4.205  27.064  0.55  9.54           C
ATOM     18  O  BASN A   1      14.004  -4.638  26.750  0.55 12.04           O
ATOM     19  HA BASN A   1      11.509  -3.543  28.530  0.55  6.85           H
ATOM     20  HB3BASN A   1      14.345  -3.366  28.720  0.55  7.32           H
ATOM     21  HB2BASN A   1      13.466  -3.076  30.012  0.55  7.32           H
ATOM     22 HD21BASN A   1      12.200  -0.384  27.892  0.55 11.12           H
ATOM     23 HD22BASN A   1      11.685  -1.776  28.021  0.55 11.12           H
ATOM     24  N   GLU A  25      17.119   0.317  19.886  1.00  4.46           N
ATOM     25  CA  GLU A  25      17.951   0.073  21.066  1.00  4.05           C
ATOM     26  C   GLU A  25      17.083  -0.029  22.303  1.00  3.92           C
ATOM     27  O   GLU A  25      16.053   0.642  22.379  1.00  4.03           O
ATOM     28  CB  GLU A  25      18.993   1.197  21.334  1.00  5.58           C
ATOM     29  CG  GLU A  25      20.050   1.356  20.259  1.00  6.70           C
ATOM     30  CD AGLU A  25      21.014   2.440  20.770  0.47  7.87           C
ATOM     31  OE1AGLU A  25      22.134   2.188  21.301  0.47  9.47           O
ATOM     32  OE2AGLU A  25      20.526   3.611  20.664  0.47  7.81           O
ATOM     33  H  AGLU A  25      17.281  -0.213  19.228  0.47  4.46           H
ATOM     34  HA AGLU A  25      18.437  -0.749  20.896  0.47  4.05           H
ATOM     35  HB3AGLU A  25      19.451   1.000  22.166  0.47  5.58           H
ATOM     36  HB2AGLU A  25      18.522   2.042  21.407  0.47  5.58           H
ATOM     37  HG2AGLU A  25      20.532   0.524  20.128  0.47  6.70           H
ATOM     38  HG3AGLU A  25      19.648   1.642  19.424  0.47  6.70           H
ATOM     39  CD BGLU A  25      20.877   2.608  20.405  0.53  8.00           C
ATOM     40  OE1BGLU A  25      21.950   2.605  19.769  0.53  8.34           O
ATOM     41  OE2BGLU A  25      20.372   3.528  21.137  0.53  7.48           O
ATOM     42  H  BGLU A  25      17.281  -0.213  19.228  0.53  4.46           H
ATOM     43  HA BGLU A  25      18.437  -0.749  20.896  0.53  4.05           H
ATOM     44  HB2BGLU A  25      18.522   2.042  21.407  0.53  5.58           H
ATOM     45  HB3BGLU A  25      19.451   1.000  22.166  0.53  5.58           H
ATOM     46  HG2BGLU A  25      20.652   0.597  20.298  0.53  6.70           H
ATOM     47  HG3BGLU A  25      19.613   1.389  19.393  0.53  6.70           H
ATOM     48  N   VAL A  26      17.559  -0.774  23.291  1.00  4.51           N
ATOM     49  CA  VAL A  26      16.931  -0.812  24.603  1.00  4.35           C
ATOM     50  C   VAL A  26      17.394   0.427  25.384  1.00  4.53           C
ATOM     51  O   VAL A  26      18.565   0.603  25.717  1.00  6.63           O
ATOM     52  CB  VAL A  26      17.285  -2.120  25.322  1.00  6.28           C
ATOM     53  CG1 VAL A  26      16.710  -2.090  26.737  1.00  8.38           C
ATOM     54  CG2 VAL A  26      16.793  -3.328  24.541  1.00  7.99           C
ATOM     55  H   VAL A  26      18.255  -1.274  23.226  1.00  4.51           H
ATOM     56  HA  VAL A  26      15.964  -0.787  24.535  1.00  4.35           H
ATOM     57  HB  VAL A  26      18.249  -2.206  25.383  1.00  6.28           H
ATOM     58 HG11 VAL A  26      16.091  -2.830  26.842  1.00  8.38           H
ATOM     59 HG13 VAL A  26      16.246  -1.249  26.872  1.00  8.38           H
ATOM     60 HG12 VAL A  26      17.436  -2.173  27.375  1.00  8.38           H
ATOM     61 HG21 VAL A  26      15.829  -3.273  24.448  1.00  7.99           H
ATOM     62 HG23 VAL A  26      17.033  -4.135  25.023  1.00  7.99           H
ATOM     63 HG22 VAL A  26      17.210  -3.329  23.665  1.00  7.99           H
ATOM     64  N   VAL A  27      16.430   1.306  25.641  1.00  4.09           N
ATOM     65  CA  VAL A  27      16.630   2.593  26.282  1.00  4.04           C
ATOM     66  C   VAL A  27      15.546   2.821  27.325  1.00  3.81           C
ATOM     67  O   VAL A  27      14.515   2.145  27.331  1.00  4.70           O
ATOM     68  CB  VAL A  27      16.645   3.734  25.238  1.00  4.32           C
ATOM     69  CG1 VAL A  27      17.867   3.681  24.325  1.00  5.54           C
ATOM     70  CG2 VAL A  27      15.361   3.732  24.400  1.00  4.82           C
ATOM     71  H   VAL A  27      15.605   1.169  25.441  1.00  4.09           H
ATOM     72  HA  VAL A  27      17.490   2.598  26.731  1.00  4.04           H
ATOM     73  HB  VAL A  27      16.695   4.566  25.735  1.00  4.32           H
ATOM     74 HG11 VAL A  27      18.669   3.757  24.866  1.00  5.54           H
ATOM     75 HG13 VAL A  27      17.822   4.417  23.695  1.00  5.54           H
ATOM     76 HG12 VAL A  27      17.868   2.836  23.848  1.00  5.54           H
ATOM     77 HG21 VAL A  27      15.288   2.883  23.937  1.00  4.82           H
ATOM     78 HG23 VAL A  27      15.402   4.457  23.757  1.00  4.82           H
ATOM     79 HG22 VAL A  27      14.600   3.854  24.988  1.00  4.82           H
ATOM     80  N   THR A  28      15.797   3.825  28.171  1.00  3.67           N
ATOM     81  CA  THR A  28      14.793   4.268  29.126  1.00  3.78           C
ATOM     82  C   THR A  28      14.774   5.796  29.120  1.00  3.63           C
ATOM     83  O   THR A  28      15.846   6.424  29.117  1.00  4.10           O
ATOM     84  CB  THR A  28      15.041   3.745  30.540  1.00  4.19           C
ATOM     85  OG1 THR A  28      16.342   4.127  31.000  1.00  4.27           O
ATOM     86  CG2 THR A  28      14.895   2.227  30.627  1.00  5.21           C
ATOM     87  H   THR A  28      16.538   4.259  28.207  1.00  3.67           H
ATOM     88  HA  THR A  28      13.919   3.927  28.880  1.00  3.78           H
ATOM     89  HB  THR A  28      14.370   4.137  31.121  1.00  4.19           H
ATOM     90  HG1 THR A  28      16.409   4.964  31.009  1.00  4.27           H
ATOM     91 HG21 THR A  28      13.998   1.964  30.366  1.00  5.21           H
ATOM     92 HG23 THR A  28      15.059   1.930  31.536  1.00  5.21           H
ATOM     93 HG22 THR A  28      15.533   1.798  30.036  1.00  5.21           H
ATOM     94  N   PRO A  29      13.578   6.399  29.122  1.00  3.53           N
ATOM     95  CA  PRO A  29      12.262   5.767  29.064  1.00  3.65           C
ATOM     96  C   PRO A  29      11.992   5.133  27.695  1.00  3.56           C
ATOM     97  O   PRO A  29      12.719   5.286  26.720  1.00  5.32           O
ATOM     98  CB  PRO A  29      11.293   6.907  29.372  1.00  4.74           C
ATOM     99  CG  PRO A  29      12.017   8.133  28.814  1.00  4.83           C
ATOM    100  CD  PRO A  29      13.477   7.862  29.201  1.00  4.02           C
ATOM    101  HA  PRO A  29      12.161   5.094  29.755  1.00  3.65           H
ATOM    102  HB3 PRO A  29      11.153   6.983  30.329  1.00  4.74           H
ATOM    103  HB2 PRO A  29      10.445   6.761  28.923  1.00  4.74           H
ATOM    104  HG2 PRO A  29      11.908   8.184  27.852  1.00  4.83           H
ATOM    105  HG3 PRO A  29      11.684   8.944  29.229  1.00  4.83           H
ATOM    106  HD3 PRO A  29      13.659   8.173  30.102  1.00  4.02           H
ATOM    107  HD2 PRO A  29      14.085   8.285  28.574  1.00  4.02           H
ATOM    108  N   MET A  30      10.875   4.415  27.628  1.00  4.82           N
ATOM    109  CA  MET A  30      10.472   3.706  26.429  1.00  3.92           C
ATOM    110  C   MET A  30       9.962   4.670  25.386  1.00  3.81           C
ATOM    111  O   MET A  30       9.205   5.595  25.690  1.00  4.84           O
ATOM    112  CB  MET A  30       9.331   2.734  26.807  1.00  4.57           C
ATOM    113  CG  MET A  30       8.873   1.861  25.654  1.00  4.88           C
ATOM    114  SD  MET A  30      10.130   0.613  25.308  1.00  5.54           S
ATOM    115  H   MET A  30      10.323   4.323  28.281  1.00  4.82           H
ATOM    116  CE AMET A  30       9.590   0.158  23.659  0.80  5.92           C
ATOM    117  HA AMET A  30      11.225   3.215  26.064  0.80  3.92           H
ATOM    118  HB2AMET A  30       9.640   2.150  27.517  0.80  4.57           H
ATOM    119  HB3AMET A  30       8.568   3.250  27.110  0.80  4.57           H
ATOM    120  HG3AMET A  30       8.744   2.405  24.861  0.80  4.88           H
ATOM    121  HG2AMET A  30       8.044   1.415  25.889  0.80  4.88           H
ATOM    122  HE1AMET A  30       8.681  -0.178  23.706  0.80  5.92           H
ATOM    123  HE3AMET A  30      10.180  -0.529  23.312  0.80  5.92           H
ATOM    124  HE2AMET A  30       9.624   0.942  23.088  0.80  5.92           H
ATOM    125  CE BMET A  30       9.842  -0.500  26.625  0.20 14.35           C
ATOM    126  HA BMET A  30      11.225   3.215  26.064  0.20  3.92           H
ATOM    127  HB2BMET A  30       9.640   2.150  27.517  0.20  4.57           H
ATOM    128  HB3BMET A  30       8.568   3.250  27.110  0.20  4.57           H
ATOM    129  HG2BMET A  30       8.044   1.415  25.889  0.20  4.88           H
ATOM    130  HG3BMET A  30       8.744   2.405  24.861  0.20  4.88           H
ATOM    131  HE1BMET A  30      10.467  -1.238  26.557  0.20 14.35           H
ATOM    132  HE3BMET A  30       8.932  -0.831  26.566  0.20 14.35           H
ATOM    133  HE2BMET A  30       9.971  -0.034  27.466  0.20 14.35           H
ATOM    134  N   GLY A  31      10.375   4.449  24.142  1.00  3.74           N
ATOM    135  CA  GLY A  31       9.846   5.184  23.024  1.00  4.22           C
ATOM    136  C   GLY A  31       8.887   4.357  22.184  1.00  3.46           C
ATOM    137  O   GLY A  31       8.350   3.328  22.642  1.00  3.89           O
ATOM    138  H   GLY A  31      10.971   3.867  23.926  1.00  3.74           H
ATOM    139  HA2 GLY A  31       9.371   5.965  23.350  1.00  4.22           H
ATOM    140  HA3 GLY A  31      10.577   5.473  22.456  1.00  4.22           H
ATOM    141  N   ILE A  32       8.655   4.794  20.950  1.00  3.51           N
ATOM    142  CA  ILE A  32       7.877   3.976  20.038  1.00  3.38           C
ATOM    143  C   ILE A  32       8.559   2.602  19.944  1.00  3.44           C
ATOM    144  O   ILE A  32       9.782   2.531  19.754  1.00  3.63           O
ATOM    145  CB  ILE A  32       7.807   4.651  18.658  1.00  3.57           C
ATOM    146  CG1 ILE A  32       7.052   5.977  18.750  1.00  4.18           C
ATOM    147  CG2 ILE A  32       7.139   3.697  17.669  1.00  4.59           C
ATOM    148  CD1 ILE A  32       7.076   6.794  17.464  1.00  4.77           C
ATOM    149  H   ILE A  32       8.932   5.542  20.629  1.00  3.51           H
ATOM    150  HA  ILE A  32       6.968   3.864  20.357  1.00  3.38           H
ATOM    151  HB  ILE A  32       8.704   4.849  18.347  1.00  3.57           H
ATOM    152 HG13 ILE A  32       7.454   6.516  19.449  1.00  4.18           H
ATOM    153 HG12 ILE A  32       6.124   5.792  18.965  1.00  4.18           H
ATOM    154 HG21 ILE A  32       7.095   4.124  16.799  1.00  4.59           H
ATOM    155 HG23 ILE A  32       6.244   3.493  17.983  1.00  4.59           H
ATOM    156 HG22 ILE A  32       7.663   2.883  17.612  1.00  4.59           H
ATOM    157 HD11 ILE A  32       6.579   7.615  17.602  1.00  4.77           H
ATOM    158 HD13 ILE A  32       6.668   6.276  16.753  1.00  4.77           H
ATOM    159 HD12 ILE A  32       7.997   6.999  17.237  1.00  4.77           H
ATOM    160  N   PRO A  33       7.817   1.492  20.034  1.00  3.75           N
ATOM    161  CA  PRO A  33       8.459   0.197  19.927  1.00  4.43           C
ATOM    162  C   PRO A  33       9.162   0.012  18.587  1.00  3.94           C
ATOM    163  O   PRO A  33       8.668   0.406  17.535  1.00  4.50           O
ATOM    164  CB  PRO A  33       7.289  -0.810  20.093  1.00  6.13           C
ATOM    165  CG  PRO A  33       6.309  -0.027  20.992  1.00  6.43           C
ATOM    166  CD  PRO A  33       6.394   1.383  20.414  1.00  4.84           C
ATOM    167  HA  PRO A  33       9.093   0.062  20.649  1.00  4.43           H
ATOM    168  HB2 PRO A  33       6.894  -1.020  19.232  1.00  6.13           H
ATOM    169  HB3 PRO A  33       7.597  -1.623  20.523  1.00  6.13           H
ATOM    170  HG3 PRO A  33       6.600  -0.051  21.917  1.00  6.43           H
ATOM    171  HG2 PRO A  33       5.412  -0.390  20.917  1.00  6.43           H
ATOM    172  HD2 PRO A  33       5.819   1.476  19.638  1.00  4.84           H
ATOM    173  HD3 PRO A  33       6.161   2.047  21.082  1.00  4.84           H
ATOM    174  N   ALA A  34      10.310  -0.693  18.648  1.00  4.36           N
ATOM    175  CA  ALA A  34      11.123  -0.882  17.448  1.00  4.68           C
ATOM    176  C   ALA A  34      10.374  -1.543  16.299  1.00  4.72           C
ATOM    177  O   ALA A  34      10.729  -1.285  15.141  1.00  5.19           O
ATOM    178  CB  ALA A  34      12.370  -1.700  17.814  1.00  7.22           C
ATOM    179  H   ALA A  34      10.625  -1.060  19.359  1.00  4.36           H
ATOM    180  HA  ALA A  34      11.383  -0.007  17.121  1.00  4.68           H
ATOM    181  HB1 ALA A  34      12.910  -1.826  17.018  1.00  7.22           H
ATOM    182  HB3 ALA A  34      12.092  -2.561  18.164  1.00  7.22           H
ATOM    183  HB2 ALA A  34      12.878  -1.218  18.486  1.00  7.22           H
HETATM  184  O   HOH A 173       9.659  -2.091  27.019  0.29  6.50           O
HETATM  185  O   HOH A 208      24.121   2.918  18.439  0.20  7.03           O
HETATM  186  O   HOH A 234      10.755   0.622  29.147  0.23  8.84           O
HETATM  187  O   HOH A 247      23.962   1.544  18.100  0.26  4.79           O
HETATM  188  O   HOH A 248      13.860   5.233  22.899  1.00  5.00           O
HETATM  189  O   HOH A 249      19.834  -2.408  23.079  1.00  5.00           O
END
"""

def frozen():
  """(hierarchy, crystal symmetry) of pdb_str_h."""
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str_h)
  return inp.construct_hierarchy(), inp.crystal_symmetry()

def fresh_model(hierarchy, symmetry):
  # A new model each run: probe2 adds Phantom Hydrogens to the waters of the one it runs on.
  return mmtbx.model.manager(model_input=None, pdb_hierarchy=hierarchy.deep_copy(),
                             crystal_symmetry=symmetry, log=null_out())

def run_probe2(model, args):
  parser = iotbx.cli_parser.CCTBXParser(program_class=probe2.Program, logger=null_out())
  parser.parse_args(args + ["output.write_files=False", "ignore_lack_of_explicit_hydrogens=True"])
  dm = DataManager()
  dm.add_model("model", model)
  p2 = probe2.Program(dm, parser.working_phil.extract(), master_phil=parser.master_phil,
                      logger=null_out())
  p2.overrideModel(model)
  # The header carries the time of the run.
  return re.sub(r"run \d{4}-\d\d-\d\d \d\d:\d\d:\d\d", "run <time>", p2.run()[1])

def two_models(hierarchy):
  h = hierarchy.deep_copy()
  h.only_model().id = "1"
  second = h.only_model().detached_copy()
  second.id = "2"
  h.append_model(second)
  return h

def fingerprint(out):
  """(line count, hash) of probe2 output, leaving out the lines with its version and command."""
  skip = ("@caption", "command:", "program:")
  body = [l for l in out.splitlines() if not l.lstrip().startswith(skip)]
  return len(body), hashlib.sha256("\n".join(body).encode("utf-8")).hexdigest()[:12]

RAW, KIN = ["output.format=raw"], ["output.format=kinemage"]

# One case per ContactDB input (set_params, atom flags, bonds) and per kinemage writer option,
# since kinemage without count_dots or contact_summary is written from ContactDB's arrays.
# A change to radii or dot output changes these; the assertion lists every new value.
DOTS_CASES = [
  (KIN, 1443, "d6036127ffb8"),
  (RAW, 1420, "d54941fa4aac"),
  (["output.format=json"], 1422, "fedd0694e019"),
  (["output.count_dots=True"], 33, "ccd04a48fee7"),
  # Waters as source take the Phantom Hydrogen path.
  (["approach=both", "source_selection=water", "target_selection=all"] + RAW, 456, "f7e87532d07c"),
  (["approach=both", "source_selection=protein", "target_selection=water"] + KIN, 482, "96c4c6c679fc"),
  (["approach=once", "source_selection=protein", "target_selection=not protein"] + RAW, 239, "7d179a454a3c"),
  (["output.contact_summary=True"] + KIN, 1448, "5b7c972e3a2b"),
  (["output.condensed=True"] + RAW, 244, "7af1dbb9b8da"),
  (["minimum_occupancy=0.5"] + RAW, 1366, "d48b67069b77"),
  (["excluded_bond_chain_length=3"] + RAW, 2452, "28c49022791d"),
  (["excluded_bond_chain_length=5"] + RAW, 1169, "f36934dd3981"),
  (["include_mainchain_mainchain=False"] + RAW, 1186, "04a667e834c2"),
  (["include_water_water=True", "source_selection=all"] + RAW, 1858, "9032c9ed59f4"),
  (["overlap_scale_factor=0.3"] + RAW, 1420, "847fa022c2cf"),
  (["probe.probe_radius=0.0"] + RAW, 322, "aa865fe4eb2a"),
  # Weak H-bond and worse overlap dots need these; the lower cutoffs give this model bad overlaps.
  (["probe.allow_weak_hydrogen_bonds=True", "output.separate_worse_clashes=True",
    "probe.clash_cutoff=0.1", "probe.worse_clash_cutoff=0.2"] + RAW, 1459, "e29c7edc7040"),
  (["output.separate_worse_clashes=True"] + KIN, 1446, "fe06808b5670"),
  (["output.report_hydrogen_bonds=False", "output.report_clashes=False"] + RAW, 1098, "b9817d8b8e50"),
  (["output.report_vdws=False"] + KIN, 338, "36ce208b2347"),
  (["output.bin_gaps=True"] + KIN, 1514, "f046ec48103b"),
  (["output.color_by_gap=False", "output.merge_contacts=False"] + KIN, 1444, "4c50995f7bab"),
  (["output.atoms_are_masters=True", "output.altid_as_pointmaster=False"] + KIN, 1443, "d9bd0ea6c1cb"),
]
# Two identical models: dots between them would be all clashes.
TWO_MODELS_CASE = (RAW, 2840, "80312111a108")

def exercise_dots():
  hierarchy, symmetry = frozen()
  cases = [(hierarchy,) + c for c in DOTS_CASES] + [(two_models(hierarchy),) + TWO_MODELS_CASE]
  changed = []
  for h, args, lines, digest in cases:
    got = fingerprint(run_probe2(fresh_model(h, symmetry), args))
    if got != (lines, digest):
      changed.append("  %r," % ((args,) + got,))
  assert not changed, "Changed dot fingerprints:\n" + "\n".join(changed)

def run_clashscore2(**kwds):
  dm = DataManager()
  dm.process_model_str("model", pdb_str)
  probe_params = iotbx.phil.parse(Helpers.probe_phil_parameters).extract()
  return clashscore2.clashscore2(probe_parameters=probe_params.probe, data_manager=dm,
                                 keep_hydrogens=False, out=null_out(), **kwds)

def exercise_contact_records():
  manager = clashscore2.probe_clashscore_manager
  assert manager.use_contact_records
  records = run_clashscore2()
  manager.use_contact_records = False
  try:
    dots = run_clashscore2()
  finally:
    manager.use_contact_records = True
  assert len(records.results) > 0
  assert records.as_JSON() == dots.as_JSON()
  assert records.get_clashscore() == dots.get_clashscore()

def exercise_keep_processed_model():
  assert run_clashscore2().processed_model is None
  kept = run_clashscore2(keep_processed_model=True)
  assert kept.processed_model is kept.processed_master
  assert kept.processed_model.get_restraints_manager() is not None

def exercise_kinemage_neutron_distances():
  from mmtbx.kinemage.validation import make_probe_dots_from_model
  kept = run_clashscore2(keep_processed_model=True, nuclear=True).processed_model
  h, symmetry = kept.get_hierarchy(), kept.crystal_symmetry()
  dots = lambda **kw: re.sub(r"run \d{4}-\d\d-\d\d \d\d:\d\d:\d\d", "run <time>",
                             make_probe_dots_from_model(**kw))
  nuclear = dots(model_manager=fresh_model(h, symmetry), use_neutron_distances=True)
  assert nuclear != dots(model_manager=fresh_model(h, symmetry))
  assert nuclear == dots(model_manager=kept, processed=True, use_neutron_distances=True)

if (__name__ == "__main__"):
  if (libtbx.env.has_module("reduce") and libtbx.env.has_module("probe")):
    t0 = time.time()
    exercise_dots()
    exercise_contact_records()
    exercise_keep_processed_model()
    exercise_kinemage_neutron_distances()
    print("OK. Time: %8.3f" % (time.time() - t0))
  else:
    print("Skipping tst_probe2_contact_records: probe or reduce not configured")
