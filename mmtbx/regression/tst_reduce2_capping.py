"""reduce2 with qmready=True capping=True: every non-H atom whose valence is
unsatisfied because a bonded non-H atom is missing gets as many H as the
missing bonds need, as the very last step (mmtbx/hydrogens/capping.py). The
capped model is written without being interpreted again.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from libtbx.utils import null_out, Sorry
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

# LYS A 1 truncated at CD (no CE, NZ; phenix_regression
# ala_lys_arg_ser_tyr_neutron_hd.pdb), a chain break, ASP A 5 - GLY A 6
# (1ywf A 128-129, the Glu cut to Gly). Charge: each fragment NH3+ +1 and
# COO- -1, ASP -1, the truncated LYS (no NZ) 0 -> -1
pdb_lys_asp_gly = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   LYS A   1      27.221  24.554  25.195  1.00  5.98           N
ATOM      2  CA  LYS A   1      25.851  24.605  24.700  1.00  6.43           C
ATOM      3  C   LYS A   1      25.466  26.033  24.327  1.00  6.30           C
ATOM      4  O   LYS A   1      25.212  26.865  25.198  1.00  6.91           O
ATOM      5  CB  LYS A   1      24.884  24.064  25.755  1.00  7.86           C
ATOM      6  CG  LYS A   1      23.430  24.044  25.312  1.00 11.05           C
ATOM      7  CD  LYS A   1      22.527  23.473  26.393  1.00 16.83           C
ATOM      8  N   ASP A   5       9.580   7.432   8.699  1.00 37.74           N
ATOM      9  CA  ASP A   5      10.692   8.369   8.630  1.00 38.33           C
ATOM     10  C   ASP A   5      10.828   9.231   9.884  1.00 37.17           C
ATOM     11  O   ASP A   5      11.938   9.511  10.322  1.00 36.50           O
ATOM     12  CB  ASP A   5      10.625   9.197   7.354  1.00 39.43           C
ATOM     13  CG  ASP A   5      11.281   8.467   6.194  1.00 43.72           C
ATOM     14  OD1 ASP A   5      10.549   7.831   5.406  1.00 49.92           O
ATOM     15  OD2 ASP A   5      12.531   8.382   6.067  1.00 51.10           O
ATOM     16  N   GLY A   6       9.691   9.552  10.506  1.00 36.16           N
ATOM     17  CA  GLY A   6       9.668  10.262  11.762  1.00 35.61           C
ATOM     18  C   GLY A   6      10.374   9.397  12.850  1.00 34.61           C
ATOM     19  O   GLY A   6      11.228   9.895  13.566  1.00 33.43           O
END
"""
expected_charge = -1

# No dictionary for UNL (three carbons, tst_capping.py), next to ASP-GLY
pdb_asp_gly_unl = "".join(l + "\n" for l in pdb_lys_asp_gly.splitlines()
  if not (l.startswith("ATOM") and l[17:20] == "LYS")).replace("END\n", """\
HETATM   20  C1  UNL B   1      40.000  40.000  40.000  1.00 20.00           C
HETATM   21  C2  UNL B   1      41.530  40.000  40.000  1.00 20.00           C
HETATM   22  C3  UNL B   1      42.103  41.419  40.000  1.00 20.00           C
END
""")

# THR A 46 (make_fake_anomalous_data.py, shifted by 10 A) without CG2, and a
# CL ion at the site of the CG2
pdb_thr_cl = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   THR A  46      14.783  17.153  14.239  1.00 19.37           N
ATOM      2  CA  THR A  46      15.952  18.027  14.226  1.00 19.58           C
ATOM      3  C   THR A  46      16.923  17.556  15.292  1.00 19.22           C
ATOM      4  O   THR A  46      16.516  16.968  16.298  1.00 16.37           O
ATOM      5  CB  THR A  46      15.588  19.477  14.547  1.00 17.37           C
ATOM      6  OG1 THR A  46      14.999  19.537  15.849  1.00 18.86           O
HETATM    7 CL    CL A  47      14.599  20.035  13.514  1.00 21.06          CL
END
"""

# TYR A 5 without CE1 (tst_capping.py): CD1 and CZ both miss it, an
# internal gap capping refuses
pdb_tyr_no_ce1 = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   TYR A   5      25.496  22.215  27.097  1.00  5.99           N
ATOM      2  CA  TYR A   5      25.049  23.419  26.405  1.00  6.90           C
ATOM      3  C   TYR A   5      23.849  23.122  25.512  1.00  6.28           C
ATOM      4  O   TYR A   5      22.884  22.491  25.943  1.00  7.24           O
ATOM      5  OXT TYR A   5      23.846  23.524  24.329  1.00  6.28           O
ATOM      6  CB  TYR A   5      26.189  24.009  25.571  1.00  6.79           C
ATOM      7  CG  TYR A   5      25.823  25.289  24.854  1.00 10.06           C
ATOM      8  CD1 TYR A   5      25.369  25.268  23.542  1.00 10.36           C
ATOM      9  CD2 TYR A   5      25.931  26.519  25.489  1.00 12.73           C
ATOM     10  CE2 TYR A   5      25.597  27.692  24.838  1.00 11.68           C
ATOM     11  CZ  TYR A   5      25.149  27.644  23.535  1.00 19.05           C
ATOM     12  OH  TYR A   5      24.816  28.809  22.883  1.00 14.32           O
END
"""

# DNA: 1bna chain A 1-3, a 5'-OH start (no P on DC 1) and a 3' end (DC 3)
pdb_dna3 = """\
CRYST1   24.870   40.390   66.200  90.00  90.00  90.00 P 21 21 21
ATOM      1  O5'  DC A   1      18.935  34.195  25.617  1.00 64.35           O
ATOM      2  C5'  DC A   1      19.130  33.921  24.219  1.00 44.69           C
ATOM      3  C4'  DC A   1      19.961  32.668  24.100  1.00 31.28           C
ATOM      4  O4'  DC A   1      19.360  31.583  24.852  1.00 37.45           O
ATOM      5  C3'  DC A   1      20.172  32.122  22.694  1.00 46.72           C
ATOM      6  O3'  DC A   1      21.350  31.325  22.681  1.00 48.89           O
ATOM      7  C2'  DC A   1      18.948  31.223  22.647  1.00 30.88           C
ATOM      8  C1'  DC A   1      19.231  30.482  23.944  1.00 36.58           C
ATOM      9  N1   DC A   1      18.070  29.661  24.380  1.00 40.51           N
ATOM     10  C2   DC A   1      18.224  28.454  25.015  1.00 16.62           C
ATOM     11  O2   DC A   1      19.360  28.014  25.214  1.00 27.75           O
ATOM     12  N3   DC A   1      17.143  27.761  25.377  1.00 20.55           N
ATOM     13  C4   DC A   1      15.917  28.226  25.120  1.00 34.72           C
ATOM     14  N4   DC A   1      14.828  27.477  25.444  1.00 40.31           N
ATOM     15  C5   DC A   1      15.719  29.442  24.471  1.00 30.78           C
ATOM     16  C6   DC A   1      16.843  30.171  24.101  1.00 25.90           C
ATOM     17  P    DG A   2      22.409  31.286  21.483  1.00 58.85           P
ATOM     18  OP1  DG A   2      23.536  32.157  21.851  1.00 57.82           O
ATOM     19  OP2  DG A   2      21.822  31.459  20.139  1.00 78.33           O
ATOM     20  O5'  DG A   2      22.840  29.751  21.498  1.00 40.36           O
ATOM     21  C5'  DG A   2      23.543  29.175  22.594  1.00 47.19           C
ATOM     22  C4'  DG A   2      23.494  27.709  22.279  1.00 47.81           C
ATOM     23  O4'  DG A   2      22.193  27.252  22.674  1.00 38.76           O
ATOM     24  C3'  DG A   2      23.693  27.325  20.807  1.00 28.58           C
ATOM     25  O3'  DG A   2      24.723  26.320  20.653  1.00 40.44           O
ATOM     26  C2'  DG A   2      22.273  26.885  20.416  1.00 21.14           C
ATOM     27  C1'  DG A   2      21.721  26.304  21.716  1.00 33.95           C
ATOM     28  N9   DG A   2      20.237  26.470  21.780  1.00 34.00           N
ATOM     29  C8   DG A   2      19.526  27.584  21.429  1.00 36.47           C
ATOM     30  N7   DG A   2      18.207  27.455  21.636  1.00 32.37           N
ATOM     31  C5   DG A   2      18.083  26.212  22.142  1.00 15.06           C
ATOM     32  C6   DG A   2      16.904  25.525  22.545  1.00 11.88           C
ATOM     33  O6   DG A   2      15.739  25.916  22.518  1.00 21.30           O
ATOM     34  N1   DG A   2      17.197  24.279  23.037  1.00 15.44           N
ATOM     35  C2   DG A   2      18.434  23.717  23.155  1.00  9.63           C
ATOM     36  N2   DG A   2      18.508  22.456  23.668  1.00 16.69           N
ATOM     37  N3   DG A   2      19.537  24.360  22.770  1.00 30.98           N
ATOM     38  C4   DG A   2      19.290  25.594  22.274  1.00 18.56           C
ATOM     39  P    DC A   3      25.064  25.621  19.252  1.00 44.67           P
ATOM     40  OP1  DC A   3      26.506  25.316  19.220  1.00 53.89           O
ATOM     41  OP2  DC A   3      24.559  26.412  18.115  1.00 57.79           O
ATOM     42  O5'  DC A   3      24.260  24.246  19.327  1.00 35.42           O
ATOM     43  C5'  DC A   3      24.584  23.285  20.335  1.00 45.75           C
ATOM     44  C4'  DC A   3      23.523  22.233  20.245  1.00 43.02           C
ATOM     45  O4'  DC A   3      22.256  22.844  20.453  1.00 36.85           O
ATOM     46  C3'  DC A   3      23.424  21.557  18.903  1.00 40.14           C
ATOM     47  O3'  DC A   3      24.121  20.309  18.928  1.00 49.62           O
ATOM     48  C2'  DC A   3      21.930  21.406  18.661  1.00 53.79           C
ATOM     49  C1'  DC A   3      21.278  21.966  19.909  1.00 22.18           C
ATOM     50  N1   DC A   3      20.196  22.889  19.521  1.00 25.44           N
ATOM     51  C2   DC A   3      18.909  22.584  19.816  1.00 19.81           C
ATOM     52  O2   DC A   3      18.685  21.512  20.382  1.00 29.92           O
ATOM     53  N3   DC A   3      17.935  23.447  19.502  1.00 21.59           N
ATOM     54  C4   DC A   3      18.217  24.603  18.897  1.00 14.01           C
ATOM     55  N4   DC A   3      17.221  25.499  18.629  1.00 26.88           N
ATOM     56  C5   DC A   3      19.526  24.945  18.571  1.00 27.59           C
ATOM     57  C6   DC A   3      20.537  24.048  18.899  1.00 27.05           C
TER
END
"""

# RNA: 1ehz chain A 2-4 without the P, OP1, OP2 of C 2: a 5'-OH start and a
# 3' end (G 4); every residue has O2'
pdb_rna3 = """\
CRYST1   54.981   33.389   61.921  90.00  90.20  90.00 P 1 21 1
ATOM     28  O5'   C A   2      55.563  49.261  54.342  1.00 98.27           O
ATOM     29  C5'   C A   2      55.925  49.246  55.742  1.00 95.40           C
ATOM     30  C4'   C A   2      56.836  48.075  56.049  1.00 93.33           C
ATOM     31  O4'   C A   2      56.122  46.828  55.830  1.00 92.18           O
ATOM     32  C3'   C A   2      58.090  47.947  55.197  1.00 92.75           C
ATOM     33  O3'   C A   2      59.174  48.753  55.651  1.00 92.89           O
ATOM     34  C2'   C A   2      58.416  46.463  55.298  1.00 91.81           C
ATOM     35  O2'   C A   2      59.140  46.136  56.466  1.00 91.36           O
ATOM     36  C1'   C A   2      57.022  45.836  55.356  1.00 90.59           C
ATOM     37  N1    C A   2      56.570  45.364  54.029  1.00 88.84           N
ATOM     38  C2    C A   2      57.094  44.157  53.520  1.00 88.64           C
ATOM     39  O2    C A   2      57.921  43.516  54.198  1.00 88.97           O
ATOM     40  N3    C A   2      56.686  43.721  52.301  1.00 87.36           N
ATOM     41  C4    C A   2      55.802  44.437  51.597  1.00 87.11           C
ATOM     42  N4    C A   2      55.430  43.972  50.397  1.00 86.30           N
ATOM     43  C5    C A   2      55.259  45.660  52.089  1.00 86.87           C
ATOM     44  C6    C A   2      55.663  46.080  53.296  1.00 88.01           C
ATOM     45  P     G A   3      60.184  49.419  54.574  1.00 92.31           P
ATOM     46  OP1   G A   3      61.015  50.422  55.295  1.00 92.97           O
ATOM     47  OP2   G A   3      59.371  49.857  53.404  1.00 91.56           O
ATOM     48  O5'   G A   3      61.137  48.219  54.105  1.00 88.57           O
ATOM     49  C5'   G A   3      62.175  47.724  54.969  1.00 83.44           C
ATOM     50  C4'   G A   3      62.769  46.443  54.422  1.00 79.87           C
ATOM     51  O4'   G A   3      61.734  45.427  54.299  1.00 78.36           O
ATOM     52  C3'   G A   3      63.405  46.499  53.040  1.00 78.97           C
ATOM     53  O3'   G A   3      64.741  47.029  53.060  1.00 79.76           O
ATOM     54  C2'   G A   3      63.359  45.032  52.608  1.00 77.19           C
ATOM     55  O2'   G A   3      64.411  44.256  53.155  1.00 77.80           O
ATOM     56  C1'   G A   3      62.018  44.572  53.194  1.00 73.98           C
ATOM     57  N9    G A   3      60.934  44.675  52.202  1.00 68.20           N
ATOM     58  C8    G A   3      60.024  45.702  52.050  1.00 65.03           C
ATOM     59  N7    G A   3      59.252  45.556  51.003  1.00 62.99           N
ATOM     60  C5    G A   3      59.655  44.348  50.447  1.00 59.95           C
ATOM     61  C6    G A   3      59.189  43.675  49.292  1.00 55.65           C
ATOM     62  O6    G A   3      58.287  44.013  48.522  1.00 53.32           O
ATOM     63  N1    G A   3      59.893  42.491  49.072  1.00 54.00           N
ATOM     64  C2    G A   3      60.906  42.006  49.876  1.00 55.46           C
ATOM     65  N2    G A   3      61.512  40.873  49.479  1.00 48.16           N
ATOM     66  N3    G A   3      61.312  42.605  50.983  1.00 56.69           N
ATOM     67  C4    G A   3      60.666  43.774  51.193  1.00 61.76           C
ATOM     68  P     G A   4      65.295  47.868  51.793  1.00 79.34           P
ATOM     69  OP1   G A   4      66.538  48.562  52.246  1.00 80.87           O
ATOM     70  OP2   G A   4      64.193  48.679  51.209  1.00 79.00           O
ATOM     71  O5'   G A   4      65.720  46.752  50.724  1.00 75.17           O
ATOM     72  C5'   G A   4      66.789  45.843  51.019  1.00 68.95           C
ATOM     73  C4'   G A   4      66.749  44.634  50.114  1.00 65.13           C
ATOM     74  O4'   G A   4      65.484  43.939  50.258  1.00 61.83           O
ATOM     75  C3'   G A   4      66.881  44.840  48.611  1.00 62.79           C
ATOM     76  O3'   G A   4      68.230  44.977  48.176  1.00 61.75           O
ATOM     77  C2'   G A   4      66.318  43.538  48.064  1.00 60.58           C
ATOM     78  O2'   G A   4      67.283  42.514  48.122  1.00 59.59           O
ATOM     79  C1'   G A   4      65.192  43.241  49.051  1.00 58.29           C
ATOM     80  N9    G A   4      63.923  43.716  48.500  1.00 53.36           N
ATOM     81  C8    G A   4      63.204  44.843  48.842  1.00 49.19           C
ATOM     82  N7    G A   4      62.140  45.013  48.107  1.00 46.88           N
ATOM     83  C5    G A   4      62.144  43.926  47.243  1.00 45.95           C
ATOM     84  C6    G A   4      61.246  43.573  46.206  1.00 43.95           C
ATOM     85  O6    G A   4      60.182  44.136  45.874  1.00 42.46           O
ATOM     86  N1    G A   4      61.672  42.428  45.540  1.00 40.34           N
ATOM     87  C2    G A   4      62.788  41.686  45.867  1.00 42.60           C
ATOM     88  N2    G A   4      63.034  40.588  45.135  1.00 38.91           N
ATOM     89  N3    G A   4      63.612  41.994  46.850  1.00 44.40           N
ATOM     90  C4    G A   4      63.239  43.117  47.480  1.00 48.65           C
TER
END
"""

# the same RNA with the deposited 5' phosphate of C 2 (P OP1 OP2, no OP3)
pdb_rna3_5p = """\
CRYST1   54.981   33.389   61.921  90.00  90.20  90.00 P 1 21 1
ATOM     25  P     C A   2      54.635  50.420  53.741  1.00100.19           P
ATOM     26  OP1   C A   2      55.145  51.726  54.238  1.00100.19           O
ATOM     27  OP2   C A   2      54.465  50.204  52.269  1.00100.19           O
ATOM     28  O5'   C A   2      55.563  49.261  54.342  1.00 98.27           O
ATOM     29  C5'   C A   2      55.925  49.246  55.742  1.00 95.40           C
ATOM     30  C4'   C A   2      56.836  48.075  56.049  1.00 93.33           C
ATOM     31  O4'   C A   2      56.122  46.828  55.830  1.00 92.18           O
ATOM     32  C3'   C A   2      58.090  47.947  55.197  1.00 92.75           C
ATOM     33  O3'   C A   2      59.174  48.753  55.651  1.00 92.89           O
ATOM     34  C2'   C A   2      58.416  46.463  55.298  1.00 91.81           C
ATOM     35  O2'   C A   2      59.140  46.136  56.466  1.00 91.36           O
ATOM     36  C1'   C A   2      57.022  45.836  55.356  1.00 90.59           C
ATOM     37  N1    C A   2      56.570  45.364  54.029  1.00 88.84           N
ATOM     38  C2    C A   2      57.094  44.157  53.520  1.00 88.64           C
ATOM     39  O2    C A   2      57.921  43.516  54.198  1.00 88.97           O
ATOM     40  N3    C A   2      56.686  43.721  52.301  1.00 87.36           N
ATOM     41  C4    C A   2      55.802  44.437  51.597  1.00 87.11           C
ATOM     42  N4    C A   2      55.430  43.972  50.397  1.00 86.30           N
ATOM     43  C5    C A   2      55.259  45.660  52.089  1.00 86.87           C
ATOM     44  C6    C A   2      55.663  46.080  53.296  1.00 88.01           C
ATOM     45  P     G A   3      60.184  49.419  54.574  1.00 92.31           P
ATOM     46  OP1   G A   3      61.015  50.422  55.295  1.00 92.97           O
ATOM     47  OP2   G A   3      59.371  49.857  53.404  1.00 91.56           O
ATOM     48  O5'   G A   3      61.137  48.219  54.105  1.00 88.57           O
ATOM     49  C5'   G A   3      62.175  47.724  54.969  1.00 83.44           C
ATOM     50  C4'   G A   3      62.769  46.443  54.422  1.00 79.87           C
ATOM     51  O4'   G A   3      61.734  45.427  54.299  1.00 78.36           O
ATOM     52  C3'   G A   3      63.405  46.499  53.040  1.00 78.97           C
ATOM     53  O3'   G A   3      64.741  47.029  53.060  1.00 79.76           O
ATOM     54  C2'   G A   3      63.359  45.032  52.608  1.00 77.19           C
ATOM     55  O2'   G A   3      64.411  44.256  53.155  1.00 77.80           O
ATOM     56  C1'   G A   3      62.018  44.572  53.194  1.00 73.98           C
ATOM     57  N9    G A   3      60.934  44.675  52.202  1.00 68.20           N
ATOM     58  C8    G A   3      60.024  45.702  52.050  1.00 65.03           C
ATOM     59  N7    G A   3      59.252  45.556  51.003  1.00 62.99           N
ATOM     60  C5    G A   3      59.655  44.348  50.447  1.00 59.95           C
ATOM     61  C6    G A   3      59.189  43.675  49.292  1.00 55.65           C
ATOM     62  O6    G A   3      58.287  44.013  48.522  1.00 53.32           O
ATOM     63  N1    G A   3      59.893  42.491  49.072  1.00 54.00           N
ATOM     64  C2    G A   3      60.906  42.006  49.876  1.00 55.46           C
ATOM     65  N2    G A   3      61.512  40.873  49.479  1.00 48.16           N
ATOM     66  N3    G A   3      61.312  42.605  50.983  1.00 56.69           N
ATOM     67  C4    G A   3      60.666  43.774  51.193  1.00 61.76           C
ATOM     68  P     G A   4      65.295  47.868  51.793  1.00 79.34           P
ATOM     69  OP1   G A   4      66.538  48.562  52.246  1.00 80.87           O
ATOM     70  OP2   G A   4      64.193  48.679  51.209  1.00 79.00           O
ATOM     71  O5'   G A   4      65.720  46.752  50.724  1.00 75.17           O
ATOM     72  C5'   G A   4      66.789  45.843  51.019  1.00 68.95           C
ATOM     73  C4'   G A   4      66.749  44.634  50.114  1.00 65.13           C
ATOM     74  O4'   G A   4      65.484  43.939  50.258  1.00 61.83           O
ATOM     75  C3'   G A   4      66.881  44.840  48.611  1.00 62.79           C
ATOM     76  O3'   G A   4      68.230  44.977  48.176  1.00 61.75           O
ATOM     77  C2'   G A   4      66.318  43.538  48.064  1.00 60.58           C
ATOM     78  O2'   G A   4      67.283  42.514  48.122  1.00 59.59           O
ATOM     79  C1'   G A   4      65.192  43.241  49.051  1.00 58.29           C
ATOM     80  N9    G A   4      63.923  43.716  48.500  1.00 53.36           N
ATOM     81  C8    G A   4      63.204  44.843  48.842  1.00 49.19           C
ATOM     82  N7    G A   4      62.140  45.013  48.107  1.00 46.88           N
ATOM     83  C5    G A   4      62.144  43.926  47.243  1.00 45.95           C
ATOM     84  C6    G A   4      61.246  43.573  46.206  1.00 43.95           C
ATOM     85  O6    G A   4      60.182  44.136  45.874  1.00 42.46           O
ATOM     86  N1    G A   4      61.672  42.428  45.540  1.00 40.34           N
ATOM     87  C2    G A   4      62.788  41.686  45.867  1.00 42.60           C
ATOM     88  N2    G A   4      63.034  40.588  45.135  1.00 38.91           N
ATOM     89  N3    G A   4      63.612  41.994  46.850  1.00 44.40           N
ATOM     90  C4    G A   4      63.239  43.117  47.480  1.00 48.65           C
TER
END
"""

def run(model_text, prefix, extra_args=(), fn_out=None, in_suffix=".pdb"):
  """fn_out None: the default output name (with the charge for add)"""
  fn_in = prefix + in_suffix
  with open(fn_in, "w") as f:
    f.write(model_text)
  if fn_out is not None and os.path.exists(fn_out): os.remove(fn_out)
  args = [fn_in, "overwrite=True"]
  if fn_out is not None: args.append("output.filename=" + fn_out)
  return run_program(program_class=reduce2.Program, logger=null_out(),
    args=args + list(extra_args))

def check_sorry(model_text, prefix, extra_args, *words):
  fn_out = prefix + "_out.pdb"
  try:
    run(model_text, prefix, extra_args, fn_out=fn_out)
  except Sorry as e:
    for word in words:
      assert word in str(e), (word, str(e))
    assert not os.path.exists(fn_out)
    return str(e)
  raise AssertionError("no Sorry, %s written" % fn_out)

def names_by_residue(hierarchy):
  result = {}
  for rg in hierarchy.residue_groups():
    result.setdefault(rg.resseq.strip(), []).extend(
      a.name.strip() for a in rg.atoms())
  return result

def check_capped(hierarchy):
  names = names_by_residue(hierarchy)
  for name in ("HCE", "HD2", "HD3"):
    assert name in names["1"], (name, names["1"])
  # the break: OXT on LYS 1, NH3+ on ASP 5
  assert "OXT" in names["1"], names["1"]
  for name in ("H1", "H2", "H3"):
    assert name in names["5"], (name, names["5"])
  assert "OXT" in names["6"], names["6"]

def test_add():
  prefix = "tst_reduce2_capping_1"
  fn_out = prefix + "H_charge%d.pdb" % expected_charge
  fn_txt = prefix + "H_charge%d.txt" % expected_charge
  for fn in (fn_out, fn_txt):
    if os.path.exists(fn): os.remove(fn)
  r = run(pdb_lys_asp_gly, prefix, ["qmready=True", "capping=True"])
  assert r.total_charge == expected_charge, r.total_charge
  assert os.path.isfile(fn_out), fn_out
  assert os.path.isfile(fn_txt), fn_txt
  # read back without processing
  check_capped(iotbx.pdb.input(fn_out).construct_hierarchy())
  # the returned model is the written one
  check_capped(r.model.get_hierarchy())
  with open(fn_txt) as f:
    description = f.read()
  assert "Capping" in description, description
  # placement leaves the truncated CD bare; capping gives it the dictionary
  # names first
  assert "chain A LYS 1: CD +3 H (HD2 HD3 HCE); missing CE NZ" in \
    description, description
  assert "Capped 1 atoms in 1 residues with 3 H" in description, description

def test_cif_output():
  prefix = "tst_reduce2_capping_2"
  fn_out = prefix + "_out.cif"
  r = run(pdb_lys_asp_gly, prefix, ["qmready=True", "capping=True"],
    fn_out=fn_out)
  assert r.total_charge == expected_charge, r.total_charge
  assert os.path.isfile(fn_out), fn_out
  check_capped(iotbx.pdb.input(fn_out).construct_hierarchy())

def test_requires_qmready():
  message = check_sorry(pdb_lys_asp_gly, "tst_reduce2_capping_3",
    ["capping=True"])
  assert message == "capping=True requires qmready=True.", message

def test_remove():
  check_sorry(pdb_lys_asp_gly, "tst_reduce2_capping_4",
    ["capping=True", "qmready=True", "approach=remove"], "capping", "remove")

def test_guess_flag_alone():
  check_sorry(pdb_lys_asp_gly, "tst_reduce2_capping_5",
    ["capping_guess_for_unknowns=True"], "capping_guess_for_unknowns",
    "capping=True")

def test_clash():
  # THR without CG2: CB has CA and OG1, its two caps have fixed directions;
  # a CL ion where CG2 was is 0.4 A from one of them
  check_sorry(pdb_thr_cl, "tst_reduce2_capping_6",
    ["qmready=True", "capping=True"], "Capping H", "THR A  46",
    "CL    CL A  47")

def test_unknown_residue():
  check_sorry(pdb_asp_gly_unl, "tst_reduce2_capping_7",
    ["qmready=True", "capping=True"], "capping=True: no dictionary for",
    "UNL B 1", "capping_guess_for_unknowns=True")
  prefix = "tst_reduce2_capping_8"
  fn_out = prefix + "_out.pdb"
  r = run(pdb_asp_gly_unl, prefix, ["qmready=True", "capping=True",
    "capping_guess_for_unknowns=True"], fn_out=fn_out)
  assert os.path.isfile(fn_out), fn_out
  # ASP-GLY -1, propane 0: the bonds guessed by distance reach electrons
  assert r.total_charge == -1, r.total_charge
  hierarchy = iotbx.pdb.input(fn_out).construct_hierarchy()
  atoms = hierarchy.atoms()
  unl = [a for a in atoms if a.parent().resname == "UNL"]
  carbons = [a for a in unl if a.element.strip() == "C"]
  h = [a for a in unl if a.element_is_hydrogen()]
  assert len(carbons) == 3, [a.name for a in unl]
  # CH3-CH2-CH3: 3 + 2 + 3
  assert len(h) == 8, [a.name for a in unl]
  with open(prefix + "_out.txt") as f:
    description = f.read()
  assert "Capped by element valence (no dictionary):" in description, \
    description
  assert "open-shell cap on" not in description, description

def test_qmready_without_capping_unchanged():
  check_sorry(pdb_lys_asp_gly, "tst_reduce2_capping_9", ["qmready=True"],
    "Model must be atom-complete for qmready=true")

def test_internal_gap():
  # the Sorry of cap_dangling_atoms() reaches the user, nothing is written
  check_sorry(pdb_tyr_no_ce1, "tst_reduce2_capping_10",
    ["qmready=True", "capping=True"], "Cannot cap chain A TYR 5",
    "atom CE1 is missing", "CD1 and CZ", "Rebuild the missing atom")

def test_optimize():
  # approach=optimize needs H: those of a plain add run (no qmready, no
  # capping), then qmready=True capping=True approach=optimize caps the
  # truncated LYS CD
  prefix = "tst_reduce2_capping_11"
  fn_h = prefix + "_h.pdb"
  run(pdb_lys_asp_gly, prefix, fn_out=fn_h)
  with open(fn_h) as f:
    with_h = f.read()
  fn_out = prefix + "_opt.pdb"
  r = run(with_h, prefix + "_in", ["qmready=True", "capping=True",
    "approach=optimize"], fn_out=fn_out)
  assert os.path.isfile(fn_out), fn_out
  # the charge is only calculated by approach=add
  assert r.total_charge is None, r.total_charge
  for hierarchy in (iotbx.pdb.input(fn_out).construct_hierarchy(),
                    r.model.get_hierarchy()):
    names = names_by_residue(hierarchy)
    for name in ("HCE", "HD2", "HD3"):
      assert name in names["1"], (name, names["1"])
  with open(prefix + "_opt.txt") as f:
    description = f.read()
  assert "\nCapping\n" in description, description
  assert "chain A LYS 1: CD +3 H (HD2 HD3 HCE); missing CE NZ" in \
    description, description

def h_names_of(hierarchy, resseq):
  return [a.name.strip() for a in hierarchy.atoms()
          if a.parent().parent().resseq.strip() == str(resseq)
          and a.element.strip() in ("H", "D")]

def description_of(prefix, charge):
  with open(prefix + "H_charge%d.txt" % charge) as f:
    return f.read()

def test_na_termini():
  """
  Nucleic-acid chain ends are not dangling with qmready: reduce2 places
  HO5' on the 5'-OH end and HO3' on the 3' end (the polymer dictionaries
  have neither), so capping adds nothing and the charge is the two
  phosphates, -2. Without qmready the H are not placed (unchanged
  behaviour).
  """
  prefix = "tst_reduce2_capping_na1"
  r = run(pdb_dna3, prefix, ["qmready=True", "capping=True"])
  h = r.model.get_hierarchy()
  assert "HO5'" in h_names_of(h, 1), h_names_of(h, 1)
  assert "HO3'" in h_names_of(h, 3), h_names_of(h, 3)
  assert "HO3'" not in h_names_of(h, 1) and "HO3'" not in h_names_of(h, 2)
  assert "HO5'" not in h_names_of(h, 2) and "HO5'" not in h_names_of(h, 3)
  assert r.total_charge == -2, r.total_charge
  assert "Capped 0 atoms" in description_of(prefix, -2)
  # qmready alone places them too
  r = run(pdb_dna3, prefix + "_q", ["qmready=True"])
  h = r.model.get_hierarchy()
  assert "HO5'" in h_names_of(h, 1) and "HO3'" in h_names_of(h, 3)
  assert r.total_charge == -2, r.total_charge
  # plain reduce2: not placed
  r = run(pdb_dna3, prefix + "_plain", [])
  h = r.model.get_hierarchy()
  assert "HO5'" not in h_names_of(h, 1) and "HO3'" not in h_names_of(h, 3)

def test_rna():
  """
  RNA 5'-OH trimer: HO5', HO3' and every HO2' are placed by reduce2, nothing
  is capped, charge -2. Base of G 3 removed: C1' gets one cap named after
  the missing N9 (HN9); O3' of G 4 removed: C3' gets HO3' (the missing
  atom's name) and no reduce2 HO3'.
  """
  prefix = "tst_reduce2_capping_rna"
  r = run(pdb_rna3, prefix, ["qmready=True", "capping=True"])
  h = r.model.get_hierarchy()
  assert "HO5'" in h_names_of(h, 2), h_names_of(h, 2)
  assert "HO3'" in h_names_of(h, 4), h_names_of(h, 4)
  for resseq in (2, 3, 4):
    assert "HO2'" in h_names_of(h, resseq), h_names_of(h, resseq)
  assert r.total_charge == -2, r.total_charge
  assert "Capped 0 atoms" in description_of(prefix, -2)
  base = ["N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"]
  text = "\n".join(l for l in pdb_rna3.splitlines()
    if not (l.startswith("ATOM") and l[22:26].strip() == "3"
            and l[12:16].strip() in base)) + "\n"
  r = run(text, prefix + "_base", ["qmready=True", "capping=True"])
  h = r.model.get_hierarchy()
  assert "HN9" in h_names_of(h, 3), h_names_of(h, 3)
  assert r.total_charge == -2, r.total_charge
  d = description_of(prefix + "_base", -2)
  assert "chain A G 3: C1' +1 H (HN9)" in d, d
  text = "\n".join(l for l in pdb_rna3.splitlines()
    if not (l.startswith("ATOM") and l[22:26].strip() == "4"
            and l[12:16].strip() == "O3'")) + "\n"
  r = run(text, prefix + "_o3", ["qmready=True", "capping=True"])
  h = r.model.get_hierarchy()
  assert h_names_of(h, 4).count("HO3'") == 1, h_names_of(h, 4)
  assert "chain A G 4: C3' +1 H (HO3')" in description_of(prefix + "_o3", -2)

def test_dna_base_removed():
  """
  DNA: base of DG 2 removed -> C1' +1 H (HN9); the chain ends still get
  their hydroxyl H; charge -2.
  """
  prefix = "tst_reduce2_capping_dna_base"
  base = ["N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"]
  text = "\n".join(l for l in pdb_dna3.splitlines()
    if not (l.startswith("ATOM") and l[22:26].strip() == "2"
            and l[12:16].strip() in base)) + "\n"
  r = run(text, prefix, ["qmready=True", "capping=True"])
  h = r.model.get_hierarchy()
  assert "HN9" in h_names_of(h, 2), h_names_of(h, 2)
  assert "HO5'" in h_names_of(h, 1) and "HO3'" in h_names_of(h, 3)
  assert r.total_charge == -2, r.total_charge
  assert "chain A DG 2: C1' +1 H (HN9)" in description_of(prefix, -2)

def split_base(text, resseq, names, shift=0.3):
  """text with the named atoms of residue resseq split into altlocs A and
  B (occupancy 0.5 each, B shifted along x)"""
  out = []
  for l in text.splitlines():
    if (l.startswith("ATOM") and l[22:26].strip() == str(resseq)
        and l[12:16].strip() in names):
      a = l[:16] + "A" + l[17:54] + "  0.50" + l[60:]
      x = float(l[30:38]) + shift
      b = l[:16] + "B" + l[17:30] + "%8.3f" % x + l[38:54] + "  0.50" + l[60:]
      out += [a, b]
    else:
      out.append(l)
  return "\n".join(out) + "\n"

def test_na_termini_altloc():
  """
  A 3' end whose base is split into altlocs A and B (backbone blank) still
  gets HO3' (looked up across the conformer), once per conformer at most;
  a 5'-OH end split the same way gets HO5'.
  """
  prefix = "tst_reduce2_capping_na_alt"
  base = ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]
  text = split_base(pdb_dna3, 3, base)
  r = run(text, prefix, ["qmready=True", "capping=True"])
  h = r.model.get_hierarchy()
  names = h_names_of(h, 3)
  assert names.count("HO3'") >= 1, names
  assert names.count("HO3'") <= 2, names
  assert "Capped 0 atoms" in description_of(prefix, r.total_charge)
  text = split_base(pdb_dna3, 1, base)
  r = run(text, prefix + "_5p", ["qmready=True"])
  names = h_names_of(r.model.get_hierarchy(), 1)
  assert 1 <= names.count("HO5'") <= 2, names

def test_5p_phosphate_keeps_charge():
  """
  qmready alone on the RNA trimer with its deposited 5' phosphate (no OP3):
  no HO5' (the residue has P), no HOP2/HOP3 (protonation is reduce2's), HO3'
  on the 3' end.
  """
  prefix = "tst_reduce2_capping_5p"
  r = run(pdb_rna3_5p, prefix, ["qmready=True"])
  h = r.model.get_hierarchy()
  assert "HO5'" not in h_names_of(h, 2), h_names_of(h, 2)
  assert not [n for n in h_names_of(h, 2) if n.startswith("HOP")]
  assert "HO3'" in h_names_of(h, 4), h_names_of(h, 4)

def test_phosphorus_refused():
  """
  A P missing an oxygen is never capped (a P-H is a phosphonate): the
  deposited 5' phosphate without OP3, and a P whose O5' is gone at a
  backbone break, both stop with a Sorry naming the atom and the remedy, and
  no file is written.
  """
  prefix = "tst_reduce2_capping_p"
  msg = check_sorry(pdb_rna3_5p, prefix, ["qmready=True", "capping=True"],
    "Cannot cap", "C A 2 P", "OP3", "phosphonate", "OP3 at a 5' end")
  text = "\n".join(l for l in pdb_dna3.splitlines()
    if not (l.startswith("ATOM") and l[22:26].strip() == "2"
            and l[12:16].strip() in ("O5'", "C5'"))) + "\n"
  check_sorry(text, prefix + "_break", ["qmready=True", "capping=True"],
    "Cannot cap", "DG A 2 P", "O5'", "phosphonate")

if __name__ == "__main__":
  test_add()
  test_cif_output()
  test_requires_qmready()
  test_remove()
  test_guess_flag_alone()
  test_clash()
  test_unknown_residue()
  test_qmready_without_capping_unchanged()
  test_internal_gap()
  test_optimize()
  test_na_termini()
  test_rna()
  test_dna_base_removed()
  test_na_termini_altloc()
  test_5p_phosphate_keeps_charge()
  test_phosphorus_refused()
  print("OK")
