from __future__ import absolute_import, division, print_function
import math
import iotbx.pdb
import mmtbx.model
from libtbx.utils import Sorry, null_out
from scitbx import matrix
from scitbx.array_family import flex
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.hydrogens import capping
from mmtbx.ligands import electrons

# ------------------------------------------------------------------------------
# Fragments: real coordinates (sources at each string), boxed in a P1 cell.
# Chain ends carry OXT, as qmready gives them: a C without OXT is dangling.

# ARG A 3 truncated at CB, SER A 4 (OXT); phenix_regression
# ala_lys_arg_ser_tyr_neutron_hd.pdb
pdb_arg = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A   3      27.911  23.777  26.214  1.00  5.99           N
ATOM      2  CA  ARG A   3      27.558  25.104  25.723  1.00  6.90           C
ATOM      3  C   ARG A   3      26.047  25.238  25.567  1.00  6.28           C
ATOM      4  O   ARG A   3      25.450  26.196  26.055  1.00  7.24           O
ATOM      5  CB  ARG A   3      28.246  25.375  24.384  1.00  6.79           C
ATOM      6  N   SER A   4      25.443  24.270  24.884  1.00  5.99           N
ATOM      7  CA  SER A   4      24.000  24.258  24.649  1.00  6.90           C
ATOM      8  C   SER A   4      23.547  25.461  23.816  1.00  6.28           C
ATOM      9  O   SER A   4      23.472  25.372  22.590  1.00  7.24           O
ATOM     10  OXT SER A   4      23.253  26.536  24.382  1.00  6.28           O
ATOM     11  CB  SER A   4      23.234  24.186  25.975  1.00  6.79           C
ATOM     12  OG  SER A   4      21.834  24.229  25.757  1.00  8.89           O
END
"""

# LYS truncated at CD (no CE, NZ); same source
pdb_lys = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   LYS A   2      27.221  24.554  25.195  1.00  5.98           N
ATOM      2  CA  LYS A   2      25.851  24.605  24.700  1.00  6.43           C
ATOM      3  C   LYS A   2      25.466  26.033  24.327  1.00  6.30           C
ATOM      4  O   LYS A   2      25.212  26.865  25.198  1.00  6.91           O
ATOM      5  OXT LYS A   2      25.408  26.362  23.122  1.00  6.30           O
ATOM      6  CB  LYS A   2      24.884  24.064  25.755  1.00  7.86           C
ATOM      7  CG  LYS A   2      23.430  24.044  25.312  1.00 11.05           C
ATOM      8  CD  LYS A   2      22.527  23.473  26.393  1.00 16.83           C
END
"""

# ALA-ALA, the first C without O; same source (ALA 1, LYS 2 cut to ALA)
pdb_ala_ala = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      27.278  23.666  27.030  1.00  7.86           N
ATOM      2  CA  ALA A   1      27.257  24.002  25.578  1.00  8.28           C
ATOM      3  C   ALA A   1      25.838  23.913  25.027  1.00  7.07           C
ATOM      4  CB  ALA A   1      28.181  23.072  24.810  1.00  7.41           C
ATOM      5  N   ALA A   2      25.121  25.032  25.071  1.00  5.98           N
ATOM      6  CA  ALA A   2      23.751  25.083  24.576  1.00  6.43           C
ATOM      7  C   ALA A   2      23.366  26.511  24.203  1.00  6.30           C
ATOM      8  O   ALA A   2      23.112  27.343  25.074  1.00  6.91           O
ATOM      9  OXT ALA A   2      23.308  26.840  22.999  1.00  6.30           O
ATOM     10  CB  ALA A   2      22.784  24.542  25.631  1.00  7.86           C
END
"""

# TYR without CE1; same source
pdb_tyr = """\
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

# SER: CB/OG in altloc A, altloc B only CB; same source
pdb_ser_altloc = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   SER A   4      26.979  24.457  25.160  1.00  5.99           N
ATOM      2  CA  SER A   4      25.536  24.445  24.925  1.00  6.90           C
ATOM      3  C   SER A   4      25.083  25.648  24.092  1.00  6.28           C
ATOM      4  O   SER A   4      25.008  25.559  22.866  1.00  7.24           O
ATOM      5  OXT SER A   4      24.789  26.723  24.658  1.00  6.28           O
ATOM      6  CB ASER A   4      24.770  24.373  26.251  0.50  6.79           C
ATOM      7  OG ASER A   4      23.370  24.416  26.033  0.50  8.89           O
ATOM      8  CB BSER A   4      24.462  24.376  26.013  0.50  6.79           C
END
"""

# LYS cut to CB, CA in altlocs A and B, CB blank; same source
pdb_lys_ca_altloc = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   LYS A   2      26.522  24.256  25.481  1.00  5.98           N
ATOM      2  CA ALYS A   2      25.302  24.307  24.986  0.50  6.43           C
ATOM      3  CA BLYS A   2      25.002  24.307  24.986  0.50  6.43           C
ATOM      4  C   LYS A   2      24.767  25.735  24.613  1.00  6.30           C
ATOM      5  O   LYS A   2      24.513  26.567  25.484  1.00  6.91           O
ATOM      6  OXT LYS A   2      24.709  26.064  23.409  1.00  6.30           O
ATOM      7  CB  LYS A   2      24.185  23.766  26.041  1.00  7.86           C
END
"""

# ASP without OD2 (1ywf A 128), CG 1.33 A from NZ of a LYS (isopeptide)
pdb_isopeptide = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A   1      16.406  17.151  22.826  1.00 37.74           N
ATOM      2  CA  ASP A   1      17.518  18.088  22.757  1.00 38.33           C
ATOM      3  C   ASP A   1      17.654  18.950  24.011  1.00 37.17           C
ATOM      4  O   ASP A   1      18.764  19.230  24.449  1.00 36.50           O
ATOM      5  OXT ASP A   1      16.630  19.371  24.592  1.00 37.17           O
ATOM      6  CB  ASP A   1      17.451  18.916  21.481  1.00 39.43           C
ATOM      7  CG  ASP A   1      18.107  18.186  20.321  1.00 43.72           C
ATOM      8  OD1 ASP A   1      17.375  17.550  19.533  1.00 49.92           O
ATOM      9  N   LYS B  10      22.389  23.110  15.559  1.00 20.00           N
ATOM     10  CA  LYS B  10      22.355  22.703  16.959  1.00 20.00           C
ATOM     11  C   LYS B  10      23.766  22.439  17.474  1.00 20.00           C
ATOM     12  O   LYS B  10      24.583  21.819  16.792  1.00 20.00           O
ATOM     13  OXT LYS B  10      24.098  22.860  18.603  1.00 20.00           O
ATOM     14  CB  LYS B  10      21.489  21.452  17.126  1.00 20.00           C
ATOM     15  CG  LYS B  10      21.295  21.012  18.569  1.00 20.00           C
ATOM     16  CD  LYS B  10      20.454  19.750  18.661  1.00 20.00           C
ATOM     17  CE  LYS B  10      20.243  19.318  20.104  1.00 20.00           C
ATOM     18  NZ  LYS B  10      19.427  18.096  20.186  1.00 20.00           N
END
"""

# CYS SG 2.4 A from a ZN, no HG (tst_add_hydrogen_1 pdb_str_005)
pdb_cys_zn = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  N   CYS A  48      21.448  18.823  19.531  1.00 99.72           N
ATOM      2  CA  CYS A  48      20.658  19.965  19.081  1.00106.82           C
ATOM      3  C   CYS A  48      20.863  21.161  20.006  1.00104.36           C
ATOM      4  O   CYS A  48      21.602  21.086  20.984  1.00104.56           O
ATOM      5  OXT CYS A  48      20.269  22.234  19.764  1.00104.36           O
ATOM      6  CB  CYS A  48      19.174  19.605  19.004  1.00105.79           C
ATOM      7  SG  CYS A  48      18.302  19.694  20.580  1.00101.00           S
HETATM    8  ZN   ZN A 102      17.687  17.434  21.052  1.00 98.67          ZN
END
"""

# HEM with FE (mmtbx/monomer_library/tst_pdb_interpretation.py)
pdb_hem = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1A HEM A 154      22.182  19.328  19.676  1.00 18.41           C
HETATM    2  C1B HEM A 154      20.152  22.097  17.158  1.00 20.12           C
HETATM    3  C1C HEM A 154      16.649  21.676  19.512  1.00 13.75           C
HETATM    4  C1D HEM A 154      18.691  18.972  22.085  1.00 14.87           C
HETATM    5  C2A HEM A 154      23.401  19.284  18.897  1.00 21.34           C
HETATM    6  C2B HEM A 154      19.612  23.065  16.229  1.00 18.28           C
HETATM    7  C2C HEM A 154      15.376  21.617  20.206  1.00 20.20           C
HETATM    8  C2D HEM A 154      19.275  18.085  23.070  1.00 18.53           C
HETATM    9  C3A HEM A 154      23.238  20.121  17.843  1.00 19.89           C
HETATM   10  C3B HEM A 154      18.333  23.333  16.617  1.00 14.10           C
HETATM   11  C3C HEM A 154      15.582  20.912  21.356  1.00 20.16           C
HETATM   12  C3D HEM A 154      20.586  17.939  22.746  1.00 17.55           C
HETATM   13  C4A HEM A 154      21.905  20.681  17.947  1.00 14.01           C
HETATM   14  C4B HEM A 154      18.074  22.505  17.786  1.00 15.74           C
HETATM   15  C4C HEM A 154      16.914  20.331  21.250  1.00 11.44           C
HETATM   16  C4D HEM A 154      20.787  18.598  21.473  1.00 12.93           C
HETATM   17  CAA HEM A 154      24.613  18.388  19.198  1.00 17.09           C
HETATM   18  CAB HEM A 154      17.293  24.171  16.134  1.00 18.06           C
HETATM   19  CAC HEM A 154      14.833  20.663  22.544  1.00 17.18           C
HETATM   20  CAD HEM A 154      21.699  17.275  23.565  1.00 19.52           C
HETATM   21  CBA HEM A 154      24.742  17.123  18.310  1.00 15.28           C
HETATM   22  CBB HEM A 154      17.501  25.428  15.492  1.00 10.31           C
HETATM   23  CBC HEM A 154      13.633  21.418  22.838  1.00 23.40           C
HETATM   24  CBD HEM A 154      22.329  18.120  24.686  1.00 26.70           C
HETATM   25  CGA HEM A 154      23.429  16.349  18.186  1.00 14.87           C
HETATM   26  CGD HEM A 154      23.490  17.473  25.414  1.00 23.81           C
HETATM   27  CHA HEM A 154      21.982  18.542  20.791  1.00 13.26           C
HETATM   28  CHB HEM A 154      21.411  21.561  17.015  1.00 14.01           C
HETATM   29  CHC HEM A 154      16.843  22.468  18.401  1.00 13.33           C
HETATM   30  CHD HEM A 154      17.403  19.447  22.183  1.00 11.65           C
HETATM   31  CMA HEM A 154      24.221  20.460  16.723  1.00 14.40           C
HETATM   32  CMB HEM A 154      20.384  23.640  15.038  1.00 11.05           C
HETATM   33  CMC HEM A 154      14.100  22.263  19.684  1.00 21.38           C
HETATM   34  CMD HEM A 154      18.529  17.481  24.263  1.00 17.96           C
HETATM   35  NA  HEM A 154      21.302  20.220  19.096  1.00 10.68           N
HETATM   36  NB  HEM A 154      19.218  21.821  18.132  1.00 18.79           N
HETATM   37  NC  HEM A 154      17.540  20.828  20.129  1.00 15.30           N
HETATM   38  ND  HEM A 154      19.627  19.252  21.110  1.00 15.13           N
HETATM   39  O1A HEM A 154      23.015  15.798  19.230  1.00 22.79           O
HETATM   40  O1D HEM A 154      23.232  16.916  26.506  1.00 47.56           O
HETATM   41  O2A HEM A 154      22.886  16.339  17.066  1.00 31.78           O
HETATM   42  O2D HEM A 154      24.607  17.548  24.863  1.00 40.84           O
HETATM   43  FE  HEM A 154      19.392  20.451  19.556  1.00 17.04          FE
END
"""

# CYS A 128 - CYS A 232 disulfide (tst_add_hydrogen_1)
pdb_disulfide = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  N   CYS A 128      20.467  17.512  18.745  1.00 11.73           N
ATOM      2  CA  CYS A 128      19.408  18.379  18.322  1.00 11.51           C
ATOM      3  C   CYS A 128      18.363  17.574  17.583  1.00 11.99           C
ATOM      4  O   CYS A 128      18.674  16.561  16.958  1.00 16.86           O
ATOM      5  OXT CYS A 128      17.169  17.942  17.608  1.00 11.99           O
ATOM      6  CB  CYS A 128      19.968  19.533  17.532  1.00 12.94           C
ATOM      7  SG  CYS A 128      21.122  20.540  18.463  1.00 12.36           S
ATOM      8  N   CYS A 232      19.014  22.068  22.852  1.00  4.75           N
ATOM      9  CA  CYS A 232      20.207  21.244  22.650  1.00  5.50           C
ATOM     10  C   CYS A 232      21.505  22.052  22.672  1.00  5.70           C
ATOM     11  O   CYS A 232      22.560  21.544  23.033  1.00  8.25           O
ATOM     12  OXT CYS A 232      21.489  23.250  22.316  1.00  5.70           O
ATOM     13  CB  CYS A 232      20.087  20.393  21.399  1.00  6.69           C
ATOM     14  SG  CYS A 232      19.969  21.413  19.863  1.00  7.70           S
END
"""

# DG C 3 - DC C 4 (phenix_regression dnatest.pdb), OP3 added on DG 3
pdb_dna = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  OP3  DG C   3      17.763  26.299  18.994  1.00 89.49           O
ATOM      2  P    DG C   3      16.959  25.260  18.230  1.00 89.49           P
ATOM      3  OP1  DG C   3      16.740  25.910  16.911  1.00 83.87           O
ATOM      4  OP2  DG C   3      15.821  25.031  19.147  1.00 79.77           O
ATOM      5  O5'  DG C   3      17.692  23.872  17.966  1.00 77.33           O
ATOM      6  C5'  DG C   3      18.928  23.856  17.277  1.00 70.14           C
ATOM      7  C4'  DG C   3      19.616  22.545  17.586  1.00 77.88           C
ATOM      8  O4'  DG C   3      20.049  22.531  18.967  1.00 79.71           O
ATOM      9  C1'  DG C   3      20.140  21.182  19.370  1.00 71.32           C
ATOM     10  N9   DG C   3      20.102  21.263  20.814  1.00 65.44           N
ATOM     11  C8   DG C   3      19.097  21.804  21.576  1.00 62.30           C
ATOM     12  N7   DG C   3      19.349  21.778  22.859  1.00 64.01           N
ATOM     13  C5   DG C   3      20.602  21.172  22.947  1.00 59.06           C
ATOM     14  C4   DG C   3      21.078  20.857  21.694  1.00 64.05           C
ATOM     15  N3   DG C   3      22.246  20.268  21.364  1.00 67.55           N
ATOM     16  C2   DG C   3      22.967  19.990  22.443  1.00 70.79           C
ATOM     17  N2   DG C   3      24.164  19.403  22.305  1.00 79.26           N
ATOM     18  N1   DG C   3      22.584  20.265  23.731  1.00 57.47           N
ATOM     19  C6   DG C   3      21.387  20.864  24.077  1.00 54.39           C
ATOM     20  O6   DG C   3      21.135  21.062  25.265  1.00 61.49           O
ATOM     21  C2'  DG C   3      18.943  20.509  18.719  1.00 76.50           C
ATOM     22  C3'  DG C   3      18.757  21.291  17.411  1.00 78.46           C
ATOM     23  O3'  DG C   3      19.261  20.615  16.263  1.00 81.31           O
ATOM     24  P    DC C   4      18.722  19.151  15.867  1.00 89.77           P
ATOM     25  OP1  DC C   4      19.149  18.872  14.477  1.00 84.77           O
ATOM     26  OP2  DC C   4      17.278  19.085  16.198  1.00 82.88           O
ATOM     27  O5'  DC C   4      19.551  18.174  16.835  1.00 76.71           O
ATOM     28  C5'  DC C   4      20.958  18.115  16.736  1.00 60.51           C
ATOM     29  C4'  DC C   4      21.509  17.297  17.877  1.00 69.72           C
ATOM     30  O4'  DC C   4      21.038  17.779  19.147  1.00 81.94           O
ATOM     31  C1'  DC C   4      21.297  16.759  20.083  1.00 80.11           C
ATOM     32  N1   DC C   4      20.616  17.186  21.348  1.00 74.12           N
ATOM     33  C2   DC C   4      21.252  17.025  22.589  1.00 70.45           C
ATOM     34  O2   DC C   4      22.364  16.493  22.632  1.00 65.16           O
ATOM     35  N3   DC C   4      20.623  17.454  23.715  1.00 70.95           N
ATOM     36  C4   DC C   4      19.411  18.028  23.626  1.00 72.43           C
ATOM     37  N4   DC C   4      18.818  18.441  24.756  1.00 72.26           N
ATOM     38  C5   DC C   4      18.745  18.209  22.370  1.00 60.90           C
ATOM     39  C6   DC C   4      19.384  17.781  21.272  1.00 67.55           C
ATOM     40  C2'  DC C   4      20.771  15.537  19.354  1.00 91.28           C
ATOM     41  C3'  DC C   4      21.040  15.867  17.877  1.00 88.62           C
ATOM     42  O3'  DC C   4      22.103  15.124  17.323  1.00106.71           O
END
"""

# ATP without PG O1G O2G O3G (2idx, phenix_regression atp.pdb)
pdb_atp = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  PB  ATP A 607      16.663  19.802  21.870  1.00 29.78           P
HETATM    2  O1B ATP A 607      15.778  19.448  20.701  1.00 30.78           O
HETATM    3  O2B ATP A 607      17.432  18.688  22.532  1.00 28.46           O
HETATM    4  O3B ATP A 607      15.706  20.641  22.868  1.00 30.38           O
HETATM    5  PA  ATP A 607      18.775  21.604  22.401  1.00 27.66           P
HETATM    6  O1A ATP A 607      18.321  22.937  22.911  1.00 28.83           O
HETATM    7  O2A ATP A 607      19.279  20.562  23.364  1.00 25.52           O
HETATM    8  O3A ATP A 607      17.652  20.991  21.427  1.00 28.97           O
HETATM    9  O5' ATP A 607      19.936  21.870  21.304  1.00 28.26           O
HETATM   10  C5' ATP A 607      20.059  23.082  20.548  1.00 28.60           C
HETATM   11  C4' ATP A 607      20.159  22.688  19.075  1.00 29.21           C
HETATM   12  O4' ATP A 607      21.269  21.813  18.914  1.00 28.34           O
HETATM   13  C3' ATP A 607      18.972  21.911  18.533  1.00 28.59           C
HETATM   14  O3' ATP A 607      17.988  22.744  17.922  1.00 29.36           O
HETATM   15  C2' ATP A 607      19.571  20.954  17.523  1.00 28.64           C
HETATM   16  O2' ATP A 607      19.701  21.534  16.233  1.00 29.52           O
HETATM   17  C1' ATP A 607      20.969  20.755  18.022  1.00 29.40           C
HETATM   18  N9  ATP A 607      21.029  19.525  18.824  1.00 30.12           N
HETATM   19  C8  ATP A 607      20.095  19.071  19.685  1.00 29.95           C
HETATM   20  N7  ATP A 607      20.484  17.914  20.252  1.00 30.22           N
HETATM   21  C5  ATP A 607      21.712  17.640  19.759  1.00 30.55           C
HETATM   22  C6  ATP A 607      22.703  16.571  19.948  1.00 30.09           C
HETATM   23  N6  ATP A 607      22.451  15.539  20.810  1.00 28.47           N
HETATM   24  N1  ATP A 607      23.840  16.670  19.212  1.00 30.02           N
HETATM   25  C2  ATP A 607      24.118  17.672  18.367  1.00 30.00           C
HETATM   26  N3  ATP A 607      23.264  18.685  18.165  1.00 31.31           N
HETATM   27  C4  ATP A 607      22.069  18.702  18.819  1.00 31.02           C
END
"""

# NAG on ASN ND2 (1a7s A 114 / A 401, no O1)
pdb_nag = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASN A 114      21.891  23.807  18.636  1.00  8.61           N
ATOM      2  CA  ASN A 114      22.241  23.467  19.985  1.00  9.29           C
ATOM      3  C   ASN A 114      23.611  23.952  20.437  1.00 10.96           C
ATOM      4  O   ASN A 114      24.057  23.663  21.528  1.00 15.40           O
ATOM      5  OXT ASN A 114      24.292  24.664  19.668  1.00 10.96           O
ATOM      6  CB  ASN A 114      22.069  21.988  20.335  1.00 11.35           C
ATOM      7  CG  ASN A 114      20.673  21.491  20.260  1.00 11.55           C
ATOM      8  OD1 ASN A 114      19.695  22.273  20.252  1.00 11.66           O
ATOM      9  ND2 ASN A 114      20.496  20.175  20.159  1.00 12.97           N
HETATM   10  C1  NAG A 401      19.225  19.596  20.144  1.00 13.79           C
HETATM   11  C2  NAG A 401      19.279  18.205  19.432  1.00 15.71           C
HETATM   12  C3  NAG A 401      17.827  17.650  19.434  1.00 17.70           C
HETATM   13  C4  NAG A 401      17.366  17.569  20.894  1.00 20.69           C
HETATM   14  C5  NAG A 401      17.395  19.019  21.500  1.00 21.37           C
HETATM   15  C6  NAG A 401      16.912  19.122  22.888  1.00 28.32           C
HETATM   16  C7  NAG A 401      20.719  18.030  17.433  1.00 17.48           C
HETATM   17  C8  NAG A 401      20.832  18.280  16.005  1.00 22.32           C
HETATM   18  N2  NAG A 401      19.575  18.463  18.023  1.00 15.00           N
HETATM   19  O3  NAG A 401      17.868  16.331  18.871  1.00 21.47           O
HETATM   20  O4  NAG A 401      16.001  17.117  20.883  1.00 25.59           O
HETATM   21  O5  NAG A 401      18.811  19.413  21.474  1.00 17.19           O
HETATM   22  O6  NAG A 401      17.687  18.343  23.765  1.00 38.39           O
HETATM   23  O7  NAG A 401      21.484  17.380  17.989  1.00 19.28           O
END
"""

# SO4, ideal tetrahedron
pdb_so4 = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  S   SO4 A   1      20.000  20.000  20.000  1.00 20.00           S
HETATM    2  O1  SO4 A   1      20.849  20.849  20.849  1.00 20.00           O
HETATM    3  O2  SO4 A   1      20.849  19.151  19.151  1.00 20.00           O
HETATM    4  O3  SO4 A   1      19.151  20.849  19.151  1.00 20.00           O
HETATM    5  O4  SO4 A   1      19.151  19.151  20.849  1.00 20.00           O
END
"""

# three carbons in a chain in a residue without a dictionary
pdb_unl = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  UNL A   1      20.000  20.000  20.000  1.00 20.00           C
HETATM    2  C2  UNL A   1      21.530  20.000  20.000  1.00 20.00           C
HETATM    3  C3  UNL A   1      22.103  21.419  20.000  1.00 20.00           C
END
"""

pdb_water_ions = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1      20.000  20.000  20.000  1.00 20.00           O
HETATM    2  ZN   ZN A   2      24.000  20.000  20.000  1.00 20.00          ZN
HETATM    3  CL   CL A   3      20.000  24.000  20.000  1.00 20.00          CL
END
"""

pdb_ala_n_ca = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      25.011  24.832  25.726  1.00  7.86           N
ATOM      2  CA  ALA A   1      24.990  25.168  24.274  1.00  8.28           C
END
"""

pdb_lone_cb = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  CB  ALA A   1      25.000  25.000  25.000  1.00  7.41           C
END
"""

# 4-chloro-1H-pyrazole (CCD JZ6, model coordinates shifted); CLAA has a
# 4-character name
pdb_jz6 = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1 CLAA JZ6 A   1      20.549  23.114  20.060  1.00 20.00          CL
HETATM    2  CAB JZ6 A   1      19.257  20.908  20.981  1.00 20.00           C
HETATM    3  CAC JZ6 A   1      18.757  21.384  18.946  1.00 20.00           C
HETATM    4  NAD JZ6 A   1      18.379  20.025  20.513  1.00 20.00           N
HETATM    5  NAE JZ6 A   1      18.090  20.307  19.340  1.00 20.00           N
HETATM    6  CAF JZ6 A   1      19.480  21.741  19.984  1.00 20.00           C
END
"""

# ------------------------------------------------------------------------------

def model_with_h(pdb_str):
  '''
  Model from pdb_str with H placed as reduce2 places them; the returned model
  carries the restraints manager built in place_hydrogens.run().
  '''
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  obj = reduce_hydrogen.place_hydrogens(
    model                 = model,
    use_neutron_distances = True,
    n_terminal_charge     = 'all_chain_starts',
    exclude_water         = False,
    stop_for_unknowns     = False,
    keep_existing_H       = True,
    mon_lib_srv           = model.get_mon_lib_srv())
  obj.run()
  return obj.get_model()

def records_by_name(result):
  '''
  {(resname, resseq, altloc, name): record}
  '''
  return dict(((r.key[4], r.key[2], r.key[5], r.key[6]), r)
              for r in result.records)

def without(pdb_str, resseq, names):
  '''
  pdb_str without the atoms of residue resseq named in names.
  '''
  keep = []
  for line in pdb_str.splitlines():
    if (line.startswith(('ATOM', 'HETATM')) and
        line[22:26].strip() == str(resseq) and line[12:16].strip() in names):
      continue
    keep.append(line)
  return '\n'.join(keep) + '\n'

def find(pdb_str, guess_for_unknowns=False):
  return capping.find_dangling_atoms(model_with_h(pdb_str),
    guess_for_unknowns=guess_for_unknowns)

def h_names_on(model, resseq, name):
  '''
  Names of the H bonded to atom name of residue resseq.
  '''
  atoms = model.get_hierarchy().atoms()
  table = capping.bond_table(model)
  for a in atoms:
    if (a.parent().parent().resseq.strip() == str(resseq) and
        a.name.strip() == name):
      return sorted(atoms[j].name.strip() for j in table.get(a.i_seq, ())
                    if atoms[j].element_is_hydrogen())
  raise AssertionError('no atom %s %s' % (resseq, name))

# ------------------------------------------------------------------------------

def test_arg_truncated_cb():
  '''
  CB: V 4 (CA CG HB2 HB3), used CA = 1 -> 3 H, tetrahedral. N (NH3+) and
  the peptide C (O and N of SER 4 in place of OXT) are not dangling.
  '''
  result = find(pdb_arg)
  assert len(result.records) == 1, result.records
  r = result.records[0]
  assert r.key == ('', 'A', '3', '', 'ARG', '', 'CB'), r.key
  assert r.n_fill == 3, r.n_fill
  assert r.coordination == 4, r.coordination
  assert r.missing_names == ['CG'], r.missing_names
  assert sorted(r.h_names) == ['HB2', 'HB3'], r.h_names
  assert r.source == 'dictionary', r.source
  assert result.unknown_residues == [] and result.skipped == []

def test_lys_truncated_cd():
  '''
  CD: V 4 (CG CE HD2 HD3), used CG = 1 -> 3. CG keeps its H, not dangling.
  '''
  by_name = records_by_name(find(pdb_lys))
  assert list(by_name) == [('LYS', '2', '', 'CD')], list(by_name)
  r = by_name[('LYS', '2', '', 'CD')]
  assert r.n_fill == 3, r.n_fill
  assert r.missing_names == ['CE'], r.missing_names
  assert r.h_names == ['HD2', 'HD3'], r.h_names

def test_backbone_c_missing_o():
  '''
  C of ALA 1 without O: V 4 (CA, O double, OXT), used CA + N(ALA 2) = 2 -> 2.
  '''
  by_name = records_by_name(find(pdb_ala_ala))
  assert list(by_name) == [('ALA', '1', '', 'C')], list(by_name)
  r = by_name[('ALA', '1', '', 'C')]
  assert r.n_fill == 2, r.n_fill
  assert r.coordination == 4, r.coordination
  assert r.missing_names == ['O'], r.missing_names

def test_tyr_missing_ce1():
  '''
  CCD Kekule ring: CG=CD1, CE1=CZ. CD1: V 4 (CG double, CE1, HD1), used
  CG 2 + HD1 = 3 -> 1. CZ: V 4 (CE1 double, CE2, OH), used 2 -> 2.
  '''
  by_name = records_by_name(find(pdb_tyr))
  assert sorted(by_name) == [('TYR', '5', '', 'CD1'), ('TYR', '5', '', 'CZ')], \
    sorted(by_name)
  cd1 = by_name[('TYR', '5', '', 'CD1')].n_fill
  cz = by_name[('TYR', '5', '', 'CZ')].n_fill
  assert set([cd1, cz]) == set([1, 2]), (cd1, cz)
  assert cd1 + cz == 3, (cd1, cz)

def test_altloc_b_only():
  '''
  Only CB B misses OG: V 4, used CA (blank) = 1 -> 3.
  '''
  result = find(pdb_ser_altloc)
  assert len(result.records) == 1, result.records
  r = result.records[0]
  assert r.key[5] == 'B' and r.key[6] == 'CB', r.key
  assert r.n_fill == 3, r.n_fill

def test_neighbour_altlocs():
  '''
  Blank CB bonded to CA A and CA B: one neighbour in one conformer, so
  coordination 1 + 3 = 4 (not 5).
  '''
  result = find(pdb_lys_ca_altloc)
  assert len(result.records) == 1, result.records
  r = result.records[0]
  assert r.key[5] == '' and r.key[6] == 'CB', r.key
  assert r.n_fill == 3, r.n_fill
  assert r.coordination == 4, r.coordination

def test_isopeptide_not_dangling():
  '''
  ASP CG without OD2, linked to LYS NZ: V 4 (CB, OD1 double, OD2), used
  CB + OD1 2 + NZ = 4 -> 0.
  '''
  model = model_with_h(pdb_isopeptide)
  table = capping.bond_table(model)
  atoms = model.get_hierarchy().atoms()
  cg = [a for a in atoms if a.name.strip() == 'CG' and
        a.parent().resname == 'ASP'][0]
  assert 'NZ' in [atoms[j].name.strip() for j in table[cg.i_seq]]
  result = capping.find_dangling_atoms(model)
  assert result.records == [], records_by_name(result)

def test_metal_not_dangling():
  '''
  CYS SG on ZN has no missing heavy partner. HEM without the methyl CMA:
  C3A: V 4 (C2A, C4A double, CMA), used 3 -> 1; FE is not a partner of the
  pyrrole N, and FE itself is never examined.
  '''
  result = find(pdb_cys_zn)
  assert result.records == [], records_by_name(result)
  result = find(without(pdb_hem, 154, ['CMA']))
  by_name = records_by_name(result)
  assert list(by_name) == [('HEM', '154', '', 'C3A')], list(by_name)
  r = by_name[('HEM', '154', '', 'C3A')]
  assert r.n_fill == 1, r.n_fill
  assert r.coordination == 3, r.coordination
  assert r.missing_names == ['CMA'], r.missing_names
  # HEM without C1A: NA V 2 (C1A, C4A; FE is not a partner), used C4A = 1
  # -> 1, coordination C4A + 1 = 2 (the coordination bond is not counted).
  # C1A is also a partner of C2A and CHA: an internal gap, not capped.
  model = model_with_h(without(pdb_hem, 154, ['C1A']))
  result = capping.find_dangling_atoms(model)
  r = records_by_name(result)[('HEM', '154', '', 'NA')]
  assert (r.n_fill, r.coordination) == (1, 2), (r.n_fill, r.coordination)
  assert r.missing_names == ['C1A'], r.missing_names
  assert [g.missing_name for g in result.internal_gaps] == ['C1A'], \
    result.internal_gaps
  assert_sorry(lambda: capping.cap_dangling_atoms(model, result), 'C1A',
    'Rebuild')
  # without the whole C1A corner (C1A C2A C3A CAA CMA CHA) no gap is left:
  # NA gets one H, oriented about C4A only, FE not counted
  model = model_with_h(without(pdb_hem, 154,
    ['C1A', 'C2A', 'C3A', 'CAA', 'CMA', 'CHA']))
  result = capping.find_dangling_atoms(model)
  assert result.internal_gaps == [], result.internal_gaps
  r = records_by_name(result)[('HEM', '154', '', 'NA')]
  assert (r.n_fill, r.coordination) == (1, 2), (r.n_fill, r.coordination)
  capped = capping.cap_dangling_atoms(model, result)
  added = dict((key[6], names) for key, names in capped.added)
  assert len(added['NA']) == 1, capped.added

def test_disulfide():
  '''
  SG-SG: no missing heavy partner. CYS 232 without SG: CB V 4 (CA SG HB2
  HB3), used CA = 1 -> 3.
  '''
  result = find(pdb_disulfide)
  assert result.records == [], records_by_name(result)
  by_name = records_by_name(find(without(pdb_disulfide, 232, ['SG'])))
  assert list(by_name) == [('CYS', '232', '', 'CB')], list(by_name)
  assert by_name[('CYS', '232', '', 'CB')].n_fill == 3

def test_nucleic_5p_missing_o3p():
  '''
  5' P of DG 3 without OP3: V 5 (OP3, OP1 double, OP2, O5'), used 4 -> 1;
  bonded OP1 OP2 O5' + 1 -> tetrahedral. P of DC 4 (O3' link) is not.
  '''
  by_name = records_by_name(find(without(pdb_dna, 3, ['OP3'])))
  assert list(by_name) == [('DG', '3', '', 'P')], list(by_name)
  r = by_name[('DG', '3', '', 'P')]
  assert r.n_fill == 1, r.n_fill
  assert r.coordination == 4, r.coordination
  assert r.missing_names == ['OP3'], r.missing_names

def test_internal_nucleotide_p():
  '''
  P of DC 4 has no OP3, O3' of DG 3 takes its place: V 5, used 5 -> 0.
  '''
  result = find(pdb_dna)
  assert result.records == [], records_by_name(result)

def test_ccd_ligand():
  '''
  ATP without the gamma phosphate. O3B: V 2 (PB, PG), used PB = 1 -> 1,
  linear. PB has all its partners.
  '''
  result = find(pdb_atp)
  assert len(result.records) == 1, records_by_name(result)
  r = result.records[0]
  assert r.key[6] == 'O3B', r.key
  assert r.n_fill == 1, r.n_fill
  assert r.coordination == 2, r.coordination
  assert r.missing_names == ['PG'], r.missing_names

def test_unknown_residue():
  '''
  No dictionary for UNL: listed; guessed, the bonds come from distances
  (no restraints): C 4, terminal used 1 -> 3, middle used 2 -> 2.
  '''
  model = model_with_h(pdb_unl)
  result = capping.find_dangling_atoms(model)
  assert result.records == [], records_by_name(result)
  assert result.unknown_residues == ['UNL A 1'], result.unknown_residues
  result = capping.find_dangling_atoms(model, guess_for_unknowns=True)
  by_name = records_by_name(result)
  assert sorted(by_name) == [('UNL', '1', '', n) for n in ('C1', 'C2', 'C3')], \
    sorted(by_name)
  for name, n_fill in (('C1', 3), ('C2', 2), ('C3', 3)):
    r = by_name[('UNL', '1', '', name)]
    assert r.source == 'guess', r.source
    assert r.n_fill == n_fill, (name, r.n_fill)
    assert r.coordination == 4, (name, r.coordination)

def test_guessed_bonds_and_extra_sigma():
  '''
  Guessed UNL propane: the two C-C bonds come from distances, not the
  restraints, so electrons gets them as sigma bonds on both ends: with the
  caps every C has 4 and propane is neutral. A dictionary residue has no
  guessed bonds.
  '''
  model = model_with_h(pdb_unl)
  result = capping.find_dangling_atoms(model, guess_for_unknowns=True)
  i1, i2, i3 = [i_seq_of(model, 1, n) for n in ('C1', 'C2', 'C3')]
  assert result.guessed_bonds == sorted(
    [tuple(sorted((i1, i2))), tuple(sorted((i2, i3)))]), result.guessed_bonds
  extra = capping.extra_sigma_bonds(result)
  assert extra == {i1: 4, i2: 4, i3: 4}, extra
  dist = electrons.run(model=model, extra_sigma_bonds=extra,
    return_formal_charges=True)
  assert dist.get_total_charge() == 0, dist.get_total_charge()
  for i_seq in (i1, i2, i3):
    assert dist[i_seq] == 0, (i_seq, dist[i_seq])
  model = model_with_h(pdb_arg)
  result = capping.find_dangling_atoms(model)
  assert result.guessed_bonds == [], result.guessed_bonds
  extra = capping.extra_sigma_bonds(result)
  assert extra == {i_seq_of(model, 3, 'CB'): 3}, extra

def test_glycan():
  '''
  NAG C1 without O1: V 4 (C2 O1 O5 H1), used C2 + O5 + ND2 (link) + H1 -> 0.
  Without C6 and O6: C5 V 4 (C4 C6 O5 H5), used C4 + O5 (+ H5 if placed).
  No second NAG on O4: O4 has no missing heavy partner and gets HO4 from
  placement.
  '''
  model = model_with_h(pdb_nag)
  result = capping.find_dangling_atoms(model)
  assert result.records == [], records_by_name(result)
  assert h_names_on(model, 401, 'O4') == ['HO4']
  model = model_with_h(without(pdb_nag, 401, ['C6', 'O6']))
  by_name = records_by_name(capping.find_dangling_atoms(model))
  assert list(by_name) == [('NAG', '401', '', 'C5')], list(by_name)
  r = by_name[('NAG', '401', '', 'C5')]
  assert r.coordination == 4, r.coordination
  assert r.n_fill + len(h_names_on(model, 401, 'C5')) == 2, r.n_fill

def test_nucleic_sugar_and_base():
  '''
  DC 4 without C2': C1' (O4' N1 H1') and C3' (C4' O3' H3') lose one partner.
  DG 3 without C8: N9 V 3 (C1' C4 C8), used 2 -> 1; N7 V 3 (C8 double,
  C5), used 1 -> 2.
  '''
  by_name = records_by_name(find(without(pdb_dna, 4, ["C2'"])))
  assert sorted(by_name) == [('DC', '4', '', "C1'"), ('DC', '4', '', "C3'")], \
    sorted(by_name)
  by_name = records_by_name(find(without(pdb_dna, 3, ['C8'])))
  assert sorted(by_name) == [('DG', '3', '', 'N7'), ('DG', '3', '', 'N9')], \
    sorted(by_name)
  n7 = by_name[('DG', '3', '', 'N7')].n_fill
  n9 = by_name[('DG', '3', '', 'N9')].n_fill
  assert (n7, n9) == (2, 1), (n7, n9)
  assert n7 + n9 == 3

def test_ala_n_ca_only():
  '''
  ALA with N and CA only. CA: V 4 (N C CB HA), HA not placed, used N = 1 ->
  3, missing C and CB. N has no missing heavy partner.
  '''
  result = find(pdb_ala_n_ca)
  assert len(result.records) == 1, records_by_name(result)
  r = result.records[0]
  assert r.key[6] == 'CA', r.key
  assert r.n_fill == 3, r.n_fill
  assert r.coordination == 4, r.coordination
  assert r.missing_names == ['C', 'CB'], r.missing_names
  assert r.h_names == ['HA'], r.h_names

def test_water_and_ions_untouched():
  result = find(pdb_water_ions)
  assert result.records == [], records_by_name(result)
  assert result.unknown_residues == [] and result.skipped == [], result.skipped

def test_lone_atom_sorry():
  '''
  CB alone: 4 H needed and no neighbour to orient them.
  '''
  model = model_with_h(pdb_lone_cb)
  try:
    capping.find_dangling_atoms(model)
  except Sorry as e:
    assert 'CB' in str(e), str(e)
  else: raise AssertionError('Sorry expected')

def test_skipped_coordination():
  '''
  SO4 (CCD: O1, O2 double, O3, O4 single), S V 6. Without O3: used 5 -> 1,
  coordination 3 + 1 = 4. Without O1: used 4 -> 2, coordination 5: skipped.
  '''
  by_name = records_by_name(find(without(pdb_so4, 1, ['O3'])))
  assert list(by_name) == [('SO4', '1', '', 'S')], list(by_name)
  r = by_name[('SO4', '1', '', 'S')]
  assert r.n_fill == 1, r.n_fill
  assert r.coordination == 4, r.coordination
  result = find(without(pdb_so4, 1, ['O1']))
  assert result.records == [], records_by_name(result)
  assert len(result.skipped) == 1, result.skipped
  label, reason = result.skipped[0]
  assert label.startswith('SO4 A 1') and label.endswith('S'), label
  assert 'coordination 5' in reason, reason

def test_sulfate_refused():
  '''
  SO4 without O3: S is a record (n_fill 1, coordination 4) but hypervalent,
  so capping refuses it. CYS without SG is a CB record, not hypervalent, and
  a thioether/thiol S (2 substituents) is not either.
  '''
  model = model_with_h(without(pdb_so4, 1, ['O3']))
  result = capping.find_dangling_atoms(model)
  assert [r.key[6] for r in result.hypervalent] == ['S'], result.hypervalent
  msg = assert_sorry(lambda: capping.cap_dangling_atoms(model, result),
    'Cannot cap', 'SO4 A 1', 'O3', 'sulfur hydride')
  assert "OP3 at a 5' end" not in msg, msg
  result = find(without(pdb_cys_zn, 48, ['SG']))
  assert result.hypervalent == [], result.hypervalent
  atom = iotbx.pdb.hierarchy.atom()
  atom.element = ' S'
  assert not capping._is_hypervalent(atom, 2)
  assert capping._is_hypervalent(atom, 3)
  atom.element = ' P'
  assert capping._is_hypervalent(atom, 2)
  atom.element = ' C'
  assert not capping._is_hypervalent(atom, 4)

# ------------------------------------------------------------------------------

def test_cut_bonds_by_distance_conformers():
  '''
  Conformer A of pdb_ser_altloc selected (with the blank atoms): the
  unselected B atoms sit at bonding distance of the blank backbone but are
  conformers of selected residue groups, not cut partners.
  '''
  model = model_with_h(pdb_ser_altloc)
  h = model.get_hierarchy()
  sel = h.atom_selection_cache().selection("altloc A or altloc ' '")
  assert sel.count(False) > 0
  assert capping.cut_bonds_by_distance(h, sel) == {}, \
    capping.cut_bonds_by_distance(h, sel)

def test_cap_model_matches_two_step():
  '''
  cap_model on pdb_arg (ARG 3 at CB): same caps and charge as find +
  extra_sigma_bonds + electrons + cap_dangling_atoms by hand.
  '''
  a, b = model_with_h(pdb_arg), model_with_h(pdb_arg)
  dangling = capping.find_dangling_atoms(a)
  charge = electrons.run(model=a,
    extra_sigma_bonds=capping.extra_sigma_bonds(dangling))
  capping.cap_dangling_atoms(a, dangling)
  r = capping.cap_model(b)
  assert r.total_charge == charge, (r.total_charge, charge)
  assert r.capped.n_added == 3 and r.warnings == [], \
    (r.capped.n_added, r.warnings)
  assert [x.name for x in a.get_hierarchy().atoms()] == \
    [x.name for x in b.get_hierarchy().atoms()]
  assert a.get_hierarchy().atoms().extract_xyz().max_distance(
    b.get_hierarchy().atoms().extract_xyz()) < 1.e-6

def test_cap_model_cached_dangling_unprocessed():
  '''
  The records of a processed model placed on an unprocessed manager of a
  copy of its hierarchy (no restraints): identical caps, no charge.
  '''
  model = model_with_h(pdb_arg)
  dangling = capping.find_dangling_atoms(model)
  copy = mmtbx.model.manager(model_input=None,
    pdb_hierarchy=model.get_hierarchy().deep_copy(),
    crystal_symmetry=model.crystal_symmetry(), log=null_out())
  r = capping.cap_model(copy, dangling=dangling, calculate_charge=False,
    mon_lib_srv=model.get_mon_lib_srv())
  ref = capping.cap_dangling_atoms(model, dangling)
  assert r.total_charge is None and r.formal_charges is None
  assert r.capped.lines == ref.lines, (r.capped.lines, ref.lines)
  a, b = model.get_hierarchy().atoms(), copy.get_hierarchy().atoms()
  assert list(a.extract_name()) == list(b.extract_name())
  assert a.extract_xyz().max_distance(b.extract_xyz()) < 1.e-6

def test_cap_model_unknown_residue_sorry():
  try:
    capping.cap_model(model_with_h(pdb_unl))
  except Sorry as e:
    assert 'capping=True: no dictionary for' in str(e), str(e)
  else:
    raise AssertionError('no Sorry')


# ------------------------------------------------------------------------------
# Geometry

ORIGIN = matrix.col((0., 0., 0.))
TOLERANCE = 1.0 # degrees

def bond_to(direction, length=1.5):
  '''
  Site of a known neighbour of the origin along direction.
  '''
  return matrix.col(direction).normalize() * length

def angle_between(a, b):
  return a.angle(b, deg=True)

def assert_angle(a, b, expected):
  angle = angle_between(a, b)
  assert abs(angle - expected) < TOLERANCE, (angle, expected)

def min_distance(directions, sites):
  '''
  Smallest distance from an H at 1 A along any of directions (parent at the
  origin) to any of sites.
  '''
  return min((site - d).length() for d in directions for site in sites)

def environment_of(*sites):
  return flex.vec3_double([tuple(s) for s in sites])

def no_environment():
  return flex.vec3_double()

def test_steric_number():
  cases = [('N', 2, 3), ('N', 3, 4), ('O', 2, 4), ('S', 2, 4), ('S', 4, 4),
           ('C', 2, 2), ('C', 3, 3), ('P', 4, 4)]
  for element, coordination, expected in cases:
    assert capping.steric_number(element, coordination) == expected, \
      (element, coordination)
  assert (capping.TETRAHEDRAL, capping.TRIGONAL, capping.LINEAR) == (4, 3, 2)

def test_tetra_three_known():
  known = [bond_to(d) for d in ((1,0,0.35), (-0.5,0.87,0.35), (-0.5,-0.87,0.35))]
  result = capping.fill_directions(ORIGIN, known, capping.TETRAHEDRAL, 1,
                                   no_environment())
  assert len(result) == 1
  assert abs(result[0].length() - 1) < 1e-6
  for k in known: assert_angle(result[0], k, 109.5)

def tetra_two_known():
  c = math.cos(math.radians(109.5/2))
  s = math.sin(math.radians(109.5/2))
  return [bond_to((c, s, 0)), bond_to((c, -s, 0))]

def test_tetra_two_known():
  known = tetra_two_known()
  result = capping.fill_directions(ORIGIN, known, capping.TETRAHEDRAL, 2,
                                   no_environment())
  assert len(result) == 2
  assert_angle(result[0], result[1], 109.5)
  for d in result:
    assert abs(d.length() - 1) < 1e-6
    for k in known: assert_angle(d, k, 109.5)

def test_tetra_two_known_one_fill():
  known = tetra_two_known()
  both = capping.fill_directions(ORIGIN, known, capping.TETRAHEDRAL, 2,
                                 no_environment())
  # an environment atom close to the first position: the other one is chosen
  for near, far in ((both[0], both[1]), (both[1], both[0])):
    env = environment_of(near * 2.5)
    result = capping.fill_directions(ORIGIN, known, capping.TETRAHEDRAL, 1,
                                     env)
    assert len(result) == 1
    assert (result[0] - far).length() < 1e-6
    for k in known: assert_angle(result[0], k, 109.5)

def test_tetra_one_known_rotation():
  known = [bond_to((0, 0, 1))]
  axis = -known[0].normalize()
  eclipsed = capping.fill_directions(ORIGIN, known, capping.TETRAHEDRAL, 3,
                                     no_environment())
  assert len(eclipsed) == 3
  for d in eclipsed: assert_angle(d, known[0], 109.5)
  # a neighbour of the known atom, exactly eclipsed with the first H
  env = environment_of(eclipsed[0] * 2.5)
  result = capping.fill_directions(ORIGIN, known, capping.TETRAHEDRAL, 3, env)
  assert len(result) == 3
  for i, d in enumerate(result):
    assert abs(d.length() - 1) < 1e-6
    assert_angle(d, known[0], 109.5)
    for e in result[i+1:]: assert_angle(d, e, 109.5)
  sites = [eclipsed[0] * 2.5]
  staggered = [d.rotate_around_origin(axis, 60, deg=True) for d in eclipsed]
  assert min_distance(result, sites) > min_distance(eclipsed, sites) + 0.1
  assert abs(min_distance(result, sites) - min_distance(staggered, sites)) \
    < 1e-3

def test_bent_one_known():
  known = [bond_to((0, 0, 1))]
  result = capping.fill_directions(ORIGIN, known, capping.TETRAHEDRAL, 1,
                                   no_environment())
  assert len(result) == 1
  assert_angle(result[0], known[0], 109.5)

def test_amine_one_known():
  known = [bond_to((0, 0, 1))]
  result = capping.fill_directions(ORIGIN, known, capping.TETRAHEDRAL, 2,
                                   environment_of(matrix.col((1., 1., -2.))))
  assert len(result) == 2
  assert_angle(result[0], result[1], 109.5)
  for d in result: assert_angle(d, known[0], 109.5)

def test_trigonal_two_known():
  known = [bond_to((1, 0, 0)), bond_to((-0.5, 0.866, 0))]
  result = capping.fill_directions(ORIGIN, known, capping.TRIGONAL, 1,
                                   no_environment())
  assert len(result) == 1
  expected = -(known[0].normalize() + known[1].normalize()).normalize()
  assert (result[0] - expected).length() < 1e-6

def test_trigonal_one_known():
  known = [bond_to((0, 0, 1))]
  result = capping.fill_directions(ORIGIN, known, capping.TRIGONAL, 2,
                                   no_environment())
  assert len(result) == 2
  assert_angle(result[0], result[1], 120)
  for d in result: assert_angle(d, known[0], 120)

def test_trigonal_one_known_one_fill():
  known = [bond_to((0, 0, 1))]
  both = capping.fill_directions(ORIGIN, known, capping.TRIGONAL, 2,
                                 no_environment())
  # the single H goes to the position farther from the environment atom
  env = environment_of(both[0] * 2.5)
  result = capping.fill_directions(ORIGIN, known, capping.TRIGONAL, 1, env)
  assert len(result) == 1
  assert_angle(result[0], known[0], 120)
  assert (result[0] - both[1]).length() < 1e-6

def test_linear():
  known = [bond_to((0.3, -0.4, 0.8))]
  result = capping.fill_directions(ORIGIN, known, capping.LINEAR, 1,
                                   no_environment())
  assert len(result) == 1
  assert_angle(result[0], known[0], 180)

def test_no_known_raises():
  try:
    capping.fill_directions(ORIGIN, [], capping.TETRAHEDRAL, 1,
                            no_environment())
  except ValueError:
    return
  raise AssertionError('no ValueError for an empty known list')

# ------------------------------------------------------------------------------
# cut bonds (rule 9): a fragment cut out of a larger model

def i_seq_of(model, resseq, name):
  for a in model.get_hierarchy().atoms():
    if (a.parent().parent().resseq.strip() == str(resseq) and
        a.name.strip() == name):
      return a.i_seq
  raise AssertionError((resseq, name))

def test_cut_bonds_backbone_n():
  '''
  ALA 2 of pdb_ala_ala cut from ALA 1: N 2 (CA, H) lost its peptide C.
  cut_bonds gives it one H, named H2 (its free dictionary H); without
  cut_bonds it is not examined (no heavy dictionary partner missing).
  '''
  model = model_with_h(pdb_ala_ala)
  model = model.select(model.selection('resseq 2'))
  n = i_seq_of(model, 2, 'N')
  assert ('ALA', '2', '', 'N') not in records_by_name(
    capping.find_dangling_atoms(model))
  result = capping.find_dangling_atoms(model, cut_bonds={n: 1})
  r = records_by_name(result)[('ALA', '2', '', 'N')]
  assert (r.n_fill, r.coordination, r.missing_names, r.h_names) == \
    (1, 3, [], ['H2']), (r.n_fill, r.coordination, r.missing_names, r.h_names)
  assert capping.extra_sigma_bonds(result)[n] == 1
  capping.cap_dangling_atoms(model, result)
  atoms = model.get_hierarchy().atoms()
  n_atom = atoms[i_seq_of(model, 2, 'N')]
  near = sorted(a.name.strip() for a in atoms
                if a.element_is_hydrogen() and a.distance(n_atom) < 1.2)
  assert near == ['H', 'H2'], near

def test_cut_bonds_keeps_dictionary_record():
  '''
  ALA 1 cut from ALA 2: C 1 (no O in this fixture) misses O and OXT
  (leaving, nothing replaces it now): n_fill 3 already counts the absent N;
  a cut bond on it changes nothing and adds no record.
  '''
  model = model_with_h(pdb_ala_ala)
  model = model.select(model.selection('resseq 1'))
  c = i_seq_of(model, 1, 'C')
  plain = records_by_name(capping.find_dangling_atoms(model))[
    ('ALA', '1', '', 'C')]
  result = capping.find_dangling_atoms(model, cut_bonds={c: 1})
  assert len([r for r in result.records if r.i_seq == c]) == 1
  assert records_by_name(result)[('ALA', '1', '', 'C')].n_fill == \
    plain.n_fill == 3

def test_cut_bonds_ignored_on_h_and_metal():
  model = model_with_h(pdb_cys_zn)
  atoms = model.get_hierarchy().atoms()
  zn = [a.i_seq for a in atoms if a.element.strip().upper() == 'ZN'][0]
  h = [a.i_seq for a in atoms if a.element_is_hydrogen()][0]
  result = capping.find_dangling_atoms(model, cut_bonds={zn: 1, h: 1})
  assert not [r for r in result.records if r.i_seq in (zn, h)], \
    records_by_name(result)

def test_cut_bonds_by_distance_disulfide():
  '''
  CYS 232 of pdb_disulfide selected alone: its SG lost SG 128.
  '''
  model = model_with_h(pdb_disulfide)
  h = model.get_hierarchy()
  sel = h.atom_selection_cache().selection('resseq 232')
  cut = capping.cut_bonds_by_distance(h, sel)
  frag = h.select(sel)
  frag.atoms().reset_i_seq()
  sg = [a.i_seq for a in frag.atoms() if a.name.strip() == 'SG'][0]
  assert cut == {sg: 1}, cut
  # the whole model: nothing is cut
  assert capping.cut_bonds_by_distance(
    h, flex.bool(h.atoms().size(), True)) == {}

def test_cut_bonds_by_distance_metal():
  '''
  pdb_cys_zn without its ZN: the coordination bond is no cut bond.
  '''
  model = model_with_h(pdb_cys_zn)
  h = model.get_hierarchy()
  sel = ~h.atom_selection_cache().selection('element ZN')
  assert capping.cut_bonds_by_distance(h, sel) == {}


# ------------------------------------------------------------------------------
# Building the caps

XH_TOLERANCE = 0.01 # A
ANGLE_TOLERANCE = 2.0 # degrees

# pdb_arg renumbered: ARG A 1, SER A 2
pdb_arg_1 = pdb_arg.replace('ARG A   3', 'ARG A   1').replace(
  'SER A   4', 'SER A   2')

def cap(pdb_str):
  '''
  (model, cap_dangling_atoms result) for pdb_str with H placed.
  '''
  model = model_with_h(pdb_str)
  result = capping.cap_dangling_atoms(model,
    capping.find_dangling_atoms(model))
  return model, result

def atom_named(model, resseq, name, chain_id='A', altloc=''):
  for a in model.get_hierarchy().atoms():
    ag = a.parent()
    if (ag.parent().parent().id.strip() == chain_id and
        ag.parent().resseq.strip() == str(resseq) and
        ag.altloc.strip() == altloc and a.name.strip() == name):
      return a
  raise AssertionError('no atom %s %s %s' % (chain_id, resseq, name))

def site(atom):
  return matrix.col(atom.xyz)

def bond_angle(a, b, c):
  '''
  Angle a-b-c in degrees.
  '''
  return (site(a) - site(b)).angle(site(c) - site(b), deg=True)

def written(model):
  '''
  {(chain, resseq, name with padding): atom} of model_as_pdb() re-read.
  '''
  pdb_inp = iotbx.pdb.input(lines=model.model_as_pdb().split('\n'),
    source_info=None)
  hierarchy = pdb_inp.construct_hierarchy()
  result = {}
  for a in hierarchy.atoms():
    rg = a.parent().parent()
    result[(rg.parent().id, rg.resseq.strip(), a.name)] = a
  return result

def assert_sorry(function, *words):
  try:
    function()
  except Sorry as e:
    for word in words: assert word in str(e), (word, str(e))
    return str(e)
  raise AssertionError('Sorry expected')

def test_cap_arg_cb():
  '''
  Truncated ARG CB: 3 H, HB2 HB3 from the dictionary, HCG for CG; C-H
  neutron length; altloc, occupancy and B of CB.
  '''
  model, result = cap(pdb_arg_1)
  key = ('', 'A', '1', '', 'ARG', '', 'CB')
  assert result.added == [(key, ['HB2', 'HB3', 'HCG'])], result.added
  assert (result.n_added, result.n_atoms, result.n_residues) == (3, 1, 1)
  assert result.warnings == []
  assert len(result.lines) == 1, result.lines
  line = result.lines[0]
  assert line.startswith('chain A ARG 1: CB +3 H'), line
  assert line == ('chain A ARG 1: CB +3 H (HB2 HB3 HCG); missing CG CD NE '
    'CZ NH1 NH2'), line
  cb = atom_named(model, 1, 'CB')
  xh = reduce_hydrogen.get_library_neutron_xh()['C']
  assert abs(xh - 1.09) < XH_TOLERANCE, xh
  hs = [atom_named(model, 1, n) for n in ('HB2', 'HB3', 'HCG')]
  ca = atom_named(model, 1, 'CA')
  for h in hs:
    assert h.element.strip() == 'H', h.element
    assert abs((site(h) - site(cb)).length() - xh) < XH_TOLERANCE
    assert abs(bond_angle(ca, cb, h) - 109.5) < ANGLE_TOLERANCE
    assert h.parent().altloc == cb.parent().altloc
    assert (h.occ, h.b, h.segid, h.hetero) == (cb.occ, cb.b, cb.segid,
                                                cb.hetero)
  # i_seq reset; the caps written and read back with padded names
  atoms = model.get_hierarchy().atoms()
  assert list(atoms.extract_i_seq()) == list(range(atoms.size()))
  out = written(model)
  for name in (' HB2', ' HB3', ' HCG'):
    assert ('A', '1', name) in out, sorted(out)
  # an altloc is one more field
  model, result = cap(pdb_ser_altloc)
  assert result.lines == ['chain A SER 4 B: CB +3 H (HB2 HB3 HOG); missing OG'],\
    result.lines

def test_cap_c_missing_o():
  '''
  Backbone C of ALA 1 without O: 2 H (HO, HC1), coordination 4, tetrahedral
  about CA and N of ALA 2.
  '''
  model, result = cap(pdb_ala_ala)
  assert len(result.added) == 1, result.added
  key, names = result.added[0]
  assert key[6] == 'C' and names == ['HO', 'HC1'], result.added
  c, ca = atom_named(model, 1, 'C'), atom_named(model, 1, 'CA')
  n = atom_named(model, 2, 'N')
  hs = [atom_named(model, 1, name) for name in names]
  for h in hs:
    assert abs((site(h) - site(c)).length() - 1.09) < XH_TOLERANCE
    # the CA-C-N angle is 116, so the angles to CA and N are 107.7
    for other in (ca, n):
      assert abs(bond_angle(other, c, h) - 109.5) < 3, bond_angle(other, c, h)
  assert abs(bond_angle(hs[0], c, hs[1]) - 109.5) < 3
  out = written(model)
  assert ('A', '1', ' HO ') in out and ('A', '1', ' HC1') in out, sorted(out)

def test_cap_5p_refused():
  '''
  5' P of DG 3 without OP3 is detected (n_fill 1) but hypervalent: a P-H
  would be a phosphonate. cap_dangling_atoms refuses it with a Sorry naming
  the atom, the missing OP3 and the remedy; nothing is added.
  '''
  model = model_with_h(without(pdb_dna, 3, ['OP3']))
  result = capping.find_dangling_atoms(model)
  assert [r.key[6] for r in result.records] == ['P'], records_by_name(result)
  assert [r.key[6] for r in result.hypervalent] == ['P']
  n_atoms = model.get_number_of_atoms()
  msg = assert_sorry(lambda: capping.cap_dangling_atoms(model, result),
    'Cannot cap', 'DG C 3', 'OP3', 'phosphonate', "OP3 at a 5' end")
  assert model.get_number_of_atoms() == n_atoms
  assert ('C', '3', 'HOP3') not in written(model)
  # library lengths by distance flag; C with 3 substituents is aromatic C-H
  expected = reduce_hydrogen._covalent_radius('P') + 0.31
  assert capping.xh_length('P', 4, False) == expected
  assert capping.xh_length('C', 3, False) == \
    reduce_hydrogen.get_library_xray_xh()['C_aromatic']
  assert capping.xh_length('C', 4, True) == \
    reduce_hydrogen.get_library_neutron_xh()['C']

def test_cap_names_parent():
  '''
  An O, N or S parent names its cap after itself (HO5' on O5', HO3' on O3',
  HNZ on NZ) when no dictionary H name is free; a C parent after the missing
  atom (HO3' on a C3' without O3', HCG on a CB); a taken name falls through
  to the next candidate.
  '''
  keep = [] # the residue groups: an atom group's parent is a weak link
  def group(parent_name, element, taken=()):
    rg = iotbx.pdb.hierarchy.residue_group()
    keep.append(rg)
    ag = iotbx.pdb.hierarchy.atom_group(resname='DC')
    rg.append_atom_group(ag)
    for name, el in [(parent_name, element)] + [(t, ' H') for t in taken]:
      a = iotbx.pdb.hierarchy.atom()
      a.name = name
      a.element = el
      ag.append_atom(a)
    return ag, ag.atoms()[0]
  ag, o5 = group(" O5'", ' O')
  assert capping.cap_names(ag, [], ['P'], 1, parent=o5) == ["HO5'"]
  assert capping.cap_names(ag, [], ['P'], 1) == ['HP']
  ag, o3 = group(" O3'", ' O')
  assert capping.cap_names(ag, [], ['P'], 1, parent=o3) == ["HO3'"]
  ag, nz = group(' NZ ', ' N')
  assert capping.cap_names(ag, [], ['CE'], 1, parent=nz) == ['HNZ']
  assert capping.cap_names(ag, ['HZ1'], ['CE'], 2, parent=nz) == ['HZ1', 'HNZ']
  ag, c3 = group(" C3'", ' C')
  assert capping.cap_names(ag, [], ["O3'"], 1, parent=c3) == ["HO3'"]
  ag, o5 = group(" O5'", ' O', taken=["HO5'"])
  assert capping.cap_names(ag, [], ['P'], 1, parent=o5) == ['HP']

def test_cap_o3b_bent():
  '''
  ATP without the gamma phosphate: the H on O3B is bent (O has 2 lone
  pairs), 109.5 from PB, not 180.
  '''
  model, result = cap(pdb_atp)
  assert len(result.added) == 1, result.added
  key, names = result.added[0]
  # an O parent is named after itself (hydroxyl convention), not after PG
  assert key[6] == 'O3B' and names == ['HO3B'], result.added
  o3b, pb = atom_named(model, 607, 'O3B'), atom_named(model, 607, 'PB')
  h = atom_named(model, 607, 'HO3B')
  xh = reduce_hydrogen.get_library_neutron_xh()['O']
  assert abs((site(h) - site(o3b)).length() - xh) < XH_TOLERANCE
  angle = bond_angle(pb, o3b, h)
  assert abs(angle - 109.5) < ANGLE_TOLERANCE, angle
  assert h.hetero, 'cap of a HETATM is HETATM'

def test_two_chains():
  '''
  Chains A and B, each with a truncated CB: both capped.
  '''
  lines = [l for l in pdb_arg_1.splitlines() if l.startswith('ATOM')]
  chain_b = []
  for l in lines:
    x = float(l[30:38]) + 10.
    chain_b.append(l[:21] + 'B' + l[22:30] + '%8.3f' % x + l[38:])
  pdb_str = '\n'.join([pdb_arg_1.splitlines()[0]] + lines + chain_b) + '\nEND\n'
  model, result = cap(pdb_str)
  keys = [k for k, names in result.added]
  assert keys == [('', 'A', '1', '', 'ARG', '', 'CB'),
                  ('', 'B', '1', '', 'ARG', '', 'CB')], keys
  assert (result.n_atoms, result.n_residues, result.n_added) == (2, 2, 6)
  assert result.lines[1].startswith('chain B ARG 1: CB +3 H'), result.lines
  for chain_id in ('A', 'B'):
    atom_named(model, 1, 'HCG', chain_id)

def test_names_partially_present():
  '''
  HB2 kept from the input on a CB missing CG: 2 H, HB3 and HCG, tetrahedral
  about CA and HB2. reduce2 drops an H it cannot ride (HB2 without CG, also
  with keep_existing_H=True), so the model is processed as read.
  '''
  model, result = cap(pdb_arg_1)
  hb2 = atom_named(model, 1, 'HB2')
  line = ('ATOM      6  HB2 ARG A   1    %8.3f%8.3f%8.3f  1.00  6.79'
          '           H' % hb2.xyz)
  pdb_str = pdb_arg_1.replace('ATOM      6  N   SER', line +
    '\nATOM      6  N   SER')
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split('\n'), source_info=None)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  found = capping.find_dangling_atoms(model)
  assert len(found.records) == 1 and found.records[0].n_fill == 2, \
    records_by_name(found)
  result = capping.cap_dangling_atoms(model, found)
  key, names = result.added[0]
  assert set(names) == set(['HB3', 'HCG']), names
  cb, ca = atom_named(model, 1, 'CB'), atom_named(model, 1, 'CA')
  hb2 = atom_named(model, 1, 'HB2')
  for name in names:
    h = atom_named(model, 1, name)
    for other in (ca, hb2):
      angle = bond_angle(other, cb, h)
      assert abs(angle - 109.5) < ANGLE_TOLERANCE, (name, angle)

def test_name_fallback():
  '''
  JZ6 without CLAA: 'HCLAA' is too long, the H on CAF is HC1; trigonal.
  '''
  model, result = cap(without(pdb_jz6, 1, ['CLAA']))
  assert len(result.added) == 1, result.added
  key, names = result.added[0]
  assert key[6] == 'CAF' and names == ['HC1'], result.added
  caf, h = atom_named(model, 1, 'CAF'), atom_named(model, 1, 'HC1')
  xh = reduce_hydrogen.get_library_neutron_xh()['C_aromatic']
  assert abs((site(h) - site(caf)).length() - xh) < XH_TOLERANCE
  for name in ('CAB', 'CAC'):
    angle = bond_angle(atom_named(model, 1, name), caf, h)
    assert abs(angle - 126) < 3, (name, angle) # in a 5-ring, 360-108 / 2
  assert result.lines[0].endswith('; missing CLAA'), result.lines
  assert ('A', '1', ' HC1') in written(model)

def jz6_with_cl(shift=0., cell=None):
  '''
  JZ6 without CLAA and a CL ion 0.8 A beyond the only position of the H on
  CAF, shifted by shift along x; in a cubic P1 cell of edge cell if given.
  '''
  pdb_str = without(pdb_jz6, 1, ['CLAA'])
  model, result = cap(pdb_str)
  caf, h = atom_named(model, 1, 'CAF'), atom_named(model, 1, 'HC1')
  d = (site(h) - site(caf)).normalize()
  cl = site(h) + d * 0.8 + matrix.col((shift, 0, 0))
  pdb_str = pdb_str.replace('END', 'HETATM    7 CL    CL A   2    '
    '%8.3f%8.3f%8.3f  1.00 20.00          CL\nEND' % tuple(cl))
  if cell is not None:
    pdb_str = pdb_str.replace('CRYST1   40.000   40.000   40.000',
      'CRYST1 %8.3f %8.3f %8.3f' % (cell, cell, cell))
  return pdb_str

def test_clash_sorry():
  '''
  The only position of the H (2 known neighbours, trigonal) is 0.8 A from a
  CL ion: Sorry naming both.
  '''
  model = model_with_h(jz6_with_cl())
  found = capping.find_dangling_atoms(model)
  assert len(found.records) == 1, records_by_name(found)
  message = assert_sorry(lambda: capping.cap_dangling_atoms(model, found),
    'Capping H', 'HC1', 'CL', '0.80 A')
  assert 'JZ6' in message, message

def test_symmetry_clash_sorry_cap_second():
  '''
  As test_symmetry_clash_sorry with the CL before JZ6 in the file, so the
  cap is the second atom of the clashing pair: the operator printed is still
  the one that brings the CL to the cap.
  '''
  cell = 8.
  lines = jz6_with_cl(shift=-cell, cell=cell).splitlines()
  cl = [l for l in lines if l[17:20] == ' CL']
  rest = [l for l in lines if l[17:20] != ' CL']
  pdb_str = '\n'.join(rest[:1] + cl + rest[1:]) + '\n'
  model = model_with_h(pdb_str)
  found = capping.find_dangling_atoms(model)
  assert len(found.records) == 1, records_by_name(found)
  message = assert_sorry(lambda: capping.cap_dangling_atoms(model, found),
    'Capping H', 'HC1', 'CL', '0.80 A', 'symmetry mate x+1,y,z')
  assert message.index('HC1') < message.index('CL '), message

def test_symmetry_clash_sorry():
  '''
  As test_clash_sorry, the CL one cell edge away: its symmetry mate is 0.8 A
  from the H.
  '''
  cell = 8.
  model = model_with_h(jz6_with_cl(shift=-cell, cell=cell))
  found = capping.find_dangling_atoms(model)
  assert len(found.records) == 1, records_by_name(found)
  assert_sorry(lambda: capping.cap_dangling_atoms(model, found),
    'Capping H', 'HC1', 'CL', '0.80 A', 'symmetry mate x+1,y,z')

def test_keys_survive_sort():
  '''
  The parents' atom_key is found again after sorting a copy of the capped
  hierarchy; the caps stay in the parents' atom groups.
  '''
  model, result = cap(pdb_arg_1)
  hierarchy = model.get_hierarchy().deep_copy()
  hierarchy.sort_atoms_in_place()
  by_key = dict((capping.atom_key(a), a) for a in hierarchy.atoms())
  for key, names in result.added:
    assert key in by_key, key
    ag_names = [a.name.strip() for a in by_key[key].parent().atoms()]
    for name in names: assert name in ag_names, (name, ag_names)

def moved(pdb_str, R, t):
  '''
  pdb_str with every atom moved by x -> R * (x - c) + c + t, c the centre of
  the 50 A cell.
  '''
  c = matrix.col((25, 25, 25))
  out = []
  for l in pdb_str.splitlines():
    if l.startswith(('ATOM', 'HETATM')):
      x = matrix.col((float(l[30:38]), float(l[38:46]), float(l[46:54])))
      l = l[:30] + '%8.3f%8.3f%8.3f' % tuple(R * (x - c) + c + t) + l[54:]
    out.append(l)
  return '\n'.join(out) + '\n'

def dihedral(model, resseq, names):
  from scitbx.math import dihedral_angle
  return dihedral_angle(sites=[atom_named(model, resseq, n).xyz
    for n in names], deg=True)

def test_rotation_frame_independent():
  '''
  ARG truncated at CB in four rigid-body frames: the CH3 rotation is the
  same (the bonded CA does not take part in the score) and staggered about
  CA, not eclipsed with N or C.
  '''
  frames = [(matrix.sqr((1,0,0, 0,1,0, 0,0,1)), (0, 0, 0))]
  for k, axis in enumerate([(1,0,0), (0,1,0), (0.3,0.5,0.8)]):
    R = matrix.col(axis).normalize().axis_and_angle_as_r3_rotation_matrix(
      37. * (k+1), deg=True)
    frames.append((R, (1.5 * k, -2., 0.7)))
  angles = []
  for R, t in frames:
    model, result = cap(moved(pdb_arg, R, matrix.col(t)))
    angles.append(dihedral(model, 3, ('N', 'CA', 'CB', 'HCG')))
  for a in angles:
    d = abs((a - angles[0] + 180.) % 360. - 180.)
    assert d < 1., angles
  staggered = min(abs((angles[0] - s + 180.) % 360. - 180.)
                  for s in (60., 180., 300.))
  assert staggered < 15., angles

def test_two_models():
  '''
  ARG truncated at CB in MODEL 1 and MODEL 2 (rotated 3 degrees): the caps
  of one MODEL are not checked against, nor oriented by, the other.
  '''
  lines = pdb_arg.splitlines()
  cryst = [l for l in lines if l.startswith('CRYST1')]
  atoms = [l for l in lines if l.startswith('ATOM')]
  R = matrix.col((0, 0, 1)).axis_and_angle_as_r3_rotation_matrix(3., deg=True)
  atoms2 = moved('\n'.join(atoms), R, matrix.col((0, 0, 0))).splitlines()
  pdb_str = '\n'.join(cryst + ['MODEL        1'] + atoms + ['ENDMDL',
    'MODEL        2'] + atoms2 + ['ENDMDL', 'END']) + '\n'
  model, result = cap(pdb_str)
  assert [(k[0], k[6], names) for k, names in result.added] == [
    ('1', 'CB', ['HB2', 'HB3', 'HCG']), ('2', 'CB', ['HB2', 'HB3', 'HCG'])], \
    result.added
  assert result.n_added == 6, result.n_added
  # the same rotation in both MODELs (it would differ if MODEL 2 counted)
  model_1 = cap(pdb_arg)[0]
  a = dihedral(model_1, 3, ('N', 'CA', 'CB', 'HCG'))
  for m in model.get_hierarchy().models():
    sites = {}
    for atom in m.atoms():
      if atom.parent().parent().resseq.strip() == '3':
        sites[atom.name.strip()] = atom.xyz
    from scitbx.math import dihedral_angle
    b = dihedral_angle(sites=[sites[n] for n in ('N', 'CA', 'CB', 'HCG')],
      deg=True)
    assert abs((a - b + 180.) % 360. - 180.) < 1., (a, b)

def gap_parents(gap):
  return sorted(k[6] for k in gap.parent_keys)

def test_internal_gap_detected():
  '''
  A missing heavy atom that two present dangling atoms of the residue miss is
  an internal gap: TYR without CE1 (CD1, CZ), ARG without CG only (CB, CD).
  The records are still reported. A truncation leaves no gap.
  '''
  result = find(pdb_tyr)
  assert len(result.records) == 2, records_by_name(result)
  assert len(result.internal_gaps) == 1, result.internal_gaps
  gap = result.internal_gaps[0]
  assert gap.missing_name == 'CE1', gap.missing_name
  assert gap_parents(gap) == ['CD1', 'CZ'], gap_parents(gap)
  assert gap.residue_key == ('', 'A', '5', '', 'TYR'), gap.residue_key
  result = find(without(pdb_arg_full, 3, ['CG']))
  assert [(g.missing_name, gap_parents(g)) for g in result.internal_gaps] == \
    [('CG', ['CB', 'CD'])], result.internal_gaps
  for pdb_str in (pdb_arg, pdb_lys, pdb_ala_ala, pdb_ser_altloc):
    result = find(pdb_str)
    assert result.records, records_by_name(result)
    assert result.internal_gaps == [], result.internal_gaps

def test_cap_internal_gap_sorry():
  '''
  cap_dangling_atoms refuses an internal gap before placing anything: one
  line per gap, at most 20.
  '''
  model = model_with_h(pdb_tyr)
  found = capping.find_dangling_atoms(model)
  n_atoms = model.get_hierarchy().atoms_size()
  message = assert_sorry(lambda: capping.cap_dangling_atoms(model, found),
    'Cannot cap chain A TYR 5: atom CE1 is missing but its neighbours CD1 '
    'and CZ are present; capping both would overlap. Rebuild the missing '
    'atom (or remove the fragment) before running capping.')
  assert model.get_hierarchy().atoms_size() == n_atoms
  found.internal_gaps = found.internal_gaps * 25
  message = assert_sorry(lambda: capping.cap_dangling_atoms(model, found),
    '... and 5 more')
  assert message.count('Cannot cap') == 20, message

def test_bond_table_reused():
  '''
  The detection result carries the bond table; capping does not build it
  again.
  '''
  model = model_with_h(pdb_arg)
  found = capping.find_dangling_atoms(model)
  assert found.bond_table == capping.bond_table(model)
  build = capping.bond_table
  def fail(model): raise AssertionError('bond table built again')
  capping.bond_table = fail
  try:
    result = capping.cap_dangling_atoms(model, found)
  finally:
    capping.bond_table = build
  assert result.n_added == 3, result.added

def test_coordination_mismatch():
  '''
  A record whose coordination disagrees with the neighbours found when
  building: RuntimeError naming the atom and the counts, not a bare assert.
  '''
  model = model_with_h(pdb_arg)
  found = capping.find_dangling_atoms(model)
  found.records[0].coordination = 5
  try:
    capping.cap_dangling_atoms(model, found)
  except RuntimeError as e:
    for word in ('ARG A 3 CB', '1 neighbours', '3 H', 'coordination 5'):
      assert word in str(e), (word, str(e))
  else:
    raise AssertionError('RuntimeError expected')

def test_no_ccd_warning():
  '''
  CCD lookup failing for a residue the restraint dictionary knows: one log
  line per residue name.
  '''
  from six.moves import cStringIO as StringIO
  import mmtbx.chemical_components
  lines = [l for l in pdb_arg.splitlines() if l.startswith('ATOM')]
  chain_b = [l[:21] + 'B' + l[22:30] + '%8.3f' % (float(l[30:38]) + 10.) +
             l[38:] for l in lines]
  model = model_with_h('\n'.join([pdb_arg.splitlines()[0]] + lines + chain_b)
    + '\nEND\n')
  get = mmtbx.chemical_components.get_cif_dictionary
  def fail(code, *args, **kwds): raise RuntimeError('no CCD')
  mmtbx.chemical_components.get_cif_dictionary = fail
  log = StringIO()
  try:
    capping.find_dangling_atoms(model, log=log)
  finally:
    mmtbx.chemical_components.get_cif_dictionary = get
  text = log.getvalue()
  warned = [l for l in text.splitlines() if 'restraint dictionary' in l]
  assert len(warned) == 2, text
  assert len([l for l in warned if ' ARG' in l]) == 1, text
  assert len([l for l in warned if ' SER' in l]) == 1, text

def test_extra_sigma_bonds_unknown_i_seq():
  model = model_with_h(pdb_lys)
  assert_sorry(lambda: charge(model, {100000: 1}), 'extra_sigma_bonds',
    '100000')

# ------------------------------------------------------------------------------
# Charge
#
# electrons on the uncapped model with extra_sigma_bonds = {i_seq: n_fill}
# gives the charge of the capped molecule (one electron per cap H, one bond
# consumed on each side).

# LYS A 2 complete (pdb_lys with CE, NZ); same source
pdb_lys_full = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   LYS A   2      27.221  24.554  25.195  1.00  5.98           N
ATOM      2  CA  LYS A   2      25.851  24.605  24.700  1.00  6.43           C
ATOM      3  C   LYS A   2      25.466  26.033  24.327  1.00  6.30           C
ATOM      4  O   LYS A   2      25.212  26.865  25.198  1.00  6.91           O
ATOM      5  OXT LYS A   2      25.408  26.362  23.122  1.00  6.30           O
ATOM      6  CB  LYS A   2      24.884  24.064  25.755  1.00  7.86           C
ATOM      7  CG  LYS A   2      23.430  24.044  25.312  1.00 11.05           C
ATOM      8  CD  LYS A   2      22.527  23.473  26.393  1.00 16.83           C
ATOM      9  CE  LYS A   2      21.073  23.454  25.950  1.00 26.78           C
ATOM     10  NZ  LYS A   2      20.176  22.897  27.000  1.00 31.27           N
END
"""

# ARG A 3 complete, OXT at the N of SER 4; same source as pdb_arg
pdb_arg_full = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A   3      27.911  23.777  26.214  1.00  5.99           N
ATOM      2  CA  ARG A   3      27.558  25.104  25.723  1.00  6.90           C
ATOM      3  C   ARG A   3      26.047  25.238  25.567  1.00  6.28           C
ATOM      4  O   ARG A   3      25.450  26.196  26.055  1.00  7.24           O
ATOM      5  OXT ARG A   3      25.443  24.270  24.884  1.00  6.28           O
ATOM      6  CB  ARG A   3      28.246  25.375  24.384  1.00  6.79           C
ATOM      7  CG  ARG A   3      29.765  25.350  24.452  1.00  8.29           C
ATOM      8  CD  ARG A   3      30.397  25.644  23.099  1.00  8.22           C
ATOM      9  NE  ARG A   3      30.051  24.639  22.094  1.00  9.42           N
ATOM     10  CZ  ARG A   3      29.061  24.750  21.208  1.00 13.50           C
ATOM     11  NH1 ARG A   3      28.285  25.829  21.172  1.00 15.22           N
ATOM     12  NH2 ARG A   3      28.845  23.768  20.345  1.00 12.90           N
END
"""

# ASP A 1 complete, OXT at the N of GLU 2
# (mmtbx/regression/tst_reduce2_qmready_charge.py)
pdb_asp_full = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A   1       9.580   7.432   8.699  1.00 37.74           N
ATOM      2  CA  ASP A   1      10.692   8.369   8.630  1.00 38.33           C
ATOM      3  C   ASP A   1      10.828   9.231   9.884  1.00 37.17           C
ATOM      4  O   ASP A   1      11.938   9.511  10.322  1.00 36.50           O
ATOM      5  OXT ASP A   1       9.691   9.552  10.506  1.00 37.17           O
ATOM      6  CB  ASP A   1      10.625   9.197   7.354  1.00 39.43           C
ATOM      7  CG  ASP A   1      11.281   8.467   6.194  1.00 43.72           C
ATOM      8  OD1 ASP A   1      10.549   7.831   5.406  1.00 49.92           O
ATOM      9  OD2 ASP A   1      12.531   8.382   6.067  1.00 51.10           O
END
"""

# SER A 4 alone (from pdb_arg)
pdb_ser_full = without(pdb_arg, 3, ['N', 'CA', 'C', 'O', 'CB'])

def i_seq_of(model, resseq, name):
  return atom_named(model, resseq, name).i_seq

def charge(model, extra_sigma_bonds=None):
  return electrons.run(model=model, extra_sigma_bonds=extra_sigma_bonds)

def check_charge(full, resseq, cut_at, cut_names, full_charge):
  '''
  Complete residue: full_charge; truncated after cut_at with its n_fill
  pre-consumed: neutral.
  '''
  model = model_with_h(full)
  assert charge(model) == full_charge, charge(model)
  model = model_with_h(without(full, resseq, cut_names))
  extra = {i_seq_of(model, resseq, cut_at): 3}
  assert charge(model, extra) == 0, charge(model, extra)

def test_charge_lys():
  check_charge(pdb_lys_full, 2, 'CD', ['CE', 'NZ'], 1)

def test_charge_asp():
  check_charge(pdb_asp_full, 1, 'CB', ['CG', 'OD1', 'OD2'], -1)

def test_charge_arg():
  check_charge(pdb_arg_full, 3, 'CB',
    ['CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'], 1)

def test_charge_ser():
  check_charge(pdb_ser_full, 4, 'CB', ['OG'], 0)

def test_charge_parity():
  '''
  Leftover electrons on the parent: 0 with all its missing bonds filled,
  1 (open shell) with one short.
  '''
  model = model_with_h(pdb_lys)
  i_cd = i_seq_of(model, 2, 'CD')
  dist = electrons.run(model=model, extra_sigma_bonds={i_cd: 3},
    return_formal_charges=True)
  assert dist[i_cd] == 0, dist[i_cd]
  assert dist.get_total_charge() == 0, dist.get_total_charge()
  dist = electrons.run(model=model, extra_sigma_bonds={i_cd: 2},
    return_formal_charges=True)
  assert dist[i_cd] == 1, dist[i_cd]

def run():
  test_arg_truncated_cb()
  test_lys_truncated_cd()
  test_backbone_c_missing_o()
  test_tyr_missing_ce1()
  test_altloc_b_only()
  test_neighbour_altlocs()
  test_isopeptide_not_dangling()
  test_metal_not_dangling()
  test_disulfide()
  test_nucleic_5p_missing_o3p()
  test_internal_nucleotide_p()
  test_ccd_ligand()
  test_unknown_residue()
  test_guessed_bonds_and_extra_sigma()
  test_glycan()
  test_nucleic_sugar_and_base()
  test_ala_n_ca_only()
  test_water_and_ions_untouched()
  test_lone_atom_sorry()
  test_skipped_coordination()
  test_cut_bonds_backbone_n()
  test_cut_bonds_keeps_dictionary_record()
  test_cut_bonds_ignored_on_h_and_metal()
  test_cut_bonds_by_distance_disulfide()
  test_cut_bonds_by_distance_metal()
  test_cut_bonds_by_distance_conformers()
  test_cap_model_matches_two_step()
  test_cap_model_cached_dangling_unprocessed()
  test_cap_model_unknown_residue_sorry()
  test_steric_number()
  test_tetra_three_known()
  test_tetra_two_known()
  test_tetra_two_known_one_fill()
  test_tetra_one_known_rotation()
  test_bent_one_known()
  test_amine_one_known()
  test_trigonal_two_known()
  test_trigonal_one_known()
  test_trigonal_one_known_one_fill()
  test_linear()
  test_no_known_raises()
  test_cap_arg_cb()
  test_cap_c_missing_o()
  test_cap_5p_refused()
  test_sulfate_refused()
  test_cap_names_parent()
  test_cap_o3b_bent()
  test_two_chains()
  test_names_partially_present()
  test_name_fallback()
  test_clash_sorry()
  test_symmetry_clash_sorry()
  test_symmetry_clash_sorry_cap_second()
  test_keys_survive_sort()
  test_rotation_frame_independent()
  test_two_models()
  test_internal_gap_detected()
  test_cap_internal_gap_sorry()
  test_bond_table_reused()
  test_coordination_mismatch()
  test_no_ccd_warning()
  test_extra_sigma_bonds_unknown_i_seq()
  test_charge_lys()
  test_charge_asp()
  test_charge_arg()
  test_charge_ser()
  test_charge_parity()

if __name__ == '__main__':
  run()
  print('OK')
