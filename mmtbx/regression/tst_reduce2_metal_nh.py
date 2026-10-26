"""reduce2: no H on a backbone amide N (with a preceding residue), a
thymine/uracil N3, a sulfonamide N or an azolium N bound to a metal (A5; a
primary sulfonamide keeps one H). Alkali and alkaline-earth metals do not count;
other N-H next to a metal are kept and listed; His ring N are left to the
Optimizer. Same machinery as thiols (tst_reduce2_metal_thiol.py): window,
symmetry, occupancy coverage, supplied H kept.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.hydrogens.reduce_hydrogen import residue_is_amino_acid
from mmtbx.monomer_library import server

# 1xmk: the ATCUN-type Ni site (Gly 288 chain start, Ser 289 and His 290
# amide N, His 290 ND1)
pdb_1xmk = """\
CRYST1   35.557   43.526   45.471  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   GLY A 288      12.334  23.080  30.892  1.00 20.00           N
ATOM      2  CA  GLY A 288      11.259  22.629  29.968  1.00 20.00           C
ATOM      3  C   GLY A 288      11.734  21.310  29.346  1.00 20.00           C
ATOM      4  O   GLY A 288      11.046  20.536  28.686  1.00 20.00           O
ATOM      5  N   SER A 289      12.992  21.082  29.673  1.00 20.00           N
ATOM      6  CA  SER A 289      13.766  20.002  29.117  1.00 20.00           C
ATOM      7  C   SER A 289      14.970  19.728  30.022  1.00 20.00           C
ATOM      8  O   SER A 289      15.660  18.700  29.826  1.00 20.00           O
ATOM      9  CB  SER A 289      14.288  20.329  27.695  1.00 20.00           C
ATOM     10  OG  SER A 289      14.983  21.589  27.736  1.00 20.00           O
ATOM     11  N   HIS A 290      15.185  20.681  30.963  1.00 20.00           N
ATOM     12  CA  HIS A 290      16.475  20.629  31.642  1.00 20.00           C
ATOM     13  C   HIS A 290      16.332  19.771  32.904  1.00 20.00           C
ATOM     14  O   HIS A 290      16.691  20.326  33.953  1.00 20.00           O
ATOM     15  CB  HIS A 290      17.030  22.019  31.962  1.00 20.00           C
ATOM     16  CG  HIS A 290      16.054  22.940  32.675  1.00 20.00           C
ATOM     17  ND1 HIS A 290      14.740  23.013  32.337  1.00 20.00           N
ATOM     18  CD2 HIS A 290      16.315  23.854  33.689  1.00 20.00           C
ATOM     19  CE1 HIS A 290      14.147  23.928  33.126  1.00 20.00           C
ATOM     20  NE2 HIS A 290      15.081  24.445  33.946  1.00 20.00           N
TER
HETATM   21 NI    NI A 398      13.849  21.950  30.990  1.00 20.00          NI
END
"""

# the same with supplied amide H on Ser 289 and His 290 (0.86 A, in plane)
pdb_1xmk_h = pdb_1xmk.replace(
  "ATOM     11  N   HIS A 290", "ATOM     30  H   SER A 289      13.367  21.603"
  "  30.245  1.00 20.00           H\nATOM     11  N   HIS A 290").replace(
  "TER\n", "ATOM     31  H   HIS A 290      14.596  21.280  31.148  1.00 20.00"
  "           H\nTER\n")

# 5gsk: thymine N3 on Hg (occupancy set to 1.00); the Hg of DT 3 is a
# symmetry mate (y,x,-z)
pdb_5gsk = """\
CRYST1   40.318   40.318   24.980  90.00  90.00  90.00 P 43 21 2
ATOM      1  P    DT A   3      10.011  19.045   4.108  1.00 20.00           P
ATOM      2  OP1  DT A   3      10.235  19.681   5.425  1.00 20.00           O
ATOM      3  OP2  DT A   3      10.906  17.959   3.643  1.00 20.00           O
ATOM      4  O5'  DT A   3       8.532  18.470   4.070  1.00 20.00           O
ATOM      5  C5'  DT A   3       7.452  19.242   4.612  1.00 20.00           C
ATOM      6  C4'  DT A   3       6.153  18.502   4.416  1.00 20.00           C
ATOM      7  O4'  DT A   3       5.899  18.337   3.006  1.00 20.00           O
ATOM      8  C3'  DT A   3       6.130  17.086   4.976  1.00 20.00           C
ATOM      9  O3'  DT A   3       5.798  17.096   6.359  1.00 20.00           O
ATOM     10  C2'  DT A   3       5.039  16.439   4.154  1.00 20.00           C
ATOM     11  C1'  DT A   3       5.228  17.096   2.794  1.00 20.00           C
ATOM     12  N1   DT A   3       6.037  16.289   1.874  1.00 20.00           N
ATOM     13  C2   DT A   3       5.382  15.302   1.178  1.00 20.00           C
ATOM     14  O2   DT A   3       4.185  15.095   1.285  1.00 20.00           O
ATOM     15  N3   DT A   3       6.183  14.567   0.341  1.00 20.00           N
ATOM     16  C4   DT A   3       7.540  14.721   0.135  1.00 20.00           C
ATOM     17  O4   DT A   3       8.121  13.998  -0.668  1.00 20.00           O
ATOM     18  C5   DT A   3       8.169  15.755   0.922  1.00 20.00           C
ATOM     19  C7   DT A   3       9.642  15.978   0.784  1.00 20.00           C
ATOM     20  C6   DT A   3       7.398  16.472   1.749  1.00 20.00           C
ATOM     21  P    DT A   6       6.997   3.938   7.325  1.00 20.00           P
ATOM     22  OP1  DT A   6       6.310   2.987   8.233  1.00 20.00           O
ATOM     23  OP2  DT A   6       8.123   4.774   7.822  1.00 20.00           O
ATOM     24  O5'  DT A   6       7.536   3.150   6.055  1.00 20.00           O
ATOM     25  C5'  DT A   6       6.736   2.142   5.437  1.00 20.00           C
ATOM     26  C4'  DT A   6       7.417   1.644   4.186  1.00 20.00           C
ATOM     27  O4'  DT A   6       7.577   2.730   3.251  1.00 20.00           O
ATOM     28  C3'  DT A   6       8.829   1.104   4.367  1.00 20.00           C
ATOM     29  O3'  DT A   6       8.835  -0.259   4.763  1.00 20.00           O
ATOM     30  C2'  DT A   6       9.378   1.210   2.965  1.00 20.00           C
ATOM     31  C1'  DT A   6       8.766   2.516   2.499  1.00 20.00           C
ATOM     32  N1   DT A   6       9.650   3.669   2.695  1.00 20.00           N
ATOM     33  C2   DT A   6      10.612   3.857   1.735  1.00 20.00           C
ATOM     34  O2   DT A   6      10.749   3.117   0.774  1.00 20.00           O
ATOM     35  N3   DT A   6      11.407   4.956   1.931  1.00 20.00           N
ATOM     36  C4   DT A   6      11.343   5.859   2.972  1.00 20.00           C
ATOM     37  O4   DT A   6      12.114   6.813   3.000  1.00 20.00           O
ATOM     38  C5   DT A   6      10.328   5.588   3.961  1.00 20.00           C
ATOM     39  C7   DT A   6      10.197   6.503   5.136  1.00 20.00           C
ATOM     40  C6   DT A   6       9.544   4.516   3.778  1.00 20.00           C
TER
HETATM   41 HG    HG A 101      12.924   5.399   0.656  1.00 20.00          HG
END
"""

# the 1xmk excerpt with the Gly 288 C (and O) in altloc A only, a shared
# Ser 289 N, and the rest of Ser 289 in A and B: H/B has no preceding C in its
# conformer
pdb_1xmk_altloc = """\
CRYST1   35.557   43.526   45.471  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   GLY A 288      12.334  23.080  30.892  1.00 20.00           N
ATOM      2  CA  GLY A 288      11.259  22.629  29.968  1.00 20.00           C
ATOM      3  C  AGLY A 288      11.734  21.310  29.346  0.50 20.00           C
ATOM      4  O  AGLY A 288      11.046  20.536  28.686  0.50 20.00           O
ATOM      5  N   SER A 289      12.992  21.082  29.673  1.00 20.00           N
ATOM      6  CA ASER A 289      13.766  20.002  29.117  0.50 20.00           C
ATOM      7  C  ASER A 289      14.970  19.728  30.022  0.50 20.00           C
ATOM      8  O  ASER A 289      15.660  18.700  29.826  0.50 20.00           O
ATOM      9  CB ASER A 289      14.288  20.329  27.695  0.50 20.00           C
ATOM     10  OG ASER A 289      14.983  21.589  27.736  0.50 20.00           O
ATOM     11  CA BSER A 289      13.766  20.002  29.117  0.50 20.00           C
ATOM     12  C  BSER A 289      14.970  19.728  30.022  0.50 20.00           C
ATOM     13  O  BSER A 289      15.660  18.700  29.826  0.50 20.00           O
ATOM     14  CB BSER A 289      14.388  20.129  27.695  0.50 20.00           C
ATOM     15  OG BSER A 289      13.583  21.089  27.036  0.50 20.00           O
ATOM     16  N   HIS A 290      15.185  20.681  30.963  1.00 20.00           N
ATOM     17  CA  HIS A 290      16.475  20.629  31.642  1.00 20.00           C
ATOM     18  C   HIS A 290      16.332  19.771  32.904  1.00 20.00           C
ATOM     19  O   HIS A 290      16.691  20.326  33.953  1.00 20.00           O
ATOM     20  CB  HIS A 290      17.030  22.019  31.962  1.00 20.00           C
ATOM     21  CG  HIS A 290      16.054  22.940  32.675  1.00 20.00           C
ATOM     22  ND1 HIS A 290      14.740  23.013  32.337  1.00 20.00           N
ATOM     23  CD2 HIS A 290      16.315  23.854  33.689  1.00 20.00           C
ATOM     24  CE1 HIS A 290      14.147  23.928  33.126  1.00 20.00           C
ATOM     25  NE2 HIS A 290      15.081  24.445  33.946  1.00 20.00           N
TER
HETATM   26 NI    NI A 398      13.849  21.950  30.990  1.00 20.00          NI
END
"""

# the same with the Gly 288 C shared and no O atoms at all: the peptide bond
# needs no observed carbonyl O
pdb_1xmk_altloc_no_o = "".join(l + "\n" for l in pdb_1xmk_altloc.splitlines()
  if not (l.startswith("ATOM") and l[76:78].strip() == "O")).replace(
  "ATOM      3  C  AGLY A 288      11.734  21.310  29.346  0.50",
  "ATOM      3  C   GLY A 288      11.734  21.310  29.346  1.00")

# a C-terminal amide cap: ligand NH2 (not an amino acid), its N bonded to the
# C of Gly 2, a fully occupied Ni 1.82 A from it
pdb_nh2_cap = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1       5.000   5.000   5.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       6.458   5.000   5.000  1.00 20.00           C
ATOM      3  C   ALA A   1       7.009   6.420   5.000  1.00 20.00           C
ATOM      4  O   ALA A   1       6.252   7.390   5.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       6.986   4.235   6.210  1.00 20.00           C
ATOM      6  N   GLY A   2       8.330   6.540   5.000  1.00 20.00           N
ATOM      7  CA  GLY A   2       8.960   7.850   5.000  1.00 20.00           C
ATOM      8  C   GLY A   2      10.470   7.740   5.000  1.00 20.00           C
ATOM      9  O   GLY A   2      11.050   6.650   5.000  1.00 20.00           O
HETATM   10  N   NH2 A   3      11.130   8.890   5.000  1.00 20.00           N
HETATM   11 NI    NI A 101      10.580  10.620   5.000  1.00 20.00          NI
END
"""

# 6BCC (neutron, carbonic anhydrase II): ethoxzolamide N1, a primary
# sulfonamide, on Zn 301 (1.99 A); one D on N1
pdb_6bcc = """\
CRYST1   42.800   41.700   72.900  90.00 104.60  90.00 P 1 21 1
HETATM 5202 ZN    ZN A 301      -6.880  -1.575  15.379  1.00 15.07          ZN
HETATM 5203  C1  EZL A 302      -5.300   2.313  16.112  1.00 15.81           C
HETATM 5204  C2  EZL A 302      -3.701   3.574  15.345  1.00 14.28           C
HETATM 5205  C3  EZL A 302      -4.808   4.300  14.866  1.00 16.55           C
HETATM 5206  C4  EZL A 302      -4.670   5.465  14.149  1.00 17.59           C
HETATM 5207  C5  EZL A 302      -3.421   5.946  13.861  1.00 20.50           C
HETATM 5208  C6  EZL A 302      -2.315   5.246  14.318  1.00 17.78           C
HETATM 5209  C7  EZL A 302      -2.425   4.065  15.050  1.00 17.32           C
HETATM 5210  C8  EZL A 302      -2.065   7.833  13.190  1.00 29.65           C
HETATM 5211  C9  EZL A 302      -1.991   8.420  14.582  1.00 30.59           C
HETATM 5213  N1  EZL A 302      -5.523  -0.404  16.246  1.00 15.25           N
HETATM 5214  N2  EZL A 302      -3.990   2.451  16.038  1.00 15.10           N
HETATM 5215  O1  EZL A 302      -7.553   1.144  16.775  1.00 16.21           O
HETATM 5216  O2  EZL A 302      -5.653   1.000  18.390  1.00 14.69           O
HETATM 5217  O3  EZL A 302      -3.284   7.129  13.109  1.00 24.39           O
HETATM 5218  S1  EZL A 302      -6.072   0.956  16.948  1.00 15.13           S
HETATM 5219  S2  EZL A 302      -6.250   3.503  15.347  1.00 16.10           S
END
"""

# 6GCY (neutron): 3QR NAZ, an N between C=O and SO2, on Zn 301 (2.02 A);
# no D on NAZ
pdb_6gcy = """\
CRYST1   42.610   41.940   72.770  90.00 104.01  90.00 P 1 21 1
HETATM 4636 ZN    ZN A 301      -2.779  -0.535  85.951  1.00  5.96          ZN2+
HETATM 4637  C1  3QR A 302      -0.792  11.401  81.717  1.00 62.07           C
HETATM 4638  C2  3QR A 302      -1.791  12.513  81.356  1.00 74.03           C
HETATM 4639  C3  3QR A 302      -1.028  13.642  80.752  1.00 85.79           C
HETATM 4640  C4  3QR A 302      -0.211  14.376  81.807  1.00 84.41           C
HETATM 4641  C5  3QR A 302       0.849  13.346  82.141  1.00 78.66           C
HETATM 4642  C6  3QR A 302       1.483  13.568  83.446  1.00 38.53           C
HETATM 4643  CAB 3QR A 302      -1.104   5.661  85.259  1.00  9.65           C
HETATM 4644  CAC 3QR A 302      -0.815   7.816  84.385  1.00 29.16           C
HETATM 4645  CAD 3QR A 302      -1.292   9.075  84.724  1.00 20.49           C
HETATM 4646  CAG 3QR A 302      -1.010  10.407  83.995  1.00 21.67           C
HETATM 4647  CAT 3QR A 302      -0.357   5.058  84.241  1.00  8.99           C
HETATM 4648  CAU 3QR A 302      -0.184   3.743  84.213  1.00 11.84           C
HETATM 4649  CAV 3QR A 302      -0.792   3.085  85.168  1.00  8.15           C
HETATM 4650  CAW 3QR A 302      -1.594   3.623  86.240  1.00  8.21           C
HETATM 4651  CAX 3QR A 302      -1.752   4.982  86.239  1.00  7.01           C
HETATM 4652  CAY 3QR A 302      -0.756   1.491  85.307  1.00  8.10           C
HETATM 4657  NAA 3QR A 302      -1.296   7.049  85.266  1.00 13.55           N
HETATM 4658  NAE 3QR A 302      -2.061   8.998  85.731  1.00 50.17           N
HETATM 4659  NAF 3QR A 302      -2.115   7.676  86.078  1.00 15.46           N
HETATM 4660  NAZ 3QR A 302      -1.621   1.064  86.382  1.00  5.14           N
HETATM 4661  O1  3QR A 302      -1.097  10.291  82.561  1.00 66.97           O
HETATM 4662  O2  3QR A 302      -2.753  12.282  80.314  1.00 68.22           O
HETATM 4663  O3  3QR A 302      -2.062  14.438  80.053  1.00 78.63           O
HETATM 4664  O4  3QR A 302       0.527  15.586  81.518  1.00 77.65           O
HETATM 4665  O5  3QR A 302       0.354  11.982  82.368  1.00 70.10           O
HETATM 4666  O6  3QR A 302       1.628  12.210  83.813  1.00 52.51           O
HETATM 4667  OBB 3QR A 302      -3.761   2.243  87.332  1.00  5.60           O
HETATM 4668  OBC 3QR A 302      -1.452   2.450  88.667  1.00  5.30           O
HETATM 4669  OBD 3QR A 302      -0.300   0.714  84.636  1.00  8.68           O
HETATM 4670  SBA 3QR A 302      -2.201   2.296  87.335  1.00  4.88           S
END
"""

# 4PE3: imidazole A 402 N1 on Zn 407 (2.02 A; Zn occupancy 0.96 in 4PE3,
# set to 1.00)
pdb_4pe3 = """\
CRYST1   79.072   97.938   35.788  90.00  90.00  90.00 P 21 21 2
HETATM 4915  C2  IMD A 402       5.879  90.895  49.538  1.00 12.43           C
HETATM 4916  C4  IMD A 402       8.044  91.155  49.329  1.00 14.02           C
HETATM 4917  C5  IMD A 402       7.529  90.793  48.091  1.00 12.30           C
HETATM 4918  N1  IMD A 402       6.185  90.618  48.249  1.00 11.41           N
HETATM 4919  N3  IMD A 402       7.022  91.225  50.200  1.00 14.77           N
HETATM 4943 ZN    ZN A 407       4.830  90.129  46.839  1.00 10.75          ZN
END
"""

# 6MFL: OPV A 402 NAM on Fe 401 (2.14 A), a neutral imidazole with H2 on
# NAM and none on NAO
pdb_6mfl = """\
CRYST1   37.762  137.055   56.208  90.00  90.02  90.00 P 1 21 1
HETATM 8826 FE    FE A 401      -1.011  -6.350 -10.509  1.00 14.46          FE
HETATM 8827  N   OPV A 402       1.335  -2.597 -10.595  1.00 15.78           N
HETATM 8828  CA  OPV A 402       0.907  -1.218 -10.699  1.00 16.68           C
HETATM 8829  C   OPV A 402      -0.504  -1.058 -10.087  1.00 13.05           C
HETATM 8830  O   OPV A 402      -0.892  -1.690  -9.083  1.00 15.83           O
HETATM 8831  CB  OPV A 402       0.517  -0.678 -12.125  1.00 19.64           C
HETATM 8832  OG1 OPV A 402      -0.486   0.294 -11.814  1.00 18.45           O
HETATM 8833  CG2 OPV A 402       1.659  -0.033 -12.954  1.00 16.15           C
HETATM 8834  CAF OPV A 402       4.189  -6.615  -9.638  1.00 16.56           C
HETATM 8835  CAG OPV A 402       3.113  -7.498  -9.853  1.00 14.03           C
HETATM 8836  CAH OPV A 402       4.032  -5.221  -9.762  1.00 19.27           C
HETATM 8837  CAI OPV A 402      -3.447  -4.494 -10.052  1.00 16.73           C
HETATM 8838  CAJ OPV A 402      -2.195  -3.635 -11.659  1.00 17.79           C
HETATM 8839  CAK OPV A 402      -3.309  -1.232 -12.004  1.00 12.10           C
HETATM 8840  CAL OPV A 402      -2.824  -0.154 -10.995  1.00 14.39           C
HETATM 8841  CAQ OPV A 402       2.562  -3.174 -10.208  1.00 14.72           C
HETATM 8842  CAS OPV A 402       1.836  -6.991 -10.167  1.00 13.34           C
HETATM 8843  CAT OPV A 402      -3.148  -2.639 -11.376  1.00 17.67           C
HETATM 8844  CAU OPV A 402       1.669  -5.580 -10.293  1.00 13.11           C
HETATM 8845  CAV OPV A 402       2.779  -4.676 -10.110  1.00 15.06           C
HETATM 8846  NAM OPV A 402      -2.395  -4.759 -10.851  1.00 13.91           N
HETATM 8847  NAO OPV A 402      -3.913  -3.204 -10.335  1.00 11.35           N
HETATM 8848  NAY OPV A 402      -1.301  -0.128 -10.770  1.00 13.65           N
HETATM 8849  OAB OPV A 402       3.454  -2.378  -9.858  1.00 17.85           O
HETATM 8850  OAD OPV A 402       0.751  -7.789 -10.367  1.00 12.77           O
HETATM 8851  OAE OPV A 402       0.459  -5.040 -10.656  1.00 10.90           O
END
"""

# Trp (from 1r27) with a Zn 2.10 A from NE1 along the N-H direction
pdb_trp = """\
ATOM      1  N   TRP A   1       3.153 138.508 109.939  1.00 20.00           N
ATOM      2  CA  TRP A   1       3.837 138.040 108.735  1.00 20.00           C
ATOM      3  C   TRP A   1       3.699 136.515 108.549  1.00 20.00           C
ATOM      4  O   TRP A   1       3.729 136.015 107.423  1.00 20.00           O
ATOM      5  CB  TRP A   1       5.325 138.409 108.784  1.00 20.00           C
ATOM      6  CG  TRP A   1       6.132 137.483 109.640  1.00 20.00           C
ATOM      7  CD1 TRP A   1       6.330 137.566 110.993  1.00 20.00           C
ATOM      8  CD2 TRP A   1       6.773 136.276 109.213  1.00 20.00           C
ATOM      9  NE1 TRP A   1       7.052 136.479 111.432  1.00 20.00           N
ATOM     10  CE2 TRP A   1       7.337 135.671 110.362  1.00 20.00           C
ATOM     11  CE3 TRP A   1       6.923 135.641 107.969  1.00 20.00           C
ATOM     12  CZ2 TRP A   1       8.040 134.462 110.304  1.00 20.00           C
ATOM     13  CZ3 TRP A   1       7.622 134.437 107.914  1.00 20.00           C
ATOM     14  CH2 TRP A   1       8.172 133.861 109.076  1.00 20.00           C
TER
HETATM   15 ZN    ZN B 101       7.624 136.117 113.420  1.00 20.00          ZN
END
"""

NOT_PLACED = ('N-H not placed: backbone amide N, thymine/uracil N3, '
  'sulfonamide or azolium N bound to a metal, taken as deprotonated')
UNRESOLVED = 'N-H placed at metal sites with unresolved protonation'
SUPPLIED = 'Supplied N-H at metal sites'
REPORTED = 'N-H next to a metal, not changed'

def run(prefix, pdb_str, args=()):
  '''reduce2 on pdb_str. :return: (output hierarchy, description, log).'''
  fn, out, txt = prefix + '_in.pdb', prefix + '_out.pdb', prefix + '.txt'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  for n in (out, txt):
    if os.path.exists(n): os.remove(n)
  log = StringIO()
  run_program(program_class=reduce2.Program, logger=log, args=[fn,
    'output.filename=' + out, 'output.description_file_name=' + txt]
    + list(args))
  with open(txt) as f:
    description = f.read()
  return iotbx.pdb.input(out).construct_hierarchy(), description, log.getvalue()

def h_on(hierarchy, resseq, name):
  '''Names of the H/D within 1.2 A of atom name in residue resseq.'''
  atoms = list(hierarchy.atoms())
  parent = [a for a in atoms if a.parent().parent().resseq_as_int() == resseq
            and a.name.strip() == name]
  assert len(parent) == 1, (resseq, name)
  return sorted(a.name.strip() for a in atoms
    if a.element.strip() in ('H', 'D') and a.distance(parent[0]) < 1.2)

def section(text, title):
  '''The lines of a report section (after its title, up to a line not
  indented or the Optimizer's timing lines, which follow it in the
  description).'''
  lines = text.splitlines()
  assert title in lines, (title, text)
  result = []
  for l in lines[lines.index(title) + 1:]:
    if not l.startswith('  ') or l.strip().startswith('Time to'): break
    result.append(l.strip())
  return result

def test_1xmk():
  h, description, log = run('tst_reduce2_metal_nh_1xmk', pdb_1xmk)
  assert h_on(h, 289, 'N') == [] and h_on(h, 290, 'N') == []
  # metal-bound chain start (metal_bound_amines_as_nh2, on by default): NH2,
  # listed there, no longer in the "not changed" list
  assert h_on(h, 288, 'N') == ['H2', 'H3']
  assert h_on(h, 290, 'ND1') == []
  for text in (description, log):
    assert section(text, NOT_PLACED) == [
      'H   SER A 289  (N   SER A 289 - NI    NI A 398 1.80 A (geometric))',
      'H   HIS A 290  (N   HIS A 290 - NI    NI A 398 1.84 A (geometric))'], text
    assert ('  N GLY A 288: metal-bound chain start, placed as NH2 '
      '(n_terminal_charge=residue_one would place 1 H, wrong for a terminus '
      'and for a chain break)') in text.splitlines(), text
    assert REPORTED not in text, text

def test_thymine_hg():
  h, description, log = run('tst_reduce2_metal_nh_5gsk', pdb_5gsk)
  assert h_on(h, 3, 'N3') == [] and h_on(h, 6, 'N3') == []
  lines = section(description, NOT_PLACED)
  assert len(lines) == 2, lines
  assert lines[0].startswith('H3   DT A   3') and 'sym y,x,-z' in lines[0], lines
  assert lines[1].startswith('H3   DT A   6') and 'sym' not in lines[1], lines

def test_alkaline_earth():
  for metal in ('MG', 'CA'):
    pdb = pdb_1xmk.replace('NI    NI A 398', '%-5s %2s A 398' % (metal, metal)
      ).replace('          NI\n', '          %2s\n' % metal)
    h, description, log = run('tst_reduce2_metal_nh_%s' % metal.lower(), pdb)
    assert h_on(h, 289, 'N') == ['H'] and h_on(h, 290, 'N') == ['H'], metal
    assert 'N-H' not in description, description

def test_trp_reported():
  h, description, log = run('tst_reduce2_metal_nh_trp', pdb_trp)
  assert h_on(h, 1, 'NE1') == ['HE1']
  assert NOT_PLACED not in description
  assert section(description, REPORTED) == [
    'NE1 TRP A   1 (HE1)  (NE1 TRP A   1 - ZN    ZN B 101 2.10 A (geometric))'
    ], description

def test_supplied_kept():
  for args in (['keep_existing_H=True'], ['approach=optimize']):
    h, description, log = run('tst_reduce2_metal_nh_supplied', pdb_1xmk_h,
      args)
    assert h_on(h, 289, 'N') == ['H'] and h_on(h, 290, 'N') == ['H'], args
    lines = description.splitlines()
    title = [l for l in lines if l.startswith(SUPPLIED)]
    assert title == [SUPPLIED + ': kept as supplied; backbone amide N, '
      'thymine/uracil N3, sulfonamide and azolium N bound to metals are '
      'usually deprotonated'], lines
    entries = section(description, title[0])
    assert entries[0].startswith('H   SER A 289  (N   SER A 289 - NI'), entries
    assert entries[1].startswith('H   HIS A 290  (N   HIS A 290 - NI'), entries
    assert NOT_PLACED not in description

def test_partial_metal():
  pdb = pdb_1xmk.replace('NI A 398      13.849  21.950  30.990  1.00',
                         'NI A 398      13.849  21.950  30.990  0.50')
  h, description, log = run('tst_reduce2_metal_nh_partial', pdb)
  assert h_on(h, 289, 'N') == ['H'] and h_on(h, 290, 'N') == ['H']
  lines = section(description, UNRESOLVED)
  assert len(lines) == 2, lines
  assert lines[0].startswith('H   SER A 289: N occupancy 1.00 exceeds '
    'coordinating NI occupancy 0.50.'), lines
  assert NOT_PLACED not in description

def test_ligand_n_kept():
  # an N bonded to another residue's C is an amide only when both residues
  # are amino acids
  h, description, log = run('tst_reduce2_metal_nh_cap', pdb_nh2_cap)
  assert h_on(h, 3, 'N') == ['HN1', 'HN2']
  assert NOT_PLACED not in description
  assert section(description, REPORTED) == [
    'N   NH2 A   3 (HN1, HN2)  (N   NH2 A   3 - NI    NI A 101 1.82 A '
    '(geometric))'], description

def test_sulfonamide():
  '''Sulfonamide N on Zn: ethoxzolamide keeps HN1A, which the Optimizer turns
  off the Zn (Zn-N1-HN1A at least 90 degrees); 3QR NAZ keeps no H'''
  from scitbx import matrix
  h, description, log = run('tst_reduce2_metal_nh_6bcc', pdb_6bcc)
  assert h_on(h, 302, 'N1') == ['HN1A'], h_on(h, 302, 'N1')
  sites = dict((a.name.strip(), matrix.col(a.xyz)) for a in h.atoms())
  angle = (sites['ZN'] - sites['N1']).angle(sites['HN1A'] - sites['N1'],
    deg=True)
  assert angle >= 90, angle
  assert section(description, NOT_PLACED) == ['HN1 EZL A 302  (N1  EZL A 302 '
    '- ZN    ZN A 301 1.99 A (geometric))'], description
  h, description, log = run('tst_reduce2_metal_nh_6gcy', pdb_6gcy)
  assert h_on(h, 302, 'NAZ') == [], h_on(h, 302, 'NAZ')

def test_azolium():
  '''Imidazolium on Zn (4PE3) loses HN1 on the bound N1 and keeps HN3; the
  neutral imidazole of OPV on Fe (6MFL) keeps H2 on the bound NAM'''
  h, description, log = run('tst_reduce2_metal_nh_4pe3', pdb_4pe3)
  assert h_on(h, 402, 'N1') == [] and h_on(h, 402, 'N3') == ['HN3'], (
    h_on(h, 402, 'N1'), h_on(h, 402, 'N3'))
  assert section(description, NOT_PLACED) == ['HN1 IMD A 402  (N1  IMD A 402 '
    '- ZN    ZN A 407 2.02 A (geometric))'], description
  h, description, log = run('tst_reduce2_metal_nh_6mfl', pdb_6mfl)
  assert h_on(h, 402, 'NAM') == ['H2'], h_on(h, 402, 'NAM')
  assert NOT_PLACED not in description, description

def test_conformer_without_peptide_c():
  # Ser 289 H/A: Gly 288 C(A) and CA(A), an amide; H/B: no C of its
  # conformer bonded to N, kept and reported
  h, description, log = run('tst_reduce2_metal_nh_altloc', pdb_1xmk_altloc)
  n289 = [a for a in h.atoms() if a.parent().parent().resseq_as_int() == 289
          and a.name.strip() == 'N'][0]
  hs = [(a.name.strip(), a.parent().altloc) for a in h.atoms()
        if a.element.strip() == 'H' and a.distance(n289) < 1.2]
  assert hs == [('H', 'B')], hs
  assert section(description, NOT_PLACED) == [
    'H  ASER A 289  (N   SER A 289 - NI    NI A 398 1.80 A (geometric))',
    'H   HIS A 290  (N   HIS A 290 - NI    NI A 398 1.84 A (geometric))'
    ], description
  assert 'N   SER A 289 (H altloc B)  (N   SER A 289 - NI    NI A 398 1.80 A '\
    '(geometric))' in section(description, REPORTED), description

def test_amino_acid_classification():
  # outside the residue-name classes: the dictionary's is_peptide()
  mon_lib_srv = server.server()
  for resname, expected in (('9FZ', True), ('9G2', True), ('NH2', False)):
    assert residue_is_amino_acid(resname, mon_lib_srv) == expected, resname

def test_peptide_without_o():
  h, description, log = run('tst_reduce2_metal_nh_no_o',
    pdb_1xmk_altloc_no_o)
  assert not [a for a in h.atoms() if a.element.strip() == 'O']
  n289 = [a for a in h.atoms() if a.parent().parent().resseq_as_int() == 289
          and a.name.strip() == 'N'][0]
  assert not [a for a in h.atoms()
              if a.element.strip() == 'H' and a.distance(n289) < 1.2]
  lines = section(description, NOT_PLACED)
  assert lines[:2] == [
    'H  ASER A 289  (N   SER A 289 - NI    NI A 398 1.80 A (geometric))',
    'H  BSER A 289  (N   SER A 289 - NI    NI A 398 1.80 A (geometric))'
    ], description

def test_rerun_stable():
  run('tst_reduce2_metal_nh_rerun0', pdb_1xmk)
  with open('tst_reduce2_metal_nh_rerun0_out.pdb') as f:
    first = f.read()
  for args in ([], ['approach=optimize']):
    h, description, log = run('tst_reduce2_metal_nh_rerun1', first, args)
    assert h_on(h, 289, 'N') == [] and h_on(h, 290, 'N') == [], args

if __name__ == '__main__':
  test_1xmk()
  test_thymine_hg()
  test_alkaline_earth()
  test_trp_reported()
  test_supplied_kept()
  test_partial_metal()
  test_ligand_n_kept()
  test_sulfonamide()
  test_azolium()
  test_conformer_without_peptide_c()
  test_amino_acid_classification()
  test_peptide_without_o()
  test_rerun_stable()
  print('OK')
