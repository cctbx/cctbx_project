"""reduce2 Optimizer: H deletions decided in an earlier altloc pass are kept.
Passes run C, B, A; each atom's deletion is decided by the last pass that
contained it, so an altloc-specific H is decided by its own pass and a blank H
by the last pass. Before, the deletion set was reset per pass and only the last
pass's deletions were applied.
Model: 1US0 His A 163 (whole residue A/B), His A 306 (ring blank, N/CA split,
so its ring H are copied into A and B) and His A 46 (blank) with residues
within 4 A of their rings; the states and pass decisions match the full model.
Pass B: His 163 "HD1Placed HE2NotPlaced", His 306 "HD1NotPlaced"; pass A: His
163 both H, His 306 "HD1NotPlaced"; His 46 "HD1NotPlaced" in both (control).
The reset also dropped the deletions of every model but the last of a
multi-model file (one Optimizer runs all models).
Shared (blank) H, both ways, each with its "Passes disagree" line: 1US0 His A
187 (blank) with residues within 4 A of its ring and Ser A 282 (CB/OG in A, B,
C). Pass C holds only blank atoms of the residues without a C conformer, so
His 187 meets no split neighbour there and keeps both ring H; passes B and A
delete HD1, which stays deleted. With the labels shifted (A to B, B to C, Ser
282 C to A) the pass without those neighbours is the last one: HD1 deleted in
passes C and B, kept.
"""
import os, re
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
import iotbx.pdb
from mmtbx.programs import reduce2

m_1us0 = """\
CRYST1   49.280   66.590   47.260  90.00  92.40  90.00 P 1 21 1
ATOM      1  N   HIS A  46       7.102  -8.519  20.780  1.00  2.96           N
ATOM      2  CA  HIS A  46       6.614  -9.791  20.290  1.00  2.95           C
ATOM      3  C   HIS A  46       7.725 -10.620  19.665  1.00  2.89           C
ATOM      4  O   HIS A  46       7.781 -11.839  19.841  1.00  3.35           O
ATOM      5  CB  HIS A  46       5.493  -9.544  19.270  1.00  3.48           C
ATOM      6  CG  HIS A  46       4.985 -10.800  18.650  1.00  3.74           C
ATOM      7  ND1 HIS A  46       4.291 -11.754  19.378  1.00  4.29           N
ATOM      8  CD2 HIS A  46       5.065 -11.272  17.378  1.00  4.70           C
ATOM      9  CE1 HIS A  46       3.968 -12.748  18.553  1.00  5.04           C
ATOM     10  NE2 HIS A  46       4.414 -12.487  17.357  1.00  5.35           N
ATOM     11  N   VAL A  47       8.623  -9.974  18.906  1.00  2.85           N
ATOM     12  CA  VAL A  47       9.657 -10.722  18.198  1.00  2.83           C
ATOM     13  C   VAL A  47      10.630 -11.397  19.157  1.00  2.78           C
ATOM     14  O   VAL A  47      11.356 -12.311  18.768  1.00  3.33           O
ATOM     15  CB  VAL A  47      10.382  -9.790  17.199  1.00  3.00           C
ATOM     16  CG1 VAL A  47      11.358  -8.864  17.897  1.00  3.36           C
ATOM     17  CG2 VAL A  47      11.071 -10.566  16.083  1.00  3.60           C
ATOM     18  N   TYR A  82       3.846  -6.932  10.615  1.00  3.67           N
ATOM     19  CA  TYR A  82       2.500  -6.530  10.978  1.00  4.09           C
ATOM     20  C   TYR A  82       2.256  -5.053  10.644  1.00  4.50           C
ATOM     21  O   TYR A  82       1.361  -4.451  11.209  1.00  5.48           O
ATOM     22  CB  TYR A  82       2.205  -6.782  12.465  1.00  4.34           C
ATOM     23  CG  TYR A  82       2.425  -8.204  12.894  1.00  4.48           C
ATOM     24  CD1 TYR A  82       1.484  -9.212  12.778  1.00  6.09           C
ATOM     25  CD2 TYR A  82       3.650  -8.597  13.450  1.00  5.43           C
ATOM     26  CE1 TYR A  82       1.739 -10.522  13.205  1.00  6.44           C
ATOM     27  CE2 TYR A  82       3.905  -9.883  13.872  1.00  6.09           C
ATOM     28  CZ  TYR A  82       2.934 -10.837  13.756  1.00  5.63           C
ATOM     29  OH  TYR A  82       3.255 -12.113  14.163  1.00  7.02           O
ATOM     30  N   THR A 113      10.848  -0.442   8.785  1.00  3.81           N
ATOM     31  CA  THR A 113      10.622  -1.289   7.644  1.00  3.64           C
ATOM     32  C   THR A 113       9.485  -2.280   7.855  1.00  3.83           C
ATOM     33  O   THR A 113       9.498  -3.061   8.801  1.00  4.78           O
ATOM     34  CB  THR A 113      11.911  -2.120   7.432  1.00  3.93           C
ATOM     35  OG1 THR A 113      12.275  -2.755   8.648  1.00  4.10           O
ATOM     36  CG2 THR A 113      13.110  -1.281   7.018  1.00  4.76           C
ATOM     37  N   GLY A 114       8.536  -2.297   6.909  1.00  3.88           N
ATOM     38  CA  GLY A 114       7.533  -3.337   6.874  1.00  4.07           C
ATOM     39  C   GLY A 114       7.988  -4.504   6.004  1.00  3.63           C
ATOM     40  O   GLY A 114       8.598  -4.300   4.948  1.00  3.84           O
ATOM     41  N   PHE A 121       9.284 -14.172  10.798  1.00  3.62           N
ATOM     42  CA  PHE A 121       9.579 -12.817  11.197  1.00  3.37           C
ATOM     43  C   PHE A 121      10.479 -12.075  10.206  1.00  3.22           C
ATOM     44  O   PHE A 121      10.433 -10.849  10.152  1.00  3.56           O
ATOM     45  CB  PHE A 121      10.222 -12.801  12.581  1.00  3.63           C
ATOM     46  CG  PHE A 121       9.277 -13.196  13.710  1.00  3.76           C
ATOM     47  CD1 PHE A 121       7.905 -12.978  13.658  1.00  4.75           C
ATOM     48  CD2 PHE A 121       9.812 -13.753  14.864  1.00  4.06           C
ATOM     49  CE1 PHE A 121       7.090 -13.270  14.762  1.00  5.74           C
ATOM     50  CE2 PHE A 121       8.996 -14.034  15.958  1.00  4.99           C
ATOM     51  CZ  PHE A 121       7.642 -13.791  15.899  1.00  5.47           C
ATOM     52  N   PRO A 132      10.024  -6.962   1.323  1.00  4.08           N
ATOM     53  CA  PRO A 132       9.297  -5.837   1.890  1.00  4.21           C
ATOM     54  C   PRO A 132       7.795  -5.987   1.776  1.00  4.30           C
ATOM     55  O   PRO A 132       7.266  -6.698   0.925  1.00  5.60           O
ATOM     56  CB  PRO A 132       9.734  -4.631   1.034  1.00  5.17           C
ATOM     57  CG  PRO A 132      11.054  -5.074   0.400  1.00  5.58           C
ATOM     58  CD  PRO A 132      10.833  -6.547   0.153  1.00  4.74           C
ATOM     59  N  AHIS A 163      16.783  10.331  11.429  0.60  4.57           N
ATOM     60  CA AHIS A 163      16.333  11.697  11.471  0.60  5.57           C
ATOM     61  C  AHIS A 163      15.236  11.990  10.437  0.60  5.73           C
ATOM     62  O  AHIS A 163      14.211  12.606  10.704  0.60  6.47           O
ATOM     63  CB AHIS A 163      17.576  12.593  11.290  0.60  6.31           C
ATOM     64  CG AHIS A 163      17.372  14.039  11.401  0.60  7.51           C
ATOM     65  ND1AHIS A 163      18.351  14.790  11.987  0.60 12.51           N
ATOM     66  CD2AHIS A 163      16.461  14.915  11.072  0.60  8.19           C
ATOM     67  CE1AHIS A 163      17.978  16.080  11.917  0.60 14.99           C
ATOM     68  NE2AHIS A 163      16.822  16.141  11.289  0.60  9.46           N
ATOM     69  N  BHIS A 163      16.653  10.488  11.436  0.40  7.46           N
ATOM     70  CA BHIS A 163      16.005  11.780  11.701  0.40  7.75           C
ATOM     71  C  BHIS A 163      15.042  12.105  10.533  0.40  8.52           C
ATOM     72  O  BHIS A 163      14.036  12.715  10.882  0.40  8.72           O
ATOM     73  CB BHIS A 163      16.986  12.891  11.987  0.40  9.69           C
ATOM     74  CG BHIS A 163      17.482  13.617  10.803  0.40 12.73           C
ATOM     75  ND1BHIS A 163      18.476  13.045  10.014  0.40 16.25           N
ATOM     76  CD2BHIS A 163      17.249  14.846  10.325  0.40 16.59           C
ATOM     77  CE1BHIS A 163      18.774  13.917   9.050  0.40 16.30           C
ATOM     78  NE2BHIS A 163      18.056  15.017   9.206  0.40 17.83           N
ATOM     79  N  ALEU A 164      15.497  11.558   9.215  0.60  5.79           N
ATOM     80  CA ALEU A 164      14.485  11.866   8.181  0.60  7.69           C
ATOM     81  C  ALEU A 164      13.189  11.126   8.408  0.60  5.71           C
ATOM     82  O  ALEU A 164      12.108  11.628   8.035  0.60  6.21           O
ATOM     83  CB ALEU A 164      15.130  11.571   6.834  0.60 10.89           C
ATOM     84  CG ALEU A 164      16.275  12.542   6.551  0.60 14.93           C
ATOM     85  CD1ALEU A 164      16.842  12.203   5.192  0.60 28.12           C
ATOM     86  CD2ALEU A 164      15.917  14.038   6.554  0.60 23.68           C
ATOM     87  N  BLEU A 164      15.340  11.718   9.302  0.40  7.83           N
ATOM     88  CA BLEU A 164      14.466  11.842   8.126  0.40  7.56           C
ATOM     89  C  BLEU A 164      13.107  11.214   8.470  0.40  6.53           C
ATOM     90  O  BLEU A 164      12.056  11.781   8.195  0.40  7.00           O
ATOM     91  CB BLEU A 164      15.055  11.262   6.839  0.40 11.14           C
ATOM     92  CG BLEU A 164      14.246  11.145   5.556  0.40 15.07           C
ATOM     93  CD1BLEU A 164      15.050  11.027   4.283  0.40 16.25           C
ATOM     94  CD2BLEU A 164      13.400   9.806   5.554  0.40 13.62           C
ATOM     95  N  ATYR A 198      18.083  16.110  17.640  0.70  4.92           N
ATOM     96  CA ATYR A 198      16.891  16.333  16.843  0.70  5.02           C
ATOM     97  C  ATYR A 198      15.637  15.900  17.588  0.70  4.74           C
ATOM     98  O  ATYR A 198      14.667  16.638  17.728  0.70  5.16           O
ATOM     99  CB ATYR A 198      17.013  15.612  15.484  0.70  5.91           C
ATOM    100  CG ATYR A 198      15.766  15.841  14.685  0.70  5.77           C
ATOM    101  CD1ATYR A 198      15.447  17.095  14.152  0.70  6.94           C
ATOM    102  CD2ATYR A 198      14.814  14.843  14.482  0.70  7.47           C
ATOM    103  CE1ATYR A 198      14.299  17.342  13.429  0.70  7.49           C
ATOM    104  CE2ATYR A 198      13.635  15.083  13.767  0.70  8.66           C
ATOM    105  CZ ATYR A 198      13.394  16.328  13.239  0.70  7.40           C
ATOM    106  OH ATYR A 198      12.229  16.482  12.569  0.70  9.50           O
ATOM    107  N  BTYR A 198      17.781  16.304  18.148  0.30  5.77           N
ATOM    108  CA BTYR A 198      16.641  16.583  17.299  0.30  6.00           C
ATOM    109  C  BTYR A 198      15.365  16.109  18.013  0.30  4.99           C
ATOM    110  O  BTYR A 198      14.378  16.812  18.164  0.30  6.53           O
ATOM    111  CB BTYR A 198      16.790  15.938  15.927  0.30  6.23           C
ATOM    112  CG BTYR A 198      15.512  16.184  15.143  0.30  6.73           C
ATOM    113  CD1BTYR A 198      15.078  17.465  14.765  0.30  9.37           C
ATOM    114  CD2BTYR A 198      14.715  15.091  14.811  0.30  8.99           C
ATOM    115  CE1BTYR A 198      13.899  17.652  14.070  0.30  9.39           C
ATOM    116  CE2BTYR A 198      13.535  15.282  14.108  0.30 10.11           C
ATOM    117  CZ BTYR A 198      13.147  16.557  13.743  0.30  9.66           C
ATOM    118  OH BTYR A 198      11.992  16.732  13.049  0.30 14.54           O
ATOM    119  C   HIS A 306      14.677  -0.747   1.918  1.00  5.64           C
ATOM    120  O   HIS A 306      15.585  -0.087   2.444  1.00  6.43           O
ATOM    121  CB  HIS A 306      14.146  -2.665   3.467  1.00  5.27           C
ATOM    122  CG  HIS A 306      12.697  -2.414   3.685  1.00  4.67           C
ATOM    123  ND1 HIS A 306      12.135  -1.146   3.697  1.00  4.78           N
ATOM    124  CD2 HIS A 306      11.710  -3.317   3.918  1.00  4.82           C
ATOM    125  CE1 HIS A 306      10.859  -1.311   3.914  1.00  4.70           C
ATOM    126  NE2 HIS A 306      10.554  -2.615   4.056  1.00  4.79           N
ATOM    127  N  AHIS A 306      16.032  -2.645   1.879  0.62  5.79           N
ATOM    128  CA AHIS A 306      14.642  -2.250   2.067  0.62  5.62           C
ATOM    129  N  BHIS A 306      15.963  -2.917   1.965  0.38  5.83           N
ATOM    130  CA BHIS A 306      14.666  -2.254   2.070  0.38  5.27           C
ATOM    131  N   LYS A 307      13.690  -0.181   1.249  1.00  5.76           N
ATOM    132  CA  LYS A 307      13.645   1.246   1.022  1.00  7.23           C
ATOM    133  C   LYS A 307      13.738   2.061   2.306  1.00  6.39           C
ATOM    134  O   LYS A 307      14.250   3.165   2.279  1.00  7.85           O
ATOM    135  CE  LYS A 307       9.813   2.137  -1.245  0.56 22.13           C
ATOM    136  NZ  LYS A 307       8.823   2.971  -1.988  0.56 28.73           N
ATOM    137  CB ALYS A 307      12.257   1.523   0.421  0.56 10.80           C
ATOM    138  CG ALYS A 307      11.868   2.899   0.041  0.56 17.07           C
ATOM    139  CD ALYS A 307      11.008   3.032  -1.150  0.56 23.99           C
ATOM    140  CB BLYS A 307      12.422   1.647   0.206  0.44  9.14           C
ATOM    141  CG BLYS A 307      11.140   1.687   0.849  0.44  9.05           C
ATOM    142  CD BLYS A 307      10.024   1.922  -0.168  0.44 12.28           C
ATOM    143  N   ASP A 308      13.190   1.513   3.392  1.00  5.19           N
ATOM    144  CA  ASP A 308      13.127   2.236   4.652  1.00  5.17           C
ATOM    145  C   ASP A 308      14.176   1.754   5.657  1.00  4.62           C
ATOM    146  O   ASP A 308      14.107   2.101   6.844  1.00  4.88           O
ATOM    147  CB  ASP A 308      11.724   2.142   5.232  1.00  5.47           C
ATOM    148  CG  ASP A 308      10.689   2.801   4.317  1.00  7.68           C
ATOM    149  OD1 ASP A 308      11.037   3.879   3.748  1.00  9.52           O
ATOM    150  OD2 ASP A 308       9.605   2.252   4.160  1.00 10.20           O
HETATM  151  C1  CIT A 400      21.468  13.190  11.677  0.38  8.99           C
HETATM  152  C2  CIT A 400      22.280  11.994  11.364  0.38  7.86           C
HETATM  153  C3  CIT A 400      21.718  10.592  11.722  0.38  7.31           C
HETATM  154  C4  CIT A 400      22.680   9.482  11.410  0.38  7.81           C
HETATM  155  C5  CIT A 400      23.997   9.419  12.127  0.38  8.12           C
HETATM  156  C6  CIT A 400      20.308  10.327  11.067  0.38  5.91           C
HETATM  157  O1  CIT A 400      20.565  13.088  12.605  0.38  7.28           O
HETATM  158  O2  CIT A 400      21.741  14.362  11.155  0.38  9.59           O
HETATM  159  O3  CIT A 400      24.138   8.452  13.038  0.38 10.61           O
HETATM  160  O4  CIT A 400      24.894  10.125  11.705  0.38 14.93           O
HETATM  161  O5  CIT A 400      19.396   9.939  11.959  0.38  7.12           O
HETATM  162  O6  CIT A 400      20.126  10.734   9.907  0.38  9.70           O
HETATM  163  O7  CIT A 400      21.618  10.670  13.172  0.38  8.00           O
HETATM  164  C1  CIT A 450       6.260 -16.763  24.168  0.51  3.55           C
HETATM  165  C2  CIT A 450       6.050 -17.259  22.751  0.51  4.50           C
HETATM  166  C3  CIT A 450       4.553 -17.365  22.401  0.51  4.83           C
HETATM  167  C4  CIT A 450       4.431 -18.117  21.078  0.51  5.20           C
HETATM  168  C5  CIT A 450       3.007 -18.440  20.727  0.51  5.80           C
HETATM  169  C6  CIT A 450       3.988 -15.955  22.237  0.51  4.77           C
HETATM  170  O1  CIT A 450       5.958 -15.573  24.450  0.51  4.07           O
HETATM  171  O2  CIT A 450       6.743 -17.602  24.951  0.51  3.60           O
HETATM  172  O3  CIT A 450       2.178 -17.525  20.672  0.51  5.95           O
HETATM  173  O4  CIT A 450       2.767 -19.628  20.386  0.51  7.53           O
HETATM  174  O5  CIT A 450       3.064 -15.660  23.074  0.51  5.33           O
HETATM  175  O6  CIT A 450       4.423 -15.171  21.396  0.51  5.42           O
HETATM  176  O7  CIT A 450       3.823 -18.019  23.401  0.51  5.58           O
HETATM  177  O   HOH A2057       0.989 -14.318  16.884  1.00 15.54           O
HETATM  178  O   HOH A2058       5.660 -16.118  17.735  0.50  6.23           O
HETATM  179  O   HOH A2121       4.119 -14.751  15.842  1.00  6.55           O
HETATM  180  O   HOH A2139       4.142 -11.800  22.044  0.51  2.78           O
HETATM  181  O   HOH A2194       4.150  -7.497  16.773  1.00  4.77           O
HETATM  182  O   HOH A2201       0.878 -11.446  16.612  1.00 11.31           O
HETATM  183  O   HOH A2227       2.037  -9.545  18.887  0.35 11.24           O
HETATM  184  O   HOH A2228       2.277  -9.390  17.798  0.65  5.86           O
HETATM  185  O   HOH A2255       3.863 -11.729  22.030  0.49 11.34           O
HETATM  186  O   HOH A2284       8.499   0.709   5.936  0.56  6.26           O
HETATM  187  O   HOH A2308       8.716  -1.097   0.845  1.00  7.88           O
HETATM  188  O   HOH A2309       7.519  -1.331   4.340  1.00  8.93           O
HETATM  189  O   HOH A2352      15.451  18.525  10.386  1.00 18.97           O
HETATM  190  O   HOH A2579      11.469  -1.743   0.391  1.00  6.20           O
HETATM  191  O   HOH A2599      18.218  10.669   8.288  1.00 12.65           O
HETATM  192  O   HOH A2600      20.973  12.054  12.732  0.63 18.43           O
HETATM  193  O   HOH A2603      19.570   9.769  11.804  0.63  9.09           O
HETATM  194  O   HOH A2604       2.724 -15.774  18.517  0.51  9.43           O
HETATM  195  O   HOH A2610       6.544 -14.408  19.687  1.00  6.62           O
"""

m_his187 = """\
CRYST1   49.280   66.590   47.260  90.00  92.40  90.00 P 1 21 1
ATOM      1  N   GLU A 185      19.650   2.564  19.498  1.00  3.24           N
ATOM      2  CA  GLU A 185      20.880   2.084  20.108  1.00  3.27           C
ATOM      3  C   GLU A 185      21.699   3.287  20.574  1.00  3.32           C
ATOM      4  O   GLU A 185      21.876   4.248  19.820  1.00  3.90           O
ATOM      5  CB  GLU A 185      21.685   1.316  19.075  1.00  3.52           C
ATOM      6  CG  GLU A 185      22.915   0.672  19.672  1.00  4.41           C
ATOM      7  CD  GLU A 185      23.700  -0.175  18.686  1.00  5.08           C
ATOM      8  OE1 GLU A 185      23.095  -0.661  17.700  1.00  5.65           O
ATOM      9  OE2 GLU A 185      24.897  -0.390  18.959  1.00  6.70           O
ATOM     10  N   CYS A 186      22.251   3.226  21.792  1.00  3.34           N
ATOM     11  CA  CYS A 186      23.063   4.352  22.231  1.00  3.42           C
ATOM     12  C   CYS A 186      23.944   3.924  23.393  1.00  3.43           C
ATOM     13  O   CYS A 186      23.437   3.339  24.357  1.00  3.70           O
ATOM     14  CB ACYS A 186      22.122   5.488  22.658  0.70  3.77           C
ATOM     15  SG ACYS A 186      23.034   6.948  23.170  0.70  5.03           S
ATOM     16  CB BCYS A 186      22.198   5.514  22.639  0.30  4.57           C
ATOM     17  SG BCYS A 186      22.935   7.156  23.001  0.30  8.84           S
ATOM     18  N   HIS A 187      25.221   4.236  23.270  1.00  3.75           N
ATOM     19  CA  HIS A 187      26.225   3.866  24.258  1.00  3.69           C
ATOM     20  C   HIS A 187      27.495   4.648  23.897  1.00  3.65           C
ATOM     21  O   HIS A 187      27.523   5.363  22.893  1.00  4.01           O
ATOM     22  CB  HIS A 187      26.444   2.336  24.274  1.00  3.89           C
ATOM     23  CG  HIS A 187      26.659   1.764  22.927  1.00  3.90           C
ATOM     24  ND1 HIS A 187      27.787   2.039  22.161  1.00  4.24           N
ATOM     25  CD2 HIS A 187      25.783   0.988  22.211  1.00  4.25           C
ATOM     26  CE1 HIS A 187      27.562   1.417  21.005  1.00  4.37           C
ATOM     27  NE2 HIS A 187      26.379   0.780  21.001  1.00  4.44           N
ATOM     28  N   PRO A 188      28.567   4.545  24.689  1.00  3.60           N
ATOM     29  CA  PRO A 188      29.749   5.371  24.377  1.00  3.80           C
ATOM     30  C   PRO A 188      30.387   5.130  23.032  1.00  3.88           C
ATOM     31  O   PRO A 188      31.148   5.992  22.596  1.00  4.70           O
ATOM     32  CB  PRO A 188      30.716   5.079  25.520  1.00  4.03           C
ATOM     33  CG  PRO A 188      29.783   4.755  26.685  1.00  4.15           C
ATOM     34  CD  PRO A 188      28.646   3.959  26.040  1.00  3.84           C
ATOM     35  CB  TYR A 189      31.161   2.318  20.873  1.00  4.91           C
ATOM     36  CG  TYR A 189      32.212   1.832  21.814  1.00  4.64           C
ATOM     37  CD1 TYR A 189      33.153   2.656  22.427  1.00  5.01           C
ATOM     38  CD2 TYR A 189      32.298   0.465  22.082  1.00  5.33           C
ATOM     39  CE1 TYR A 189      34.117   2.143  23.284  1.00  5.53           C
ATOM     40  CE2 TYR A 189      33.257  -0.053  22.913  1.00  5.78           C
ATOM     41  CZ  TYR A 189      34.192   0.781  23.527  1.00  5.71           C
ATOM     42  OH  TYR A 189      35.102   0.214  24.322  1.00  7.60           O
ATOM     43  N  ATYR A 189      30.179   3.957  22.382  0.70  3.81           N
ATOM     44  CA ATYR A 189      30.751   3.760  21.088  0.70  4.26           C
ATOM     45  C  ATYR A 189      29.859   4.180  19.894  0.70  4.04           C
ATOM     46  O  ATYR A 189      30.312   4.189  18.749  0.70  4.60           O
ATOM     47  N  BTYR A 189      30.072   4.022  22.367  0.30  4.47           N
ATOM     48  CA BTYR A 189      30.638   3.770  21.034  0.30  4.24           C
ATOM     49  C  BTYR A 189      29.628   4.269  19.967  0.30  4.19           C
ATOM     50  O  BTYR A 189      30.000   4.201  18.798  0.30  5.01           O
ATOM     51  N  ALEU A 190      28.628   4.548  20.230  0.70  3.78           N
ATOM     52  CA ALEU A 190      27.606   5.007  19.274  0.70  3.71           C
ATOM     53  C  ALEU A 190      26.781   6.055  20.023  0.70  3.36           C
ATOM     54  O  ALEU A 190      25.693   5.796  20.563  0.70  3.68           O
ATOM     55  CB ALEU A 190      26.715   3.852  18.800  0.70  4.15           C
ATOM     56  CG ALEU A 190      25.758   4.277  17.672  0.70  4.34           C
ATOM     57  CD1ALEU A 190      26.459   4.437  16.339  0.70  5.82           C
ATOM     58  CD2ALEU A 190      24.674   3.225  17.536  0.70  5.31           C
ATOM     59  N  BLEU A 190      28.428   4.746  20.343  0.30  5.13           N
ATOM     60  CA BLEU A 190      27.378   5.229  19.418  0.30  4.89           C
ATOM     61  C  BLEU A 190      26.506   6.219  20.135  0.30  4.99           C
ATOM     62  O  BLEU A 190      25.418   5.939  20.669  0.30  5.91           O
ATOM     63  CB BLEU A 190      26.539   4.062  18.892  0.30  4.88           C
ATOM     64  CG BLEU A 190      25.427   4.359  17.878  0.30  5.95           C
ATOM     65  CD1BLEU A 190      26.014   4.711  16.521  0.30  6.21           C
ATOM     66  CD2BLEU A 190      24.487   3.169  17.835  0.30  6.61           C
ATOM     67  N   SER A 210      20.679  -1.334  24.144  1.00  3.19           N
ATOM     68  CA  SER A 210      22.006  -1.906  24.172  1.00  3.32           C
ATOM     69  C   SER A 210      23.019  -0.993  24.876  1.00  3.35           C
ATOM     70  O   SER A 210      24.130  -0.788  24.387  1.00  3.92           O
ATOM     71  CB  SER A 210      22.518  -2.307  22.785  1.00  3.79           C
ATOM     72  OG  SER A 210      21.607  -3.164  22.119  1.00  3.96           O
ATOM     73  N   LEU A 227      32.574  -2.417  25.682  1.00  6.82           N
ATOM     74  CA  LEU A 227      32.380  -1.365  26.640  1.00  6.14           C
ATOM     75  C   LEU A 227      32.705  -1.800  28.062  1.00  6.48           C
ATOM     76  O   LEU A 227      33.351  -1.077  28.810  1.00  7.14           O
ATOM     77  CB  LEU A 227      30.932  -0.862  26.570  1.00  5.61           C
ATOM     78  CG  LEU A 227      30.585  -0.081  25.326  1.00  5.35           C
ATOM     79  CD1 LEU A 227      29.083  -0.001  25.174  1.00  5.85           C
ATOM     80  CD2 LEU A 227      31.216   1.298  25.409  1.00  5.53           C
ATOM     81  N  ATRP A 295      29.755   2.637  16.151  0.50  7.09           N
ATOM     82  CA ATRP A 295      29.679   1.177  16.205  0.50  6.69           C
ATOM     83  C  ATRP A 295      28.311   0.694  16.667  0.50  5.65           C
ATOM     84  O  ATRP A 295      27.849   1.037  17.755  0.50  6.52           O
ATOM     85  CB ATRP A 295      30.758   0.591  17.124  0.50  7.61           C
ATOM     86  CG ATRP A 295      30.675  -0.918  17.200  0.50  7.75           C
ATOM     87  CD1ATRP A 295      31.075  -1.817  16.250  0.50  9.69           C
ATOM     88  CD2ATRP A 295      30.143  -1.706  18.258  0.50  7.16           C
ATOM     89  NE1ATRP A 295      30.837  -3.088  16.679  0.50 10.02           N
ATOM     90  CE2ATRP A 295      30.261  -3.077  17.912  0.50  8.94           C
ATOM     91  CE3ATRP A 295      29.594  -1.450  19.501  0.50  8.20           C
ATOM     92  CZ2ATRP A 295      29.834  -4.130  18.725  0.50  8.89           C
ATOM     93  CZ3ATRP A 295      29.175  -2.492  20.314  0.50  8.59           C
ATOM     94  CH2ATRP A 295      29.290  -3.827  19.949  0.50  9.44           C
ATOM     95  N  BTRP A 295      30.180   2.779  16.151  0.50  6.46           N
ATOM     96  CA BTRP A 295      30.180   1.312  16.148  0.50  6.55           C
ATOM     97  C  BTRP A 295      28.766   0.862  16.524  0.50  5.03           C
ATOM     98  O  BTRP A 295      28.302   1.013  17.664  0.50  5.11           O
ATOM     99  CB BTRP A 295      31.195   0.806  17.116  0.50  7.68           C
ATOM    100  CG BTRP A 295      31.230  -0.697  17.155  0.50  9.55           C
ATOM    101  CD1BTRP A 295      31.642  -1.486  16.145  0.50 11.08           C
ATOM    102  CD2BTRP A 295      30.889  -1.555  18.224  0.50 10.12           C
ATOM    103  NE1BTRP A 295      31.558  -2.810  16.525  0.50 14.46           N
ATOM    104  CE2BTRP A 295      31.099  -2.873  17.813  0.50 12.93           C
ATOM    105  CE3BTRP A 295      30.417  -1.394  19.512  0.50 10.92           C
ATOM    106  CZ2BTRP A 295      30.869  -3.993  18.580  0.50 17.07           C
ATOM    107  CZ3BTRP A 295      30.180  -2.495  20.327  0.50 13.58           C
ATOM    108  CH2BTRP A 295      30.413  -3.769  19.842  0.50 18.10           C
ATOM   2640  N   SER A 282      30.487  16.632  34.463  1.00  4.20           N
ATOM   2641  CA  SER A 282      30.421  17.016  33.067  1.00  4.34           C
ATOM   2642  C   SER A 282      31.649  16.510  32.295  1.00  3.91           C
ATOM   2643  O   SER A 282      31.527  16.107  31.139  1.00  4.24           O
ATOM   2644  CB ASER A 282      30.380  18.557  33.061  0.60  5.87           C
ATOM   2645  OG ASER A 282      30.223  18.988  31.756  0.60  6.95           O
ATOM   2646  CB BSER A 282      29.963  18.420  32.693  0.20  4.90           C
ATOM   2647  OG BSER A 282      30.941  19.278  33.244  0.20  5.80           O
ATOM   2648  CB CSER A 282      30.406  18.536  32.891  0.20  4.80           C
ATOM   2649  OG CSER A 282      29.221  18.983  33.515  0.20  4.95           O
"""

def shift_labels(model_text):
  '''A to B, B to C; Ser 282: C to A, its A and B dropped'''
  out = []
  for l in model_text.split('\n'):
    if l.startswith(('ATOM', 'HETATM')):
      alt = l[16]
      if l[22:26] == ' 282':
        if alt in ('A', 'B'):
          continue
        if alt == 'C':
          l = l[:16] + 'A' + l[17:]
      else:
        l = l[:16] + {'A': 'B', 'B': 'C'}.get(alt, alt) + l[17:]
    out.append(l)
  return '\n'.join(out)

def two_models(model_text):
  lines = model_text.split('\n')
  atoms = '\n'.join(l for l in lines if l.startswith(('ATOM', 'HETATM')))
  return '\n'.join([lines[0], 'MODEL        1', atoms, 'ENDMDL',
    'MODEL        2', atoms, 'ENDMDL', 'END', ''])

def run(prefix, extra=(), model_text=m_1us0, model_index=0,
    residues=(46, 163, 306)):
  '''reduce2 with flips; returns (description, {(resseq, altloc): atom names})
  of the model with this index'''
  fn_in = prefix + '_in.pdb'
  fn_out = prefix + '_out.pdb'
  fn_desc = prefix + '_desc.txt'
  with open(fn_in, "w") as f:
    f.write(model_text)
  for fn in (fn_out, fn_desc):
    if os.path.exists(fn): os.remove(fn)
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "approach=add", "add_flip_movers=True",
      "output.filename=" + fn_out, "output.description_file_name=" + fn_desc,
      "overwrite=True"] + list(extra))
  with open(fn_desc) as f:
    desc = f.read()
  names = {}
  h = iotbx.pdb.input(fn_out).construct_hierarchy()
  for rg in h.models()[model_index].residue_groups():
    if rg.resseq_as_int() in residues:
      for ag in rg.atom_groups():
        key = (rg.resseq_as_int(), ag.altloc.strip())
        if key not in names:
          names[key] = set()
        for a in ag.atoms():
          names[key].add(a.name.strip())
  return desc, names

def his_poses(desc):
  '''{(pass altloc, His line label): pose}'''
  out, alt = {}, None
  for l in desc.split('\n'):
    m = re.search(r"Running Reduce optimization on model index \d+, "
      r"alternate '(.*)'", l)
    if m:
      alt = m.group(1)
      continue
    m = re.match(r"\s*HisFlip at (chain A .*?) NE2 .* pose (.*)$", l)
    if m:
      out[(alt, m.group(1))] = m.group(2)
  return out

def deleting(desc):
  return [l.strip() for l in desc.split('\n')
    if l.strip().startswith('Deleting chain')]

def check(desc, names, n_models=1):
  p = his_poses(desc)
  # the fragment still makes the decisions the test is about
  assert 'HE2NotPlaced' in p[('B', 'chain A BHIS 163')], p
  assert 'HD1Placed HE2Placed' in p[('A', 'chain A AHIS 163')], p
  for alt in ('B', 'A'):
    assert 'HD1NotPlaced HE2Placed' in p[(alt, 'chain A HIS 306')], p
    assert 'HD1NotPlaced HE2Placed' in p[(alt, 'chain A HIS 46')], p
  # ring H as H or D (hd_and_h writes exchange pairs in new altlocs)
  ring = {}
  for k, v in names.items():
    ring[k] = set(n.replace('DD1', 'HD1').replace('DE2', 'HE2') for n in v)
  def anywhere(resseq):
    out = set()
    for k, v in ring.items():
      if k[0] == resseq:
        out |= v
    return out
  # His 163: HE2 deleted in its own pass B, both ring H kept in A
  assert 'HD1' in ring[(163, 'B')] and 'HE2' not in ring[(163, 'B')], names
  assert {'HD1', 'HE2'} <= ring[(163, 'A')], names
  # His 306: HD1 deleted in both passes, HE2 kept
  assert 'HD1' not in anywhere(306) and 'HE2' in anywhere(306), names
  # control, blank His 46: HD1 deleted, HE2 kept
  assert 'HD1' not in anywhere(46) and 'HE2' in anywhere(46), names
  d = deleting(desc)
  assert sorted(d) == sorted(n_models * ['Deleting chain A HIS 46 HD1',
    'Deleting chain A BHIS 163 HE2', 'Deleting chain A AHIS 306 HD1',
    'Deleting chain A BHIS 306 HD1']), d

def test_xray():
  desc, names = run('tst_reduce2_altloc_deletions_x')
  check(desc, names)

def test_neutron():
  '''hd_and_h: the H are deleted before the ring H become D'''
  desc, names = run('tst_reduce2_altloc_deletions_n',
    ["neutron_option=hd_and_h", "use_neutron_distances=True"])
  check(desc, names)

def test_two_models():
  '''the first of two models keeps its deletions'''
  desc, names = run('tst_reduce2_altloc_deletions_m',
    model_text=two_models(m_1us0))
  check(desc, names, n_models=2)

def check_his187(desc, names, passes, last):
  '''passes: (pass altloc, 'HD1Placed' or 'HD1NotPlaced') in pass order'''
  p = his_poses(desc)
  for alt, pose in passes:
    assert p[(alt, 'chain A HIS 187')].startswith('Unflipped ' + pose), p
  deleted = [alt for alt, pose in passes if pose == 'HD1NotPlaced']
  kept = [alt for alt, pose in passes if pose == 'HD1Placed']
  line = ("Passes disagree about deleting chain A HIS 187 HD1 (deleted in %s;"
    " kept in %s); the last pass decides: %s" % (
    " ".join("'%s'" % a for a in deleted), " ".join("'%s'" % a for a in kept),
    last))
  assert [l.strip() for l in desc.split('\n') if 'Passes disagree' in l] \
    == [line], desc
  assert 'HE2' in names[(187, '')], names
  if last == 'kept':
    assert 'HD1' in names[(187, '')], names
    assert deleting(desc) == [], deleting(desc)
  else:
    assert 'HD1' not in names[(187, '')], names
    assert deleting(desc) == ['Deleting chain A HIS 187 HD1'], deleting(desc)

def test_shared_h_deleted_by_last_pass():
  '''blank HD1 kept in pass C, deleted in B and A: deleted'''
  desc, names = run('tst_reduce2_altloc_deletions_187d', model_text=m_his187,
    residues=(187,))
  check_his187(desc, names, (('C', 'HD1Placed'), ('B', 'HD1NotPlaced'),
    ('A', 'HD1NotPlaced')), 'deleted')

def test_shared_h_kept_by_last_pass():
  '''blank HD1 deleted in passes C and B, kept in A: kept'''
  desc, names = run('tst_reduce2_altloc_deletions_187k',
    model_text=shift_labels(m_his187), residues=(187,))
  check_his187(desc, names, (('C', 'HD1NotPlaced'), ('B', 'HD1NotPlaced'),
    ('A', 'HD1Placed')), 'kept')

if __name__ == "__main__":
  test_xray()
  test_neutron()
  test_two_models()
  test_shared_h_deleted_by_last_pass()
  test_shared_h_kept_by_last_pass()
  print("OK")
