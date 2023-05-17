
from __future__ import absolute_import, division, print_function
from mmtbx.command_line import sort_hetatms
import iotbx.pdb
from six.moves import cStringIO as StringIO
import os.path as op
import os

def exercise():
  pdb_raw = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
ATOM     21  N   GLN A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     22  CA  GLN A   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     23  C   GLN A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     24  O   GLN A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     25  CB  GLN A   4       0.656   2.148   1.711  1.00  9.80           C
ATOM     26  CG  GLN A   4       1.944   1.458   1.213  1.00 10.25           C
ATOM     27  CD  GLN A   4       2.504   2.044  -0.089  1.00 12.43           C
ATOM     28  OE1 GLN A   4       2.744   3.268  -0.190  1.00 14.62           O
ATOM     29  NE2 GLN A   4       2.750   1.161  -1.091  1.00  9.05           N
TER
ATOM     30  N   GLN B   1       2.154   1.821   4.871  1.00 10.38           N
ATOM     31  CA  GLN B   1       3.270   2.361   5.640  1.00 11.39           C
ATOM     32  C   GLN B   1       4.594   1.768   5.172  1.00 11.52           C
ATOM     33  O   GLN B   1       4.768   0.546   5.054  1.00 12.05           O
ATOM     34  CB  GLN B   1       3.056   2.183   7.147  1.00 11.96           C
ATOM     35  CG  GLN B   1       1.829   2.950   7.647  1.00 10.81           C
ATOM     36  CD  GLN B   1       1.344   2.414   8.954  1.00 13.10           C
ATOM     37  OE1 GLN B   1       0.774   1.325   9.002  1.00 10.65           O
ATOM     38  NE2 GLN B   1       1.549   3.187  10.039  1.00 12.30           N
ATOM     39  N   ASN B   2       5.514   2.664   4.856  1.00 11.99           N
ATOM     40  CA  ASN B   2       6.831   2.310   4.318  1.00 12.30           C
ATOM     41  C   ASN B   2       7.854   2.761   5.324  1.00 13.40           C
ATOM     42  O   ASN B   2       8.219   3.943   5.374  1.00 13.92           O
ATOM     43  CB  ASN B   2       7.065   3.016   2.993  1.00 12.13           C
ATOM     44  CG  ASN B   2       5.961   2.735   2.003  1.00 12.77           C
ATOM     45  OD1 ASN B   2       5.798   1.604   1.551  1.00 14.27           O
ATOM     46  ND2 ASN B   2       5.195   3.747   1.679  1.00 10.07           N
ATOM     47  N   TYR B   3       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR B   3       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR B   3      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR B   3      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR B   3       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR B   3       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR B   3       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR B   3       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR B   3       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR B   3       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR B   3       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR B   3       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR B   3      11.358   2.999   7.612  1.00 17.49           O
TER
HETATM   61  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM   62  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM   63  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM   64  O   HOH A  11      11.808   4.179   9.970  1.00 23.99           O
HETATM   65  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
HETATM   66  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM   67  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
END
"""
  pdb_file = "unsorted.pdb"
  with open(pdb_file, "w") as f:
    f.write(pdb_raw)
  if (op.exists("unsorted_sorted.pdb")):
    os.remove("unsorted_sorted.pdb")
  out = StringIO()
  sort_hetatms.run(
    args=[pdb_file, "--verbose"],
    out=out)
  assert ("""Residue group pdb=" O   HOH    10 " added to chain A (distance = 2.814, symop = -x-1,y-1/2,-z)""" in out.getvalue())
  hierarchy = iotbx.pdb.input("unsorted_sorted.pdb").construct_hierarchy()
  chains = hierarchy.models()[0].chains()
  assert (len(chains) == 4)
  rgsA = chains[-2].residue_groups()
  rgsB = chains[-1].residue_groups()
  assert (len(rgsA) == 3) and (len(rgsB) == 4)
  sort_hetatms.run(
    args=[pdb_file, "preserve_chain_id=True", "renumber=False"],
    out=StringIO())
  hierarchy = iotbx.pdb.input("unsorted_sorted.pdb").construct_hierarchy()
  chains = hierarchy.models()[0].chains()
  assert (len(chains) == 3)
  assert (chains[-1].id == "A")
  pdb_raw2 = """\
CRYST1   21.937   10.000   23.477  90.00 107.08  90.00 P 1 21 1      2
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
TER      21      ASN A   3
HETATM   61  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM   62  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM   63  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM   64  O   HOH A  11       6.127   3.965  10.368  1.00 23.99           O
HETATM   65  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
HETATM   66  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM   67  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
END
"""
  pdb_file = "unsorted2.pdb"
  with open(pdb_file, "w") as f:
    f.write(pdb_raw2)
  if (op.exists("unsorted2_sorted.pdb")):
    os.remove("unsorted2_sorted.pdb")
  out = StringIO()
  sort_hetatms.run(
    args=[pdb_file, "--verbose", "--remove_waters_outside_radius"],
    out=out)
  assert ("Water     11  is not near any polymer chain, will delete" in
          out.getvalue())
  out = StringIO()
  sort_hetatms.run(
    args=[pdb_file, "--verbose",],
    out=out)
  assert ("Residue group A   11  is not near any macromolecule chain" in
          out.getvalue())
  pdb_raw = """\
HEADER    LIGASE                                  31-AUG-98   1BS2
CRYST1  100.374  100.374  204.341  90.00  90.00  90.00 P 43 21 2     8
ATOM   1094  N   LYS A 143      22.992  81.325 -18.373  1.00 63.31           N
ATOM   1095  CA  LYS A 143      24.003  82.214 -18.922  1.00 53.83           C
ATOM   1096  C   LYS A 143      24.725  81.433 -20.002  1.00 49.46           C
ATOM   1097  O   LYS A 143      25.026  80.249 -19.833  1.00 48.59           O
ATOM   1098  CB  LYS A 143      24.980  82.658 -17.835  1.00 58.17           C
ATOM   1099  CG  LYS A 143      24.831  84.114 -17.432  1.00 65.92           C
ATOM   1100  CD  LYS A 143      23.422  84.434 -16.967  1.00 69.21           C
ATOM   1101  CE  LYS A 143      23.253  85.929 -16.742  1.00 67.51           C
ATOM   1102  NZ  LYS A 143      21.905  86.254 -16.193  1.00 74.33           N
ATOM   1103  N   LYS A 144      24.985  82.095 -21.120  1.00 46.52           N
ATOM   1104  CA  LYS A 144      25.641  81.453 -22.246  1.00 48.43           C
ATOM   1105  C   LYS A 144      26.946  82.162 -22.583  1.00 49.19           C
ATOM   1106  O   LYS A 144      27.014  83.395 -22.623  1.00 47.13           O
ATOM   1107  CB  LYS A 144      24.704  81.475 -23.458  1.00 55.29           C
ATOM   1108  CG  LYS A 144      25.183  80.676 -24.659  1.00 63.36           C
ATOM   1109  CD  LYS A 144      24.334  80.964 -25.899  1.00 61.42           C
ATOM   1110  CE  LYS A 144      22.874  80.579 -25.705  1.00 58.72           C
ATOM   1111  NZ  LYS A 144      22.078  80.793 -26.946  1.00 61.59           N
ATOM   1112  N   VAL A 145      27.986  81.373 -22.819  1.00 42.70           N
ATOM   1113  CA  VAL A 145      29.280  81.932 -23.158  1.00 38.93           C
ATOM   1114  C   VAL A 145      29.869  81.233 -24.374  1.00 40.05           C
ATOM   1115  O   VAL A 145      29.715  80.021 -24.553  1.00 45.68           O
ATOM   1116  CB  VAL A 145      30.274  81.815 -21.972  1.00 30.33           C
ATOM   1117  CG1 VAL A 145      30.431  80.356 -21.569  1.00 35.25           C
ATOM   1118  CG2 VAL A 145      31.627  82.424 -22.351  1.00 15.32           C
ATOM   1119  N   ILE A 146      30.519  82.013 -25.225  1.00 36.19           N
ATOM   1120  CA  ILE A 146      31.157  81.457 -26.399  1.00 42.19           C
ATOM   1121  C   ILE A 146      32.640  81.522 -26.106  1.00 46.89           C
ATOM   1122  O   ILE A 146      33.161  82.577 -25.735  1.00 42.54           O
ATOM   1123  CB  ILE A 146      30.886  82.283 -27.666  1.00 43.82           C
ATOM   1124  CG1 ILE A 146      29.406  82.226 -28.030  1.00 46.47           C
ATOM   1125  CG2 ILE A 146      31.723  81.741 -28.816  1.00 39.50           C
ATOM   1126  CD1 ILE A 146      29.048  83.115 -29.197  1.00 46.73           C
ATOM   1127  N   ILE A 147      33.311  80.387 -26.250  1.00 47.64           N
ATOM   1128  CA  ILE A 147      34.740  80.323 -26.019  1.00 47.32           C
ATOM   1129  C   ILE A 147      35.406  79.905 -27.326  1.00 49.89           C
ATOM   1130  O   ILE A 147      35.276  78.766 -27.777  1.00 46.22           O
ATOM   1131  CB  ILE A 147      35.078  79.325 -24.887  1.00 43.61           C
ATOM   1132  CG1 ILE A 147      34.402  79.775 -23.586  1.00 40.57           C
ATOM   1133  CG2 ILE A 147      36.585  79.257 -24.688  1.00 40.81           C
ATOM   1134  CD1 ILE A 147      34.584  78.825 -22.424  1.00 34.81           C
ATOM   1135  N   GLU A 148      36.102  80.855 -27.940  1.00 49.93           N
ATOM   1136  CA  GLU A 148      36.788  80.613 -29.199  1.00 50.03           C
ATOM   1137  C   GLU A 148      38.263  80.328 -28.913  1.00 44.28           C
ATOM   1138  O   GLU A 148      38.943  81.145 -28.301  1.00 43.55           O
ATOM   1139  CB  GLU A 148      36.646  81.846 -30.106  1.00 47.09           C
ATOM   1140  CG  GLU A 148      37.243  81.679 -31.489  1.00 40.60           C
ATOM   1141  CD  GLU A 148      37.982  82.916 -31.943  1.00 48.19           C
ATOM   1142  OE1 GLU A 148      37.324  83.935 -32.246  1.00 52.87           O
ATOM   1143  OE2 GLU A 148      39.228  82.873 -31.982  1.00 45.83           O
ATOM   1144  N   PHE A 149      38.754  79.172 -29.347  1.00 37.67           N
ATOM   1145  CA  PHE A 149      40.147  78.828 -29.106  1.00 39.79           C
ATOM   1146  C   PHE A 149      40.775  77.925 -30.159  1.00 42.70           C
ATOM   1147  O   PHE A 149      40.079  77.314 -30.966  1.00 45.66           O
ATOM   1148  CB  PHE A 149      40.305  78.199 -27.718  1.00 37.53           C
ATOM   1149  CG  PHE A 149      39.387  77.039 -27.461  1.00 33.29           C
ATOM   1150  CD1 PHE A 149      38.009  77.226 -27.387  1.00 29.22           C
ATOM   1151  CD2 PHE A 149      39.903  75.760 -27.270  1.00 28.26           C
ATOM   1152  CE1 PHE A 149      37.156  76.154 -27.124  1.00 26.87           C
ATOM   1153  CE2 PHE A 149      39.062  74.681 -27.006  1.00 29.29           C
ATOM   1154  CZ  PHE A 149      37.683  74.879 -26.934  1.00 30.14           C
ATOM   1155  N   SER A 150      42.106  77.852 -30.120  1.00 44.78           N
ATOM   1156  CA  SER A 150      42.924  77.077 -31.057  1.00 38.94           C
ATOM   1157  C   SER A 150      42.996  77.839 -32.378  1.00 35.66           C
ATOM   1158  O   SER A 150      44.039  78.374 -32.718  1.00 38.02           O
ATOM   1159  CB  SER A 150      42.353  75.676 -31.283  1.00 32.93           C
ATOM   1160  OG  SER A 150      43.261  74.895 -32.047  1.00 31.61           O
ATOM   1161  N   SER A 151      41.879  77.879 -33.104  1.00 36.86           N
ATOM   1162  CA  SER A 151      41.745  78.598 -34.377  1.00 32.32           C
ATOM   1163  C   SER A 151      42.958  78.629 -35.304  1.00 36.45           C
ATOM   1164  O   SER A 151      43.531  79.692 -35.565  1.00 35.55           O
ATOM   1165  CB  SER A 151      41.305  80.038 -34.107  1.00 27.79           C
ATOM   1166  OG  SER A 151      40.065  80.069 -33.429  1.00 35.13           O
ATOM   1167  N   PRO A 152      43.356  77.468 -35.836  1.00 40.61           N
ATOM   1168  CA  PRO A 152      44.512  77.444 -36.736  1.00 42.49           C
ATOM   1169  C   PRO A 152      44.150  77.973 -38.125  1.00 42.14           C
ATOM   1170  O   PRO A 152      42.991  78.288 -38.397  1.00 40.80           O
ATOM   1171  CB  PRO A 152      44.883  75.968 -36.760  1.00 42.76           C
ATOM   1172  CG  PRO A 152      43.534  75.304 -36.693  1.00 42.42           C
ATOM   1173  CD  PRO A 152      42.812  76.112 -35.631  1.00 39.19           C
ATOM   1174  N   ASN A 153      45.148  78.086 -38.993  1.00 44.89           N
ATOM   1175  CA  ASN A 153      44.921  78.542 -40.362  1.00 45.11           C
ATOM   1176  C   ASN A 153      45.092  77.307 -41.232  1.00 46.15           C
ATOM   1177  O   ASN A 153      45.996  76.502 -41.003  1.00 50.16           O
ATOM   1178  CB  ASN A 153      45.951  79.592 -40.770  1.00 45.30           C
ATOM   1179  CG  ASN A 153      46.044  80.729 -39.781  1.00 43.47           C
ATOM   1180  OD1 ASN A 153      45.093  81.487 -39.590  1.00 46.50           O
ATOM   1181  ND2 ASN A 153      47.196  80.853 -39.143  1.00 33.35           N
TER    4875      MET A 607
HETATM 4876  N   ARG X 900      43.972  82.052 -36.884  1.00 37.31           N
HETATM 4877  CA  ARG X 900      44.639  83.205 -36.209  1.00 43.51           C
HETATM 4878  C   ARG X 900      46.150  83.005 -36.147  1.00 44.70           C
HETATM 4879  O   ARG X 900      46.842  83.928 -35.666  1.00 51.15           O
HETATM 4880  CB  ARG X 900      44.080  83.392 -34.788  1.00 31.47           C
HETATM 4881  CG  ARG X 900      42.619  83.803 -34.755  1.00 31.50           C
HETATM 4882  CD  ARG X 900      42.116  83.982 -33.337  1.00 38.70           C
HETATM 4883  NE  ARG X 900      40.700  84.343 -33.288  1.00 42.91           N
HETATM 4884  CZ  ARG X 900      40.206  85.521 -33.660  1.00 48.35           C
HETATM 4885  NH1 ARG X 900      41.014  86.471 -34.114  1.00 45.64           N
HETATM 4886  NH2 ARG X 900      38.899  85.751 -33.584  1.00 46.85           N
HETATM 4887  OXT ARG X 900      46.623  81.932 -36.581  1.00 42.67           O
HETATM 4963  O   HOH X 976      48.036  82.421 -33.012  1.00 30.31           O
HETATM 4964  O   HOH X 977      44.190  81.253 -31.435  1.00 48.82           O
HETATM 4965  O   HOH X 978      41.198  81.209 -30.407  1.00 48.77           O
"""
  pdb_file = "unsorted3.pdb"
  with open(pdb_file, "w") as f:
    f.write(pdb_raw)
  if (op.exists("unsorted3_sorted.pdb")):
    os.remove("unsorted3_sorted.pdb")
  out = StringIO()
  sort_hetatms.run(
    args=[pdb_file, "--verbose",],
    out=out)
  assert op.isfile("unsorted3_sorted.pdb")
  hierarchy = iotbx.pdb.input("unsorted3_sorted.pdb").construct_hierarchy()
  for atom in hierarchy.atoms():
    assert atom.fetch_labels().chain_id == "A"
  # now with ARG ligand minus waters, but still flagged as HETATM
  pdb_in = iotbx.pdb.input("unsorted3.pdb")
  hierarchy = pdb_in.construct_hierarchy()
  sel_str = "not resname HOH"
  sel_cache = hierarchy.atom_selection_cache()
  sel = sel_cache.selection(sel_str)
  hierarchy_new = hierarchy.select(sel)
  with open("unsorted4.pdb", "w") as f:
    f.write(hierarchy_new.as_pdb_string(
      crystal_symmetry=pdb_in.crystal_symmetry()))
  if (op.exists("unsorted4_sorted.pdb")):
    os.remove("unsorted4_sorted.pdb")
  out = StringIO()
  sort_hetatms.run(
    args=["unsorted4.pdb", "--verbose",],
    out=out)
  assert op.isfile("unsorted4_sorted.pdb")
  hierarchy = iotbx.pdb.input("unsorted4_sorted.pdb").construct_hierarchy()
  for atom in hierarchy.atoms():
    assert atom.fetch_labels().chain_id == "A"
  print("OK")

if (__name__ == "__main__"):
  exercise()
