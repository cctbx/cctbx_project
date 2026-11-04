"""Test flipping symmetric residues"""
from __future__ import absolute_import, division, print_function
import sys
from six.moves import cStringIO as StringIO
from libtbx import easy_run
import iotbx.pdb
import iotbx.pdb.hierarchy

pdbs = {
  "ASP" : """
ATOM      1  N   ASP L   1      55.410  -0.618  -1.595  1.00 33.13           N
ATOM      2  CA  ASP L   1      55.727  -2.042  -1.664  1.00 50.61           C
ATOM      3  C   ASP L   1      56.021  -2.621  -0.277  1.00 42.48           C
ATOM      4  O   ASP L   1      56.306  -1.887   0.661  1.00 36.07           O
ATOM      5  CB  ASP L   1      56.856  -2.319  -2.677  1.00 50.51           C
ATOM      6  CG  ASP L   1      58.182  -1.647  -2.310  1.00 51.45           C
ATOM      7  OD1 ASP L   1      59.243  -2.181  -2.696  1.00 52.89           O
ATOM      8  OD2 ASP L   1      58.188  -0.585  -1.659  1.00 51.81           O
ATOM      9  H1  ASP L   1      55.244  -0.311  -2.414  1.00 33.13           H
ATOM     10  H2  ASP L   1      54.692  -0.498  -1.082  1.00 33.13           H
ATOM     11  H3  ASP L   1      56.100  -0.177  -1.247  1.00 33.13           H
ATOM     12  HA  ASP L   1      54.941  -2.506  -1.993  1.00 50.61           H
ATOM     13  HB2 ASP L   1      56.994  -3.277  -2.743  1.00 50.51           H
ATOM     14  HB3 ASP L   1      56.577  -2.012  -3.554  1.00 50.51           H
""",
  "GLU" : """
ATOM   1691  N   GLU B   1     -45.832-106.231  -7.917  1.00 56.97           N
ATOM   1692  CA  GLU B   1     -45.794-105.198  -8.995  1.00 52.47           C
ATOM   1693  C   GLU B   1     -45.700-103.803  -8.381  1.00 46.38           C
ATOM   1694  O   GLU B   1     -46.454-103.456  -7.470  1.00 44.22           O
ATOM   1695  CB  GLU B   1     -47.033-105.306  -9.893  1.00 63.38           C
ATOM   1696  CG  GLU B   1     -46.783-104.967 -11.365  1.00 61.17           C
ATOM   1697  CD  GLU B   1     -46.547-103.486 -11.609  1.00 56.72           C
ATOM   1698  OE1 GLU B   1     -46.005-103.137 -12.680  1.00 64.27           O
ATOM   1699  OE2 GLU B   1     -46.904-102.668 -10.735  1.00 57.61           O
ATOM   1700  H1  GLU B   1     -46.103-107.006  -8.258  1.00 68.38           H
ATOM   1701  H2  GLU B   1     -45.018-106.327  -7.570  1.00 68.38           H
ATOM   1702  H3  GLU B   1     -46.397-105.973  -7.279  1.00 68.38           H
ATOM   1703  HA  GLU B   1     -45.009-105.340  -9.546  1.00 62.98           H
ATOM   1704  HB2 GLU B   1     -47.366-106.216  -9.853  1.00 76.07           H
ATOM   1705  HB3 GLU B   1     -47.710-104.695  -9.561  1.00 76.07           H
ATOM   1706  HG2 GLU B   1     -45.998-105.449 -11.670  1.00 73.41           H
ATOM   1707  HG3 GLU B   1     -47.557-105.235 -11.885  1.00 73.41           H
""",
  "PHE" : """
ATOM   1206  N   PHE L  83      77.038 -23.401   3.393  1.00 55.55           N
ATOM   1207  CA  PHE L  83      77.649 -22.073   3.419  1.00 53.48           C
ATOM   1208  C   PHE L  83      77.083 -21.233   2.280  1.00 48.61           C
ATOM   1209  O   PHE L  83      77.325 -21.514   1.106  1.00 54.11           O
ATOM   1210  CB  PHE L  83      79.177 -22.165   3.306  1.00 55.82           C
ATOM   1211  CG  PHE L  83      79.833 -20.881   2.882  1.00 60.50           C
ATOM   1212  CD1 PHE L  83      80.317 -20.728   1.593  1.00 64.71           C
ATOM   1213  CD2 PHE L  83      79.954 -19.821   3.767  1.00 65.07           C
ATOM   1214  CE1 PHE L  83      80.912 -19.543   1.196  1.00 68.67           C
ATOM   1215  CE2 PHE L  83      80.549 -18.634   3.373  1.00 68.71           C
ATOM   1216  CZ  PHE L  83      81.028 -18.497   2.086  1.00 69.85           C
ATOM   1217  H   PHE L  83      77.591 -24.047   3.262  1.00 55.55           H
ATOM   1218  HA  PHE L  83      77.440 -21.651   4.267  1.00 53.48           H
ATOM   1219  HB2 PHE L  83      79.540 -22.436   4.163  1.00 55.82           H
ATOM   1220  HB3 PHE L  83      79.404 -22.860   2.669  1.00 55.82           H
ATOM   1221  HD1 PHE L  83      80.241 -21.430   0.987  1.00 64.71           H
ATOM   1222  HD2 PHE L  83      79.632 -19.908   4.635  1.00 65.07           H
ATOM   1223  HE1 PHE L  83      81.234 -19.453   0.328  1.00 68.67           H
ATOM   1224  HE2 PHE L  83      80.626 -17.930   3.975  1.00 68.71           H
ATOM   1225  HZ  PHE L  83      81.428 -17.701   1.820  1.00 69.85           H
""",
  "TYR" : """
ATOM   1271  N   TYR L  87      67.948 -15.639   0.990  1.00 39.05           N
ATOM   1272  CA  TYR L  87      66.992 -14.569   0.769  1.00 35.23           C
ATOM   1273  C   TYR L  87      65.567 -15.103   0.846  1.00 35.65           C
ATOM   1274  O   TYR L  87      65.270 -16.177   0.327  1.00 40.62           O
ATOM   1275  CB  TYR L  87      67.220 -13.932  -0.604  1.00 33.62           C
ATOM   1276  CG  TYR L  87      68.586 -13.316  -0.803  1.00 33.61           C
ATOM   1277  CD1 TYR L  87      68.792 -11.955  -0.608  1.00 30.41           C
ATOM   1278  CD2 TYR L  87      69.669 -14.092  -1.206  1.00 36.13           C
ATOM   1279  CE1 TYR L  87      70.042 -11.388  -0.799  1.00 36.04           C
ATOM   1280  CE2 TYR L  87      70.920 -13.529  -1.399  1.00 38.98           C
ATOM   1281  CZ  TYR L  87      71.099 -12.181  -1.194  1.00 40.09           C
ATOM   1282  OH  TYR L  87      72.345 -11.628  -1.387  1.00 41.66           O
ATOM   1283  H   TYR L  87      67.946 -16.244   0.379  1.00 39.05           H
ATOM   1284  HA  TYR L  87      67.119 -13.901   1.460  1.00 35.23           H
ATOM   1285  HB2 TYR L  87      67.082 -14.608  -1.286  1.00 33.62           H
ATOM   1286  HB3 TYR L  87      66.548 -13.247  -0.744  1.00 33.62           H
ATOM   1287  HD1 TYR L  87      68.080 -11.417  -0.345  1.00 30.41           H
ATOM   1288  HD2 TYR L  87      69.551 -15.003  -1.348  1.00 36.13           H
ATOM   1289  HE1 TYR L  87      70.168 -10.477  -0.661  1.00 36.04           H
ATOM   1290  HE2 TYR L  87      71.635 -14.060  -1.666  1.00 38.98           H
ATOM   1291  HH  TYR L  87      72.884 -12.227  -1.623  1.00 41.66           H
""",
  "ARG" : """
ATOM   6724  N   ARG A 465      -8.114  20.359  -8.988  1.00 49.56           N
ATOM   6725  CA  ARG A 465      -8.113  21.468  -8.039  1.00 52.45           C
ATOM   6726  C   ARG A 465      -9.335  21.427  -7.136  1.00 54.28           C
ATOM   6727  O   ARG A 465      -9.283  21.899  -5.996  1.00 54.38           O
ATOM   6728  CB  ARG A 465      -8.072  22.802  -8.784  1.00 53.89           C
ATOM   6729  CG  ARG A 465      -6.707  23.162  -9.350  1.00 55.53           C
ATOM   6730  CD  ARG A 465      -6.579  24.655  -9.586  1.00 57.77           C
ATOM   6731  NE  ARG A 465      -6.803  25.405  -8.354  1.00 59.92           N
ATOM   6732  CZ  ARG A 465      -6.632  26.716  -8.223  1.00 61.69           C
ATOM   6733  NH1 ARG A 465      -6.849  27.289  -7.047  1.00 62.84           N
ATOM   6734  NH2 ARG A 465      -6.242  27.453  -9.254  1.00 62.15           N
ATOM   6735  H   ARG A 465      -8.015  20.588  -9.811  1.00 49.56           H
ATOM   6736  HA  ARG A 465      -7.321  21.380  -7.486  1.00 52.45           H
ATOM   6737  HB2 ARG A 465      -8.714  22.775  -9.510  1.00 53.89           H
ATOM   6738  HB3 ARG A 465      -8.356  23.506  -8.180  1.00 53.89           H
ATOM   6739  HG2 ARG A 465      -6.014  22.870  -8.737  1.00 55.53           H
ATOM   6740  HG3 ARG A 465      -6.566  22.688 -10.185  1.00 55.53           H
ATOM   6741  HD2 ARG A 465      -5.696  24.856  -9.934  1.00 57.77           H
ATOM   6742  HD3 ARG A 465      -7.219  24.935 -10.259  1.00 57.77           H
ATOM   6743  HE  ARG A 465      -7.065  24.966  -7.662  1.00 59.92           H
ATOM   6744 HH11 ARG A 465      -7.099  26.812  -6.376  1.00 62.84           H
ATOM   6745 HH12 ARG A 465      -6.740  28.137  -6.955  1.00 62.84           H
ATOM   6746 HH21 ARG A 465      -6.097  27.083 -10.017  1.00 62.15           H
ATOM   6747 HH22 ARG A 465      -6.134  28.301  -9.160  1.00 62.15           H
""",
  "PHE_missing_CE2" : """
ATOM    810  N   PHE A 105      27.518 -30.761  62.509  1.00 24.66           N
ANISOU  810  N   PHE A 105     3396   2930   3044     46    490     24       N
ATOM    811  CA  PHE A 105      27.910 -31.927  61.750  1.00 22.63           C
ANISOU  811  CA  PHE A 105     3141   2671   2786     43    489     25       C
ATOM    812  C   PHE A 105      26.878 -32.113  60.659  1.00 22.52           C
ANISOU  812  C   PHE A 105     3123   2661   2773     39    491     25       C
ATOM    813  O   PHE A 105      25.759 -31.615  60.765  1.00 23.79           O
ANISOU  813  O   PHE A 105     3280   2824   2934     38    493     26       O
ATOM    814  CB  PHE A 105      27.958 -33.165  62.645  1.00 23.99           C
ANISOU  814  CB  PHE A 105     3321   2841   2955     42    491     24       C
ATOM    815  CG  PHE A 105      26.609 -33.618  63.120  1.00 28.62           C
ANISOU  815  CG  PHE A 105     3906   3429   3538     38    495     25       C
ATOM    816  CD1 PHE A 105      26.093 -33.177  64.320  1.00 30.02           C
ANISOU  816  CD1 PHE A 105     4085   3607   3715     40    497     24       C
ATOM    817  CD2 PHE A 105      25.841 -34.478  62.353  1.00 30.68           C
ANISOU  817  CD2 PHE A 105     4166   3691   3798     34    497     25       C
ATOM    818  CE1 PHE A 105      24.849 -33.596  64.733  1.00 32.03           C
ANISOU  818  CE1 PHE A 105     4339   3863   3968     36    500     24       C
""",
  "PHE_2" : """
ATOM   3264  N   PHE A 217      30.776 -28.100 -17.162  1.00 15.42           N
ATOM   3265  CA  PHE A 217      31.363 -27.624 -18.384  1.00 14.70           C
ATOM   3266  C   PHE A 217      32.212 -28.677 -19.046  1.00 15.01           C
ATOM   3267  O   PHE A 217      33.021 -29.350 -18.413  1.00 16.31           O
ATOM   3268  CB  PHE A 217      32.177 -26.335 -18.136  1.00 16.76           C
ATOM   3269  CG  PHE A 217      31.302 -25.242 -17.626  1.00 16.55           C
ATOM   3270  CD1 PHE A 217      31.215 -24.983 -16.284  1.00 16.83           C
ATOM   3271  CD2 PHE A 217      30.503 -24.516 -18.502  1.00 18.31           C
ATOM   3272  CE1 PHE A 217      30.337 -24.029 -15.786  1.00 16.65           C
ATOM   3273  CE2 PHE A 217      29.624 -23.532 -18.011  1.00 15.24           C
ATOM   3274  CZ  PHE A 217      29.561 -23.283 -16.645  1.00 17.98           C
ATOM   3275  H   PHE A 217      31.327 -28.470 -16.616  1.00 18.51           H
ATOM   3276  HA  PHE A 217      30.648 -27.401 -18.999  1.00 17.64           H
ATOM   3277  HB2 PHE A 217      32.864 -26.510 -17.475  1.00 20.11           H
ATOM   3278  HB3 PHE A 217      32.577 -26.042 -18.970  1.00 20.11           H
ATOM   3279  HD1 PHE A 217      31.726 -25.485 -15.690  1.00 20.20           H
ATOM   3280  HD2 PHE A 217      30.540 -24.690 -19.415  1.00 21.98           H
ATOM   3281  HE1 PHE A 217      30.310 -23.859 -14.872  1.00 19.99           H
ATOM   3282  HE2 PHE A 217      29.101 -23.038 -18.600  1.00 18.29           H
ATOM   3283  HZ  PHE A 217      28.992 -22.626 -16.313  1.00 21.58           H
ATOM   4166  CD1 TYR A 277      41.891 -14.126 -14.217  1.00 26.09           C
""",
  'LEU' : '''
ATOM   1142  N   LEU B 158      12.672  -5.202  21.045  1.00 12.10      B    N
ATOM   1143  CA  LEU B 158      12.719  -4.364  19.855  1.00  9.96      B    C
ATOM   1144  C   LEU B 158      14.074  -3.681  19.767  1.00 13.93      B    C
ATOM   1145  O   LEU B 158      14.542  -3.091  20.750  1.00 10.55      B    O
ATOM   1146  CB  LEU B 158      11.594  -3.321  19.891  1.00  9.63      B    C
ATOM   1147  CG  LEU B 158      11.380  -2.482  18.634  1.00 11.09      B    C
ATOM   1148  CD1 LEU B 158       9.897  -2.151  18.478  1.00 11.88      B    C
ATOM   1149  CD2 LEU B 158      11.900  -3.216  17.391  1.00 15.14      B    C
''',
  'VAL' : '''
HETATM    1  N   VAL A   1      -0.458  -1.800  -0.548  1.00 20.00      A    N
HETATM    2  CA  VAL A   1      -0.329  -0.702   0.393  1.00 20.00      A    C
HETATM    3  C   VAL A   1      -1.671  -0.463   1.082  1.00 20.00      A    C
HETATM    4  O   VAL A   1      -2.646  -0.015   0.423  1.00 20.00      A    O
HETATM    5  CB  VAL A   1       0.092   0.561  -0.354  1.00 20.00      A    C
HETATM    6  CG2 VAL A   1       0.584   1.604   0.647  1.00 20.00      A    C
HETATM    7  CG1 VAL A   1       1.215   0.224  -1.332  1.00 20.00      A    C
HETATM    8  OXT VAL A   1      -1.802  -0.713   2.309  1.00 20.00      A    O-1
HETATM    9  H   VAL A   1      -1.391  -2.165  -0.511  1.00 20.00      A    H
HETATM   10  H2  VAL A   1       0.193  -2.524  -0.308  1.00 20.00      A    H
HETATM   11  HA  VAL A   1       0.419  -0.949   1.136  1.00 20.00      A    H
HETATM   12  HB  VAL A   1      -0.755   0.957  -0.900  1.00 20.00      A    H
HETATM   13 HG21 VAL A   1      -0.219   1.854   1.334  1.00 20.00      A    H
HETATM   14 HG22 VAL A   1       0.897   2.498   0.115  1.00 20.00      A    H
HETATM   15 HG23 VAL A   1       1.426   1.202   1.204  1.00 20.00      A    H
HETATM   16 HG11 VAL A   1       0.813  -0.349  -2.162  1.00 20.00      A    H
HETATM   17 HG12 VAL A   1       1.975  -0.362  -0.822  1.00 20.00      A    H
HETATM   18 HG13 VAL A   1       1.657   1.143  -1.707  1.00 20.00      A    H
''',
}

geo_strs = {
  "ASP" : """dihedral pdb=" CA  ASP L   1 "
           pdb=" CB  ASP L   1 "
           pdb=" CG  ASP L   1 "
           pdb=" OD1 ASP L   1 "
""",
  "GLU" : """  dihedral pdb=" CB  GLU B   1 "
           pdb=" CG  GLU B   1 "
           pdb=" CD  GLU B   1 "
           pdb=" OE1 GLU B   1 "
""",
  "PHE" : """dihedral pdb=" CA  PHE L  83 "
           pdb=" CB  PHE L  83 "
           pdb=" CG  PHE L  83 "
           pdb=" CD1 PHE L  83 "
""",
  "TYR" : """dihedral pdb=" CA  TYR L  87 "
           pdb=" CB  TYR L  87 "
           pdb=" CG  TYR L  87 "
           pdb=" CD1 TYR L  87 "
""",
  "ARG" : """dihedral pdb=" CD  ARG A 465 "
           pdb=" NE  ARG A 465 "
           pdb=" CZ  ARG A 465 "
           pdb=" NH1 ARG A 465 "
""",
  "PHE_missing_CE2" : 'Residue "A PHE  105": not complete - not flipped',
  "PHE_2" : """dihedral pdb=" CA  PHE A 217 "
           pdb=" CB  PHE A 217 "
           pdb=" CG  PHE A 217 "
           pdb=" CD1 PHE A 217 "
""",
  'LEU' : '''chirality pdb=" CG  LEU B 158 " segid="B   "
            pdb=" CB  LEU B 158 " segid="B   "
            pdb=" CD1 LEU B 158 " segid="B   "
            pdb=" CD2 LEU B 158 " segid="B   "
''',
  'VAL' : '''chirality pdb=" CB  VAL A   1 " segid="A   "
            pdb=" CA  VAL A   1 " segid="A   "
            pdb=" CG1 VAL A   1 " segid="A   "
            pdb=" CG2 VAL A   1 " segid="A   "
''',
}
values = {
  "ASP" : [ 152,  -82, -999, -999],
  "GLU" : [ 162, -162,  -18,   18],
  "PHE" : [-105,   15,   74,   16],
  "TYR" : [ -96,    6,   85,    5],
  "ARG" : [-179,  179, -999, -999],
  "PHE_missing_CE2" : [-999,-999,-999,-999],
  "PHE_2" : [97,   -7,  -79,  -11],
  # chiral
  'LEU' : [2.65, -5.24, -2.65, 0.06],
  'VAL' : [2.65, -5.37, -2.75, 0.11],
}

def largest_torsion(s, f, i):
  t = s.split(f)
  if len(t)==1: return -999
  s=t[1].splitlines()
  f="ideal   model   delta sinusoidal    sigma   weight residual"
  reading = False
  for line in s:
    print(line)
    if reading: break
    if line.find(f)>-1:
      reading=True
  d=line.split()
  print(d)
  return float(d[i])

def largest_chiral(s, f, i):
  t = s.split(f)
  if len(t)==1: return -999
  s=t[1].splitlines()
  f="both_signs  ideal   model   delta    sigma   weight residual"
  reading = False
  for line in s:
    print(line)
    if reading: break
    if line.find(f)>-1:
      reading=True
  d=line.split()
  print(d)
  return float(d[i])

def run():
  for code, lines in pdbs.items():
    #if code!='VAL': continue
    print(('-%s- ' % code)*10)
    with open("tst_symmetric_flips_%s.pdb" % code, "w") as f:
      f.write(lines)
    cmd = "phenix.pdb_interpretation tst_symmetric_flips_%s.pdb" % code
    cmd += ' flip_symmetric_amino_acids=None'
    print(cmd)
    ero = easy_run.fully_buffered(command=cmd)
    out = StringIO()
    ero.show_stdout(out=out)
    i=0
    if code in ['LEU', 'VAL']:
      d = largest_chiral(out.getvalue(),
                         geo_strs[code],
                         i+2)
      print('chiral',d)
      assert abs(d-values[code][i])<.1, '%s %s' % (d, values[code][i])
    else:
      d = largest_torsion(out.getvalue(),
                          geo_strs[code],
                          i+1)
      print('torsion',d)
      assert abs(d-values[code][i])<1, '%s %s' % (d, values[code][i])
    i+=1
    if code in ['LEU', 'VAL']:
      d = largest_chiral(out.getvalue(),
                         geo_strs[code],
                         i+2)
      print('delta',d)
      assert abs(d-values[code][i])<.1, '%s %s' % (d, values[code][i])
    else:
      d = largest_torsion(out.getvalue(),
                          geo_strs[code],
                          i+1)
      print('delta',d)
      assert abs(d-values[code][i])<1, '%s %s' % (d, values[code][i])

    pdb_inp = iotbx.pdb.input('tst_symmetric_flips_%s.pdb' % code)
    hierarchy = pdb_inp.construct_hierarchy()
    rc = hierarchy.flip_symmetric_amino_acids()
    print(rc)
    hierarchy.write_pdb_file("tst_symmetric_flips_%s.pdb" % code)
    cmd = "phenix.pdb_interpretation tst_symmetric_flips_%s.pdb" % code
    cmd += ' flip_symmetric_amino_acids=None'
    print(cmd)
    ero = easy_run.fully_buffered(command=cmd)
    out = StringIO()
    ero.show_stdout(out=out)

    i+=1
    if code in ['LEU', 'VAL']:
      d = largest_chiral(out.getvalue(),
                         geo_strs[code],
                         i)
      print('chiral',d)
      assert abs(d-values[code][i])<.1, '%s %s' % (d, values[code][i])
    else:
      d = largest_torsion(out.getvalue(),
                          geo_strs[code],
                          i-1)
      print('torsion',d)
      assert abs(d-values[code][i])<1, '%s %s' % (d, values[code][i])
    i+=1
    if code in ['LEU', 'VAL']:
      d = largest_chiral(out.getvalue(),
                         geo_strs[code],
                         i)
      print('delta',d)
      assert abs(d-values[code][i])<.1, '%s %s' % (d, values[code][i])
    else:
      d = largest_torsion(out.getvalue(),
                          geo_strs[code],
                          i-1)
      print('delta',d)
      assert abs(d-values[code][i])<1, '%s %s' % (d, values[code][i])

def rows_of(code, keep=()):
  """Heavy atoms of pdbs[code] plus the H named in keep:
  [[name, altloc, xyz], ...]"""
  rows = []
  for line in pdbs[code].splitlines():
    if line[:6] not in ("ATOM  ", "HETATM"): continue
    name = line[12:16].strip()
    if line[76:78].strip() == "H" and name not in keep: continue
    xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    rows.append([name, "", xyz])
  return rows

def split(rows, names):
  """The atoms in names become conformers A and B, B 0.1 A away."""
  result = []
  for name, altloc, xyz in rows:
    if name in names:
      result.append([name, "A", xyz])
      result.append([name, "B", (xyz[0] + 0.1, xyz[1], xyz[2])])
    else:
      result.append([name, altloc, xyz])
  return result

def exchange(rows, altloc, name1, name2):
  """Exchange two names in one conformer: that turns its test around."""
  for row in rows:
    if row[1] != altloc: continue
    if row[0] == name1: row[0] = name2
    elif row[0] == name2: row[0] = name1

def turn(rows, altloc, name, axis1, axis2):
  """Turn one atom of a conformer 180 deg about axis1-axis2; the dihedral
  through that axis changes by 180 deg."""
  from scitbx import matrix
  sites = {}
  for row in rows:
    if row[1] in ("", altloc): sites[row[0]] = matrix.col(row[2])
  origin = sites[axis1]
  axis = (sites[axis2] - origin).normalize()
  r = axis.axis_and_angle_as_r3_rotation_matrix(180, deg=True)
  for row in rows:
    if row[0] == name and row[1] == altloc:
      row[2] = (r * (sites[name] - origin) + origin).elems

def add_h(rows, name, altloc, parent):
  """An H or D 1 A from parent, in that conformer."""
  for row in rows:
    if row[0] == parent:
      rows.append([name, altloc, (row[2][0], row[2][1], row[2][2] + 1)])
      return
  assert 0, parent

def make(code, rows, resname=None):
  """Hierarchy of one residue."""
  resname = resname or code
  lines = []
  for i, (name, altloc, xyz) in enumerate(rows):
    label = name if len(name) == 4 else " %-3s" % name
    occ = 1.0 if altloc == "" else 0.5
    lines.append(
      "ATOM  %5d %s%1s%3s A  11    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s" % (
      i + 1, label, altloc or " ", resname, xyz[0], xyz[1], xyz[2], occ,
      20 + i, name[0]))
  return iotbx.pdb.input(source_info=None,
    lines="\n".join(lines)).construct_hierarchy()

def flip(hierarchy):
  """Sites by (altloc, name) before and after, and the returned info."""
  def sites():
    result = {}
    for ag in hierarchy.only_residue_group().atom_groups():
      for atom in ag.atoms():
        result[ag.altloc.strip(), atom.name.strip()] = atom.xyz
    return result
  before = sites()
  info = hierarchy.flip_symmetric_amino_acids()
  return before, sites(), info

def assert_swapped(before, after, altloc, name1, name2):
  assert before[altloc, name1] == after[altloc, name2], (altloc, name1, name2)
  assert before[altloc, name2] == after[altloc, name1], (altloc, name1, name2)

def assert_only_swapped(before, after, swapped):
  """swapped: [(altloc, name1, name2), ...]; every other atom stays."""
  moved = set()
  for altloc, name1, name2 in swapped:
    assert_swapped(before, after, altloc, name1, name2)
    moved.add((altloc, name1))
    moved.add((altloc, name2))
  for key in before:
    if key not in moved:
      assert before[key] == after[key], key

def exercise_acid_h():
  # one H on an acid: the O must not swap away from it
  for code, parent, h, d, pair in (
      ("GLU", "OE2", "HE2", "DE2", ("OE1", "OE2")),
      ("ASP", "OD2", "HD2", "DD2", ("OD1", "OD2"))):
    for name in (h, d):
      rows = rows_of(code)
      add_h(rows, name, "", parent)
      before, after, info = flip(make(code, rows))
      assert_only_swapped(before, after, [])
      assert ('"A %s   11": protonated acid (one H on the O pair) - not flipped'
        % code) in info, info
    # both H: a complete pair goes with the O
    other = h[:-1] + "1"
    rows = rows_of(code)
    add_h(rows, h, "", parent)
    add_h(rows, other, "", pair[0])
    before, after, info = flip(make(code, rows))
    assert_only_swapped(before, after,
      [("",) + pair, ("", other, h)])
    assert "Time to flip 1 residue(s)" in info

def exercise_blank_o_acid_h_in_altlocs():
  # 6KK8 Glu A 167: O without altloc, HE2 in A, DE2 in B
  rows = rows_of("GLU")
  add_h(rows, "HE2", "A", "OE2")
  add_h(rows, "DE2", "B", "OE2")
  before, after, info = flip(make("GLU", rows))
  assert_only_swapped(before, after, [])
  assert '"A GLU   11": protonated acid (one H on the O pair) - not flipped' \
    in info, info

def exercise_a_flips_b_not():
  # each conformer on its own test; one residue in the count
  rows = split(rows_of("GLU"), ["CB", "CG", "CD", "OE1", "OE2"])
  exchange(rows, "B", "OE1", "OE2")
  hierarchy = make("GLU", rows)
  before, after, info = flip(hierarchy)
  assert_only_swapped(before, after, [("A", "OE1", "OE2")])
  assert "Time to flip 1 residue(s)" in info, info
  # the same with a ring and two swapped pairs
  ring = ["CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
  rows = split(rows_of("PHE"), ring)
  exchange(rows, "B", "CD1", "CD2")
  exchange(rows, "B", "CE1", "CE2")
  before, after, info = flip(make("PHE", rows))
  assert_only_swapped(before, after,
    [("A", "CD1", "CD2"), ("A", "CE1", "CE2")])
  assert "Time to flip 1 residue(s)" in info, info

def exercise_disagree_on_blank_pair():
  # a blank O belongs to both conformers: A says no, B says yes
  rows = split(rows_of("GLU"), ["CB"])
  turn(rows, "A", "CB", "CG", "CD")
  before, after, info = flip(make("GLU", rows))
  assert_only_swapped(before, after, [])
  assert '"A GLU   11": conformers disagree - not flipped' in info, info
  # the other way around: the report does not depend on the conformer order
  rows = split(rows_of("GLU"), ["CB"])
  turn(rows, "B", "CB", "CG", "CD")
  before, after, info = flip(make("GLU", rows))
  assert_only_swapped(before, after, [])
  assert '"A GLU   11": conformers disagree - not flipped' in info, info
  # both agree: the blank pair swaps once
  rows = split(rows_of("GLU"), ["CB"])
  before, after, info = flip(make("GLU", rows))
  assert_only_swapped(before, after, [("", "OE1", "OE2")])
  assert "Time to flip 1 residue(s)" in info, info

def exercise_test_spans_blank_and_altloc():
  names = ["CB", "CG", "OD1", "OD2"]
  # complete: CA blank, the rest in A and B
  before, after, info = flip(make("ASP", split(rows_of("ASP"), names)))
  assert_only_swapped(before, after,
    [("A", "OD1", "OD2"), ("B", "OD1", "OD2")])
  assert "Time to flip 1 residue(s)" in info, info
  # the same for a chiral test
  rows = split(rows_of("VAL"), ["CB", "CG1", "CG2"])
  before, after, info = flip(make("VAL", rows))
  assert_only_swapped(before, after,
    [("A", "CG1", "CG2"), ("B", "CG1", "CG2")])
  # an atom of the test missing in B: B cannot be tested, A flips
  rows = [r for r in split(rows_of("ASP"), names) if r[:2] != ["CG", "B"]]
  before, after, info = flip(make("ASP", rows))
  assert_only_swapped(before, after, [("A", "OD1", "OD2")])
  assert '"A ASP   11": test atom missing in altloc B - not flipped' in info, info
  # A has all of its test atoms in its own group, B has no CA: master renamed B
  # on A's decision
  rows = split(rows_of("ASP"), names)
  for row in rows:
    if row[0] == "CA": row[1] = "A"
  before, after, info = flip(make("ASP", rows))
  assert_only_swapped(before, after, [("A", "OD1", "OD2")])
  assert '"A ASP   11": test atom missing in altloc B - not flipped' in info, info
  # missing for every conformer: no decision, nothing to report
  rows = [r for r in split(rows_of("ASP"), names) if r[0] != "CG"]
  before, after, info = flip(make("ASP", rows))
  assert_only_swapped(before, after, [])
  assert info == "", info
  # a blank atom of the test missing with a blank pair: B not tested, none flips
  rows = split(rows_of("GLU"), ["CB", "CG"])
  rows = [r for r in rows if r[:2] != ["CB", "B"]]
  before, after, info = flip(make("GLU", rows))
  assert_only_swapped(before, after, [])
  assert '"A GLU   11": test atom missing in altloc B - not flipped' in info, info

def exercise_no_dihedral():
  # master: TypeError in abs(None)
  def on_ca(rows, altloc):
    ca = [r for r in rows if r[0] == "CA"][0]
    for row in rows:
      if row[0] == "CB" and row[1] == altloc: row[2] = ca[2]
    return rows
  # one group, CB on CA
  before, after, info = flip(make("ASP", on_ca(rows_of("ASP"), "")))
  assert_only_swapped(before, after, [])
  assert '"A ASP   11": test atoms coincident or collinear - not flipped' in \
    info, info
  assert "test atom missing" not in info, info
  # blank N, CA, C, O; CB of B on CA; A decides for itself
  rows = split(rows_of("ASP"), ["CB", "CG", "OD1", "OD2"])
  before, after, info = flip(make("ASP", on_ca(rows, "B")))
  assert_only_swapped(before, after, [("A", "OD1", "OD2")])
  assert ('"A ASP   11": test atoms coincident or collinear in altloc B - '
    'not flipped') in info, info
  # the same with a blank pair: nothing flips, as with a missing test atom
  rows = split(rows_of("ASP"), ["CB", "CG"])
  before, after, info = flip(make("ASP", on_ca(rows, "B")))
  assert_only_swapped(before, after, [])
  assert ('"A ASP   11": test atoms coincident or collinear in altloc B - '
    'not flipped') in info, info
  # collinear
  rows = rows_of("GLU")
  line = {"CB": (0, 0, 0), "CG": (1, 0, 0), "CD": (2, 0, 0)}
  for row in rows:
    if row[0] in line: row[2] = line[row[0]]
  before, after, info = flip(make("GLU", rows))
  assert_only_swapped(before, after, [])
  assert "test atoms coincident or collinear - not flipped" in info, info

# phenix_regression pdb/pdb1jxt.ent: B beyond 90 deg through CD1 and CD2
pdb_1jxt_phe_13 = """\
ATOM    230  CA  PHE A  13       5.831   6.291   5.009  1.00  2.47           C
ATOM    233  CB APHE A  13       7.031   6.945   5.693  0.65  3.27           C
ATOM    234  CB BPHE A  13       7.149   6.917   5.744  0.35  2.82           C
ATOM    235  CG APHE A  13       7.711   8.011   4.898  0.65  2.57           C
ATOM    236  CG BPHE A  13       7.863   7.833   4.788  0.35  1.62           C
ATOM    237  CD1APHE A  13       8.669   7.676   3.937  0.65  3.32           C
ATOM    238  CD1BPHE A  13       8.961   7.361   4.045  0.35  3.15           C
ATOM    239  CD2APHE A  13       7.366   9.356   5.067  0.65  3.05           C
ATOM    240  CD2BPHE A  13       7.578   9.214   4.755  0.35  2.83           C
ATOM    241  CE1APHE A  13       9.289   8.650   3.166  0.65  3.40           C
ATOM    242  CE1BPHE A  13       9.672   8.207   3.200  0.35  3.54           C
ATOM    243  CE2APHE A  13       8.015  10.359   4.340  0.65  2.96           C
ATOM    244  CE2BPHE A  13       8.246  10.091   3.879  0.35  2.36           C
ATOM    245  CZ APHE A  13       8.962   9.987   3.373  0.65  3.58           C
ATOM    246  CZ BPHE A  13       9.266   9.542   3.081  0.35  4.87           C
"""

def exercise_second_call_renames_nothing():
  hierarchy = iotbx.pdb.input(source_info=None,
    lines=pdb_1jxt_phe_13).construct_hierarchy()
  hierarchy.flip_symmetric_amino_acids()
  after_first = hierarchy.atoms().extract_xyz().deep_copy()
  for i in range(3):
    pairs = []
    hierarchy.flip_symmetric_amino_acids(swapped_pairs=pairs)
    assert pairs == [], pairs
    assert list(hierarchy.atoms().extract_xyz()) == list(after_first)
  # Val: CG2 moved along CB-CG2; volumes 0.27, 0, -0.27, -4.93 (ideal -2.5). The
  # delta exceeds 2 before and after a swap, so a rule on the delta alone swaps
  # them on every call.
  for shift in (-0.9, -1.0, -1.1, -2.8):
    rows = rows_of("VAL")
    sites = dict((r[0], r[2]) for r in rows)
    for row in rows:
      if row[0] == "CG2":
        row[2] = tuple(x + shift*(x - c)
          for x, c in zip(sites["CG2"], sites["CB"]))
    hierarchy = make("VAL", rows)
    hierarchy.flip_symmetric_amino_acids()
    pairs = []
    hierarchy.flip_symmetric_amino_acids(swapped_pairs=pairs)
    assert pairs == [], (shift, pairs)

# 1us0: Phe A 276 A and B are copies straddling chi2 90 deg; Glu A 120 A and B
# differ (OE1 A is nearer OE2 B)
pdb_1us0_phe_276 = """\
ATOM   2560  CA APHE A 276      20.326   6.924  41.525  0.62  8.37           C
ATOM   2563  CB APHE A 276      18.811   6.689  41.733  0.62  9.44           C
ATOM   2564  CG APHE A 276      17.952   6.093  40.580  0.62  8.00           C
ATOM   2565  CD1APHE A 276      17.228   6.769  39.591  0.62  9.11           C
ATOM   2566  CD2APHE A 276      17.853   4.736  40.499  0.62  7.92           C
ATOM   2567  CE1APHE A 276      16.443   6.107  38.652  0.62  9.20           C
ATOM   2568  CE2APHE A 276      17.071   4.049  39.586  0.62  9.24           C
ATOM   2569  CZ APHE A 276      16.383   4.723  38.648  0.62  9.98           C
ATOM   2571  CA BPHE A 276      20.134   6.965  41.065  0.38  6.42           C
ATOM   2574  CB BPHE A 276      18.733   6.722  41.699  0.38  8.87           C
ATOM   2575  CG BPHE A 276      17.888   6.137  40.558  0.38  8.63           C
ATOM   2576  CD1BPHE A 276      17.238   7.030  39.714  0.38 11.38           C
ATOM   2577  CD2BPHE A 276      17.874   4.778  40.439  0.38  8.44           C
ATOM   2578  CE1BPHE A 276      16.511   6.453  38.679  0.38 12.54           C
ATOM   2579  CE2BPHE A 276      17.126   4.222  39.418  0.38  9.68           C
ATOM   2580  CZ BPHE A 276      16.511   5.083  38.577  0.38 12.17           C
"""
pdb_1us0_glu_120 = """\
ATOM   1068  CA  GLU A 120       8.284 -15.895   9.413  1.00  4.56           C
ATOM   1071  CB AGLU A 120       9.062 -16.155   8.127  0.61  5.80           C
ATOM   1072  CG AGLU A 120       9.155 -17.675   7.889  0.61  9.15           C
ATOM   1073  CD AGLU A 120       9.937 -17.891   6.590  0.61 14.69           C
ATOM   1074  OE1AGLU A 120      10.365 -18.995   6.339  0.61 23.54           O
ATOM   1075  OE2AGLU A 120      10.261 -17.020   5.784  0.61 27.15           O
ATOM   1076  CB BGLU A 120       9.077 -16.132   8.114  0.39  6.45           C
ATOM   1077  CG BGLU A 120       9.036 -17.601   7.659  0.39  6.23           C
ATOM   1078  CD BGLU A 120       9.768 -18.525   8.594  0.39  7.11           C
ATOM   1079  OE1BGLU A 120      10.478 -18.062   9.485  0.39 10.53           O
ATOM   1080  OE2BGLU A 120       9.603 -19.769   8.429  0.39  8.41           O
"""

def exercise_corresponding_conformers():
  def nearest(hierarchy, name1, name2):
    # distance of name1 in A to name1 and to name2 in B
    ags = hierarchy.only_residue_group().atom_groups()
    a = dict((x.name.strip(), x) for x in ags[0].atoms() if ags[0].altloc == "A")
    b = dict((x.name.strip(), x) for x in ags[1].atoms())
    return a[name1].distance(b[name1]), a[name1].distance(b[name2])
  def run_flip(pdb_str, occ=None):
    hierarchy = iotbx.pdb.input(source_info=None,
      lines=pdb_str).construct_hierarchy()
    if occ is not None:
      for ag in hierarchy.only_residue_group().atom_groups():
        for atom in ag.atoms():
          if ag.altloc in occ: atom.occ = occ[ag.altloc]
    pairs = []
    info = hierarchy.flip_symmetric_amino_acids(swapped_pairs=pairs)
    # a second call renames nothing
    again = []
    hierarchy.flip_symmetric_amino_acids(swapped_pairs=again)
    assert again == [], again
    return hierarchy, pairs, info
  # A wants a flip, B not: 0.62, 0.38 weigh for flipping both
  hierarchy, pairs, info = run_flip(pdb_1us0_phe_276)
  assert len(pairs) == 4, len(pairs) # CD and CE in A and B
  assert ('"A PHE  276": altlocs A, B correspond, decided together - flipped'
    ) in info, info
  same, partner = nearest(hierarchy, "CD1", "CD2")
  assert same < 0.1 and partner > 2, (same, partner)
  # other weights: neither flips, the names stay matched too
  hierarchy, pairs, info = run_flip(pdb_1us0_phe_276, occ={"A": 0.38, "B": 0.62})
  assert pairs == []
  assert ('"A PHE  276": altlocs A, B correspond, decided together - not '
    'flipped') in info, info
  same, partner = nearest(hierarchy, "CD1", "CD2")
  assert same < 0.5 and partner > 2, (same, partner)
  # a tie keeps the input names
  hierarchy, pairs, info = run_flip(pdb_1us0_phe_276, occ={"A": 0., "B": 0.})
  assert pairs == []
  assert "decided together - not flipped" in info, info
  # conformers that differ keep their own decisions: only A flips
  hierarchy = iotbx.pdb.input(source_info=None,
    lines=pdb_1us0_glu_120).construct_hierarchy()
  sites = dict(((ag.altloc, a.name.strip()), a.xyz)
    for ag in hierarchy.only_residue_group().atom_groups() for a in ag.atoms())
  info = hierarchy.flip_symmetric_amino_acids()
  assert "together" not in info, info
  after = dict(((ag.altloc, a.name.strip()), a.xyz)
    for ag in hierarchy.only_residue_group().atom_groups() for a in ag.atoms())
  assert_only_swapped(sites, after, [("A", "OE1", "OE2")])

def three_conformers():
  # 1us0 Phe A 276 A and B, with C from A: CD and CE names exchanged, then CD1,
  # CE1 and CA turned about CB-CG. A swap in A and B alone used to make C join
  # them, and the next call chose differently.
  from scitbx import matrix
  rows = []
  occ = {"A": 0.12, "B": 0.08, "C": 0.80}
  for line in pdb_1us0_phe_276.splitlines():
    rows.append([line[12:16].strip(), line[16],
      matrix.col((float(line[30:38]), float(line[38:46]), float(line[46:54])))])
  sites = dict((r[0], r[2]) for r in rows if r[1] == "A")
  origin = sites["CB"]
  axis = (sites["CG"] - origin).normalize()
  def turned(site, angle):
    r = axis.axis_and_angle_as_r3_rotation_matrix(angle, deg=True)
    return r * (site - origin) + origin
  exchange = {"CD1": "CD2", "CD2": "CD1", "CE1": "CE2", "CE2": "CE1"}
  angles = {"CD1": -2., "CE1": -2., "CA": -3.3}
  for name, altloc, site in list(rows):
    if altloc != "A": continue
    name = exchange.get(name, name)
    if name in angles: site = turned(site, angles[name])
    rows.append([name, "C", site])
  lines = []
  for name, altloc, site in rows:
    lines.append("ATOM  %5d  %-3s%1sPHE A 276    %8.3f%8.3f%8.3f%6.2f%6.2f"
      "           %1s" % (len(lines) + 1, name, altloc, site[0], site[1],
      site[2], occ[altloc], 10., name[0]))
  return "\n".join(lines)

def exercise_three_conformers():
  hierarchy = iotbx.pdb.input(source_info=None,
    lines=three_conformers()).construct_hierarchy()
  # after a PDB round trip
  hierarchy = iotbx.pdb.input(source_info=None,
    lines=hierarchy.as_pdb_string()).construct_hierarchy()
  first = []
  info = hierarchy.flip_symmetric_amino_acids(swapped_pairs=first)
  second = []
  hierarchy.flip_symmetric_amino_acids(swapped_pairs=second)
  assert second == [], (len(first), len(second))
  # all three match: each CD1 nearest the CD1 of the others
  ags = hierarchy.only_residue_group().atom_groups()
  for i, ag1 in enumerate(ags):
    for ag2 in ags[i+1:]:
      cd1, cd2 = ag1.get_atom("CD1"), ag2.get_atom("CD1")
      assert cd1.distance(cd2) < cd1.distance(ag2.get_atom("CD2")), \
        (ag1.altloc, ag2.altloc)
  assert "altlocs A, B, C correspond, decided together" in info, info

def exercise_crossed_input():
  # A's names crossed against B: each conformer's own test says no; the names
  # still end matched, the same two namings whatever the input names
  lines = []
  exchange = {" CD1": " CD2", " CD2": " CD1", " CE1": " CE2", " CE2": " CE1"}
  for line in pdb_1us0_phe_276.splitlines():
    if line[16] == "A" and line[12:16] in exchange:
      line = line[:12] + exchange[line[12:16]] + line[16:]
    lines.append(line)
  hierarchy = iotbx.pdb.input(source_info=None,
    lines=lines).construct_hierarchy()
  pairs = []
  info = hierarchy.flip_symmetric_amino_acids(swapped_pairs=pairs)
  assert "altlocs A, B correspond, decided together" in info, info
  ags = hierarchy.only_residue_group().atom_groups()
  cd1, cd2 = ags[0].get_atom("CD1"), ags[1].get_atom("CD1")
  assert cd1.distance(cd2) < cd1.distance(ags[1].get_atom("CD2"))
  second = []
  hierarchy.flip_symmetric_amino_acids(swapped_pairs=second)
  assert second == []

def exercise_contradiction_gives_no_set():
  # OD1/OD2 of A, B, C turned 60 deg apart about CB-CG: A-B and B-C by the same
  # names, A-C crossed. CA turned with them, so that |chi2| is near 90 in all
  # three: B alone wants a flip. Each pair is decided together, the three are
  # not.
  from scitbx import matrix
  from scitbx.math import dihedral_angle
  rows = rows_of("ASP")
  sites = dict((r[0], matrix.col(r[2])) for r in rows)
  chi = dihedral_angle(sites=[sites[n] for n in ("CA", "CB", "CG", "OD1")],
    deg=True)
  axis = (sites["CG"] - sites["CB"]).normalize()
  def turned(name, origin, angle):
    r = axis.axis_and_angle_as_r3_rotation_matrix(angle, deg=True)
    return (r * (sites[name] - sites[origin]) + sites[origin]).elems
  result = [r for r in rows if r[0] not in ("CA", "OD1", "OD2")]
  for altloc, angle, target in (("A", 80, 85), ("B", 140, 95), ("C", 200, 85)):
    result.append(["CA", altloc, turned("CA", "CG", angle + chi - target)])
    for name in ("OD1", "OD2"):
      result.append([name, altloc, turned(name, "CG", angle)])
  def of(altlocs):
    return [r for r in result if r[1] == "" or r[1] in altlocs]
  own = []
  for altloc in "ABC":
    pairs = []
    make("ASP", of(altloc)).flip_symmetric_amino_acids(swapped_pairs=pairs)
    if pairs: own.append(altloc)
  assert own == ["B"], own
  for altlocs in ("AB", "BC", "AC"):
    info = make("ASP", of(altlocs)).flip_symmetric_amino_acids()
    assert "altlocs %s, %s correspond" % tuple(altlocs) in info, info
  hierarchy = make("ASP", result)
  pairs = []
  info = hierarchy.flip_symmetric_amino_acids(swapped_pairs=pairs)
  assert "correspond" not in info, info
  assert sorted(a.parent().altloc for a, b in pairs) == own, own
  second = []
  hierarchy.flip_symmetric_amino_acids(swapped_pairs=second)
  assert second == []

# Different rotamers that happen to be crossed: each conformer's naming is clear
pdb_1us0_val_67 = """\
ATOM    553  CA AVAL A  67       2.821  -3.084  42.712  0.71  5.76           C
ATOM    556  CB AVAL A  67       3.744  -2.990  41.455  0.71  7.80           C
ATOM    557  CG1AVAL A  67       5.136  -2.638  41.824  0.71  9.87           C
ATOM    558  CG2AVAL A  67       3.708  -4.288  40.640  0.71  7.30           C
ATOM    560  CA BVAL A  67       3.013  -3.012  42.739  0.29  5.33           C
ATOM    563  CB BVAL A  67       4.135  -3.043  41.702  0.29  5.18           C
ATOM    564  CG1BVAL A  67       3.818  -4.057  40.589  0.29  6.80           C
ATOM    565  CG2BVAL A  67       4.231  -1.691  41.007  0.29  5.86           C
"""
pdb_3cp5_leu_13 = """\
ATOM     34  CB  LEU A  13      12.971   0.595  34.455  1.00 13.45           C
ATOM     35  CG ALEU A  13      13.660   1.308  35.615  0.79 15.19           C
ATOM     36  CD1ALEU A  13      14.058   0.299  36.710  0.79 18.23           C
ATOM     37  CD2ALEU A  13      14.935   1.905  35.031  0.79 26.92           C
ATOM     38  CG BLEU A  13      14.278  -0.047  34.007  0.21 10.72           C
ATOM     39  CD1BLEU A  13      15.322   0.991  33.614  0.21  7.47           C
ATOM     40  CD2BLEU A  13      14.810  -0.985  35.077  0.21  8.27           C
"""
pdb_2erl_glu_21 = """\
ATOM    270  CB AGLU    21       4.980  -4.649  -1.566  0.32 13.49           C
ATOM    271  CG AGLU    21       3.942  -5.699  -1.498  0.32 13.60           C
ATOM    272  CD AGLU    21       4.064  -6.621  -2.739  0.32 15.27           C
ATOM    273  OE1AGLU    21       5.043  -6.411  -3.478  0.32 16.46           O
ATOM    274  OE2AGLU    21       3.193  -7.474  -2.870  0.32 19.03           O
ATOM    279  CB BGLU    21       4.769  -4.484  -1.766  0.68 14.37           C
ATOM    280  CG BGLU    21       5.737  -5.573  -1.952  0.68 19.43           C
ATOM    281  CD BGLU    21       5.251  -6.534  -3.036  0.68 19.91           C
ATOM    282  OE1BGLU    21       4.044  -6.877  -3.026  0.68 22.58           O
ATOM    283  OE2BGLU    21       6.037  -6.819  -3.934  0.68 23.06           O
"""
pdb_7kcu_arg_10 = """\
ATOM    159  CD AARG A  10     -21.462   8.478 -11.276  0.28 31.85           C
ATOM    160  NE AARG A  10     -20.646   8.199 -10.103  0.28 23.92           N
ATOM    161  CZ AARG A  10     -20.219   6.989  -9.765  0.28 38.75           C
ATOM    162  NH1AARG A  10     -19.483   6.826  -8.679  0.28 39.86           N
ATOM    163  NH2AARG A  10     -20.525   5.941 -10.515  0.28 25.97           N
ATOM    183  CD BARG A  10     -21.261   8.504 -11.282  0.27 31.19           C
ATOM    184  NE BARG A  10     -20.649   7.542 -10.394  0.27 24.02           N
ATOM    185  CZ BARG A  10     -20.197   7.811  -9.184  0.27 32.39           C
ATOM    186  NH1BARG A  10     -19.638   6.849  -8.473  0.27 40.16           N
ATOM    187  NH2BARG A  10     -20.311   9.024  -8.681  0.27 23.10           N
ATOM    207  CD CARG A  10     -21.195   7.771 -11.486  0.45 28.11           C
ATOM    208  NE CARG A  10     -20.890   8.477 -10.269  0.45 25.80           N
ATOM    209  CZ CARG A  10     -20.224   7.966  -9.257  0.45 31.91           C
ATOM    210  NH1CARG A  10     -19.992   8.709  -8.200  0.45 19.86           N
ATOM    211  NH2CARG A  10     -19.792   6.717  -9.304  0.45 34.19           N
"""

def renamed_on_own_decisions(hierarchy):
  pairs = []
  info = hierarchy.flip_symmetric_amino_acids(swapped_pairs=pairs)
  assert "correspond" not in info, info
  again = []
  hierarchy.flip_symmetric_amino_acids(swapped_pairs=again)
  assert again == [], again
  return sorted((a.parent().altloc, a.name.strip()) for a, b in pairs)

def exercise_clear_conformers_keep_own_decisions():
  for pdb_str, expected in (
      (pdb_1us0_val_67, []),
      (pdb_3cp5_leu_13, []),
      (pdb_2erl_glu_21, []),
      (pdb_7kcu_arg_10, [("A", "NH1"), ("B", "NH1"), ("C", "NH1")])):
    hierarchy = iotbx.pdb.input(source_info=None,
      lines=pdb_str).construct_hierarchy()
    result = renamed_on_own_decisions(hierarchy)
    assert result == expected, (pdb_str[17:26], result)

def exercise_chiral_namings_close_keep_own_decisions():
  # CA moved to volume -0.9 in both conformers of 1us0 Val A 67: each alone
  # keeps its names (delta 1.6, exchanged 3.4), the two namings differ by less
  # than 2.
  from scitbx import matrix
  hierarchy = iotbx.pdb.input(source_info=None,
    lines=pdb_1us0_val_67).construct_hierarchy()
  for ag in hierarchy.only_residue_group().atom_groups():
    def site(name):
      return matrix.col(ag.get_atom(name).xyz)
    normal = (site("CG1") - site("CB")).cross(site("CG2") - site("CB"))
    volume = (site("CA") - site("CB")).dot(normal)
    ag.get_atom("CA").xyz = (
      site("CA") + normal * ((-0.9 - volume) / normal.dot(normal))).elems
  assert renamed_on_own_decisions(hierarchy) == []

def exercise_namings_beyond_limit_keep_own_decisions():
  # 1us0 Phe A 276 with CA turned about CB-CG to chi2 60 in A, 120 in B: the two
  # namings differ by 59 and 58 deg, so each conformer decides alone.
  from scitbx import matrix
  from scitbx.math import dihedral_angle
  hierarchy = iotbx.pdb.input(source_info=None,
    lines=pdb_1us0_phe_276).construct_hierarchy()
  atom_groups = hierarchy.only_residue_group().atom_groups()
  for ag, target in zip(atom_groups, (60, 120)):
    def chi():
      return dihedral_angle(sites=[ag.get_atom(name).xyz
        for name in ("CA", "CB", "CG", "CD1")], deg=True)
    ca, cb, cg = [matrix.col(ag.get_atom(name).xyz)
      for name in ("CA", "CB", "CG")]
    r = (cg - cb).normalize().axis_and_angle_as_r3_rotation_matrix(
      chi() - target, deg=True)
    ag.get_atom("CA").xyz = (r * (ca - cg) + cg).elems
    assert abs(chi() - target) < 0.01, chi()
  assert renamed_on_own_decisions(hierarchy) == [("B", "CD1"), ("B", "CE1")]

# phenix_regression pdb/2izq_a_9_13.pdb: three residues on one position
pdb_2izq_a_11 = """\
ATOM    221  CA ATRP A  11      21.094  12.032  35.602  0.50  5.24           C
ATOM    224  CB ATRP A  11      20.690  10.588  35.288  0.50  6.15           C
ATOM    225  CG ATRP A  11      19.252  10.269  35.140  0.50  5.91           C
ATOM    226  CD1ATRP A  11      18.524  10.178  33.986  0.50  7.01           C
ATOM    227  CD2ATRP A  11      18.371   9.973  36.236  0.50  5.97           C
ATOM    229  CE2ATRP A  11      17.132   9.708  35.665  0.50  7.37           C
ATOM    245  CA CPHE A  11      20.950  12.135  35.430  0.15  5.92           C
ATOM    248  CB CPHE A  11      20.768  10.667  34.994  0.15  6.01           C
ATOM    249  CG CPHE A  11      19.330  10.235  34.845  0.15  7.05           C
ATOM    250  CD1CPHE A  11      18.847   9.877  33.587  0.15  8.78           C
ATOM    251  CD2CPHE A  11      18.533  10.174  35.995  0.15  7.70           C
ATOM    252  CE1CPHE A  11      17.551   9.436  33.473  0.15 10.43           C
ATOM    253  CE2CPHE A  11      17.230   9.752  35.854  0.15  9.27           C
ATOM    254  CZ CPHE A  11      16.789   9.396  34.594  0.15 10.98           C
ATOM    256  CA BTYR A  11      21.106  11.838  35.524  0.35  5.51           C
ATOM    259  CB BTYR A  11      20.694  10.354  35.327  0.35  5.65           C
ATOM    260  CG BTYR A  11      19.188  10.175  35.507  0.35  7.68           C
ATOM    261  CD1BTYR A  11      18.548  10.134  34.268  0.35  9.45           C
ATOM    263  CD2BTYR A  11      18.463  10.012  36.681  0.35  9.08           C
ATOM    265  CE1BTYR A  11      17.195   9.960  34.223  0.35 10.76           C
ATOM    267  CE2BTYR A  11      17.100   9.826  36.693  0.35 11.29           C
ATOM    269  CZ BTYR A  11      16.546   9.812  35.432  0.35 11.90           C
ATOM    271  OH BTYR A  11      15.178   9.650  35.313  0.35 19.29           O
"""

def exercise_microheterogeneity():
  # the first resname (TRP) has no flip; PHE and TYR each decide for themselves
  hierarchy = iotbx.pdb.input(source_info=None,
    lines=pdb_2izq_a_11).construct_hierarchy()
  pairs = []
  info = hierarchy.flip_symmetric_amino_acids(swapped_pairs=pairs)
  assert "Time to flip 2 residue(s)" in info, info
  renamed = sorted((atom1.parent().resname, atom1.parent().altloc, atom1.name)
    for atom1, atom2 in pairs)
  assert renamed == [("PHE", "C", " CD1"), ("PHE", "C", " CE1"),
                     ("TYR", "B", " CD1"), ("TYR", "B", " CE1")], renamed

def exercise_pair_spans_blank_and_altloc():
  # 8EYP Tyr A 102: the two atoms of a pair cannot be exchanged
  rows = split(rows_of("TYR"), ["CD2"])
  before, after, info = flip(make("TYR", rows))
  assert_only_swapped(before, after, [])
  assert '"A TYR   11": pair spans blank and altloc - not flipped' in info, info

def exercise_incomplete_h_in_altloc():
  # 7KCU Arg A 69: blank NH, HH22 missing in A
  names = ["HH11", "HH12", "HH21", "HH22"]
  rows = [r for r in split(rows_of("ARG", keep=names), names)
    if r[:2] != ["HH22", "A"]]
  before, after, info = flip(make("ARG", rows))
  assert_only_swapped(before, after, [])
  assert '"A ARG   11": not complete - not flipped' in info, info
  # the H in altlocs only, nothing blank in a pair: only B swaps
  rows = [r for r in split(rows_of("ARG", keep=names),
    names + ["NH1", "NH2"]) if r[:2] != ["HH22", "A"]]
  before, after, info = flip(make("ARG", rows))
  assert_only_swapped(before, after, [("B", "NH1", "NH2"),
    ("B", "HH11", "HH21"), ("B", "HH12", "HH22")])
  assert '"A ARG   11": not complete in altloc A - not flipped' in info, info

def exercise_d_ring():
  # Tyr shared Phe's pair list: the D pairs were listed twice and swapped back
  for code in ("PHE", "TYR"):
    names = ["HD1", "HD2", "HE1", "HE2"]
    rows = rows_of(code, keep=names)
    for row in rows:
      if row[0] in names: row[0] = "D" + row[0][1:]
    before, after, info = flip(make(code, rows))
    assert_only_swapped(before, after, [("", "CD1", "CD2"), ("", "CE1", "CE2"),
      ("", "DD1", "DD2"), ("", "DE1", "DE2")])
    rows = rows_of(code, keep=names)
    before, after, info = flip(make(code, rows))
    assert_only_swapped(before, after, [("", "CD1", "CD2"), ("", "CE1", "CE2"),
      ("", "HD1", "HD2"), ("", "HE1", "HE2")])

def second_group(hierarchy, altloc, names):
  """Move the atoms in names to a second atom group of that altloc."""
  rg = hierarchy.only_residue_group()
  first = [ag for ag in rg.atom_groups() if ag.altloc == altloc][0]
  group = iotbx.pdb.hierarchy.atom_group(altloc=altloc, resname=first.resname)
  for name in names:
    atom = first.get_atom(name)
    first.remove_atom(atom)
    group.append_atom(atom)
  rg.append_atom_group(group)
  assert len([c for c in rg.conformers() if c.altloc == altloc]) == 1
  return hierarchy

def glu_with_he2_in_second_group():
  rows = rows_of("GLU")
  add_h(rows, "HE2", "", "OE2")
  for row in rows: row[1] = "A"
  return second_group(make("GLU", rows), "A", ["HE2"])

def distance(atoms, name1, name2):
  sites = dict([(atom.name.strip(), atom) for atom in atoms])
  return sites[name1].distance(sites[name2])

def exercise_conformer_in_several_atom_groups():
  # The reader merges an altloc; code that builds a hierarchy may not.
  hierarchy = glu_with_he2_in_second_group()
  before, after, info = flip(hierarchy)
  assert_only_swapped(before, after, [])
  assert ('"A GLU   11": protonated acid (one H on the O pair) in altloc A - '
    'not flipped') in info, info
  assert abs(distance(hierarchy.atoms(), "OE2", "HE2") - 1) < 1.e-3

def exercise_conformer_in_several_atom_groups_no_acid_h():
  # flipped completely: a pair across the groups, and a pair in the second group
  for code, names, swapped in (
      ("GLU", ["OE2"], [("A", "OE1", "OE2")]),
      ("PHE", ["CE1", "CE2"], [("A", "CD1", "CD2"), ("A", "CE1", "CE2")])):
    rows = rows_of(code)
    for row in rows: row[1] = "A"
    before, after, info = flip(second_group(make(code, rows), "A", names))
    assert_only_swapped(before, after, swapped)
    assert "not flipped" not in info, info
    assert "Time to flip 1 residue(s)" in info

def exercise_conformer_in_several_atom_groups_process():
  import mmtbx.model
  from libtbx.utils import null_out
  model = mmtbx.model.manager(model_input=None,
    pdb_hierarchy=glu_with_he2_in_second_group(), log=null_out())
  params = model.get_default_pdb_interpretation_params()
  params.pdb_interpretation.sort_atoms = False
  params.pdb_interpretation.flip_symmetric_amino_acids = True
  model.process(pdb_interpretation_params=params)
  atoms = model.get_hierarchy().atoms()
  assert abs(distance(atoms, "OE2", "HE2") - 1) < 1.e-3

def exercise_model_id_in_info():
  rows = rows_of("GLU")
  add_h(rows, "HE2", "", "OE2")
  text = make("GLU", rows).as_pdb_string()
  hierarchy = iotbx.pdb.input(source_info=None, lines=(
    "MODEL        1\n%sENDMDL\nMODEL        2\n%sENDMDL\n" % (text, text)
    ).splitlines()).construct_hierarchy()
  info = hierarchy.flip_symmetric_amino_acids()
  for model_id in ("1", "2"):
    assert ('"MODEL %s A GLU   11": protonated acid (one H on the O pair) - not '
      'flipped' % model_id) in info, info

def exercise_selected_hierarchy_and_residue_identity():
  # select() shares atoms: their parent can be gone. The report has the icode.
  rows = split(rows_of("GLU"), ["CB"])
  turn(rows, "A", "CB", "CG", "CD")
  hierarchy = make("GLU", rows)
  hierarchy.only_residue_group().icode = "X"
  selected = hierarchy.select(hierarchy.atom_selection_cache().selection("all"))
  del hierarchy
  info = selected.flip_symmetric_amino_acids()
  assert '"A GLU   11X": conformers disagree - not flipped' in info, info

def exercise_site_values_follow_the_site():
  # all values of a site go with it; what names the atom stays
  site = ["xyz", "sigxyz", "occ", "sigocc", "b", "sigb", "uij", "fp", "fdp",
          "resolution"]
  name = ["name", "serial", "element", "charge", "segid", "hetero", "i_seq",
          "tmp"]
  if iotbx.pdb.hierarchy.atom.has_siguij(): site.append("siguij")
  def distinct(i):
    return {"sigxyz": (0.01 * i, 0.02, 0.03), "occ": 0.7 - 0.4 * i,
      "sigocc": 0.01 + 0.01 * i, "b": 20. + i, "sigb": 1. + i,
      "uij": (0.1 + i, 0.2, 0.3, 0.01, 0.02, 0.03),
      "fp": -1. - i, "fdp": 2. + i, "resolution": 1.5 + i,
      "siguij": (0.001 * (i + 1), 0.002, 0.003, 0.004, 0.005, 0.006),
      "charge": "%d-" % (i + 1), "segid": "S%d  " % i, "tmp": 7 + i,
      "hetero": bool(i)}
  hierarchy = make("GLU", split(rows_of("GLU"), ["OE1", "OE2"]))
  hierarchy.atoms().reset_i_seq()
  for ag in hierarchy.only_residue_group().atom_groups():
    if ag.altloc == "": continue
    for i, atom_name in enumerate(("OE1", "OE2")):
      atom = ag.get_atom(atom_name)
      for attr, value in distinct(i).items():
        if attr == "siguij" and "siguij" not in site: continue
        setattr(atom, attr, value)
  def values(ag, atom_name, attrs):
    atom = ag.get_atom(atom_name)
    return [getattr(atom, attr) for attr in attrs]
  before = {}
  for ag in hierarchy.only_residue_group().atom_groups():
    if ag.altloc == "": continue
    for atom_name in ("OE1", "OE2"):
      before[ag.altloc, atom_name] = (values(ag, atom_name, site),
                                      values(ag, atom_name, name))
  assert hierarchy.flip_symmetric_amino_acids() != ""
  for ag in hierarchy.only_residue_group().atom_groups():
    if ag.altloc == "": continue
    for atom_name, partner in (("OE1", "OE2"), ("OE2", "OE1")):
      assert values(ag, atom_name, site) == before[ag.altloc, partner][0]
      assert values(ag, atom_name, name) == before[ag.altloc, atom_name][1]

def run_in_process():
  exercise_acid_h()
  exercise_blank_o_acid_h_in_altlocs()
  exercise_a_flips_b_not()
  exercise_disagree_on_blank_pair()
  exercise_test_spans_blank_and_altloc()
  exercise_no_dihedral()
  exercise_second_call_renames_nothing()
  exercise_corresponding_conformers()
  exercise_three_conformers()
  exercise_crossed_input()
  exercise_contradiction_gives_no_set()
  exercise_clear_conformers_keep_own_decisions()
  exercise_chiral_namings_close_keep_own_decisions()
  exercise_namings_beyond_limit_keep_own_decisions()
  exercise_microheterogeneity()
  exercise_pair_spans_blank_and_altloc()
  exercise_incomplete_h_in_altloc()
  exercise_d_ring()
  exercise_conformer_in_several_atom_groups()
  exercise_conformer_in_several_atom_groups_no_acid_h()
  exercise_model_id_in_info()
  exercise_selected_hierarchy_and_residue_identity()
  exercise_site_values_follow_the_site()

if __name__=="__main__":
  import libtbx.load_env
  run_in_process()
  if libtbx.env.find_in_repositories('chem_data') is not None:
    exercise_conformer_in_several_atom_groups_process()
    args = sys.argv[1:]
    del sys.argv[1:]
    run(*tuple(args))
  else:
    print('chem_data is missing, skipping test')
