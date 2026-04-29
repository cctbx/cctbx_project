"""Test flipping symmetric residues"""
from __future__ import absolute_import, division, print_function
import sys
from six.moves import cStringIO as StringIO
from libtbx import easy_run
import iotbx.pdb

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

if __name__=="__main__":
  import libtbx.load_env
  if libtbx.env.find_in_repositories('chem_data') is not None:
    args = sys.argv[1:]
    del sys.argv[1:]
    run(*tuple(args))
  else:
    print('chem_data is missing, skipping test')
