from __future__ import division
import sys
from StringIO import StringIO
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
}
values = {
  "ASP" : [ 152,  -82, -999, -999],
  "GLU" : [ 162, -162,  -18,   18],
  "PHE" : [-105,   15, -999, -999],
  "TYR" : [ -96,    6, -999, -999],
  "ARG" : [-179,  179, -999, -999],
}

def largest_torsion(s, f, i):
  t = s.split(f)
  if len(t)==1: return -999
  s=t[1].splitlines()
  f="ideal   model   delta sinusoidal    sigma   weight residual"
  reading = False
  for line in s:
    print line
    if reading: break
    if line.find(f)>-1:
      reading=True
  d=line.split()
  print d
  return float(d[i])

def run():
  for code, lines in pdbs.items():
    print ('-%s- ' % code)*10
    f=file("tst_symmetric_flips_%s.pdb" % code, "wb").write(lines)
    cmd = "phenix.pdb_interpretation tst_symmetric_flips_%s.pdb" % code
    print cmd
    ero = easy_run.fully_buffered(command=cmd)
    out = StringIO()
    ero.show_stdout(out=out)
    i=0
    d = largest_torsion(out.getvalue(),
                        geo_strs[code],
                        i+1)
    print 'torsion',d
    assert abs(d-values[code][i])<1
    i+=1
    d = largest_torsion(out.getvalue(),
                        geo_strs[code],
                        i+1)
    print 'delta',d
    assert abs(d-values[code][i])<1

    pdb_inp = iotbx.pdb.input('tst_symmetric_flips_%s.pdb' % code)
    hierarchy = pdb_inp.construct_hierarchy()
    rc = hierarchy.flip_symmetric_amino_acids()
    print rc
    hierarchy.write_pdb_file("tst_symmetric_flips_%s.pdb" % code)
    cmd = "phenix.pdb_interpretation tst_symmetric_flips_%s.pdb" % code
    print cmd
    ero = easy_run.fully_buffered(command=cmd)
    out = StringIO()
    ero.show_stdout(out=out)

    i+=1
    d = largest_torsion(out.getvalue(),
                        geo_strs[code],
                        i-1)
    print 'torsion',d
    assert abs(d-values[code][i])<1
    i+=1
    d = largest_torsion(out.getvalue(),
                        geo_strs[code],
                        i-1)
    print 'delta',d
    assert abs(d-values[code][i])<1

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
