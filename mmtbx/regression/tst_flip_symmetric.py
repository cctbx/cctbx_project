import os, sys
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
}
values = {
  "ASP" : [ 152,  -82, -999, -999],
  "GLU" : [ 162, -162,  -18,   18],
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
    print code
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
