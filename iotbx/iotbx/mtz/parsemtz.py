from iotbx import mtz

if __name__=="__main__":
 for c in xrange(1):
  p = mtz.Mtz("test.mtz")
  print "spacegroup:",p.SpaceGroup()
  print "Unitcell:", tuple(p.UnitCell(0))

  print "Columns:", tuple(p.columns())

  H = p.getColumn("H")
  K = p.getColumn("K")
  L = p.getColumn("L")
  FC = p.getColumn("I")
  PHIC = p.getColumn("SIGI")

  for i in xrange(p.size()):
    if H(i)==K(i)==L(i):
     print H(i),K(i),L(i),FC(i),PHIC(i)
