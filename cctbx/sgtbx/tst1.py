# $Id$

import sgtbx
print sgtbx.STBF
print sgtbx.CRBF
print sgtbx.CTBF
s = sgtbx.parse_string("x+1/2,y,z")
print s.string()
print s.where()
m = sgtbx.RTMx(s)
print m.as_xyz(0, 0, "xyz", ",")
print m.as_xyz(0, 1, "xyz", ",")
print m.as_tuple()
print m.as_tuple(12, 72)
print m.modPositive()
mm = m * m
print mm.as_xyz()
mi = m.inverse()
print mi.as_xyz()
print mi.as_xyz(0, 1, "xyz", ",")
print (m * mi).as_xyz()
m = sgtbx.RTMx("2*x,2*y,2*z", "", 72, 12)
mi = m.inverse()
print mi.isValid()
print mi.as_xyz()
m = sgtbx.RTMx("-x,-y,z")
print (m * m).as_xyz()
m = sgtbx.RTMx("y+z-1/3,x+z-1/3,x+y", "", 2, 3)
print m
print m.inverse()
assert (m * m.inverse()).isUnit()
assert (m.multiply(m.inverse())).isUnit()
mi = m.inverse_with_cancel()
print mi
mmi = m.multiply(mi)
print mmi
assert mmi.isUnit()

SgOps = sgtbx.SpaceGroup("P 61")
for i in xrange(SgOps.OrderZ()):
  print SgOps(i).as_xyz()
  print SgOps(i).cancel().TBF()
SgOps.makeTidy()
for i in xrange(SgOps.OrderZ()):
  print SgOps(i).as_xyz()

SgOps1 = sgtbx.SpaceGroup("P 3")
SgOps2 = sgtbx.SpaceGroup("P 4")
print SgOps1 == SgOps2, "should be 0"
print SgOps1 != SgOps2, "should be 1"
SgOps2 = sgtbx.SpaceGroup("P 3")
SgOps2.makeTidy()
print SgOps1 == SgOps2, "should be 1"
print SgOps1 != SgOps2, "should be 0"

SgOps = sgtbx.SpaceGroup("F 4 2 3")
Z2POp = SgOps.getZ2POp()
print Z2POp.M().as_xyz()
print Z2POp.InvM().as_xyz()
P2ZOp = Z2POp.swap()
print P2ZOp.M().as_xyz()
print P2ZOp.InvM().as_xyz()
SgOpsP = SgOps.ChangeBasis(Z2POp)
SgOpsZ = SgOpsP.ChangeBasis(P2ZOp)
assert SgOps == SgOpsZ
x = (0.123, 0.456, 0.789)
px = Z2POp(x)
zx = P2ZOp(px)
print "x (%.6g, %.6g, %.6g)" % x
print "px (%.6g, %.6g, %.6g)" % px
print "zx (%.6g, %.6g, %.6g)" % zx
assert ("(%.6g, %.6g, %.6g)" % x) == ("(%.6g, %.6g, %.6g)" % zx)

SgOps = sgtbx.SpaceGroup("P 2")
SgOps.ParseHallSymbol("P 2x")
print SgOps.OrderZ()
assert SgOps == sgtbx.SpaceGroup("P 2 2")

SgOps = sgtbx.SpaceGroup("P 4c -2ab -1b -1ac -1a -1bc")
Z2POp = SgOps.getZ2POp()
print Z2POp.M().as_xyz(), Z2POp.InvM().as_xyz()

import pickle
pstr = pickle.dumps(SgOps)
print pstr
up = pickle.loads(pstr)
print sgtbx.SpaceGroupInfo(SgOps).BuildHallSymbol(1)
print sgtbx.SpaceGroupInfo(up).BuildHallSymbol(1)

enantiomorphic_pairs = {}
cbop = sgtbx.ChOfBasisOp(sgtbx.RTMx("x+1/12,y+1/12,z+1/12"))
iter = sgtbx.SpaceGroupSymbolIterator()
#print iter.next().ExtendedHermann_Mauguin()
for s in iter:
  SgOps = sgtbx.SpaceGroup(s.Hall())#ChangeBasis(cbop)
  SgInfo = sgtbx.SpaceGroupInfo(SgOps)
  ch1 = SgInfo.getChangeOfHandOp()
  print s.ExtendedHermann_Mauguin(), ch1.M()
  try:
    b = sgtbx.Brick(SgInfo)
  except RuntimeError, e:
    print e
  print b
  for i in xrange(3):
    for j in xrange(2):
      p = b(i, j)
      print p.Point(), p.Off()
  e = SgOps.ChangeBasis(ch1)
  eInfo = sgtbx.SpaceGroupInfo(e)
  if (e != SgOps):
    assert SgInfo.isEnantiomorphic()
    enantiomorphic_pairs[s.SgNumber()] = eInfo.SgNumber()
  else:
    assert not SgInfo.isEnantiomorphic()
  ch2 = eInfo.getChangeOfHandOp()
  ch1ch2 = ch1.M().multiply(ch2.M()).modPositive()
  assert ch1ch2.isUnit() # this is sufficient, but not necessary.
  SECg = sgtbx.SymEquivCoordinates(SgOps, (0.123,0.456,0.789))
  flippedX = ch1(SECg(0))
  SECf = sgtbx.SymEquivCoordinates(e, flippedX)
  Fg = SECg.StructureFactor((3,5,7))
  Ff = SECf.StructureFactor((3,5,7))
  d = abs(Fg) - abs(Ff)
  assert d < 1.e-5
for p in enantiomorphic_pairs.keys():
  print sgtbx.SpaceGroupSymbols(p).Hermann_Mauguin(),
  print sgtbx.SpaceGroupSymbols(enantiomorphic_pairs[p]).Hermann_Mauguin()
assert len(enantiomorphic_pairs) == 22
