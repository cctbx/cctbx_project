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

sgo = sgtbx.SgOps("P 3")
for i in xrange(sgo.OrderZ()):
  print sgo(i).as_xyz()
sgo.makeTidy()
for i in xrange(sgo.OrderZ()):
  print sgo(i).as_xyz()

sgo1 = sgtbx.SgOps("P 3")
sgo2 = sgtbx.SgOps("P 4")
print sgo1 == sgo2, "should be 0"
print sgo1 != sgo2, "should be 1"
sgo2 = sgtbx.SgOps("P 3")
sgo2.makeTidy()
print sgo1 == sgo2, "should be 1"
print sgo1 != sgo2, "should be 0"

sgo = sgtbx.SgOps("F 4 2 3")
Z2POp = sgo.getZ2POp()
print Z2POp.M().as_xyz()
print Z2POp.InvM().as_xyz()
P2ZOp = Z2POp.swap()
print P2ZOp.M().as_xyz()
print P2ZOp.InvM().as_xyz()
sgop = sgo.ChangeBasis(Z2POp)
sgoz = sgop.ChangeBasis(P2ZOp)
assert sgo == sgoz

sgo = sgtbx.SgOps("P 2")
sgo.ParseHallSymbol("P 2x")
print sgo.OrderZ()
assert sgo == sgtbx.SgOps("P 2 2")
