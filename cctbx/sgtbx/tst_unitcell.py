# $Id$

from cctbx import uctbx
from cctbx import sgtbx

uc = uctbx.UnitCell((1, 1, 1, 90, 90, 120))
G = uc.getMetricalMatrix()
SgOps = sgtbx.SpaceGroup(sgtbx.parse_string("P 3"))
print "OK"
SgOps.CheckMetricalMatrix(G)
SgOps = sgtbx.SpaceGroup(sgtbx.parse_string("P 4 3*"))
try:
  SgOps.CheckMetricalMatrix(G)
except RuntimeError, e:
  print "Expected:", e.args[0]
else:
  raise SystemError
uc = uctbx.UnitCell((1, 1, 1, 90, 90, 90))
G = uc.getMetricalMatrix()
SgOps = sgtbx.SpaceGroup(sgtbx.parse_string("P 4 3*"))
print "OK"
