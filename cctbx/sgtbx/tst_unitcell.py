# $Id$

import uctbx
import sgtbx

uc = uctbx.UnitCell((1, 1, 1, 90, 90, 120))
G = uc.get_MetricalMatrix()
sgo = sgtbx.SgOps(sgtbx.parse_string("P 3"))
print "OK"
sgo.CheckMetricalMatrix(G)
sgo = sgtbx.SgOps(sgtbx.parse_string("P 4 3*"))
try:
  sgo.CheckMetricalMatrix(G)
except RuntimeError, e:
  print "Expected:", e.args[0]
else:
  raise SystemError
uc = uctbx.UnitCell((1, 1, 1, 90, 90, 90))
G = uc.get_MetricalMatrix()
sgo = sgtbx.SgOps(sgtbx.parse_string("P 4 3*"))
print "OK"
