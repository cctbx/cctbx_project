# $Id$

import sgtbx
import uctbx

UnitCell = uctbx.UnitCell((10, 10, 15, 90, 90, 120))
print UnitCell
Symbols = sgtbx.SpaceGroupSymbols("P 62 2 2")
SgOps = sgtbx.SgOps(Symbols.Hall())
SgOps.CheckUnitCell(UnitCell)
for i in xrange(SgOps.OrderZ()):
  print SgOps(i).as_xyz()
