# $Id$

import sgtbx
import uctbx

UnitCell = uctbx.UnitCell((11, 12, 13, 90, 100, 90))
print UnitCell
Symbols = sgtbx.SpaceGroupSymbols("C 2")
SgOps = sgtbx.SgOps(Symbols.Hall())
SgOps.CheckUnitCell(UnitCell)
for i in xrange(SgOps.OrderZ()):
  print SgOps(i).as_xyz()
