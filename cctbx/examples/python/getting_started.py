from cctbx import sgtbx
from cctbx import uctbx

UnitCell = uctbx.UnitCell((11, 12, 13, 90, 100, 90))
print UnitCell
Symbols = sgtbx.SpaceGroupSymbols("C 2")
SgOps = sgtbx.SpaceGroup(Symbols)
SgOps.CheckUnitCell(UnitCell)
for M in SgOps: print M
