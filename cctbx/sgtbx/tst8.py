# $Id$

from cctbx import uctbx
from cctbx import sgtbx

uc = uctbx.UnitCell([])
sg = sgtbx.SpaceGroup("F 4 2")
sginfo = sg.Info()
print sginfo.BuildLookupSymbol()
WTab = sgtbx.WyckoffTable(sginfo)
WTab.expand(sg)
for p in WTab:
  print p.M(), p.Letter(), p.SpecialOp()
  for o in p: print "   ", o
p = WTab("a")
print p.M(), p.Letter(), p.SpecialOp()
X = (-1.45, 2.45, 0.2)
MinMateDistance = 0.05
MinMateDistance = 9.21
MinMateDistance = 0.3
SnapParameters = sgtbx.SpecialPositionSnapParameters(uc, sg, 0, MinMateDistance)
SS = sgtbx.SiteSymmetry(SnapParameters, X, 1)
print SS.OriginalPosition()
print SS.SnapPosition()
print SS.DistanceMoved()
print SS.ShortestDistance()
print SS.isWellBehaved()
print SS.M()
print SS.SpecialOp()
print SS.PointGroupType()
for M in SS:
  print M
SpecialPositionTolerances = sgtbx.SpecialPositionTolerances(
                                      uc, sg, MinMateDistance, MinMateDistance)
SE = sgtbx.SymEquivCoordinates(SpecialPositionTolerances, X)
print SE.M()
for sx in SE:
  print sx
SE = sgtbx.SymEquivCoordinates(SS)
print SE.M()
for sx in SE:
  print sx
p = WTab.getWyckoffMapping(SS).WP()
print p.M(), p.Letter(), p.SpecialOp()
