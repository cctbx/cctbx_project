# $Id$

import sys
import random
import uctbx
import sgtbx

ShortCut = "--ShortCut" in sys.argv
StandardOnly = "--StandardOnly" in sys.argv
RandomSeed0 = "--RandomSeed0" in sys.argv
Endless = "--Endless" in sys.argv

if (RandomSeed0):
  random.seed(0)

def get_unitcell(SgType):
  if (143 <= SgType.SgNumber() < 195):
    RefUnitCell = uctbx.UnitCell((10, 10, 10, 90, 90, 120))
  else:
    RefUnitCell = uctbx.UnitCell((10, 10, 10, 90, 90, 90))
  return RefUnitCell.ChangeBasis(SgType.CBOp().M().as_tuple()[0])

if (ShortCut):
  settings = ("B 2 2 2",)
else:
  from settings import * # see examples/python/make_settings.py


def OneCycle():
  print "# Starting over"
  for LookupSymbol in settings:
    if (StandardOnly and LookupSymbol[:5] == "Hall:"): continue
    SgSymbols = sgtbx.SpaceGroupSymbols(LookupSymbol)
    HSym = SgSymbols.Hall()
    SgOps = sgtbx.SgOps(HSym)
    SgType = SgOps.getSpaceGroupType()
    print "SpaceGroup %s (%d) %s" % (
      SgOps.BuildLookupSymbol(SgType),
      SgType.SgNumber(),
      SgOps.BuildHallSymbol(SgType))
    sys.stdout.flush()
    UnitCell = get_unitcell(SgType)
    SgOps.CheckUnitCell(UnitCell)
    SnapParameters20 = \
    sgtbx.SpecialPositionSnapParameters(UnitCell, SgOps, 0, 2.0)
    SnapParameters01 = \
    sgtbx.SpecialPositionSnapParameters(UnitCell, SgOps, 0, 0.1)
    WTab = sgtbx.WyckoffTable(SgOps, SgType)
    for i in xrange(20):
      RandomX = (random.uniform(-2,2),
                 random.uniform(-2,2),
                 random.uniform(-2,2))
      print "RandomX ", RandomX
      #
      SP = sgtbx.SpecialPosition(SnapParameters20, RandomX, 1, 1)
      SWMap = WTab.getWyckoffMapping(SP)
      print SWMap.WP().M(), SWMap.WP().Letter(), SWMap.WP().SpecialOp(),
      print SP.getPointGroupType()
      WWMap = WTab.getWyckoffMapping(UnitCell, SgOps, SP.SnapPosition(), 0.1)
      assert SWMap.WP().Letter() == WWMap.WP().Letter()
      assert UnitCell.Distance2(SWMap.snap(RandomX), SP.SnapPosition()) < 1.e-5
      assert UnitCell.Distance2(WWMap.snap(RandomX), SP.SnapPosition()) < 1.e-5
      SES = sgtbx.SymEquivCoordinates(SP)
      x = SWMap.snap_to_representative(RandomX)
      d = SES.getShortestDistance2(UnitCell, x)
      assert d < 1.e-5
      x = WWMap.snap_to_representative(RandomX)
      d = SES.getShortestDistance2(UnitCell, x)
      assert d < 1.e-5
      #
      WWMap = WTab.getWyckoffMapping(UnitCell, SgOps, RandomX, 2.0)
      SP = sgtbx.SpecialPosition(SnapParameters01, WWMap.snap(RandomX), 0, 1)
      print WWMap.WP().M(), WWMap.WP().Letter(), WWMap.WP().SpecialOp(),
      print SP.getPointGroupType()
      assert SP.DistanceMoved2() < 1.e-5
      SWMap = WTab.getWyckoffMapping(SP)
      assert SWMap.WP().Letter() == WWMap.WP().Letter()
      assert UnitCell.Distance2(SWMap.snap(RandomX), SP.SnapPosition()) < 1.e-5
      assert UnitCell.Distance2(WWMap.snap(RandomX), SP.SnapPosition()) < 1.e-5
      SEW = sgtbx.SymEquivCoordinates(WWMap, RandomX)
      x = WWMap.snap_to_representative(RandomX)
      d = SEW.getShortestDistance2(UnitCell, x)
      assert d < 1.e-5
      x = SWMap.snap_to_representative(WWMap.snap(RandomX))
      d = SEW.getShortestDistance2(UnitCell, x)
      assert d < 1.e-5
    print

print "# RandomSeed0 =", RandomSeed0
while 1:
  OneCycle()
  if (not Endless): break
