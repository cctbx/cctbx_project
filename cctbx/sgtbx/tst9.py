# $Id$

import sys
import random
from cctbx_boost import uctbx
from cctbx_boost import sgtbx

ShortCut = "--ShortCut" in sys.argv
StandardOnly = "--StandardOnly" in sys.argv
CheckLettersOnly = "--CheckLettersOnly" in sys.argv
RandomSeed0 = "--RandomSeed0" in sys.argv
Endless = "--Endless" in sys.argv

if (RandomSeed0):
  random.seed(0)

def get_unitcell(SgInfo):
  if (143 <= SgInfo.SgNumber() < 195):
    RefUnitCell = uctbx.UnitCell((10, 10, 10, 90, 90, 120))
  else:
    RefUnitCell = uctbx.UnitCell((10, 10, 10, 90, 90, 90))
  return RefUnitCell.ChangeBasis(SgInfo.CBOp().M().as_tuple()[0])

if (ShortCut):
  settings = ("P 61 2 2",)
else:
  from settings import * # see examples/python/make_settings.py


def OneCycle():
  nUndetermined = 0
  nMismatches = 0
  print "# Starting over"
  for LookupSymbol in settings:
    if (StandardOnly and LookupSymbol[:5] == "Hall:"): continue
    SgSymbols = sgtbx.SpaceGroupSymbols(LookupSymbol)
    HSym = SgSymbols.Hall()
    SgOps = sgtbx.SpaceGroup(HSym)
    SgInfo = SgOps.Info()
    print "SpaceGroup %s (%d) %s" % (
      SgInfo.BuildLookupSymbol(),
      SgInfo.SgNumber(),
      SgInfo.BuildHallSymbol())
    sys.stdout.flush()
    UnitCell = get_unitcell(SgInfo)
    SgOps.CheckUnitCell(UnitCell)
    SpecialPositionTolerances = \
    sgtbx.SpecialPositionTolerances(UnitCell, SgOps, 0.1, 0.1)
    SnapParameters = \
    sgtbx.SpecialPositionSnapParameters(UnitCell, SgOps, 1, 0.05)
    WTab = sgtbx.WyckoffTable(SgInfo)
    for WP in WTab:
      print "WyckoffPosition", WP.Letter(), WP.M(), WP.SpecialOp()
      while 1:
        RandomX = (random.uniform(-2,2),
                   random.uniform(-2,2),
                   random.uniform(-2,2))
        WX = WP.SpecialOp().multiply(RandomX)
        try:
          # make sure that we are not too close to a "more special" position
          # at the same time, probe for numerical instability
          SE = sgtbx.SymEquivCoordinates(SpecialPositionTolerances, WX)
        except RuntimeError, e:
          print e
          print "ERROR: WX =", WX
        else:
          if (SE.M() == WP.M()): break
      print "      WX =", WX
      print "  Ref WX =", SgInfo.CBOp()(WX)
      for SymOp in [random.choice(SgOps)]:
        UM = sgtbx.RTMx("x-2,y+3,z-5", "", 1, 1).multiply(SymOp)
        UMWX = UM.multiply(WX)
        SS = sgtbx.SiteSymmetry(SnapParameters, UMWX, 1)
        assert SS.DistanceMoved2() < 1.e-5
        try:
          WMap = WTab.getWyckoffMapping(SS)
        except RuntimeError, e:
          print e
          nUndetermined = nUndetermined + 1
        else:
          WPX = WMap.WP()
          if (WP.Letter() != WPX.Letter()):
            print WP.Letter(), WP.M(), WP.SpecialOp(), "input"
            print WPX.Letter(),WPX.M(),WPX.SpecialOp(),"getWyckoffMapping"
            print "ERROR: Wyckoff letter mismatch!"
            print "ERROR: Ref WX =", SgInfo.CBOp()(WX)
            nMismatches = nMismatches + 1
          WWMap = WTab.getWyckoffMapping(UnitCell, SgOps, UMWX, 0.04)
          assert WMap.WP().Letter() == WWMap.WP().Letter()
          assert UnitCell.Distance2(WWMap.snap(UMWX),SS.SnapPosition()) < 1.e-5
          if (not CheckLettersOnly):
            SES = sgtbx.SymEquivCoordinates(SS)
            print "  Mapping", WMap.Mapping()
            x = WMap.snap_to_representative(UMWX)
            d = SES.getShortestDistance2(UnitCell, x)
            assert d < 1.e-5
            x = WWMap.snap_to_representative(UMWX)
            d = SES.getShortestDistance2(UnitCell, x)
            assert d < 1.e-5
    print
  assert nUndetermined == 0
  assert nMismatches == 0

print "# CheckLettersOnly =", CheckLettersOnly
print "# RandomSeed0 =", RandomSeed0
while 1:
  OneCycle()
  if (not Endless): break
