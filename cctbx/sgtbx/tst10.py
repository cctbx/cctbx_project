import sys
import random
from cctbx_boost import sgtbx
from cctbx.development import debug_utils

def OneCycle(settings):
  print "Testing SiteSymmetry and WyckoffMapping 2"
  for LookupSymbol in settings:
    SgSymbols = sgtbx.SpaceGroupSymbols(LookupSymbol)
    HSym = SgSymbols.Hall()
    SgOps = sgtbx.SpaceGroup(HSym)
    SgInfo = SgOps.Info()
    print "SpaceGroup %s (%d) %s" % (
      SgInfo.BuildLookupSymbol(),
      SgInfo.SgNumber(),
      SgInfo.BuildHallSymbol())
    sys.stdout.flush()
    UnitCell = debug_utils.get_compatible_unit_cell(SgInfo, 1000).UnitCell
    SnapParametersLarge = \
    sgtbx.SpecialPositionSnapParameters(UnitCell, SgOps, 0, 2.0)
    SmallSnapDist2 = 1.e-5
    SnapParametersSmall = \
    sgtbx.SpecialPositionSnapParameters(UnitCell, SgOps, 0, SmallSnapDist2)
    WTab = sgtbx.WyckoffTable(SgInfo, 1)
    for i in xrange(20):
      RandomX = (random.uniform(-2,2),
                 random.uniform(-2,2),
                 random.uniform(-2,2))
      print "RandomX", debug_utils.format_round_scaled_list(RandomX)
      #
      SS = sgtbx.SiteSymmetry(SnapParametersLarge, RandomX, 1)
      SWMap = WTab.getWyckoffMapping(SS)
      print SWMap.WP().M(), SWMap.WP().Letter(), SWMap.WP().SpecialOp(),
      print SS.PointGroupType()
      WWMap = WTab.getWyckoffMapping(UnitCell, SgOps, SS.SnapPosition(), 0.1)
      assert SWMap.WP().Letter() == WWMap.WP().Letter()
      assert UnitCell.Distance2(SWMap.snap(RandomX), SS.SnapPosition()) < 1.e-5
      assert UnitCell.Distance2(WWMap.snap(RandomX), SS.SnapPosition()) < 1.e-5
      SES = sgtbx.SymEquivCoordinates(SS)
      x = SWMap.snap_to_representative(RandomX)
      d = SES.getShortestDistance2(UnitCell, x)
      assert d < 1.e-5
      x = WWMap.snap_to_representative(RandomX)
      d = SES.getShortestDistance2(UnitCell, x)
      assert d < 1.e-5
      #
      WWMap = WTab.getWyckoffMapping(UnitCell, SgOps, RandomX, 2.0)
      SS = sgtbx.SiteSymmetry(SnapParametersSmall, WWMap.snap(RandomX), 0)
      print WWMap.WP().M(), WWMap.WP().Letter(), WWMap.WP().SpecialOp(),
      print SS.PointGroupType()
      assert SS.DistanceMoved2() < SmallSnapDist2
      SWMap = WTab.getWyckoffMapping(SS)
      assert SWMap.WP().Letter() == WWMap.WP().Letter()
      assert UnitCell.Distance2(SWMap.snap(RandomX), SS.SnapPosition()) < 1.e-5
      assert UnitCell.Distance2(WWMap.snap(RandomX), SS.SnapPosition()) < 1.e-5
      SEW = sgtbx.SymEquivCoordinates(WWMap, RandomX)
      x = WWMap.snap_to_representative(RandomX)
      d = SEW.getShortestDistance2(UnitCell, x)
      assert d < 1.e-5
      x = SWMap.snap_to_representative(WWMap.snap(RandomX))
      d = SEW.getShortestDistance2(UnitCell, x)
      assert d < 1.e-5
    print

def run(timing=1):
  from tst import run_other
  run_other(sys.argv[1:], timing, OneCycle,
    ("Hall: -P 6 2c (z,1/2*x,1/2*y)",))

if (__name__ == "__main__"):
  run()
