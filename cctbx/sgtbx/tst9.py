import sys
import random
from cctbx_boost import sgtbx
from cctbx.development import debug_utils

def OneCycle(settings, CheckLettersOnly=0, ignore_exceptions=0):
  print "Testing SiteSymmetry and WyckoffMapping 1"
  nUndetermined = 0
  nMismatches = 0
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
          assert not ignore_exceptions
        else:
          if (SE.M() == WP.M()): break
      print "      WX =", debug_utils.format_round_scaled_list(WX)
      print "  Ref WX =", debug_utils.format_round_scaled_list(
        SgInfo.CBOp()(WX))
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
          assert not ignore_exceptions
        else:
          WPX = WMap.WP()
          if (WP.Letter() != WPX.Letter()):
            print WP.Letter(), WP.M(), WP.SpecialOp(), "input"
            print WPX.Letter(),WPX.M(),WPX.SpecialOp(),"getWyckoffMapping"
            print "ERROR: Wyckoff letter mismatch!"
            print "ERROR: Ref WX =", SgInfo.CBOp()(WX)
            nMismatches = nMismatches + 1
            assert not ignore_exceptions
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

def run(timing=1):
  from tst import run_other
  run_other(sys.argv[1:], timing, OneCycle,
    ("P 61 2 2",))

if (__name__ == "__main__"):
  run()
