import sys
import random
from cctbx_boost import sgtbx
from cctbx.development import debug_utils

def OneCycle(settings):
  print "Testing shortest_difference_info"
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
    coor_ref = [random.random() for i in xrange(3)]
    for m in SgOps:
      coor_cmp = [x + random.randrange(5) for x in m.multiply(coor_ref)]
      dist_info = sgtbx.shortest_difference_info(
        UnitCell, SgOps, coor_ref, coor_cmp)
      assert dist_info.dist() < 1.e-5
      coor_apply = dist_info.apply(coor_cmp)
      for i in xrange(3):
        assert (coor_ref[i] - coor_apply[i]) < 1.e-5

def run(timing=1):
  from tst import run_other
  run_other(sys.argv[1:], timing, OneCycle,
    ("P 3 1 m",))

if (__name__ == "__main__"):
  run()
