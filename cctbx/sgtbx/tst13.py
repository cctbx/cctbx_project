import sys
import random
from cctbx_boost import sgtbx
from cctbx_boost import sftbx
from cctbx.development import debug_utils

def OneCycle(settings):
  print "Testing StructureFactor"
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
    SnapParameters = \
    sgtbx.SpecialPositionSnapParameters(UnitCell, SgOps, 0, 2.)
    for i in xrange(1):
      RandomX = (random.uniform(-2,2),
                 random.uniform(-2,2),
                 random.uniform(-2,2))
      print "RandomX ", debug_utils.format_round_scaled_list(RandomX)
      SS = sgtbx.SiteSymmetry(SnapParameters, RandomX, 1)
      X = SS.SnapPosition()
      n = float(SS.M()) / SgOps.OrderZ()
      print "X ", debug_utils.format_round_scaled_list(list(X) + [n]),
      if (SS.M() != SgOps.OrderZ()): print "special",
      print
      SEC = sgtbx.SymEquivCoordinates(SS)
      RandomH = (1, 2, 3)
      #print "RandomH ", RandomH
      FC = SEC.StructureFactor(RandomH)
      FH = n * sftbx.StructureFactor(SgOps, RandomH, X, (0,0,0,0,0,0))
      #print FC, FH
      assert abs(FC - FH) < 1.e-5

def run(timing=1):
  from tst import run_other
  run_other(sys.argv[1:], timing, OneCycle,
    ("Hall: C 2y",))

if (__name__ == "__main__"):
  run()
