import sys
from cctbx_boost import sgtbx

def OneCycle(settings):
  print "Testing StructureSeminvariant"
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
    ss = sgtbx.StructureSeminvariant(SgOps)
    for i in xrange(ss.size()):
      print ss.V(i), ss.M(i)

def run(timing=1):
  from tst import run_other
  run_other(sys.argv[1:], timing, OneCycle,
    ("Hall: -P 1 (-1/2*x+1/2*y+1/2*z,1/2*x-1/2*y+1/2*z,1/2*x+1/2*y-1/2*z)",))

if (__name__ == "__main__"):
  run()
