# test StructureSeminvariant

import sys, os
from cctbx_boost import sgtbx

ShortCut = "--ShortCut" in sys.argv
StandardOnly = "--StandardOnly" in sys.argv
Endless = "--Endless" in sys.argv

if (ShortCut):
  settings = ("Hall: -P 1 (-1/2*x+1/2*y+1/2*z,1/2*x-1/2*y+1/2*z,1/2*x+1/2*y-1/2*z)",)
else:
  from settings import * # see examples/python/make_settings.py

def OneCycle():

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
    ss = sgtbx.StructureSeminvariant(SgOps)
    for i in xrange(ss.size()):
      print ss.V(i), ss.M(i)
  t = os.times()
  if (not ShortCut): print "u+s,u,s:", t[0] + t[1], t[0], t[1]

while 1:
  OneCycle()
  if (not Endless): break
