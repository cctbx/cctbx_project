# $Id$

from cctbx_boost import sgtbx
from symbols import table_hall_std530
from tst2 import parse, hkl

for hsym in table_hall_std530:
  print hsym
  SgOps = parse(hsym[6:])
  if (SgOps):
    SgOps.makeTidy()
    SgInfo = SgOps.Info()
    print "OrderZ:", SgOps.OrderZ()
    print "isChiral:", SgOps.isChiral()
    print "isEnantiomorphic:", SgInfo.isEnantiomorphic()
    ts = SgOps.MatchTabulatedSettings()
    print ts.ExtendedHermann_Mauguin()
    symbols = sgtbx.SpaceGroupSymbols(ts.ExtendedHermann_Mauguin())
    assert hsym[6:] == symbols.Hall()
    hkl(SgOps)
