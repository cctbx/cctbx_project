# $Id$

import sgtbx
from symbols import table_hall_std530
from tst2 import parse, hkl

for hsym in table_hall_std530:
  print hsym
  sgo = parse(hsym[6:])
  if (sgo):
    sgo.makeTidy()
    print "OrderZ:", sgo.OrderZ()
    print "isChiral:", sgo.isChiral()
    print "isEnantiomorphic:", sgo.isEnantiomorphic()
    ts = sgo.MatchTabulatedSettings()
    print ts.ExtendedHermann_Mauguin()
    symbols = sgtbx.SpaceGroupSymbols(ts.ExtendedHermann_Mauguin())
    assert hsym[6:] == symbols.Hall()
    hkl(sgo)
