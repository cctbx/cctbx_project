# $Id$

from symbols import table_hall_std530
from tst2 import parse, hkl

for hsym in table_hall_std530:
  print hsym
  sgo = parse(hsym[6:])
  if (sgo):
    sgo.makeTidy()
    print sgo.OrderZ()
    print sgo.isChiral()
    hkl(sgo, 1)
