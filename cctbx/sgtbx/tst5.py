# $Id$

import sgtbx
from symbols import table_hall_std530

nLTrDict = {}
for HallSymbol in table_hall_std530:
  for Z in "PABCIRHF":
    hsym = HallSymbol[6] + Z + HallSymbol[8:]
    SgOps = sgtbx.SpaceGroup(hsym)
    n = SgOps.nLTr()
    print n, hsym
    try: nLTrDict[n] = nLTrDict[n] + 1
    except: nLTrDict[n] = 1
    Z2POp = SgOps.getZ2POp()
    PSgOps = SgOps.ChangeBasis(Z2POp)
    if (SgOps != PSgOps.ChangeBasis(Z2POp.swap())):
      import sys, tst2
      print Z2POp.M().as_xyz()
      print Z2POp.InvM().as_xyz()
      tst2.show(SgOps)
      tst2.show(PSgOps.ChangeBasis(Z2POp.swap()))
      sys.exit(1)
k = nLTrDict.keys()
k.sort()
print "Distribution of number of lattice centring vectors"
print "#vectors occurrences"
for n in k: print n, nLTrDict[n]
