# $Id$

from neutron import *
for symbol in ("D", "Cu", "Cd"):
  r = NeutronNews1992Record(symbol)
  print r.Symbol(), r.BoundCohScattLength(), r.AbsCrossSect()
