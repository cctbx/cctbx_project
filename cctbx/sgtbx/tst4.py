# $Id$

from symbols import table_hall_std530
from tst2 import parse
from cctbx_boost import sgtbx

def parse_xyz(xyz):
  s = sgtbx.parse_string(xyz)
  try:
    return sgtbx.RTMx(s)
  except RuntimeError, e:
    print "-->" + s.string() + "<--"
    print ("-" * (s.where() + 3)) + "^"
    print e.args[0]
    return None

def recycle(SgOps):
  Decimal = 0
  TrFirst = 0
  LettersXYZ = "xyz"
  Separator = ","
  newSgOps = sgtbx.SpaceGroup()
  for iLIS in xrange(SgOps.OrderZ()):
    xyz = SgOps(iLIS).as_xyz(Decimal, TrFirst, LettersXYZ, Separator)
    rtmx = parse_xyz(xyz)
    if (not rtmx): return
    if (xyz != rtmx.as_xyz(Decimal, TrFirst, LettersXYZ, Separator)):
      print "Error: xyz mismatch."
      print xyz
      print rtmx.as_xyz(Decimal, TrFirst, LettersXYZ, Separator)
      raise
    newSgOps.expandSMx(rtmx)
  if (   SgOps.nLTr() != newSgOps.nLTr()
      or SgOps.fInv() != newSgOps.fInv()
      or SgOps.nSMx() != newSgOps.nSMx()):
    print "Error: newSgOps mismatch."
    print    SgOps.nLTr(),    SgOps.fInv(),    SgOps.nSMx()
    print newSgOps.nLTr(), newSgOps.fInv(), newSgOps.nSMx()
    raise

for hsym in table_hall_std530:
  print hsym
  SgOps = parse(hsym[6:])
  if (SgOps):
    recycle(SgOps)
