# $Id$

from symbols import table_hall_std530
from tst2 import parse
import sgtbx

def parse_xyz(xyz):
  s = sgtbx.parse_string(xyz)
  try:
    return sgtbx.RTMx(s)
  except RuntimeError, e:
    print "-->" + s.string() + "<--"
    print ("-" * (s.where() + 3)) + "^"
    print e.args[0]
    return None

def recycle(sgo):
  Decimal = 0
  TrFirst = 0
  LettersXYZ = "xyz"
  Separator = ","
  newsgo = sgtbx.SgOps()
  for iLIS in xrange(sgo.OrderZ()):
    xyz = sgo(iLIS).as_xyz(Decimal, TrFirst, LettersXYZ, Separator)
    rtmx = parse_xyz(xyz)
    if (not rtmx): return
    if (xyz != rtmx.as_xyz(Decimal, TrFirst, LettersXYZ, Separator)):
      print "Error: xyz mismatch."
      print xyz
      print rtmx.as_xyz(Decimal, TrFirst, LettersXYZ, Separator)
    newsgo.expandSMx(rtmx)
  print xyz
  if (   sgo.nLTr() != newsgo.nLTr()
      or sgo.fInv() != newsgo.fInv()
      or sgo.nSMx() != newsgo.nSMx()):
    print "Error: newsgo mismatch."
    print    sgo.nLTr(),    sgo.fInv(),    sgo.nSMx()
    print newsgo.nLTr(), newsgo.fInv(), newsgo.nSMx()

for hsym in table_hall_std530:
  print hsym
  sgo = parse(hsym[6:])
  if (sgo):
    recycle(sgo)
