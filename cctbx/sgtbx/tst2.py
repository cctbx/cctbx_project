# $Id$

import sys
import sgtbx

def parse(hall_symbol):
  s = sgtbx.parse_string(hall_symbol)
  try:
    return sgtbx.SgOps(s)
  except RuntimeError, e:
    print "-->" + s.string() + "<--"
    print ("-" * (s.where() + 3)) + "^"
    print e.args[0]
    return None

def show(sgo):
  print "nLTr:", sgo.nLTr();
  print "fInv:", sgo.fInv();
  print "nSMx:", sgo.nSMx();
  for iLIS in xrange(sgo.OrderZ()):
    print sgo(iLIS).as_xyz()

def show_RotMxInfo(sgo):
  for iLIS in xrange(sgo.OrderZ()):
    info = sgo(iLIS).getRotMxInfo()
    print "Rtype =", info.Rtype(),
    print "EV =", info.EV(),
    print "SenseOfRotation =", info.SenseOfRotation()

def hkl(sgo, FriedelSym = 0):
  H = (1,2,-3)
  semi = sgo.getEquivMillerIndices(H)
  for iList in xrange(semi.M(FriedelSym)):
    print semi(iList)
  for iList in xrange(semi.N()):
    print semi[iList].HR(), semi[iList].HT()
  CutP = sgo.getCutParameters(FriedelSym)
  print 'CutParameters =', CutP
  Master = sgo.getMasterIndex(H, CutP, 1)
  print Master.H()
  print Master.iMate()
  for Pretty in xrange(2):
    Master = sgo.getMasterIndex(H, CutP, Pretty)
    for iList in xrange(semi.M(FriedelSym)):
      assert(Master.H() == sgo.getMasterIndex(semi(iList), CutP, Pretty).H())

if (__name__ == '__main__'):
  sgo = parse(sys.argv[1])
  if (sgo):
    show(sgo)
    show_RotMxInfo(sgo)
    hkl(sgo)
