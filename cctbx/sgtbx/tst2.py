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

def hkl(sgo):
  H = (1,2,-3)
  semi = sgo.getEquivMillerIndices(H)
  for iList in xrange(semi.M(0)):
    print semi(iList)
  for iList in xrange(semi.N()):
    print semi[iList].HR(), semi[iList].HT()
  CutP = sgo.getCutParameters()
  print 'CutParameters =', CutP
  Master = sgo.getMasterIndex(H, CutP, 1)
  print Master.H()
  print Master.iMate()
  for Pretty in xrange(2):
    Master = sgo.getMasterIndex(H, CutP, Pretty)
    for iList in xrange(semi.M(0)):
      assert(Master.H() == sgo.getMasterIndex(semi(iList), CutP, Pretty).H())

def BuildIndices(SgOps):
  import uctbx
  UnitCell = uctbx.UnitCell((10, 10, 10, 90, 90, 90))
  MIG = sgtbx.MillerIndexGenerator(UnitCell, SgOps, 3)
  for H in MIG: print H

if (__name__ == '__main__'):
  sgo = parse(sys.argv[1])
  if (sgo):
    show(sgo)
    #show_RotMxInfo(sgo)
    #hkl(sgo)
    BuildIndices(sgo)
