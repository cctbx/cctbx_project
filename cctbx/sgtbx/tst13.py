# test StructureFactor

import sys, os
from cctbx_boost import sgtbx
from cctbx_boost import uctbx
import random

ShortCut = "--ShortCut" in sys.argv
StandardOnly = "--StandardOnly" in sys.argv
RandomSeed0 = "--RandomSeed0" in sys.argv
Endless = "--Endless" in sys.argv

if (RandomSeed0):
  random.seed(0)

def get_unitcell(SgInfo):
  if (143 <= SgInfo.SgNumber() < 195):
    RefUnitCell = uctbx.UnitCell((10, 10, 10, 90, 90, 120))
  else:
    RefUnitCell = uctbx.UnitCell((10, 10, 10, 90, 90, 90))
  return RefUnitCell.ChangeBasis(SgInfo.CBOp().M().as_tuple()[0])

if (ShortCut):
  settings = ("Hall: C 2y",)
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
    UnitCell = get_unitcell(SgInfo)
    SgOps.CheckUnitCell(UnitCell)
    SnapParameters = \
    sgtbx.SpecialPositionSnapParameters(UnitCell, SgOps, 0, 2.)
    for i in xrange(1):
      RandomX = (random.uniform(-2,2),
                 random.uniform(-2,2),
                 random.uniform(-2,2))
      print "RandomX ", RandomX
      SS = sgtbx.SiteSymmetry(SnapParameters, RandomX, 1)
      X = SS.SnapPosition()
      n = float(SS.M()) / SgOps.OrderZ()
      print "X ", X, n,
      if (SS.M() != SgOps.OrderZ()): print "special",
      print
      SEC = sgtbx.SymEquivCoordinates(SS)
      RandomH = (1, 2, 3)
      print "RandomH ", RandomH
      FC = SEC.StructureFactor(RandomH)
      FH = n * SgOps.StructureFactor(RandomH, X, (0,0,0,0,0,0))
      print FC, FH
      assert abs(FC - FH) < 1.e-5
  t = os.times()
  if (not ShortCut): print "u+s,u,s:", t[0] + t[1], t[0], t[1]

while 1:
  OneCycle()
  if (not Endless): break
