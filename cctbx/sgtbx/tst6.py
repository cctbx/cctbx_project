# $Id$

import sys, string

ShortCut = "--ShortCut" in sys.argv
StandardOnly = "--StandardOnly" in sys.argv
TidyCBOp = "--TidyCBOp" in sys.argv
QuickMode = "--Quick" in sys.argv
Endless = "--Endless" in sys.argv

from cctbx_boost import sgtbx

#Shifts = (0, 1, 2, 3, 5, 7, 11)
#Shifts = (1, 3)
Shifts = (1,)

RBF = 36
TBF = 144

if (ShortCut):
  settings = ("Hall:  P 2c 2",)
else:
  from settings import * # see examples/python/make_settings.py

def OneCycle():
  print "# Starting over"
  for LookupSymbol in settings:
    if (StandardOnly and LookupSymbol[:5] == "Hall:"): continue
    SgSymbols = sgtbx.SpaceGroupSymbols(LookupSymbol)
    HSym = SgSymbols.Hall()
    SgOps = sgtbx.SpaceGroup(HSym)
    SgNumber = SgOps.Info().SgNumber()
    RefSgOps = sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols(SgNumber))
    if (SgNumber < 75):
      RotOps = sgtbx.SpaceGroup('P 1')
    else:
      RotOps = sgtbx.SpaceGroup('P 3*')
    for Rot in RotOps:
      CBOpRot = sgtbx.ChOfBasisOp(Rot)
      SgOpsRot = SgOps.ChangeBasis(CBOpRot)
      for sx in (Shifts):
        for sy in (Shifts):
          for sz in (Shifts):
            xyz = "x+%d/12,y+%d/12,z+%d/12" % (sx, sy, sz)
            CBOp = sgtbx.ChOfBasisOp(sgtbx.RTMx(xyz))
            print HSym, CBOpRot.M().as_xyz(), CBOp.M().as_xyz()
            s = SgOpsRot.ChangeBasis(CBOp)
            print s.OrderZ()
            t = sgtbx.SpaceGroupInfo(s, TidyCBOp, RBF, TBF)
            print "Space group number:", t.SgNumber()
            print "CBMx:", t.CBOp().M().as_xyz()
            print "InvCBMx:", t.CBOp().InvM().as_xyz()
            assert t.SgNumber() == SgNumber
            if (not QuickMode):
              assert s.ChangeBasis(t.CBOp()) == RefSgOps
              assert s == RefSgOps.ChangeBasis(t.CBOp().swap())
            h = s.Info().BuildHallSymbol(TidyCBOp)
            print "BuildHallSymbol:", h
            if (not QuickMode):
              assert s == sgtbx.SpaceGroup(h)
            print

while 1:
  OneCycle()
  if (not Endless): break
