# $Id$

import sys, string

TidyCBOp = "--TidyCBOp" in sys.argv
QuickMode = "--Quick" in sys.argv
ShortCut = "--ShortCut" in sys.argv

import sgtbx

#Shifts = (0, 1, 2, 3, 5, 7, 11)
#Shifts = (1, 3)
Shifts = (1,)

RBF = 36
TBF = 144

if (ShortCut):
  table_hall = (" P 2c 2",)
else:
  table_hall = []
  for i in sgtbx.SpaceGroupSymbolIterator():
    table_hall.append(i.Hall())

for HallSymbol in table_hall:
  for Z in "PABCIRHF":
    HSym = HallSymbol[0] + Z + HallSymbol[2:]
    SgOps = sgtbx.SgOps(HSym)
    SgNumber = SgOps.getSpaceGroupType().SgNumber()
    RefSgOps = sgtbx.SgOps(sgtbx.SpaceGroupSymbols(SgNumber).Hall())
    if (SgNumber < 75):
      RotOps = sgtbx.SgOps('P 1')
    else:
      RotOps = sgtbx.SgOps('P 3*')
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
            t = s.getSpaceGroupType(TidyCBOp, RBF, TBF)
            print "Space group number:", t.SgNumber()
            print "CBMx:", t.CBOp().M().as_xyz()
            print "InvCBMx:", t.CBOp().InvM().as_xyz()
            assert t.SgNumber() == SgNumber
            if (not QuickMode):
              assert s.ChangeBasis(t.CBOp()) == RefSgOps
              assert s == RefSgOps.ChangeBasis(t.CBOp().swap())
            h = s.BuildHallSymbol(t, TidyCBOp)
            print "BuildHallSymbol:", h
            if (not QuickMode):
              assert s == sgtbx.SgOps(h)
            print
