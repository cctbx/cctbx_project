import sys
from cctbx_boost import sgtbx
from cctbx.development import debug_utils

RBFerr = "out of rotation-base-factor range"
TBFerr = "out of translation-base-factor range"

def random_expand(SgOps):
  matrix_indices = range(SgOps.OrderZ())
  debug_utils.random.shuffle(matrix_indices)
  s = sgtbx.SpaceGroup()
  for i in matrix_indices:
    s.expandSMx(SgOps(i))
  assert s.OrderZ() == SgOps.OrderZ()
  return s

def OneCycle(settings, TidyCBOp = 0):
  print "Testing SpaceGroupInfo: modified Hall symbols"
  for LookupSymbol in settings:
    SgSymbols = sgtbx.SpaceGroupSymbols(LookupSymbol)
    HallSymbol = SgSymbols.Hall().strip()
    if (HallSymbol[0] != "-"): HallSymbol = " " + HallSymbol
    for Z in "PABCIRHF":
      HSym = HallSymbol[0] + Z + HallSymbol[2:]
      SgOps = sgtbx.SpaceGroup(HSym)
      SgNumber = SgOps.Info().SgNumber()
      RefSgOps = sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols(SgNumber))
      if (SgNumber < 75):
        RotOps = sgtbx.SpaceGroup('P 1')
      else:
        RotOps = sgtbx.SpaceGroup('P 3*')
      for Rot in RotOps:
        CBOp = sgtbx.ChOfBasisOp(Rot)
        SgOpsRot = SgOps.ChangeBasis(CBOp)
        print HSym, CBOp.M().as_xyz()
        s = SgOpsRot.ChangeBasis(CBOp)
        s = random_expand(s)
        print s.OrderZ()
        try:
          t = sgtbx.SpaceGroupInfo(s, TidyCBOp)
        except RuntimeError, e:
          e = str(e)
          if (e.find(RBFerr) >= 0 or e.find(TBFerr) >= 0):
            print e
          else:
            raise
        else:
          print "Space group number:", t.SgNumber()
          print "CBMx:", t.CBOp().M().as_xyz()
          print "InvCBMx:", t.CBOp().InvM().as_xyz()
          assert t.SgNumber() == SgNumber
          assert s.ChangeBasis(t.CBOp()) == RefSgOps
          assert s == RefSgOps.ChangeBasis(t.CBOp().swap())
          try:
            l = t.BuildLookupSymbol()
          except RuntimeError, e:
            e = str(e)
            if (e.find(RBFerr) >= 0 or e.find(TBFerr) >= 0):
              print e
            else:
              raise
          else:
            print "LookupSymbol:", l
            assert s == sgtbx.SpaceGroup(sgtbx.SpaceGroupSymbols(l))
        print

def run(timing=1):
  from tst import run_other
  run_other(sys.argv[1:], timing, OneCycle,
    ("Hall: P 2c 2",))

if (__name__ == "__main__"):
  run()
