import sys
from cctbx_boost import sgtbx
from cctbx_boost import miller

def OneCycle(settings):
  print "Testing ReciprocalSpaceASU"
  for LookupSymbol in settings:
    SgSymbols = sgtbx.SpaceGroupSymbols(LookupSymbol)
    HSym = SgSymbols.Hall()
    SgOps = sgtbx.SpaceGroup(HSym)
    SgInfo = SgOps.Info()
    print "SpaceGroup %s (%d) %s" % (
      SgInfo.BuildLookupSymbol(),
      SgInfo.SgNumber(),
      SgInfo.BuildHallSymbol())
    sys.stdout.flush()
    asu = sgtbx.ReciprocalSpaceASU(SgInfo)
    print asu.CBOp().M()
    print asu.CBOp().InvM()
    print asu.ReferenceASU().LaueGroupCode()
    print asu.ReferenceASU().representation()
    print asu.ReferenceASU().getCutParameters()
    H = [0,0,0]
    for H[0] in range(-2, 3):
      for H[1] in range(-2, 3):
        for H[2] in range(-2, 3):
          SEMI = miller.SymEquivIndices(SgOps, H)
          K = None
          for S in SEMI:
            HR = S.HR()
            if (asu.isInASU(HR)):
              K = HR
              break
            mHR = (-HR[0], -HR[1], -HR[2])
            if (asu.isInASU(mHR)):
              K = mHR
              break
          assert K != None
          for S in SEMI:
            HR = S.HR()
            assert (HR == K) == asu.isInASU(HR)
            mHR = (-HR[0], -HR[1], -HR[2])
            assert (mHR == K) == asu.isInASU(mHR)

def run(timing=1):
  from tst import run_other
  run_other(sys.argv[1:], timing, OneCycle,
    ("P 3 1 2",))

if (__name__ == "__main__"):
  run()
