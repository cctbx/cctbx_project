# test CCP4_ReciprocalSpaceASU

import sys
import sgtbx

ShortCut = "--ShortCut" in sys.argv
StandardOnly = "--StandardOnly" in sys.argv
Endless = "--Endless" in sys.argv

if (ShortCut):
  settings = ("Hall: -P 6 2c (z,1/2*x,1/2*y)",)
else:
  from settings import * # see examples/python/make_settings.py

def OneCycle():

  for LookupSymbol in settings:
    if (StandardOnly and LookupSymbol[:5] == "Hall:"): continue
    SgSymbols = sgtbx.SpaceGroupSymbols(LookupSymbol)
    HSym = SgSymbols.Hall()
    SgOps = sgtbx.SgOps(HSym)
    SgType = SgOps.getSpaceGroupType()
    print "SpaceGroup %s (%d) %s" % (
      SgOps.BuildLookupSymbol(SgType),
      SgType.SgNumber(),
      SgOps.BuildHallSymbol(SgType))
    sys.stdout.flush()
    asu = sgtbx.CCP4_ReciprocalSpaceASU(SgOps)
    H = [0,0,0]
    for H[0] in range(-2, 3):
      for H[1] in range(-2, 3):
        for H[2] in range(-2, 3):
          SEMI = SgOps.getEquivMillerIndices(H)
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

while 1:
  OneCycle()
  if (not Endless): break
