# This script generates a list of non-standard space group settings.
# The settings are used for testing.
#
# usage: python make_settings.py > settings.py

import sgtbx

settings = [0]
for i in xrange(1, 231): settings.append({})

ListCBOp = []
for xyz in ("x,y,z", "z,x,y", "y,z,x"):
  ListCBOp.append(sgtbx.ChOfBasisOp(sgtbx.RTMx(xyz)))

nBuilt = 0
for i in sgtbx.SpaceGroupSymbolIterator():
  HallSymbol = i.Hall()
  for Z in "PABCIRHF":
    HSym = HallSymbol[0] + Z + HallSymbol[2:]
    for CBOp in ListCBOp:
      SgOps = sgtbx.SpaceGroup(HSym).ChangeBasis(CBOp)
      SgType = SgOps.getSpaceGroupType()
      settings[SgType.SgNumber()][SgOps.BuildLookupSymbol(SgType)] = 0
      nBuilt = nBuilt + 1
print "# nBuilt =", nBuilt

nNonRedundant = 0
print "settings = ("
for i in xrange(1, 231):
  print "#", i
  symbols = settings[i].keys()
  symbols.sort()
  for s in symbols:
    print "'" + s + "',"
    nNonRedundant = nNonRedundant + 1
print ")"
print "# nNonRedundant =", nNonRedundant
