# This script generates a list of non-standard space group settings.
# The settings are used for testing.
#
# usage: python make_settings.py > settings.py

import sgtbx

settings = [0]
for i in xrange(1, 231): settings.append({})

nBuilt = 0
for i in sgtbx.SpaceGroupSymbolIterator():
  HallSymbol = i.Hall()
  for Z in "PABCIRHF":
    HSym = HallSymbol[0] + Z + HallSymbol[2:]
    SgOps = sgtbx.SgOps(HSym)
    SgType = SgOps.getSpaceGroupType()
    settings[SgType.SgNumber()][SgOps.BuildLookupSymbol(SgType)] = 0
    nBuilt = nBuilt + 1

print "# nBuilt =", nBuilt
print "settings = ("
for i in xrange(1, 231):
  print "#", i
  symbols = settings[i].keys()
  symbols.sort()
  for s in symbols:
    print "'" + s + "',"
print ")"
