# $Id$

import math
from cctbx_boost.eltbx.henke import *
henke = Henke("SI")
print henke.Label(), henke.Z()
for i in xrange(16):
  energy = math.pow(2., i)
  fpfdp = henke(energy)
  print "%5.3f %.6g %.6g" % (energy, fpfdp.fp(), fpfdp.fdp())
