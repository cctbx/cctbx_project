# $Id$

import math
from cctbx.eltbx.sasaki import *
sasaki = Sasaki("SI")
print sasaki.Label(), sasaki.Z()
for i in xrange(16):
  energy = math.pow(2., i)
  fpfdp = sasaki(energy)
  print "%5.3f %.6g %.6g" % (energy, fpfdp.fp(), fpfdp.fdp())
