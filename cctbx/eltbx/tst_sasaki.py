# $Id$

factor_eV_Angstrom = 6626.0755 * 2.99792458 / 1.60217733

def compare(a, b, tolerance = 1.e-6):
  assert abs(a-b) < tolerance, "a=%.6g, b=%.6g" % (a,b)

def verify(sasaki_table, wave_length, fp, fdp):
  fpfdp = sasaki_table(factor_eV_Angstrom / wave_length)
  compare(fpfdp.fp(), fp)
  compare(fpfdp.fdp(), fdp)

import math
from cctbx_boost.eltbx.sasaki import *
sasaki = Sasaki("SI")
assert sasaki.Label() == "Si"
assert sasaki.Z() == 14
verify(sasaki, 0.1, -0.0226, 0.0010) # wide first
verify(sasaki, 2.89, 0.3824, 1.0517) # wide last
verify(sasaki, 6.7289, -6.9495, 4.1042) # K edge
sasaki = Sasaki("ag")
assert sasaki.Label() == "Ag"
assert sasaki.Z() == 47
verify(sasaki, 3.2426, -8.6870, 14.1062) # L1 edge
verify(sasaki, 3.5169, -16.5123, 13.7723) # L2 edge
verify(sasaki, 3.6995, -29.9836, 10.5137) # L3 edge, right before edge
verify(sasaki, 3.6996, -32.2967, 3.1759) # L3 edge, right after edge
print "OK"
