# $Id$

from cctbx_boost.eltbx.wavelengths import *

print keV_as_Angstrom(1)
print Angstrom_as_keV(10)

i = 0
while 1:
  wl = WaveLength(i)
  if (wl() == 0.): break
  print wl.Label(), wl(), wl.Energy()
  i = i + 1

for lbl in ("CrA1", "Cu", "MoA1", "AgA2"):
  wl = WaveLength(lbl)
  print wl.Label(), wl(), wl.Energy()
