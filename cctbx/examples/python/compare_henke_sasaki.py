# This script was used for debugging.
# The f' and f" obtained from the Henke tables and the Sasaki tables are
# printed side by side.

import sys
import math

from eltbx.tinypse import TinyPSE
from eltbx.henke import Henke
from eltbx.sasaki import Sasaki
from eltbx.wavelengths import keV_as_Angstrom

# determine logarithmic step size for loop over energies
steps = 20
log_energy_start = math.log(1000)
log_energy_stop = math.log(100000)
log_engergy_step = (log_energy_stop - log_energy_start) / steps

for Z in xrange(1, 93): # loop over periodic system of elements
  element_symbol = TinyPSE(Z).Symbol()
  try:
    h = Henke(element_symbol, 1)
    s = Sasaki(element_symbol, 1)
  except:
    print "Skipping:", Z, element_symbol
  else:
    assert(h.Z() == Z and s.Z() == Z)
    print Z, h.Label()
    print "Energy[eV] Wavelength[A] Henke Sasaki"
    for i in xrange(steps + 1): # loop over energies
      energy = math.exp(log_energy_start + i * log_engergy_step)
      print energy, keV_as_Angstrom(energy / 1000), h(energy)(), s(energy)()
