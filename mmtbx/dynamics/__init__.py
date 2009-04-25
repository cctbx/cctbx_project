import cctbx.array_family.flex

import boost.python
ext = boost.python.import_ext("mmtbx_dynamics_ext")
from mmtbx_dynamics_ext import *

class kinetic_energy_and_temperature(object):

  def __init__(O, velocities, masses):
    O.kinetic_energy = kinetic_energy(velocities=velocities, masses=masses)
    k_boltz = 1.380662e-03 # XXX incorrect
    dof = 3 * velocities.size()
    if (dof == 0):
      O.temperature = 0
    else:
      O.temperature = 2 * O.kinetic_energy / (dof * k_boltz)
