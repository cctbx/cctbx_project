import cctbx.array_family.flex

import boost.python
ext = boost.python.import_ext("mmtbx_dynamics_ext")
from mmtbx_dynamics_ext import *

class kinetic_energy_and_temperature(object):

  def __init__(O, velocities, masses):
    from mmtbx.dynamics.constants import boltzmann_constant_akma
    O.kinetic_energy = kinetic_energy(velocities=velocities, masses=masses)
    dof = 3 * velocities.size()
    if (dof == 0):
      O.temperature = 0
    else:
      O.temperature = 2 * O.kinetic_energy / (dof * boltzmann_constant_akma)
