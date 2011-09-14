import cctbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("mmtbx_dynamics_ext")
from mmtbx_dynamics_ext import *

def kinetic_energy_as_temperature(dof, e):
  from mmtbx.dynamics.constants import boltzmann_constant_akma as k
  return e / (0.5 * k * dof)

def temperature_as_kinetic_energy(dof, t):
  from mmtbx.dynamics.constants import boltzmann_constant_akma as k
  return t * (0.5 * k * dof)

class kinetic_energy_and_temperature(object):

  def __init__(O, velocities, masses):
    O.kinetic_energy = kinetic_energy(velocities=velocities, masses=masses)
    dof = 3 * velocities.size()
    if (dof == 0):
      O.temperature = 0
    else:
      O.temperature = kinetic_energy_as_temperature(
        dof=dof, e=O.kinetic_energy)
