from mmtbx import dynamics
from mmtbx.dynamics.constants import boltzmann_constant_akma
from scitbx.rigid_body.essence import tst_molecules
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import numbers_as_str
import random
import math
import sys

def show_velocity_scales(args):
  assert len(args) > 0
  for arg in args:
    simulation_index = int(arg)
    sim = tst_molecules.get_test_simulation_by_index(
      simulation_index,
      use_random_wells=False)
    print "q shifts (as used in refinement)"
    shift = 1
    packed_q = sim.pack_q()
    e_pot_0 = sim.e_pot
    for iq in xrange(len(packed_q)):
      q_orig = packed_q[iq]
      packed_q[iq] += shift
      sim.unpack_q(packed_q=packed_q)
      sim.e_pot_update()
      print iq, sim.e_pot
      packed_q[iq] = q_orig
    sim.unpack_q(packed_q=packed_q)
    sim.energies_and_accelerations_update()
    print "end e_pot:", sim.e_pot, sim.e_kin
    print
    print "time step 1"
    qd_e_kin_scales = flex.double()
    for B in sim.bodies:
      BJ0 = B.J
      qd0 = B.J.qd_zero
      qd = list(qd0)
      for iqd in xrange(len(qd)):
        shift = 1
        for i_pass in [0,1,2]:
          qd[iqd] = qd0[iqd] + shift
          B.qd = matrix.col(qd)
          qd[iqd] = qd0[iqd]
          B.J = BJ0.time_step_position(qd=B.qd, delta_t=1)
          sim.e_kin_update()
          sim.e_pot_update()
          print len(qd_e_kin_scales), sim.e_pot, sim.e_kin
          if (sim.e_kin == 0):
            break
          if (i_pass == 0):
            assert sim.e_kin != 0
            shift = qd_e_kin_scale = 1 / sim.e_kin**0.5
          else:
            shift *= 2**0.5
        qd_e_kin_scales.append(qd_e_kin_scale)
      B.J = BJ0
      B.qd = B.J.qd_zero
    sim.energies_and_accelerations_update()
    print "end e_pot:", sim.e_pot, sim.e_kin
    print "qd_e_kin_scales:", numbers_as_str(values=qd_e_kin_scales)
    print
    temp = 300
    qd_global = (0.5 * temp * boltzmann_constant_akma)**0.5
      # simplification of temperature_as_kinetic_energy() / dof
    rg = random.gauss
    i_qd = 0
    for B in sim.bodies:
      qd_new = []
      for qd in B.J.qd_zero:
        qd_new.append(qd + rg(mu=0, sigma=qd_e_kin_scales[i_qd]*qd_global))
        i_qd += 1
      B.qd = matrix.col(qd_new)
    assert i_qd == qd_e_kin_scales.size()
    sim.e_kin_update()
    dof = qd_e_kin_scales.size()
    print "e_kin, temp:", sim.e_kin, \
      dynamics.kinetic_energy_as_temperature(dof=dof, e=sim.e_kin)
    assert sim.e_kin != 0
    e_kin_target = dynamics.temperature_as_kinetic_energy(dof=dof, t=temp)
    factor = (e_kin_target / sim.e_kin)**0.5
    for B in sim.bodies:
      B.qd *= factor
    sim.energies_and_accelerations_update()
    print "e_kin, temp:", sim.e_kin, \
      dynamics.kinetic_energy_as_temperature(dof=dof, e=sim.e_kin)
    print

def assign_cartesian_velocities(temperature=300):
  masses = flex.double(5, 1)
  velocities = flex.vec3_double()
  velocities.reserve(masses.size())
  for j,mass in enumerate(masses):
    sigma = (boltzmann_constant_akma / mass * temperature)**0.5
    velocities.append([random.gauss(0, sigma) for i in (1,2,3)])
  kt = dynamics.kinetic_energy_and_temperature(velocities, masses)
  return kt.temperature

def exercise_cartesian_velocities():
  temps = flex.double()
  for i in xrange(100):
    temps.append(assign_cartesian_velocities())
  temps.min_max_mean().show()

def run(args):
  show_velocity_scales(args)
  exercise_cartesian_velocities()

if (__name__ == "__main__"):
  run(sys.argv[1:])
