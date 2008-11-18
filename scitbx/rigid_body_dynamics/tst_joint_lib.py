from scitbx.rigid_body_dynamics import featherstone
from scitbx.rigid_body_dynamics import joint_lib
from scitbx.rigid_body_dynamics.utils import \
  spatial_inertia_from_sites, \
  kinetic_energy
from scitbx.rigid_body_dynamics.tst_free_motion_hard import \
  potential_energy, \
  potential_f_ext_pivot_at_origin
from scitbx.array_family import flex
from scitbx import matrix
import math
import sys

class featherstone_system_model(object):

  def __init__(model, A, I, J):
    model.NB = 1
    model.pitch = [J]
    model.parent =[-1]
    model.Xtree = [A.Xtree]
    model.I = [I]

class revolute_simulation(object):

  def __init__(O, mersenne_twister):
    def random_vector():
      return matrix.col(mersenne_twister.random_double(size=3)*2-1)
    def random_angle():
      return (mersenne_twister.random_double()*2-1)*math.pi
    #
    O.sites = [random_vector()]
    O.A = joint_lib.revolute_alignment(
      pivot=random_vector(),
      normal=random_vector().normalize())
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T)
    #
    O.wells = [random_vector()]
    #
    O.J = joint_lib.revolute(qE=matrix.col([random_angle()]))
    O.qd = matrix.col([random_angle()])
    #
    O.energies_and_accelerations_update()

  def sites_moved(O):
    T = O.A.T_inv * O.J.T * O.A.T
    return [T * site for site in O.sites]

  def energies_and_accelerations_update(O):
    O.e_kin = kinetic_energy(I_spatial=O.I, v_spatial=O.J.S*O.qd)
    O.e_pot = potential_energy(
      sites=O.sites, wells=O.wells, A_T=O.A.T, J_T_inv=O.J.T_inv)
    O.f_ext = potential_f_ext_pivot_at_origin(
      sites=O.sites, wells=O.wells, A_T=O.A.T, J_T_inv=O.J.T_inv)
    O.e_tot = O.e_kin + O.e_pot
    #
    model = featherstone_system_model(A=O.A, I=O.I, J=O.J)
    q = [None]
    tau = None
    grav_accn = [0,0,0]
    O.qdd = featherstone.FDab(model, q, [O.qd], tau, [O.f_ext], grav_accn)[0]

  def dynamics_step(O, delta_t):
    O.qd = O.J.time_step_velocity(qd=O.qd, qdd=O.qdd, delta_t=delta_t)
    O.J = O.J.time_step_position(qd=O.qd, delta_t=delta_t)
    O.energies_and_accelerations_update()

def exercise_revolute_sim(out, mersenne_twister, n_dynamics_steps, delta_t):
  sim = revolute_simulation(mersenne_twister=mersenne_twister)
  e_pots = flex.double([sim.e_pot])
  e_kins = flex.double([sim.e_kin])
  for i_step in xrange(n_dynamics_steps):
    sim.dynamics_step(delta_t=delta_t)
    e_pots.append(sim.e_pot)
    e_kins.append(sim.e_kin)
  e_tots = e_pots + e_kins
  out = sys.stdout
  print >> out, "energy samples:", e_tots.size()
  print >> out, "e_pot min, max:", min(e_pots), max(e_pots)
  print >> out, "e_kin min, max:", min(e_kins), max(e_kins)
  print >> out, "e_tot min, max:", min(e_tots), max(e_tots)
  print >> out, "start e_tot:", e_tots[0]
  print >> out, "final e_tot:", e_tots[-1]
  ave = flex.sum(e_tots) / e_tots.size()
  range = flex.max(e_tots) - flex.min(e_tots)
  relative_range = range / ave
  print >> out, "ave:", ave
  print >> out, "range:", range
  print >> out, "relative range:", relative_range
  print >> out
  out.flush()
  return relative_range

def exercise_revolute(
      out=sys.stdout,
      n_trials=10,
      n_dynamics_steps=100,
      delta_t=0.01):
  mersenne_twister = flex.mersenne_twister(seed=0)
  relative_ranges = flex.double()
  for i in xrange(n_trials):
    relative_ranges.append(exercise_revolute_sim(
      out=out,
      mersenne_twister=mersenne_twister,
      n_dynamics_steps=n_dynamics_steps,
      delta_t=delta_t))
  print >> out, "relative ranges:"
  relative_ranges.min_max_mean().show(out=out, prefix="  ")
  print >> out

def run(args):
  assert len(args) == 0
  exercise_revolute()
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
