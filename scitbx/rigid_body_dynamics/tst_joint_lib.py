from scitbx.rigid_body_dynamics import featherstone
from scitbx.rigid_body_dynamics import joint_lib
from scitbx.rigid_body_dynamics.utils import \
  spatial_inertia_from_sites, \
  kinetic_energy
from scitbx.rigid_body_dynamics.test_utils import \
  potential_energy, \
  potential_f_ext_pivot_at_origin
from scitbx.array_family import flex
from scitbx import matrix
import math
import sys

def accumulate_AJA(bodies):
  result = []
  for B in bodies:
    AJA_tree = B.A.T_inv * B.J.T * B.A.T
    if (B.parent != -1):
      AJA_tree = AJA_tree * result[B.parent]
    result.append(AJA_tree)
  return result

def accumulate_Ttree(bodies):
  result = []
  for B in bodies:
    if (B.parent == -1):
      Ttree = B.A.T_inv
    else:
      Ttree = bodies[B.parent].A.T * B.A.T_inv
    result.append(Ttree)
  return result

class featherstone_system_model(object):

  def __init__(model, bodies):
    model.NB = len(bodies)
    model.pitch = []
    model.parent =[]
    model.Xtree = []
    model.I = []
    for B,Ttree in zip(bodies, accumulate_Ttree(bodies)):
      model.pitch.append(B.J)
      model.parent.append(B.parent)
      model.Xtree.append(joint_lib.T_as_X(Ttree))
      model.I.append(B.I)

class random_revolute(object):

  def __init__(O, mersenne_twister):
    def random_vector():
      return matrix.col(mersenne_twister.random_double(size=3)*2-1)
    def random_angle():
      return (mersenne_twister.random_double()*2-1)*math.pi
    #
    O.sites = [random_vector()]
    for i_trial in xrange(100): # guard against unlikely singularity
      O.A = joint_lib.revolute_alignment(
        pivot=random_vector(),
        normal=random_vector().normalize())
      if (abs(O.A.normal.cos_angle(O.sites[0] - O.A.pivot)) > 1.e-3):
        break
    else:
      raise RuntimeError
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T)
    #
    O.wells = [random_vector()]
    #
    O.J = joint_lib.revolute(qE=matrix.col([random_angle()]))
    O.qd = matrix.col([random_angle()])

class revolute_simulation(object):

  def __init__(O, mersenne_twister, NB):
    O.bodies = []
    for ib in xrange(NB):
      B = random_revolute(mersenne_twister=mersenne_twister)
      B.parent = -1
      O.bodies.append(B)
    O.energies_and_accelerations_update()

  def energies_and_accelerations_update(O):
    O.e_kin = 0
    O.e_pot = 0
    f_ext = []
    AJA_tree_list = accumulate_AJA(bodies=O.bodies)
    for B in O.bodies:
      O.e_kin += kinetic_energy(I_spatial=B.I, v_spatial=B.J.S*B.qd)
      if (B.parent == -1): AJA_tree = None
      else:                AJA_tree = AJA_tree_list[B.parent]
      O.e_pot += potential_energy(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree)
      f_ext.append(potential_f_ext_pivot_at_origin(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree))
    O.e_tot = O.e_kin + O.e_pot
    #
    model = featherstone_system_model(bodies=O.bodies)
    q = [None]*len(O.bodies)
    qd = [B.qd for B in O.bodies]
    tau = None
    grav_accn = [0,0,0]
    O.qdd = featherstone.FDab(
      model, q, qd, tau, f_ext, grav_accn, f_ext_in_ff=True)

  def dynamics_step(O, delta_t):
    for B,qdd in zip(O.bodies, O.qdd):
      B.qd = B.J.time_step_velocity(qd=B.qd, qdd=qdd, delta_t=delta_t)
      B.J = B.J.time_step_position(qd=B.qd, delta_t=delta_t)
    O.energies_and_accelerations_update()

def exercise_revolute_sim(out, mersenne_twister, n_dynamics_steps, delta_t):
  sim = revolute_simulation(mersenne_twister=mersenne_twister, NB=2)
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
