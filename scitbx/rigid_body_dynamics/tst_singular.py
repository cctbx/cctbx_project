from scitbx.rigid_body_dynamics import featherstone
from scitbx.rigid_body_dynamics import joint_lib
from scitbx.rigid_body_dynamics.utils import \
  spatial_inertia_from_sites, \
  featherstone_system_model, \
  e_kin_from_model
from scitbx.rigid_body_dynamics import test_utils
from scitbx.rigid_body_dynamics.tst_joint_lib import exercise_sim
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import null_out, show_times_at_exit
import sys

class simulation(object):

  def __init__(O, bodies):
    O.bodies = bodies
    O.energies_and_accelerations_update()

  def energies_and_accelerations_update(O):
    model = featherstone_system_model(bodies=O.bodies)
    q = [None]*len(O.bodies)
    qd = [B.qd for B in O.bodies]
    #
    O.e_kin = e_kin_from_model(model, q, qd)
    O.e_pot_and_f_ext_update()
    #
    tau = None
    grav_accn = [0,0,0]
    O.qdd = featherstone.FDab(model, q, qd, tau, O.f_ext_bf, grav_accn)

  def e_pot_and_f_ext_update(O):
    AJA_accu = []
    O.e_pot = 0
    O.f_ext_bf = []
    for B in O.bodies:
      AJA = B.A.Tb0 * B.J.Tsp * B.A.T0b
      if (B.parent == -1):
        AJA_tree = None
      else:
        AJA_tree = AJA_accu[B.parent]
        AJA = AJA_tree * AJA
      AJA_accu.append(AJA)
      e_pot_bf = test_utils.potential_energy_bf(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree)
      f_ext_using_bf = test_utils.potential_f_ext_bf(
        sites=B.sites, wells=B.wells, A=B.A, J=B.J, AJA_tree=AJA_tree)
      O.f_ext_bf.append(f_ext_using_bf)
      O.e_pot += e_pot_bf
    O.e_tot = O.e_kin + O.e_pot

  def dynamics_step(O, delta_t):
    for B,qdd in zip(O.bodies, O.qdd):
      B.qd = B.J.time_step_velocity(qd=B.qd, qdd=qdd, delta_t=delta_t)
      B.J = B.J.time_step_position(qd=B.qd, delta_t=delta_t)
    O.energies_and_accelerations_update()

  def check_d_pot_d_q(O):
    pass

class six_dof_body(object):

  def __init__(O, mersenne_twister, n_sites):
    if (n_sites > 0):
      O.sites = [matrix.col(mersenne_twister.random_double_point_on_sphere())]
    while (len(O.sites) != n_sites):
      O.sites.append(O.sites[0]
        + matrix.col(mersenne_twister.random_double_point_on_sphere()))
    O.A = joint_lib.six_dof_alignment(sites=O.sites)
    O.I = spatial_inertia_from_sites(
      sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = test_utils.create_wells(
      sites=O.sites, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    qr = matrix.col(mersenne_twister.random_double(size=3)-0.5)
    O.J = joint_lib.six_dof(type="euler_params", qE=qE, qr=qr, r_is_qr=True)
    O.qd = matrix.col(mersenne_twister.random_double(size=6)*2-1)

def exercise_six_dof(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  for n_sites in xrange(1,4):
    for i_trial in xrange(n_trials):
      body = six_dof_body(mersenne_twister=mersenne_twister, n_sites=n_sites)
      body.parent = -1
      sim = simulation(bodies=[body])
      print >> out, "six_dof number of sites:", n_sites
      relative_range = exercise_sim(
        out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
      if (out is sys.stdout):
        assert relative_range < 1.e-4

def run(args):
  assert len(args) in [0,2]
  if (len(args) == 0):
    n_trials = 3
    n_dynamics_steps = 30
    out = null_out()
  else:
    n_trials = max(1, int(args[0]))
    n_dynamics_steps = max(1, int(args[1]))
    out = sys.stdout
  show_times_at_exit()
  exercise_six_dof(
    out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
