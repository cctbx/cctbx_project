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
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, show_times_at_exit
import math
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

  def d_pot_d_q(O):
    model = featherstone_system_model(bodies=O.bodies)
    q = [None]*len(O.bodies)
    qd = [B.J.qd_zero for B in O.bodies]
    qdd = [B.J.qdd_zero for B in O.bodies]
    grav_accn = [0,0,0]
    taus = featherstone.ID(model, q, qd, qdd, O.f_ext_bf, grav_accn)
    result = []
    for B,tau in zip(O.bodies, taus):
      tau_as_d_pot_d_q = getattr(B.J, "tau_as_d_pot_d_q", None)
      if (tau_as_d_pot_d_q is None):
        result.append(tau)
      else:
        result.append(tau_as_d_pot_d_q(tau=tau))
    return result

  def d_pot_d_q_via_finite_differences(O, eps=1.e-6):
    result = []
    for B in O.bodies:
      gs = []
      J_orig = B.J
      for iq in xrange(J_orig.q_size):
        fs = []
        for signed_eps in [eps, -eps]:
          B.J = J_orig.add_finite_difference(iq=iq, signed_eps=signed_eps)
          O.e_pot_and_f_ext_update()
          fs.append(O.e_pot)
        gs.append((fs[0]-fs[1])/(2*eps))
      B.J = J_orig
      result.append(matrix.col(gs))
    O.energies_and_accelerations_update()
    return result

  def check_d_pot_d_q(O):
    qdd_orig = O.qdd
    ana = O.d_pot_d_q()
    fin = O.d_pot_d_q_via_finite_differences()
    if (0):
      for a,f in zip(ana, fin):
        print "fin:", f.elems
        print "ana:", a.elems
      print
    assert approx_equal(ana, fin)
    assert approx_equal(O.qdd, qdd_orig)

class six_dof_body(object):

  def __init__(O, mersenne_twister, n_sites):
    if (n_sites > 0):
      O.sites = [matrix.col(mersenne_twister.random_double_point_on_sphere())]
    while (len(O.sites) != n_sites):
      O.sites.append(O.sites[0]
        + matrix.col(mersenne_twister.random_double_point_on_sphere()))
    O.A = joint_lib.six_dof_alignment(sites=O.sites)
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = test_utils.create_wells(
      sites=O.sites, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    qr = matrix.col(mersenne_twister.random_double(size=3)-0.5)
    O.J = joint_lib.six_dof(type="euler_params", qE=qE, qr=qr, r_is_qr=True)
    O.qd = matrix.col(mersenne_twister.random_double(size=6)*2-1)

class spherical_body(object):

  def __init__(O, mersenne_twister, n_sites):
    if (n_sites > 0):
      O.sites = [matrix.col(mersenne_twister.random_double_point_on_sphere())]
    while (len(O.sites) != n_sites):
      O.sites.append(O.sites[0]
        + matrix.col(mersenne_twister.random_double_point_on_sphere()))
    O.A = joint_lib.spherical_alignment(sites=O.sites)
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = test_utils.create_wells(
      sites=O.sites, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    O.J = joint_lib.spherical(type="euler_params", qE=qE)
    O.qd = matrix.col(mersenne_twister.random_double(size=3)*2-1)

class revolute_body(object):

  def __init__(O, mersenne_twister, prev=None):
    if (prev is None):
      normal = matrix.col(mersenne_twister.random_double_point_on_sphere())
      pivot = matrix.col(mersenne_twister.random_double_point_on_sphere()) * 0.6
    else:
      normal = prev.A.normal
      pivot = prev.A.pivot + normal * 0.3
    O.sites = [pivot + normal * (-0.8)]
    O.A = joint_lib.revolute_alignment(pivot=pivot, normal=normal)
    O.I = spatial_inertia_from_sites(sites=O.sites, alignment_T=O.A.T0b)
    #
    O.wells = test_utils.create_wells(
      sites=O.sites, mersenne_twister=mersenne_twister)
    #
    O.J = joint_lib.revolute(qE=matrix.col([math.pi/8]))
    O.qd = matrix.col([math.pi/12])

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
      if (out is not sys.stdout):
        assert relative_range < 1.e-4

def exercise_six_dof2(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  for n_sites in xrange(1,4):
    for i_trial in xrange(n_trials):
      body1 = six_dof_body(mersenne_twister=mersenne_twister, n_sites=n_sites)
      body1.parent = -1
      body2 = six_dof_body(mersenne_twister=mersenne_twister, n_sites=n_sites)
      body2.parent = 0
      sim = simulation(bodies=[body1, body2])
      print >> out, "six_dof2 number of sites:", n_sites
      relative_range = exercise_sim(
        out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
      if (out is not sys.stdout):
        assert relative_range < 1.e-4

def exercise_revolute(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  for i_trial in xrange(n_trials):
    body = revolute_body(mersenne_twister=mersenne_twister)
    body.parent = -1
    sim = simulation(bodies=[body])
    print >> out, "revolute:"
    relative_range = exercise_sim(
      out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
    if (out is not sys.stdout):
      assert relative_range < 1.e-4

def exercise_revolute2(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  for i_trial in xrange(n_trials):
    body1 = revolute_body(mersenne_twister=mersenne_twister)
    body1.parent = -1
    body2 = revolute_body(mersenne_twister=mersenne_twister, prev=body1)
    body2.parent = 0
    sim = simulation(bodies=[body1, body2])
    print >> out, "revolute2:"
    relative_range = exercise_sim(
      out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
    if (out is not sys.stdout):
      assert relative_range < 1.e-4

def exercise_spherical(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  for n_sites in xrange(1,3):
    for i_trial in xrange(n_trials):
      body = spherical_body(mersenne_twister=mersenne_twister, n_sites=n_sites)
      body.parent = -1
      sim = simulation(bodies=[body])
      print >> out, "spherical number of sites:", n_sites
      relative_range = exercise_sim(
        out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
      if (out is not sys.stdout):
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
  if (1):
    exercise_six_dof(
      out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  if (1):
    exercise_six_dof2(
      out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  if (1):
    exercise_spherical(
      out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  if (1):
    exercise_revolute(
      out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  if (1):
    exercise_revolute2(
      out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
