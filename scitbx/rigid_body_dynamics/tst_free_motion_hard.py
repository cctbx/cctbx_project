from scitbx.rigid_body_dynamics import featherstone
from scitbx.rigid_body_dynamics import joint_lib
from scitbx.rigid_body_dynamics.utils import \
  center_of_mass_from_sites, \
  spatial_inertia_from_sites, \
  kinetic_energy
from scitbx.rigid_body_dynamics import test_utils
from free_motion_reference_impl import \
  body_inertia, \
  create_triangle_with_center_of_mass_at_origin
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, show_times_at_exit
import sys

def exercise_euler_params_qE_as_euler_angles_xyz_qE(mersenne_twister):
  for i_trial in xrange(30):
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    qr = matrix.col(mersenne_twister.random_double(size=3)-0.5)
    J = joint_lib.six_dof_euler_params(qE=qE, qr=qr)
    Jea = joint_lib.six_dof_euler_angles_xyz(qE=qE, qr=qr)
    assert approx_equal(Jea.E, J.E)
    assert approx_equal(Jea.r, J.r)
    Jep = joint_lib.six_dof_euler_params(qE=Jea.qE, qr=qr)
    assert approx_equal(Jep.E, J.E)
    assert approx_equal(Jep.r, J.r)

def exercise_T_as_X(mersenne_twister):
  for i_trial in xrange(10):
    T1 = matrix.rt((
      mersenne_twister.random_double_r3_rotation_matrix(),
      mersenne_twister.random_double(size=3)-0.5))
    T2 = matrix.rt((
      mersenne_twister.random_double_r3_rotation_matrix(),
      mersenne_twister.random_double(size=3)-0.5))
    T12 = T1 * T2
    T21 = T2 * T1
    T1i = T1.inverse_assuming_orthogonal_r()
    T2i = T2.inverse_assuming_orthogonal_r()
    X1 = joint_lib.T_as_X(T1)
    X2 = joint_lib.T_as_X(T2)
    X12 = joint_lib.T_as_X(T12)
    X21 = joint_lib.T_as_X(T21)
    X1i = joint_lib.T_as_X(T1i)
    X2i = joint_lib.T_as_X(T2i)
    assert approx_equal(X1*X2, X12)
    assert approx_equal(X2*X1, X21)
    assert approx_equal(X1.inverse(), X1i)
    assert approx_equal(X2.inverse(), X2i)

def create_triangle_with_random_center_of_mass(mersenne_twister):
  sites = create_triangle_with_center_of_mass_at_origin()
  t = matrix.col(mersenne_twister.random_double(size=3)*2-1)
  return [site+t for site in sites]

def create_wells(sites, mersenne_twister):
  "overall random rotation and translation + noise"
  r = matrix.sqr(mersenne_twister.random_double_r3_rotation_matrix())
  t = matrix.col(mersenne_twister.random_double(size=3)-0.5)
  wells = []
  for site in sites:
    t_noise = t + matrix.col(mersenne_twister.random_double(size=3)-0.5)*0.2
    wells.append(r * site + t_noise)
  return wells

class six_dof_alignment(object):

  def __init__(O, sites):
    c = center_of_mass_from_sites(sites=sites)
    O.T0b = matrix.rt(((1,0,0,0,1,0,0,0,1), -c))
    O.Tb0 = matrix.rt(((1,0,0,0,1,0,0,0,1), c))
    O.Xtree = featherstone.Xtrans(c)

class featherstone_system_model(object):

  def __init__(model, I, A, J):
    model.NB = 1
    model.pitch = [J]
    model.parent =[-1]
    model.Xtree = [A.Xtree]
    model.I = [I]

class simulation(object):

  def __init__(O, six_dof_joint, six_dof_r_is_qr, mersenne_twister):
    O.sites_F0 = create_triangle_with_random_center_of_mass(
      mersenne_twister=mersenne_twister)
    O.A = six_dof_alignment(sites=O.sites_F0)
    O.I_spatial = spatial_inertia_from_sites(
      sites=O.sites_F0, alignment_T=O.A.T0b)
    #
    O.wells = create_wells(sites=O.sites_F0, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    qr = matrix.col(mersenne_twister.random_double(size=3)-0.5)
    if (six_dof_r_is_qr):
      qr = joint_lib.RBDA_Eq_4_12(qE).transpose() * qr
    O.J = six_dof_joint(qE=qE, qr=qr, r_is_qr=six_dof_r_is_qr)
    O.v_spatial = matrix.col(mersenne_twister.random_double(size=6)*2-1)
    #
    O.energies_and_accelerations_update()

  def sites_moved(O):
    AJA = O.A.Tb0 * O.J.Tsp * O.A.T0b
    return [AJA * site for site in O.sites_F0]

  def energies_and_accelerations_update(O):
    O.e_kin = kinetic_energy(I_spatial=O.I_spatial, v_spatial=O.v_spatial)
    O.e_pot = test_utils.potential_energy(
      sites=O.sites_F0, wells=O.wells, A=O.A, J=O.J)
    f_ext_ff = test_utils.potential_f_ext_ff(
      sites=O.sites_F0, wells=O.wells, A=O.A, J=O.J)
    O.e_tot = O.e_kin + O.e_pot
    #
    e_pot_bf = test_utils.potential_energy_bf(
      sites=O.sites_F0, wells=O.wells, A=O.A, J=O.J)
    assert approx_equal(e_pot_bf, O.e_pot)
    f_ext_bf = test_utils.potential_f_ext_bf(
      sites=O.sites_F0, wells=O.wells, A=O.A, J=O.J)
    #
    model = featherstone_system_model(I=O.I_spatial, A=O.A, J=O.J)
    #
    q = [None] # already stored in joint as qE and qr
    qd = [O.v_spatial]
    tau = None
    grav_accn = [0,0,0]
    qdd_using_f_ext_ff = featherstone.FDab(
      model, q, qd, tau, [f_ext_ff], grav_accn, f_ext_in_ff=True)
    qdd_using_f_ext_bf = featherstone.FDab(
      model, q, qd, tau, [f_ext_bf], grav_accn, f_ext_in_ff=False)
    assert approx_equal(qdd_using_f_ext_bf, qdd_using_f_ext_ff)
    O.a_spatial = qdd_using_f_ext_bf[0]

  def dynamics_step(O, delta_t):
    O.v_spatial = O.J.time_step_velocity(
      v_spatial=O.v_spatial, a_spatial=O.a_spatial, delta_t=delta_t)
    O.J = O.J.time_step_position(v_spatial=O.v_spatial, delta_t=delta_t)
    O.energies_and_accelerations_update()

def run_simulation(
      out,
      six_dof_joint,
      six_dof_r_is_qr,
      mersenne_twister,
      n_dynamics_steps,
      delta_t):
  sim = simulation(
    six_dof_joint=six_dof_joint,
    six_dof_r_is_qr=six_dof_r_is_qr,
    mersenne_twister=mersenne_twister)
  sites_moved = [sim.sites_moved()]
  e_pots = flex.double([sim.e_pot])
  e_kins = flex.double([sim.e_kin])
  for i_step in xrange(n_dynamics_steps):
    sim.dynamics_step(delta_t=delta_t)
    sites_moved.append(sim.sites_moved())
    e_pots.append(sim.e_pot)
    e_kins.append(sim.e_kin)
  e_tots = e_pots + e_kins
  print >> out, six_dof_joint.__name__, "r_is_qr=%s" % str(six_dof_r_is_qr)
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
  return sites_moved, relative_range

def run_simulations(out, mersenne_twister, n_dynamics_steps, delta_t):
  mt_state = mersenne_twister.getstate()
  relative_ranges = []
  sites_moved_accu = []
  for six_dof_joint in [
        joint_lib.six_dof_euler_params,
        joint_lib.six_dof_euler_angles_xyz]:
    for six_dof_r_is_qr in [False, True]:
      mersenne_twister.setstate(mt_state)
      sites_moved, relative_range = run_simulation(
        out=out,
        six_dof_joint=six_dof_joint,
        six_dof_r_is_qr=six_dof_r_is_qr,
        mersenne_twister=mersenne_twister,
        n_dynamics_steps=n_dynamics_steps,
        delta_t=delta_t)
      sites_moved_accu.append(sites_moved)
      relative_ranges.append(relative_range)
  print >> out, "rms joints:"
  rms_max_list = flex.double()
  for other in sites_moved_accu[1:]:
    rms = flex.double()
    for sites_ref,sites_other in zip(sites_moved_accu[0], other):
      sites_ref = flex.vec3_double(sites_ref)
      rms.append(sites_ref.rms_difference(flex.vec3_double(sites_other)))
    rms.min_max_mean().show(out=out, prefix="  ")
    rms_max_list.append(flex.max(rms))
    print >> out
  sys.stdout.flush()
  return relative_ranges, rms_max_list

def exercise(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  exercise_euler_params_qE_as_euler_angles_xyz_qE(
    mersenne_twister=mersenne_twister)
  exercise_T_as_X(
    mersenne_twister=mersenne_twister)
  relative_ranges_accu = [flex.double() for i in xrange(4)]
  rms_max_list_accu = [flex.double() for i in xrange(3)]
  for i in xrange(n_trials):
    relative_ranges, rms_max_list = run_simulations(
      out=out,
      mersenne_twister=mersenne_twister,
      n_dynamics_steps=n_dynamics_steps,
      delta_t=delta_t)
    for r,a in zip(relative_ranges, relative_ranges_accu):
      a.append(r)
    for r,a in zip(rms_max_list, rms_max_list_accu):
      a.append(r)
  print >> out, "Accumulated results:"
  print >> out
  for i,accu in enumerate(relative_ranges_accu):
    print >> out, "relative ranges (%d):" % i
    accu.min_max_mean().show(out=out, prefix="  ")
    print >> out
  for i,accu in enumerate(rms_max_list_accu):
    print >> out, "rms max (0-%d):" % (i+1)
    accu.min_max_mean().show(out=out, prefix="  ")
    print >> out
  if (out is not sys.stdout):
    for accu in relative_ranges_accu:
      assert flex.max(accu) < 1.e-4
    for accu in rms_max_list_accu:
      assert flex.max(accu) < 1.e-4

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
  exercise(out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
