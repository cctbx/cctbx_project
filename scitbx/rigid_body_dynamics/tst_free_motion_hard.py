from scitbx.rigid_body_dynamics import featherstone
from scitbx.rigid_body_dynamics import joint_lib
from scitbx.rigid_body_dynamics.utils import \
  kinetic_energy
from scitbx.rigid_body_dynamics.tst_free_motion import \
  featherstone_system_model
from scitbx.rigid_body_dynamics.free_motion_reference_impl import \
  body_inertia, \
  create_triangle_with_center_of_mass_at_origin
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
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

def create_wells(sites, mersenne_twister):
  "overall random rotation and translation + noise"
  r = matrix.sqr(mersenne_twister.random_double_r3_rotation_matrix())
  t = matrix.col(mersenne_twister.random_double(size=3)-0.5)
  wells = []
  for site in sites:
    t_noise = t + matrix.col(mersenne_twister.random_double(size=3)-0.5)*0.2
    wells.append(r * site + t_noise)
  return wells

def potential_energy(sites, wells, A_T, J_T_inv):
  result = 0
  for s, w in zip(sites, wells):
    result += (A_T * s - J_T_inv * A_T * w).dot()
  return result

def potential_f_ext_pivot_at_origin(sites, wells, A_T, J_T_inv):
  f_cart = [2 * (A_T * s - J_T_inv * A_T * w) for s, w in zip(sites, wells)]
  f = matrix.col((0,0,0))
  nc = matrix.col((0,0,0))
  for s,force in zip(sites, f_cart):
    f += force
    nc += (A_T * s).cross(force)
  return matrix.col((nc, f)).resolve_partitions()

def potential_energy_no_align(sites, wells, J_T_inv):
  result = 0
  for s, w in zip(sites, wells):
    result += (s - J_T_inv * w).dot()
  return result

def potential_f_ext_no_align_pivot_at_origin(sites, wells, J_T_inv):
  f_cart = [2 * (s - J_T_inv * w) for s, w in zip(sites, wells)]
  f = matrix.col((0,0,0))
  nc = matrix.col((0,0,0))
  for s,force in zip(sites, f_cart):
    f += force
    nc += s.cross(force)
  return matrix.col((nc, f)).resolve_partitions()

class simulation(object):

  def __init__(O, six_dof_joint, mersenne_twister):
    O.sites = create_triangle_with_center_of_mass_at_origin()
    O.m = len(O.sites) # unit masses
    O.I = body_inertia(sites_cart=O.sites)
    #
    O.wells = create_wells(sites=O.sites, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    qr = matrix.col(mersenne_twister.random_double(size=3)-0.5)
    O.J = six_dof_joint(qE=qE, qr=qr)
    O.v_spatial = matrix.col(mersenne_twister.random_double(size=6)*2-1)
    #
    O.energies_and_accelerations_update()

  def sites_moved(O):
    return [O.J.T * site for site in O.sites]

  def energies_and_accelerations_update(O):
    O.e_kin = kinetic_energy(
      I_spatial=featherstone.mcI(m=O.m, c=matrix.col((0,0,0)), I=O.I),
      v_spatial=O.v_spatial)
    O.e_pot = potential_energy_no_align(
      sites=O.sites, wells=O.wells, J_T_inv=O.J.T_inv)
    O.f_ext = potential_f_ext_no_align_pivot_at_origin(
      sites=O.sites, wells=O.wells, J_T_inv=O.J.T_inv)
    O.e_tot = O.e_kin + O.e_pot
    #
    model = featherstone_system_model(m=O.m, I=O.I, J=O.J)
    q = [None] # already stored in joint as qE and qr
    qd = [O.v_spatial]
    tau = None
    grav_accn = [0,0,0]
    qdd = featherstone.FDab(model, q, qd, tau, [O.f_ext], grav_accn)
    O.a_spatial = qdd[0]

  def dynamics_step(O, delta_t):
    O.v_spatial = O.J.time_step_velocity(
      v_spatial=O.v_spatial, a_spatial=O.a_spatial, delta_t=delta_t)
    O.J = O.J.time_step_position(v_spatial=O.v_spatial, delta_t=delta_t)
    O.energies_and_accelerations_update()

def run_simulation(
      out, six_dof_joint, mersenne_twister, n_dynamics_steps, delta_t):
  sim = simulation(
    six_dof_joint=six_dof_joint,
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
  print >> out, six_dof_joint
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
    mersenne_twister.setstate(mt_state)
    sites_moved, relative_range = run_simulation(
      out=out,
      six_dof_joint=six_dof_joint,
      mersenne_twister=mersenne_twister,
      n_dynamics_steps=n_dynamics_steps,
      delta_t=delta_t)
    sites_moved_accu.append(sites_moved)
    relative_ranges.append(relative_range)
  print >> out, "rms joints:"
  rms = flex.double()
  for sites_ep,sites_ea in zip(*sites_moved_accu):
    sites_ep = flex.vec3_double(sites_ep)
    rms.append(sites_ep.rms_difference(flex.vec3_double(sites_ea)))
  rms.min_max_mean().show(out=out, prefix="  ")
  print >> out
  sys.stdout.flush()
  return relative_ranges, flex.max(rms)

def exercise(out, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  exercise_euler_params_qE_as_euler_angles_xyz_qE(
    mersenne_twister=mersenne_twister)
  relative_ranges_accu = [flex.double(), flex.double()]
  rms_max_accu = flex.double()
  for i in xrange(n_trials):
    relative_ranges, rms_max = run_simulations(
      out=out,
      mersenne_twister=mersenne_twister,
      n_dynamics_steps=n_dynamics_steps,
      delta_t=delta_t)
    for r,a in zip(relative_ranges, relative_ranges_accu):
      a.append(r)
    rms_max_accu.append(rms_max)
  for accu in relative_ranges_accu:
    print >> out, "relative ranges:"
    accu.min_max_mean().show(out=out, prefix="  ")
    print >> out
  print >> out, "rms max:"
  rms_max_accu.min_max_mean().show(out=out, prefix="  ")
  print >> out
  if (out is not sys.stdout):
    for accu in relative_ranges_accu:
      assert flex.max(accu) < 1.e-4
    assert flex.max(rms_max_accu) < 1.e-5

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
  exercise(out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
