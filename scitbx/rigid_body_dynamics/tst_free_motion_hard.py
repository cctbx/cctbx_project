from scitbx.rigid_body_dynamics import featherstone
from scitbx.rigid_body_dynamics import joint_lib
from scitbx.rigid_body_dynamics.utils import \
  spatial_inertia_from_sites, \
  kinetic_energy
from scitbx.rigid_body_dynamics import test_utils
from free_motion_reference_impl import \
  create_triangle_with_center_of_mass_at_origin
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, show_times_at_exit
import math
import sys

def exercise_euler_params_qE_as_euler_angles_xyz_qE(mersenne_twister):
  for i_trial in xrange(30):
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    qr = matrix.col(mersenne_twister.random_double(size=3)-0.5)
    J = joint_lib.six_dof(type="euler_params", qE=qE, qr=qr)
    Jea = joint_lib.six_dof(type="euler_angles_xyz", qE=qE, qr=qr)
    assert approx_equal(Jea.E, J.E)
    assert approx_equal(Jea.r, J.r)
    Jep = joint_lib.six_dof(type="euler_params", qE=Jea.qE, qr=qr)
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

class featherstone_system_model(object):

  def __init__(model, I, A, J):
    model.NB = 1
    model.pitch = [J]
    model.parent =[-1]
    model.Xtree = [joint_lib.T_as_X(A.T0b)]
    model.I = [I]

class simulation_mixin(object):

  def sites_moved(O):
    AJA = O.A.Tb0 * O.J.Tsp * O.A.T0b
    return [AJA * site for site in O.sites_F0]

  def energies_and_accelerations_update(O):
    if (O.J.S is None):
      v_spatial = O.qd
    else:
      v_spatial = O.J.S * O.qd
    O.e_kin = kinetic_energy(I_spatial=O.I_spatial, v_spatial=v_spatial)
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
    qd = [O.qd]
    tau = None
    grav_accn = [0,0,0]
    qdd_using_f_ext_ff = featherstone.FDab(
      model, q, qd, tau, [f_ext_ff], grav_accn, f_ext_in_ff=True)
    qdd_using_f_ext_bf = featherstone.FDab(
      model, q, qd, tau, [f_ext_bf], grav_accn, f_ext_in_ff=False)
    assert approx_equal(qdd_using_f_ext_bf, qdd_using_f_ext_ff)
    O.qdd = qdd_using_f_ext_bf[0]

  def dynamics_step(O, delta_t):
    O.qd = O.J.time_step_velocity(qd=O.qd, qdd=O.qdd, delta_t=delta_t)
    O.J = O.J.time_step_position(qd=O.qd, delta_t=delta_t)
    O.energies_and_accelerations_update()

class six_dof_simulation(simulation_mixin):

  def __init__(O, six_dof_type, r_is_qr, mersenne_twister):
    O.sites_F0 = create_triangle_with_random_center_of_mass(
      mersenne_twister=mersenne_twister)
    O.A = joint_lib.six_dof_alignment(sites=O.sites_F0)
    O.I_spatial = spatial_inertia_from_sites(
      sites=O.sites_F0, alignment_T=O.A.T0b)
    #
    O.wells = create_wells(sites=O.sites_F0, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    qr = matrix.col(mersenne_twister.random_double(size=3)-0.5)
    if (r_is_qr):
      qr = joint_lib.RBDA_Eq_4_12(qE).transpose() * qr
    O.J = joint_lib.six_dof(type=six_dof_type, qE=qE, qr=qr, r_is_qr=r_is_qr)
    O.qd = matrix.col(mersenne_twister.random_double(size=6)*2-1)
    #
    O.energies_and_accelerations_update()

class five_dof_simulation(simulation_mixin):

  def __init__(O, r_is_qr, mersenne_twister):
    O.sites_F0 = [matrix.col(mersenne_twister.random_double_point_on_sphere())]
    O.sites_F0.append(O.sites_F0[0]
      + matrix.col(mersenne_twister.random_double_point_on_sphere()))
    O.A = joint_lib.five_dof_alignment(sites=O.sites_F0)
    O.I_spatial = spatial_inertia_from_sites(
      sites=O.sites_F0, alignment_T=O.A.T0b)
    #
    O.wells = create_wells(sites=O.sites_F0, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col((mersenne_twister.random_double(size=3)*2-1)*math.pi/4)
    qr = matrix.col(mersenne_twister.random_double(size=3)-0.5)
    if (r_is_qr):
      qr = joint_lib.RBDA_Eq_4_7(q=qE).transpose() * qr
    O.J = joint_lib.five_dof(qE=qE, qr=qr, r_is_qr=r_is_qr)
    O.qd = matrix.col(mersenne_twister.random_double(size=5)*2-1)
    #
    O.energies_and_accelerations_update()

class five_six_dof_simulation(simulation_mixin):

  def __init__(O, six_dof_type, sim5):
    O.sites_F0 = sim5.sites_F0
    O.A = sim5.A
    O.I_spatial = sim5.I_spatial
    #
    O.wells = sim5.wells
    #
    O.J = joint_lib.six_dof(
      type=six_dof_type, qE=sim5.J.qE, qr=sim5.J.qr, r_is_qr=sim5.J.r_is_qr)
    #
    O.qd = sim5.J.S * sim5.qd
    #
    O.energies_and_accelerations_update()
    #
    assert approx_equal(O.sites_moved(), sim5.sites_moved())
    assert approx_equal(O.e_pot, sim5.e_pot)
    assert approx_equal(O.e_kin, sim5.e_kin)
    assert abs(O.qdd[2]) < 1.e-10
    assert approx_equal(O.qdd.elems[:2], sim5.qdd.elems[:2])
    assert approx_equal(O.qdd.elems[3:], sim5.qdd.elems[2:])

plot_number = [0]

def run_simulation(
      out,
      dof,
      six_dof_type,
      r_is_qr,
      mersenne_twister,
      n_dynamics_steps,
      delta_t):
  assert dof in [5,6]
  if (dof == 5):
    sim = five_dof_simulation(
      r_is_qr=r_is_qr,
      mersenne_twister=mersenne_twister)
    if (six_dof_type is None):
      sim_label = "five_dof(r_is_qr=%s)"
    else:
      sim = five_six_dof_simulation(six_dof_type=six_dof_type, sim5=sim)
      sim_label = 'five_six_dof(type="%s", r_is_qr=%%s)' % six_dof_type
  else:
    sim = six_dof_simulation(
      six_dof_type=six_dof_type,
      r_is_qr=r_is_qr,
      mersenne_twister=mersenne_twister)
    sim_label = 'six_dof(type="%s", r_is_qr=%%s)' % six_dof_type
  sim_label %= str(sim.J.r_is_qr)
  sites_moved = [sim.sites_moved()]
  e_pots = flex.double([sim.e_pot])
  e_kins = flex.double([sim.e_kin])
  for i_step in xrange(n_dynamics_steps):
    sim.dynamics_step(delta_t=delta_t)
    sites_moved.append(sim.sites_moved())
    e_pots.append(sim.e_pot)
    e_kins.append(sim.e_kin)
  e_tots = e_pots + e_kins
  print >> out, sim_label
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
  if (out is sys.stdout):
    l = sim_label \
      .replace(' ', "") \
      .replace('"', "") \
      .replace("(", "_") \
      .replace(")", "_") \
      .replace(",", "_")
    f = open("tmp%02d_%s.xy" % (plot_number[0], l), "w")
    for es in [e_pots, e_kins, e_tots]:
      for e in es:
        print >> f, e
      print >> f, "&"
    f.close()
    plot_number[0] += 1
  return sim, sim_label, sites_moved, relative_range

def run_simulations(
      out,
      dof,
      mersenne_twister,
      n_dynamics_steps,
      delta_t):
  mt_state = mersenne_twister.getstate()
  sim_labels = []
  sites_moved_accu = []
  relative_ranges = []
  for r_is_qr in [True, False]:
    for six_dof_type in [None, "euler_params", "euler_angles_xyz"]:
      if (dof == 6 and six_dof_type is None): continue
      mersenne_twister.setstate(mt_state)
      sim, sim_label, sites_moved, relative_range = run_simulation(
        out=out,
        dof=dof,
        six_dof_type=six_dof_type,
        r_is_qr=r_is_qr,
        mersenne_twister=mersenne_twister,
        n_dynamics_steps=n_dynamics_steps,
        delta_t=delta_t)
      sim_labels.append(sim_label)
      sites_moved_accu.append(sites_moved)
      relative_ranges.append(relative_range)
  rms_max_list = flex.double()
  for sim_label,other in zip(sim_labels[1:], sites_moved_accu[1:]):
    print >> out, "rms joints %s" % sim_labels[0]
    print >> out, "       vs. %s:" % sim_label
    rms = flex.double()
    for sites_ref,sites_other in zip(sites_moved_accu[0], other):
      sites_ref = flex.vec3_double(sites_ref)
      rms.append(sites_ref.rms_difference(flex.vec3_double(sites_other)))
    rms.min_max_mean().show(out=out, prefix="  ")
    rms_max_list.append(flex.max(rms))
    print >> out
  sys.stdout.flush()
  return sim_labels, relative_ranges, rms_max_list

def exercise_simulation(out, dof, n_trials, n_dynamics_steps, delta_t=0.001):
  mersenne_twister = flex.mersenne_twister(seed=0)
  sim_labels = None
  relative_ranges_accu = None
  rms_max_list_accu = None
  for i in xrange(n_trials):
    sim_labels_new, relative_ranges, rms_max_list = run_simulations(
      out=out,
      dof=dof,
      mersenne_twister=mersenne_twister,
      n_dynamics_steps=n_dynamics_steps,
      delta_t=delta_t)
    if (sim_labels is None):
      sim_labels = sim_labels_new
    else:
      assert sim_labels == sim_labels_new
    if (relative_ranges_accu is None):
      relative_ranges_accu=[flex.double() for i in xrange(len(relative_ranges))]
    else:
      assert len(relative_ranges) == len(relative_ranges_accu)
    for r,a in zip(relative_ranges, relative_ranges_accu):
      a.append(r)
    if (rms_max_list_accu is None):
      rms_max_list_accu = [flex.double() for i in xrange(len(rms_max_list))]
    else:
      assert len(rms_max_list) == len(rms_max_list_accu)
    for r,a in zip(rms_max_list, rms_max_list_accu):
      a.append(r)
  print >> out, "Accumulated results:"
  print >> out
  for sim_label,accu in zip(sim_labels, relative_ranges_accu):
    print >> out, "relative ranges %s:" % sim_label
    accu.min_max_mean().show(out=out, prefix="  ")
    print >> out
  for sim_label,accu in zip(sim_labels[1:], rms_max_list_accu):
    print >> out, "rms max %s" % sim_labels[0]
    print >> out, "    vs. %s:" % sim_label
    accu.min_max_mean().show(out=out, prefix="  ")
    print >> out
  if (out is not sys.stdout):
    for accu in relative_ranges_accu:
      assert flex.max(accu) < 1.e-4
    for i,accu in enumerate(rms_max_list_accu):
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
  mersenne_twister = flex.mersenne_twister(seed=0)
  exercise_euler_params_qE_as_euler_angles_xyz_qE(
    mersenne_twister=mersenne_twister)
  exercise_T_as_X(mersenne_twister=mersenne_twister)
  for dof in [6, 5]:
    exercise_simulation(
      out=out, dof=dof, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
