from __future__ import absolute_import, division, print_function
from scitbx.rigid_body.proto import featherstone
from scitbx.rigid_body.proto import joint_lib
from scitbx.rigid_body.proto.utils import \
  spatial_inertia_from_sites, \
  kinetic_energy, \
  T_as_X
from scitbx.rigid_body.proto import test_utils
from scitbx.rigid_body.proto.free_motion_reference_impl import \
  create_triangle_with_center_of_mass_at_origin
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, show_times_at_exit
import sys
from six.moves import range
from six.moves import zip

def exercise_euler_params_qE_as_euler_angles_xyz_qE(mersenne_twister):
  for i_trial in range(30):
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
  for i_trial in range(10):
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
    X1 = T_as_X(T1)
    X2 = T_as_X(T2)
    X12 = T_as_X(T12)
    X21 = T_as_X(T21)
    X1i = T_as_X(T1i)
    X2i = T_as_X(T2i)
    assert approx_equal(X1*X2, X12)
    assert approx_equal(X2*X1, X21)
    assert approx_equal(X1.inverse(), X1i)
    assert approx_equal(X2.inverse(), X2i)

def create_triangle_with_random_center_of_mass(mersenne_twister):
  sites = create_triangle_with_center_of_mass_at_origin()
  t = matrix.col(mersenne_twister.random_double(size=3)*2-1)
  return [site+t for site in sites]

class featherstone_system_model(object):

  def __init__(model, I, A, J):
    model.NB = 1
    model.pitch = [J]
    model.parent =[-1]
    model.Xtree = [T_as_X(A.T0b)]
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
    O.f_ext_bf = f_ext_bf
    O.qdd = qdd_using_f_ext_bf[0]

  def dynamics_step(O, delta_t):
    O.qd = O.J.time_step_velocity(qd=O.qd, qdd=O.qdd, delta_t=delta_t)
    O.J = O.J.time_step_position(qd=O.qd, delta_t=delta_t)
    O.energies_and_accelerations_update()

  def d_pot_d_q(O):
    model = featherstone_system_model(I=O.I_spatial, A=O.A, J=O.J)
    q = [None] # already stored in joint as qE and qr
    qd = [matrix.zeros(n=6)]
    qdd = [matrix.zeros(n=6)]
    grav_accn = [0,0,0]
    tau = featherstone.ID(model, q, qd, qdd, [O.f_ext_bf], grav_accn)[0]
    assert approx_equal(tau, -O.f_ext_bf)
    return O.J.tau_as_d_pot_d_q(tau=tau)

  def d_pot_d_q_via_finite_differences(O, eps=1.e-6):
    result = []
    for q in [O.J.qE, O.J.qr]:
      for i in range(len(q)):
        fs = []
        for signed_eps in [eps, -eps]:
          q_eps = list(q)
          q_eps[i] += signed_eps
          q_eps = matrix.col(q_eps)
          if (q is O.J.qE): qE = q_eps; qr=O.J.qr
          else:             qE = O.J.qE; qr=q_eps
          J = joint_lib.six_dof(
            type=O.J.type, qE=qE, qr=qr, r_is_qr=O.J.r_is_qr)
          e_pot = test_utils.potential_energy(
            sites=O.sites_F0, wells=O.wells, A=O.A, J=J)
          fs.append(e_pot)
        result.append((fs[0]-fs[1])/(2*eps))
    return matrix.col(result)

  def check_d_pot_d_q(O):
    ana = O.d_pot_d_q()
    fin = O.d_pot_d_q_via_finite_differences()
    assert approx_equal(ana, fin)

class six_dof_simulation(simulation_mixin):

  def __init__(O, six_dof_type, r_is_qr, mersenne_twister):
    O.sites_F0 = create_triangle_with_random_center_of_mass(
      mersenne_twister=mersenne_twister)
    O.A = joint_lib.six_dof_alignment(sites=O.sites_F0)
    O.I_spatial = spatial_inertia_from_sites(
      sites=O.sites_F0, alignment_T=O.A.T0b)
    #
    O.wells = test_utils.create_wells(
      sites=O.sites_F0, mersenne_twister=mersenne_twister)
    #
    qE = matrix.col(mersenne_twister.random_double(size=4)).normalize()
    qr = matrix.col(mersenne_twister.random_double(size=3)-0.5)
    if (r_is_qr):
      qr = joint_lib.RBDA_Eq_4_12(qE).transpose() * qr
    O.J = joint_lib.six_dof(type=six_dof_type, qE=qE, qr=qr, r_is_qr=r_is_qr)
    O.qd = matrix.col(mersenne_twister.random_double(size=6)*2-1)
    #
    O.energies_and_accelerations_update()

plot_prefix = 0
plot_number = [0]

def run_simulation(
      out,
      six_dof_type,
      r_is_qr,
      mersenne_twister,
      n_dynamics_steps,
      delta_t):
  sim = six_dof_simulation(
    six_dof_type=six_dof_type,
    r_is_qr=r_is_qr,
    mersenne_twister=mersenne_twister)
  sim_label = 'six_dof(type="%s", r_is_qr=%s)' % (
    six_dof_type, str(sim.J.r_is_qr))
  sim.check_d_pot_d_q()
  sites_moved = [sim.sites_moved()]
  e_pots = flex.double([sim.e_pot])
  e_kins = flex.double([sim.e_kin])
  for i_step in range(n_dynamics_steps):
    sim.dynamics_step(delta_t=delta_t)
    sites_moved.append(sim.sites_moved())
    e_pots.append(sim.e_pot)
    e_kins.append(sim.e_kin)
  e_tots = e_pots + e_kins
  sim.check_d_pot_d_q()
  print(sim_label, file=out)
  print("e_pot min, max:", min(e_pots), max(e_pots), file=out)
  print("e_kin min, max:", min(e_kins), max(e_kins), file=out)
  print("e_tot min, max:", min(e_tots), max(e_tots), file=out)
  print("start e_tot:", e_tots[0], file=out)
  print("final e_tot:", e_tots[-1], file=out)
  ave = flex.sum(e_tots) / e_tots.size()
  range_ = flex.max(e_tots) - flex.min(e_tots)
  relative_range = range_ / ave
  print("ave:", ave, file=out)
  print("range:", range_, file=out)
  print("relative range:", relative_range, file=out)
  print(file=out)
  out.flush()
  if (out is sys.stdout):
    l = sim_label \
      .replace(' ', "") \
      .replace('"', "") \
      .replace("(", "_") \
      .replace(")", "_") \
      .replace(",", "_")
    f = open("tmp_%02d_%02d_%s.xy" % (plot_prefix, plot_number[0], l), "w")
    for es in [e_pots, e_kins, e_tots]:
      for e in es:
        print(e, file=f)
      print("&", file=f)
    f.close()
    plot_number[0] += 1
  return sim, sim_label, sites_moved, e_tots, relative_range

def run_simulations(
      out,
      mersenne_twister,
      n_dynamics_steps,
      delta_t):
  mt_state = mersenne_twister.getstate()
  sim_labels = []
  sites_moved_accu = []
  e_tots_list = []
  relative_ranges = []
  for r_is_qr in [True, False]:
    for six_dof_type in ["euler_params", "euler_angles_xyz"]:
      mersenne_twister.setstate(mt_state)
      sim, sim_label, sites_moved, e_tots, relative_range = run_simulation(
        out=out,
        six_dof_type=six_dof_type,
        r_is_qr=r_is_qr,
        mersenne_twister=mersenne_twister,
        n_dynamics_steps=n_dynamics_steps,
        delta_t=delta_t)
      sim_labels.append(sim_label)
      e_tots_list.append(e_tots)
      sites_moved_accu.append(sites_moved)
      relative_ranges.append(relative_range)
  rms_max_list = flex.double()
  for sim_label,other in zip(sim_labels[1:], sites_moved_accu[1:]):
    print("rms joints %s" % sim_labels[0], file=out)
    print("       vs. %s:" % sim_label, file=out)
    rms = flex.double()
    for sites_ref,sites_other in zip(sites_moved_accu[0], other):
      sites_ref = flex.vec3_double(sites_ref)
      rms.append(sites_ref.rms_difference(flex.vec3_double(sites_other)))
    rms.min_max_mean().show(out=out, prefix="  ")
    rms_max_list.append(flex.max(rms))
    print(file=out)
  out.flush()
  return sim_labels, e_tots_list, relative_ranges, rms_max_list

def exercise_simulation(
      out, n_trials, n_dynamics_steps, delta_t=0.0001, random_seed=0):
  mersenne_twister = flex.mersenne_twister(seed=random_seed)
  sim_labels = None
  relative_ranges_accu = None
  rms_max_list_accu = None
  for i_trial in range(n_trials):
    sim_labels_new, e_tots_list, \
    relative_ranges, rms_max_list = run_simulations(
      out=out,
      mersenne_twister=mersenne_twister,
      n_dynamics_steps=n_dynamics_steps,
      delta_t=delta_t)
    if (sim_labels is None):
      sim_labels = sim_labels_new
    else:
      assert sim_labels == sim_labels_new
    if (relative_ranges_accu is None):
      relative_ranges_accu=[flex.double() for i in range(len(relative_ranges))]
    else:
      assert len(relative_ranges) == len(relative_ranges_accu)
    for r,a in zip(relative_ranges, relative_ranges_accu):
      a.append(r)
    if (rms_max_list_accu is None):
      rms_max_list_accu = [flex.double() for i in range(len(rms_max_list))]
    else:
      assert len(rms_max_list) == len(rms_max_list_accu)
    for r,a in zip(rms_max_list, rms_max_list_accu):
      a.append(r)
    if (out is sys.stdout):
      f = open("tmp_e_tots_%02d_%02d.xy" % (plot_prefix, i_trial), "w")
      print("@with g0", file=f)
      for i,l in enumerate(sim_labels):
        l = l[l.find('"')+1:].replace('"','')[:-1]
        print('@ s%d legend "%s"' % (i, l), file=f)
      for es in e_tots_list:
        for e in es:
          print(e, file=f)
        print("&", file=f)
      f.close()
  print("Accumulated results:", file=out)
  print(file=out)
  for sim_label,accu in zip(sim_labels, relative_ranges_accu):
    print("relative ranges %s:" % sim_label, file=out)
    accu.min_max_mean().show(out=out, prefix="  ")
    print(file=out)
  for sim_label,accu in zip(sim_labels[1:], rms_max_list_accu):
    print("rms max %s" % sim_labels[0], file=out)
    print("    vs. %s:" % sim_label, file=out)
    accu.min_max_mean().show(out=out, prefix="  ")
    print(file=out)
  if (out is not sys.stdout):
    for accu in relative_ranges_accu:
      assert flex.max(accu) < 1.e-4
    for i,accu in enumerate(rms_max_list_accu):
      assert flex.max(accu) < 1.e-4

def run(args):
  assert len(args) in [0,3]
  if (len(args) == 0):
    n_trials = 3
    n_dynamics_steps = 30
    random_seed = 0
    out = null_out()
  else:
    n_trials = max(1, int(args[0]))
    n_dynamics_steps = max(1, int(args[1]))
    random_seed = int(args[2])
    out = sys.stdout
  show_times_at_exit()
  mersenne_twister = flex.mersenne_twister(seed=0)
  exercise_euler_params_qE_as_euler_angles_xyz_qE(
    mersenne_twister=mersenne_twister)
  exercise_T_as_X(mersenne_twister=mersenne_twister)
  global plot_prefix
  plot_prefix = random_seed
  exercise_simulation(
    out=out, n_trials=n_trials, n_dynamics_steps=n_dynamics_steps,
    random_seed=random_seed)
  print("OK")

if (__name__ == "__main__"):
  run(sys.argv[1:])
