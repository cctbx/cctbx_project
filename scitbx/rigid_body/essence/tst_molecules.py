from __future__ import division
from scitbx.rigid_body.essence import featherstone
from scitbx.rigid_body.essence import joint_lib
from scitbx.rigid_body.essence import utils
from scitbx.graph import tst_tardy_pdb
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import null_out, show_times_at_exit
from libtbx.test_utils import approx_equal
import random
import sys

class potential_object(object):

  def __init__(O,
        sites,
        wells,
        restraint_edges,
        restraint_edge_weight=1/0.1**2,
        epsilon=1.e-100):
    O.wells = wells
    O.restraints = []
    for edge in restraint_edges:
      s = [sites[i] for i in edge]
      O.restraints.append((edge, abs(s[0]-s[1]), restraint_edge_weight))
    O.epsilon = epsilon

  def e_pot(O, sites_moved):
    result = 0
    for s, w in zip(sites_moved, O.wells):
      result += (s - w).dot()
    for edge,d_ideal,w in O.restraints:
      s = [sites_moved[i] for i in edge]
      d_model = abs(s[0]-s[1])
      if (d_model < O.epsilon): continue
      delta = d_ideal - d_model
      result += w * delta**2
    return result

  def d_e_pot_d_sites(O, sites_moved):
    result = []
    for s, w in zip(sites_moved, O.wells):
      result.append(2 * (s - w))
    for edge,d_ideal,w in O.restraints:
      s = [sites_moved[i] for i in edge]
      d_model = abs(s[0]-s[1])
      if (d_model < O.epsilon): continue
      delta = d_ideal - d_model
      g0 = -w * 2 * delta / d_model * (s[0] - s[1])
      result[edge[0]] += g0
      result[edge[1]] -= g0
    return result

class simulation(object):

  def __init__(O,
        labels, sites, bonds, cluster_manager, potential_obj, bodies):
    O.labels = labels
    O.sites = sites
    O.bonds = bonds
    O.cluster_manager = cluster_manager
    O.potential_obj = potential_obj
    O.bodies = bodies
    O.degrees_of_freedom = sum([B.J.degrees_of_freedom for B in O.bodies])
    O.flag_positions_as_changed()

  def flag_positions_as_changed(O):
    O.__featherstone_system_model = None
    O.__AJA = None
    O.__JAr = None
    O.__sites_moved = None
    O.__e_pot = None
    O.__d_e_pot_d_sites = None
    O.__f_ext_bf = None
    O.flag_velocities_as_changed()

  def flag_velocities_as_changed(O):
    O.__qdd = None
    O.__e_kin = None

  def featherstone_system_model(O):
    if (O.__featherstone_system_model is None):
      O.__featherstone_system_model = featherstone.system_model(
        bodies=O.bodies)
    return O.__featherstone_system_model

  def AJA(O):
    if (O.__AJA is None):
      O.__AJA = []
      for B in O.bodies:
        AJA = B.A.Tb0 * B.J.Tsp * B.A.T0b
        if (B.parent != -1):
          AJA = O.__AJA[B.parent] * AJA
        O.__AJA.append(AJA)
    return O.__AJA

  def JAr(O):
    if (O.__JAr is None):
      O_AJA = O.AJA()
      O.__JAr = []
      for B in O.bodies:
        JAr = B.J.Tps.r * B.A.T0b.r
        if (B.parent != -1):
          JAr *= O_AJA[B.parent].r.transpose()
        O.__JAr.append(JAr)
    return O.__JAr

  def sites_moved(O):
    if (O.__sites_moved is None):
      O_AJA = O.AJA()
      O.__sites_moved = [None] * len(O.sites)
      n_done = 0
      for iB,B in enumerate(O.bodies):
        AJA = O_AJA[iB]
        for i_seq in O.cluster_manager.clusters[iB]:
          assert O.__sites_moved[i_seq] is None
          O.__sites_moved[i_seq] = AJA * O.sites[i_seq]
          n_done += 1
      assert n_done == len(O.sites)
    return O.__sites_moved

  def e_pot(O):
    if (O.__e_pot is None):
      O.__e_pot = O.potential_obj.e_pot(
        sites_moved=O.sites_moved())
    return O.__e_pot

  def d_e_pot_d_sites(O):
    if (O.__d_e_pot_d_sites is None):
      O.__d_e_pot_d_sites = O.potential_obj.d_e_pot_d_sites(
        sites_moved=O.sites_moved())
    return O.__d_e_pot_d_sites

  def f_ext_bf(O):
    O_JAr = O.JAr()
    O_d_e_pot_d_sites = O.d_e_pot_d_sites()
    O.__f_ext_bf = []
    for iB,B in enumerate(O.bodies):
      f = matrix.col((0,0,0))
      nc = matrix.col((0,0,0))
      for i_seq in O.cluster_manager.clusters[iB]:
        s = O.sites[i_seq]
        force_bf = -(O_JAr[iB] * O_d_e_pot_d_sites[i_seq])
        f += force_bf
        nc += (B.A.T0b * s).cross(force_bf)
      O.__f_ext_bf.append(matrix.col((nc, f)).resolve_partitions())
    return O.__f_ext_bf

  def qdd(O):
    if (O.__qdd is None):
      O.__qdd = O.featherstone_system_model().FDab(
        tau=None, f_ext=O.f_ext_bf())
    return O.__qdd

  def e_kin(O):
    if (O.__e_kin is None):
      O.__e_kin = O.featherstone_system_model().e_kin()
    return O.__e_kin

  def e_tot(O):
    return O.e_kin() + O.e_pot()

  def reset_e_kin(O, e_kin_target, e_kin_epsilon=1.e-12):
    assert e_kin_target >= 0
    O_e_kin = O.e_kin()
    if (O_e_kin >= e_kin_epsilon):
      factor = (e_kin_target / O_e_kin)**0.5
      for B in O.bodies:
        B.qd *= factor
    O.flag_velocities_as_changed()

  def assign_zero_velocities(O):
    for B in O.bodies:
      B.qd = B.J.qd_zero
    O.flag_velocities_as_changed()

  def assign_random_velocities(O,
        e_kin_target=None,
        e_kin_epsilon=1.e-12,
        random_gauss=None):
    work_e_kin_target = e_kin_target
    if (e_kin_target is None):
      work_e_kin_target = 1
    elif (e_kin_target == 0):
      O.assign_zero_velocities()
      return
    else:
      assert e_kin_target >= 0
    qd_e_kin_scales = flex.double(
      O.featherstone_system_model().qd_e_kin_scales(
        e_kin_epsilon=e_kin_epsilon))
    if (O.degrees_of_freedom != 0):
      qd_e_kin_scales *= (work_e_kin_target / O.degrees_of_freedom)**0.5
    if (random_gauss is None):
      random_gauss = random.gauss
    i_qd = 0
    for B in O.bodies:
      qd_new = []
      for qd in B.J.qd_zero:
        qd_new.append(qd + random_gauss(mu=0, sigma=qd_e_kin_scales[i_qd]))
        i_qd += 1
      B.qd = matrix.col(qd_new)
    assert i_qd == O.degrees_of_freedom
    O.flag_velocities_as_changed()
    if (e_kin_target is not None):
      O.reset_e_kin(e_kin_target=e_kin_target, e_kin_epsilon=e_kin_epsilon)
    return qd_e_kin_scales

  def dynamics_step(O, delta_t):
    O_qdd = O.qdd()
    for B in O.bodies:
      B.J = B.J.time_step_position(qd=B.qd, delta_t=delta_t)
    for B,qdd in zip(O.bodies, O_qdd):
      B.qd = B.J.time_step_velocity(qd=B.qd, qdd=qdd, delta_t=delta_t)
    O.flag_positions_as_changed()

  def d_pot_d_q(O):
    return O.featherstone_system_model().d_pot_d_q(f_ext=O.f_ext_bf())

  def d_pot_d_q_via_finite_differences(O, eps=1.e-6):
    result = []
    for B in O.bodies:
      gs = []
      J_orig = B.J
      q_orig = list(J_orig.get_q())
      for iq in xrange(J_orig.q_size):
        fs = []
        for signed_eps in [eps, -eps]:
          q_eps = list(q_orig)
          q_eps[iq] += signed_eps
          B.J = J_orig.new_q(q=q_eps)
          O.flag_positions_as_changed()
          fs.append(O.e_pot())
        gs.append((fs[0]-fs[1])/(2*eps))
      B.J = J_orig
      O.flag_positions_as_changed()
      result.append(matrix.col(gs))
    return result

  def check_d_pot_d_q(O, verbose=0):
    qdd_orig = O.qdd()
    ana = O.d_pot_d_q()
    fin = O.d_pot_d_q_via_finite_differences()
    if (verbose):
      for a,f in zip(ana, fin):
        print "fin:", f.elems
        print "ana:", a.elems
      print
    assert approx_equal(ana, fin)
    assert approx_equal(O.qdd(), qdd_orig)

  def pack_q(O):
    result = flex.double()
    for B in O.bodies:
      result.extend(flex.double(B.J.get_q()))
    return result

  def unpack_q(O, packed_q):
    i = 0
    for B in O.bodies:
      n = B.J.q_size
      B.J = B.J.new_q(q=packed_q[i:i+n])
      i += n
    assert i == packed_q.size()
    O.flag_positions_as_changed()

  def minimization(O, max_iterations=None, callback_after_step=None):
    refinery(
      sim=O,
      max_iterations=max_iterations,
      callback_after_step=callback_after_step)

class refinery(object):

  def __init__(O, sim, max_iterations=None, callback_after_step=None):
    O.sim = sim
    O.callback_after_step = callback_after_step
    O.x = sim.pack_q()
    import scitbx.lbfgs
    scitbx.lbfgs.run(
      target_evaluator=O,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound=True))

  def compute_functional_and_gradients(O):
    O.sim.unpack_q(packed_q=O.x)
    f = O.sim.e_pot()
    g = flex.double()
    for d in O.sim.d_pot_d_q():
      g.extend(flex.double(d))
    return f, g

class six_dof_body(object):

  def __init__(O, sites, masses):
    mass_points = utils.mass_points(sites=sites, masses=masses)
    O.A = joint_lib.six_dof_alignment(
      center_of_mass=mass_points.center_of_mass())
    O.I = mass_points.spatial_inertia(alignment_T=O.A.T0b)
    #
    qE = matrix.col((1,0,0,0))
    qr = matrix.col((0,0,0))
    O.J = joint_lib.six_dof(qE=qE, qr=qr)
    O.qd = O.J.qd_zero

class revolute_body(object):

  def __init__(O, sites, masses, pivot, normal):
    mass_points = utils.mass_points(sites=sites, masses=masses)
    O.A = joint_lib.revolute_alignment(pivot=pivot, normal=normal)
    O.I = mass_points.spatial_inertia(alignment_T=O.A.T0b)
    #
    O.J = joint_lib.revolute(qE=matrix.col([0]))
    O.qd = O.J.qd_zero

class translational_body(object):

  def __init__(O, sites, masses):
    mass_points = utils.mass_points(sites=sites, masses=masses)
    O.A = joint_lib.translational_alignment(
      center_of_mass=mass_points.center_of_mass())
    O.I = mass_points.spatial_inertia(alignment_T=O.A.T0b)
    #
    qr = matrix.col((0,0,0))
    O.J = joint_lib.translational(qr=qr)
    O.qd = O.J.qd_zero

def construct_bodies(sites, masses, cluster_manager):
  assert len(sites) == len(masses)
  result = []
  cm = cluster_manager
  for ic,cluster in enumerate(cm.clusters):
    body_sites = [matrix.col(sites[i]) for i in cluster]
    body_masses = [masses[i] for i in cluster]
    he = cm.hinge_edges[ic]
    if (he[0] == -1):
      if (len(sites) == 1):
        body = translational_body(sites=body_sites, masses=body_masses)
      else:
        body = six_dof_body(sites=body_sites, masses=body_masses)
      body.parent = -1
    else:
      normal_sites = [matrix.col(sites[i]) for i in he]
      body = revolute_body(
        sites=body_sites,
        masses=body_masses,
        pivot=normal_sites[1],
        normal=(normal_sites[1]-normal_sites[0]).normalize())
      body.parent = cm.cluster_indices[he[1]]
    body.i_seqs = cluster
    result.append(body)
  return result

def exercise_qd_e_kin_scales(sim):
  def slow():
    result = flex.double()
    for B in sim.bodies:
      BJ0 = B.J
      qd0 = B.J.qd_zero
      qd = list(qd0)
      for iqd in xrange(len(qd)):
        qd[iqd] = qd0[iqd] + 1
        B.qd = matrix.col(qd)
        qd[iqd] = qd0[iqd]
        sim.flag_velocities_as_changed()
        B.J = BJ0.time_step_position(qd=B.qd, delta_t=1)
        e_kin = sim.e_kin()
        if (e_kin < 1.e-12):
          result.append(1)
        else:
          result.append(1 / e_kin**0.5)
      B.J = BJ0
      B.qd = B.J.qd_zero
      sim.flag_positions_as_changed()
    assert sim.e_kin() == 0
    assert len(result) == sim.degrees_of_freedom
    return result
  scales_slow = slow()
  model = sim.featherstone_system_model()
  scales_fast = model.qd_e_kin_scales()
  assert approx_equal(scales_fast, scales_slow)

def exercise_random_velocities(sim):
  for e_kin_target in [1, 1.3, 13, 0]:
    sim.assign_random_velocities(e_kin_target=e_kin_target)
    assert approx_equal(sim.e_kin(), e_kin_target)

def exercise_sim(out, n_dynamics_steps, delta_t, sim):
  sim.check_d_pot_d_q()
  e_pots = flex.double([sim.e_pot()])
  e_kins = flex.double([sim.e_kin()])
  for i_step in xrange(n_dynamics_steps):
    sim.dynamics_step(delta_t=delta_t)
    e_pots.append(sim.e_pot())
    e_kins.append(sim.e_kin())
  e_tots = e_pots + e_kins
  sim.check_d_pot_d_q()
  print >> out, "degrees of freedom:", sim.degrees_of_freedom
  print >> out, "energy samples:", e_tots.size()
  print >> out, "e_pot min, max:", min(e_pots), max(e_pots)
  print >> out, "e_kin min, max:", min(e_kins), max(e_kins)
  print >> out, "e_tot min, max:", min(e_tots), max(e_tots)
  print >> out, "start e_tot:", e_tots[0]
  print >> out, "final e_tot:", e_tots[-1]
  ave = flex.sum(e_tots) / e_tots.size()
  range = flex.max(e_tots) - flex.min(e_tots)
  if (ave == 0): relative_range = 0
  else:          relative_range = range / ave
  print >> out, "ave:", ave
  print >> out, "range:", range
  print >> out, "relative range:", relative_range
  print >> out
  out.flush()
  return relative_range

def exercise_dynamics_quick(out, sim, n_dynamics_steps, delta_t=0.0001):
  relative_range = exercise_sim(
    out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
  if (out is not sys.stdout):
    assert relative_range < 1.e-4
  print >> out

def exercise_minimization_quick(out, sim, max_iterations=3):
  print >> out, "Minimization:"
  print >> out, "  start e_pot:", sim.e_pot()
  e_pot_start = sim.e_pot()
  sim.minimization(max_iterations=max_iterations)
  print >> out, "  final e_pot:", sim.e_pot()
  e_pot_final = sim.e_pot()
  if (out is not sys.stdout):
    assert e_pot_final < e_pot_start * 0.5
  print >> out

def construct_simulation(
      labels,
      sites,
      masses,
      tardy_tree):
  cm = tardy_tree.cluster_manager
  return simulation(
    labels=labels,
    sites=sites,
    bonds=tardy_tree.edge_list,
    cluster_manager=cm,
    potential_obj=potential_object(
      sites=sites,
      wells=sites,
      restraint_edges=cm.loop_edges+cm.loop_edge_bendings),
    bodies=construct_bodies(sites=sites, masses=masses, cluster_manager=cm))

n_test_simulations = len(tst_tardy_pdb.test_cases)

def get_test_simulation_by_index(i):
  tc = tst_tardy_pdb.test_cases[i]
  tt = tc.tardy_tree_construct()
  return construct_simulation(
    labels=tc.labels,
    sites=tc.sites,
    masses=[1.0]*len(tc.sites),
    tardy_tree=tt)

def run(args):
  assert len(args) in [0,1]
  if (len(args) == 0):
    n_dynamics_steps = 100
    out = null_out()
  else:
    n_dynamics_steps = max(1, int(args[0]))
    out = sys.stdout
  show_times_at_exit()
  #
  if (1):
    random.seed(0)
    for i in xrange(n_test_simulations):
      print >> out, "test_simulation index:", i
      sim = get_test_simulation_by_index(i=i)
      exercise_qd_e_kin_scales(sim=sim)
      exercise_random_velocities(sim=sim)
      sim.assign_random_velocities(e_kin_target=1)
      assert approx_equal(sim.e_kin(), 1)
      exercise_dynamics_quick(
        out=out, sim=sim, n_dynamics_steps=n_dynamics_steps)
      exercise_minimization_quick(out=out, sim=sim)
  #
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
