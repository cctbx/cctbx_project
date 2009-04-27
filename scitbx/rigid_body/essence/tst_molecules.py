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

def random_wells(sites):
  result = []
  for site in sites:
    result.append(site+matrix.col.random(n=3, a=-1, b=1))
  return result

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

  def e_pot_and_normalization_factor(O, sites_moved):
    result = 0
    for s, w in zip(sites_moved, O.wells):
      result += (s - w).dot()
    for edge,d_ideal,w in O.restraints:
      s = [sites_moved[i] for i in edge]
      d_model = abs(s[0]-s[1])
      if (d_model < O.epsilon): continue
      delta = d_ideal - d_model
      result += w * delta**2
    return result, 1.0

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
    O.energies_and_accelerations_update()

  def AJA_update(O):
    O.AJA = []
    for B in O.bodies:
      AJA = B.A.Tb0 * B.J.Tsp * B.A.T0b
      if (B.parent != -1):
        AJA = O.AJA[B.parent] * AJA
      O.AJA.append(AJA)

  def JAr_update(O):
    O.JAr = []
    for B in O.bodies:
      JAr = B.J.Tps.r * B.A.T0b.r
      if (B.parent != -1):
        JAr *= O.AJA[B.parent].r.transpose()
      O.JAr.append(JAr)

  def sites_moved_update(O):
    O.sites_moved = [None] * len(O.sites)
    n_done = 0
    for iB,B in enumerate(O.bodies):
      AJA = O.AJA[iB]
      for i_seq in O.cluster_manager.clusters[iB]:
        assert O.sites_moved[i_seq] is None
        O.sites_moved[i_seq] = AJA * O.sites[i_seq]
        n_done += 1
    assert n_done == len(O.sites)

  def f_ext_bf_update(O, d_e_pot_d_sites):
    O.f_ext_bf = []
    for iB,B in enumerate(O.bodies):
      f = matrix.col((0,0,0))
      nc = matrix.col((0,0,0))
      for i_seq in O.cluster_manager.clusters[iB]:
        s = O.sites[i_seq]
        force_bf = -(O.JAr[iB] * d_e_pot_d_sites[i_seq])
        f += force_bf
        nc += (B.A.T0b * s).cross(force_bf)
      O.f_ext_bf.append(matrix.col((nc, f)).resolve_partitions())

  def e_kin_update(O):
    model = featherstone.system_model(bodies=O.bodies)
    O.e_kin = model.e_kin()
    return model

  def energies_and_accelerations_update(O):
    model = O.e_kin_update()
    O.e_pot_and_f_ext_update()
    O.qdd = model.FDab(tau=None, f_ext=O.f_ext_bf)

  def e_pot_update(O):
    O.AJA_update()
    O.JAr_update()
    O.sites_moved_update()
    O.e_pot, O.e_pot_normalization_factor = \
      O.potential_obj.e_pot_and_normalization_factor(
        sites_moved=O.sites_moved)

  def e_pot_and_f_ext_update(O):
    O.e_pot_update()
    O.f_ext_bf_update(
      d_e_pot_d_sites=O.potential_obj.d_e_pot_d_sites(
        sites_moved=O.sites_moved))
    O.e_tot = O.e_kin + O.e_pot

  def dynamics_step(O, delta_t, e_kin_cap=None):
    delta_t_norm = delta_t / O.e_pot_normalization_factor
    for B,qdd in zip(O.bodies, O.qdd):
      B.qd = B.J.time_step_velocity(qd=B.qd, qdd=qdd, delta_t=delta_t_norm)
    O.apply_velocity_scaling(e_kin_cap=e_kin_cap)
    for B,qdd in zip(O.bodies, O.qdd):
      B.J = B.J.time_step_position(qd=B.qd, delta_t=delta_t)
    O.energies_and_accelerations_update()

  def apply_velocity_scaling(O, e_kin_cap):
    if (e_kin_cap is None):
      O.e_kin_before_velocity_scaling = None
      return
    assert e_kin_cap > 0
    O.e_kin_before_velocity_scaling = featherstone.system_model(
      bodies=O.bodies).e_kin()
    if (O.e_kin_before_velocity_scaling > e_kin_cap):
      factor = (e_kin_cap / O.e_kin_before_velocity_scaling)**0.5
      for B in O.bodies:
        B.qd *= factor

  def assign_random_velocities(O, e_kin_target):
    assert O.degrees_of_freedom != 0
    qd_e_kin_scales = flex.double()
    for B in O.bodies:
      BJ0 = B.J
      qd0 = B.J.qd_zero
      qd = list(qd0)
      for iqd in xrange(len(qd)):
        qd[iqd] = qd0[iqd] + 1
        B.qd = matrix.col(qd)
        qd[iqd] = qd0[iqd]
        B.J = BJ0.time_step_position(qd=B.qd, delta_t=1)
        O.e_kin_update()
        assert O.e_kin != 0
        qd_e_kin_scale = 1 / O.e_kin**0.5
        qd_e_kin_scales.append(qd_e_kin_scale)
      B.J = BJ0
      B.qd = B.J.qd_zero
    O.e_kin_update()
    assert O.e_kin == 0
    assert qd_e_kin_scales.size() == O.degrees_of_freedom
    qd_e_kin_scales *= (e_kin_target / O.degrees_of_freedom)**0.5
    rg = random.gauss
    i_qd = 0
    for B in O.bodies:
      qd_new = []
      for qd in B.J.qd_zero:
        qd_new.append(qd + rg(mu=0, sigma=qd_e_kin_scales[i_qd]))
        i_qd += 1
      B.qd = matrix.col(qd_new)
    assert i_qd == O.degrees_of_freedom
    O.e_kin_update()
    assert O.e_kin != 0
    factor = (e_kin_target / O.e_kin)**0.5
    for B in O.bodies:
      B.qd *= factor
    O.energies_and_accelerations_update()

  def d_pot_d_q(O):
    return featherstone.system_model(bodies=O.bodies).d_pot_d_q(
      f_ext=O.f_ext_bf)

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
          O.e_pot_and_f_ext_update()
          fs.append(O.e_pot)
        gs.append((fs[0]-fs[1])/(2*eps))
      B.J = J_orig
      result.append(matrix.col(gs))
    O.energies_and_accelerations_update()
    return result

  def check_d_pot_d_q(O, verbose=0):
    qdd_orig = O.qdd
    ana = O.d_pot_d_q()
    fin = O.d_pot_d_q_via_finite_differences()
    if (verbose):
      for a,f in zip(ana, fin):
        print "fin:", f.elems
        print "ana:", a.elems
      print
    assert approx_equal(ana, fin)
    assert approx_equal(O.qdd, qdd_orig)

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
    O.sim.energies_and_accelerations_update()

  def compute_functional_and_gradients(O):
    O.sim.unpack_q(packed_q=O.x)
    O.sim.e_pot_and_f_ext_update()
    f = O.sim.e_pot
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

def exercise_sim(out, n_dynamics_steps, delta_t, sim):
  sim.check_d_pot_d_q()
  e_pots = flex.double([sim.e_pot])
  e_kins = flex.double([sim.e_kin])
  for i_step in xrange(n_dynamics_steps):
    sim.dynamics_step(delta_t=delta_t)
    e_pots.append(sim.e_pot)
    e_kins.append(sim.e_kin)
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

def exercise_dynamics_quick(out, sim, n_dynamics_steps, delta_t=0.001):
  relative_range = exercise_sim(
    out=out, n_dynamics_steps=n_dynamics_steps, delta_t=delta_t, sim=sim)
  if (out is not sys.stdout):
    assert relative_range < 1.e-4
  print >> out

def exercise_minimization_quick(out, sim, max_iterations=3):
  print >> out, "Minimization:"
  print >> out, "  start e_pot:", sim.e_pot
  e_pot_start = sim.e_pot
  sim.minimization(max_iterations=max_iterations)
  print >> out, "  final e_pot:", sim.e_pot
  e_pot_final = sim.e_pot
  if (out is not sys.stdout):
    assert e_pot_final < e_pot_start * 0.98
  print >> out

def exercise_apply_velocity_scaling(out, sim):
  print >> out, "exercise_apply_velocity_scaling():"
  e_kins = flex.double([sim.e_kin])
  for i_step in xrange(10):
    sim.dynamics_step(delta_t=0.1, e_kin_cap=2)
    e_kins.append(sim.e_kin)
    print >> out, "e_kin:", sim.e_kin
  assert approx_equal(e_kins, [
    0.0, 0.10596967199225664, 0.41471669284176654, 0.90125773606867454,
    1.5269900929646234, 2.0446638368687799, 2.0528390856843681,
    2.0592456511353947, 2.0643773120065561, 2.0685121042338253,
    2.0718756155399802])

def construct_simulation(
      labels,
      sites,
      masses,
      tardy_tree,
      use_random_wells=True):
  if (use_random_wells):
    wells = random_wells(sites)
  else:
    wells = sites
  cm = tardy_tree.cluster_manager
  return simulation(
    labels=labels,
    sites=sites,
    bonds=tardy_tree.edge_list,
    cluster_manager=cm,
    potential_obj=potential_object(
      sites=sites,
      wells=wells,
      restraint_edges=cm.loop_edges+cm.loop_edge_bendings),
    bodies=construct_bodies(sites=sites, masses=masses, cluster_manager=cm))

n_test_simulations = len(tst_tardy_pdb.test_cases)

def get_test_simulation_by_index(i, use_random_wells=True):
  tc = tst_tardy_pdb.test_cases[i]
  tt = tc.tardy_tree_construct()
  return construct_simulation(
    labels=tc.labels,
    sites=tc.sites,
    masses=[1.0]*len(tc.sites),
    tardy_tree=tt,
    use_random_wells=use_random_wells)

def run(args):
  assert len(args) in [0,1]
  if (len(args) == 0):
    n_dynamics_steps = 30
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
      exercise_dynamics_quick(
        out=out, sim=sim, n_dynamics_steps=n_dynamics_steps)
      exercise_minimization_quick(out=out, sim=sim)
  #
  if (1):
    random.seed(0)
    sim = get_test_simulation_by_index(i=13)
    assert sim.degrees_of_freedom == 12
    exercise_apply_velocity_scaling(out=out, sim=sim)
  #
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
