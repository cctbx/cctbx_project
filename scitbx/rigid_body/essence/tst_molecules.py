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

  def energies_and_accelerations_update(O):
    model = featherstone.system_model(bodies=O.bodies)
    O.e_kin = model.e_kin()
    O.e_pot_and_f_ext_update()
    O.qdd = model.FDab(tau=None, f_ext=O.f_ext_bf)

  def e_pot_and_f_ext_update(O):
    O.AJA_update()
    O.JAr_update()
    O.sites_moved_update()
    O.e_pot = O.potential_obj.e_pot(sites_moved=O.sites_moved)
    O.f_ext_bf_update(
      d_e_pot_d_sites=O.potential_obj.d_e_pot_d_sites(
        sites_moved=O.sites_moved))
    O.e_tot = O.e_kin + O.e_pot

  def dynamics_step(O, delta_t):
    for B,qdd in zip(O.bodies, O.qdd):
      B.qd = B.J.time_step_velocity(qd=B.qd, qdd=qdd, delta_t=delta_t)
      B.J = B.J.time_step_position(qd=B.qd, delta_t=delta_t)
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

  def minimization(O, max_iterations=None, callback_after_step=None):
    refinery(
      sim=O,
      max_iterations=max_iterations,
      callback_after_step=callback_after_step)

class refinery(object):

  def __init__(O, sim, max_iterations=None, callback_after_step=None):
    O.sim = sim
    O.callback_after_step = callback_after_step
    O.x = flex.double()
    for B in sim.bodies:
      O.x.extend(flex.double(B.J.get_q()))
    import scitbx.lbfgs
    scitbx.lbfgs.run(
      target_evaluator=O,
      termination_params=scitbx.lbfgs.termination_parameters(
       max_iterations=max_iterations))
    O.sim.energies_and_accelerations_update()

  def unpack_x(O):
    x = O.x
    i = 0
    for B in O.sim.bodies:
      n = B.J.q_size
      B.J = B.J.new_q(q=x[i:i+n])
      i += n
    assert i == x.size()
    O.sim.e_pot_and_f_ext_update()

  def compute_functional_and_gradients(O):
    O.unpack_x()
    f = O.sim.e_pot
    g = flex.double()
    for d in O.sim.d_pot_d_q():
      g.extend(flex.double(d))
    return f, g

class six_dof_body(object):

  def __init__(O, sites):
    O.A = joint_lib.six_dof_alignment(
      center_of_mass=utils.center_of_mass_from_sites(sites=sites))
    O.I = utils.spatial_inertia_from_sites(sites=sites, alignment_T=O.A.T0b)
    #
    qE = matrix.col((1,0,0,0))
    qr = matrix.col((0,0,0))
    O.J = joint_lib.six_dof(qE=qE, qr=qr)
    O.qd = O.J.qd_zero

class revolute_body(object):

  def __init__(O, sites, pivot, normal):
    O.A = joint_lib.revolute_alignment(pivot=pivot, normal=normal)
    O.I = utils.spatial_inertia_from_sites(sites=sites, alignment_T=O.A.T0b)
    #
    O.J = joint_lib.revolute(qE=matrix.col([0]))
    O.qd = O.J.qd_zero

def construct_bodies(sites, cluster_manager):
  result = []
  cm = cluster_manager
  for ic,cluster in enumerate(cm.clusters):
    he = cm.hinge_edges[ic]
    if (he[0] == -1):
      body = six_dof_body(
        sites=[matrix.col(sites[i]) for i in cluster])
      body.parent = -1
    else:
      normal_sites = [matrix.col(sites[i]) for i in he]
      body = revolute_body(
        sites=[matrix.col(sites[i]) for i in cluster],
        pivot=normal_sites[1],
        normal=(normal_sites[1]-normal_sites[0]).normalize())
      body.parent = cm.cluster_indices[he[0]]
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

n_test_simulations = len(tst_tardy_pdb.test_cases)

def get_test_simulation_by_index(i):
  tc = tst_tardy_pdb.test_cases[i]
  tt = tc.tardy_tree_construct()
  cm = tt.cluster_manager
  return simulation(
    labels=tc.labels,
    sites=tc.sites,
    bonds=tc.bonds,
    cluster_manager=cm,
    potential_obj=potential_object(
      sites=tc.sites,
      wells=random_wells(tc.sites),
      restraint_edges=cm.loop_edges+cm.loop_edge_bendings),
    bodies=construct_bodies(sites=tc.sites, cluster_manager=cm))

def run(args):
  assert len(args) in [0,1]
  if (len(args) == 0):
    n_dynamics_steps = 30
    out = null_out()
  else:
    n_dynamics_steps = max(1, int(args[0]))
    out = sys.stdout
  show_times_at_exit()
  random.seed(0)
  for i in xrange(n_test_simulations):
    sim = get_test_simulation_by_index(i=i)
    exercise_dynamics_quick(
      out=out, sim=sim, n_dynamics_steps=n_dynamics_steps)
    exercise_minimization_quick(out=out, sim=sim)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
