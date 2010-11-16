from __future__ import division
from scitbx.rigid_body.essence import tardy
from scitbx.graph import tardy_tree
from scitbx.graph import test_cases_tardy_pdb
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import null_out, show_times_at_exit
from libtbx.test_utils import approx_equal, is_above_limit
import random
import sys

if (1):
  random.seed(0)
  flex.set_random_seed(0)

class tardy_model(tardy.model):

  def d_e_pot_d_q_via_finite_differences(O, eps=1.e-6):
    result = []
    for body in O.bodies:
      gs = []
      j_orig = body.joint
      q_orig = list(j_orig.get_q())
      for iq in xrange(j_orig.q_size):
        fs = []
        for signed_eps in [eps, -eps]:
          q_eps = list(q_orig)
          q_eps[iq] += signed_eps
          body.joint = j_orig.new_q(q=q_eps)
          O.flag_positions_as_changed()
          fs.append(O.e_pot())
        gs.append((fs[0]-fs[1])/(2*eps))
      body.joint = j_orig
      O.flag_positions_as_changed()
      result.append(matrix.col(gs))
    return result

  def check_d_e_pot_d_q(O, verbose=0):
    qdd_orig = O.qdd_array()
    ana = O.d_e_pot_d_q()
    fin = O.d_e_pot_d_q_via_finite_differences()
    if (verbose):
      for a,f in zip(ana, fin):
        print "fin:", f.elems
        print "ana:", a.elems
      print
    assert approx_equal(ana, fin)
    assert approx_equal(O.qdd_array(), qdd_orig)

class potential_object(object):

  def __init__(O,
        sites,
        wells,
        restraint_edges,
        restraint_edge_weight=1/0.1**2,
        epsilon=1.e-100):
    if (isinstance(wells, flex.vec3_double)):
      wells = [matrix.col(w) for w in wells]
    O.wells = wells
    O.restraints = []
    for edge in restraint_edges:
      s = [matrix.col(sites[i]) for i in edge]
      O.restraints.append((edge, abs(s[0]-s[1]), restraint_edge_weight))
    O.epsilon = epsilon

  def e_pot(O, sites_moved):
    if (isinstance(sites_moved, flex.vec3_double)):
      sites_moved = [matrix.col(s) for s in sites_moved]
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
    if (isinstance(sites_moved, flex.vec3_double)):
      sites_moved = [matrix.col(s) for s in sites_moved]
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

def exercise_qd_e_kin_scales(tardy_model):
  def slow():
    result = flex.double()
    for body in tardy_model.bodies:
      bj0 = body.joint
      qd0 = body.joint.qd_zero
      qd = list(qd0)
      for iqd in xrange(len(qd)):
        qd[iqd] = qd0[iqd] + 1
        body.qd = matrix.col(qd)
        qd[iqd] = qd0[iqd]
        tardy_model.flag_velocities_as_changed()
        body.joint = bj0.time_step_position(qd=body.qd, delta_t=1)
        e_kin = tardy_model.e_kin()
        if (e_kin < 1.e-12):
          result.append(1)
        else:
          result.append(1 / e_kin**0.5)
      body.joint = bj0
      body.qd = body.joint.qd_zero
      tardy_model.flag_positions_as_changed()
    assert tardy_model.e_kin() == 0
    assert len(result) == tardy_model.degrees_of_freedom
    return result
  scales_slow = slow()
  scales_fast = tardy_model.qd_e_kin_scales()
  assert approx_equal(scales_fast, scales_slow)

def exercise_random_velocities(tardy_model):
  prev_qd = tardy_model.pack_qd()
  for e_kin_target in [1, 1.3, 0, 13]:
    tardy_model.assign_random_velocities(e_kin_target=e_kin_target)
    assert approx_equal(tardy_model.e_kin(), e_kin_target)
  assert not approx_equal(tardy_model.pack_qd(), prev_qd, out=None)
  tardy_model.unpack_qd(qd_packed=prev_qd)
  assert approx_equal(tardy_model.pack_qd(), prev_qd)

def exercise_accumulate_in_each_tree():
  def check(n_vertices, edge_list, clusters, nosiet, somiet):
    sites = matrix.col_list([(i,i%2,0) for i in xrange(n_vertices)])
    labels = [str(i) for i in xrange(n_vertices)]
    masses = [13, 7, 23, 19, 29, 11, 17][:n_vertices]
    assert len(masses) == n_vertices
    tt = tardy_tree.construct(sites=sites, edge_list=edge_list)
    assert len(tt.cluster_manager.fixed_hinges) == 0
    assert tt.cluster_manager.clusters == clusters
    tm = construct_tardy_model(
      labels=labels, sites=sites, masses=masses, tardy_tree=tt)
    assert tm.root_indices() == list(reversed([i for i,n in nosiet]))
    assert tm.number_of_sites_in_each_tree() == nosiet
    assert tm.sum_of_masses_in_each_tree() == somiet
  #
  check(
    n_vertices=5,
    edge_list=[],
    clusters=[[0],[1],[2],[3],[4]],
    nosiet=[(4, 1), (3, 1), (2, 1), (1, 1), (0, 1)],
    somiet=[(4, 29), (3, 19), (2, 23), (1, 7), (0, 13)])
  check(
    n_vertices=5,
    edge_list=[(0,1), (1,2), (3,4)],
    clusters=[[0,1,2], [3,4]],
    nosiet=[(1, 2), (0, 3)],
    somiet=[(1, 48), (0, 43)])
  check(
    n_vertices=5,
    edge_list=[(0,1), (1,2), (2,3), (3,4)],
    clusters=[[0,1,2], [3], [4]],
    nosiet=[(0, 5)],
    somiet=[(0, 91)])
  check(
    n_vertices=7,
    edge_list=[(0,1), (1,2), (2,3), (2,5), (3,4), (5,6)],
    clusters=[[1,2,3,5], [0], [4], [6]],
    nosiet=[(0, 7)],
    somiet=[(0, 119)])

def exercise_near_singular_hinges():
  # similar to scitbx/graph/test_cases_tardy_pdb.py, "collinear" test case
  """
    0                6
    | \             /
    |  2---3---4---5
    | /
    1
  """
  x = -0.5*3**0.5
  y = 0.5
  def build_sites(eps):
    return matrix.col_list([
      (x,-y,0),
      (x,y,0),
      (0,0,0),
      (1,0,0),
      (2,0,eps)])
  edge_list = [(0,1),(0,2),(1,2),(2,3),(3,4)]
  sites = build_sites(eps=0)
  for i,j in edge_list:
    assert approx_equal(abs(sites[i]-sites[j]), 1)
  sites = build_sites(eps=1e-5)
  labels = [str(i) for i in xrange(len(sites))]
  masses = [1] * len(sites)
  tt = tardy_tree.construct(n_vertices=len(sites), edge_list=edge_list)
  tt.build_tree()
  assert tt.cluster_manager.clusters == [[0,1,2,3], [4]]
  def e_kin_1():
    tm = construct_tardy_model(
      labels=labels,
      sites=sites,
      masses=masses,
      tardy_tree=tt)
    rnd = random.Random(0)
    tm.assign_random_velocities(e_kin_target=1, random_gauss=rnd.gauss)
    assert approx_equal(tm.e_kin(), 1, eps=1e-10)
    tm.dynamics_step(delta_t=0.01)
    return tm.e_kin()
  assert approx_equal(e_kin_1(), 60.9875715394)
  tt.fix_near_singular_hinges(sites=sites)
  assert tt.cluster_manager.fixed_hinges == [(2,3)]
  assert tt.cluster_manager.clusters == [[0,1,2,3,4]]
  assert approx_equal(e_kin_1(), 1.00004830172, eps=1e-10)
  #
  sites.append(matrix.col((3,0,0)))
  labels.append("5")
  masses.append(1)
  edge_list.append((4,5))
  tt = tardy_tree.construct(n_vertices=len(sites), edge_list=edge_list)
  tt.build_tree()
  assert tt.cluster_manager.clusters == [[0,1,2,3], [4], [5]]
  assert approx_equal(e_kin_1(), 9.55508653428)
  tt.fix_near_singular_hinges(sites=sites)
  assert tt.cluster_manager.fixed_hinges == [(2,3), (3,4)]
  assert tt.cluster_manager.clusters == [[0,1,2,3,4,5]]
  assert approx_equal(e_kin_1(), 1.00005333167, eps=1e-10)
  #
  sites.append(matrix.col((3-x,-y,0)))
  assert approx_equal(abs(sites[5] - sites[6]), 1)
  labels.append("6")
  masses.append(1)
  edge_list.append((5,6))
  tt = tardy_tree.construct(n_vertices=len(sites), edge_list=edge_list)
  tt.build_tree()
  assert tt.cluster_manager.clusters == [[0,1,2,3], [4], [5], [6]]
  assert approx_equal(e_kin_1(), 0.99994375467)
  tt.fix_near_singular_hinges(sites=sites)
  assert tt.cluster_manager.fixed_hinges == [(2,3), (3,4)]
  assert tt.cluster_manager.clusters == [[0,1,2,3,4,5], [6]]
  assert approx_equal(e_kin_1(), 1.0000438095, eps=1e-10)

def exercise_linear_velocity_manipulations(tardy_model):
  for nosiet in [None, tardy_model.number_of_sites_in_each_tree()]:
    tardy_model.assign_random_velocities(e_kin_target=17)
    if (tardy_model.degrees_of_freedom == 0):
      assert approx_equal(tardy_model.e_kin(), 0)
      assert tardy_model.mean_linear_velocity(
        number_of_sites_in_each_tree=None) is None
    else:
      assert approx_equal(tardy_model.e_kin(), 17)
      mlv = tardy_model.mean_linear_velocity(
        number_of_sites_in_each_tree=nosiet)
      if (mlv is not None):
        assert is_above_limit(value=abs(mlv), limit=1e-3)
        tardy_model.subtract_from_linear_velocities(
          number_of_sites_in_each_tree=nosiet, value=mlv)
        mlv = tardy_model.mean_linear_velocity(
          number_of_sites_in_each_tree=nosiet)
        assert approx_equal(mlv, (0,0,0))

def exercise_fixed_vertices_special_cases():
  tardy_models = []
  """
          2
         /
    0---1
  """
  x = 0.5*3**0.5
  y = 0.5
  sites = matrix.col_list([
    (0,0,0),
    (1,0,0),
    (1+x,y,0)])
  edge_list = [(0,1),(1,2)]
  for i,j in edge_list:
    assert approx_equal(abs(sites[i]-sites[j]), 1)
  labels = [str(i) for i in xrange(len(sites))]
  masses = [1] * len(sites)
  #
  tt = tardy_tree.construct(
    sites=sites,
    edge_list=edge_list,
    fixed_vertex_lists=[])
  assert tt.cluster_manager.clusters == [[0,1,2]]
  tm = construct_tardy_model(
    labels=labels, sites=sites, masses=masses, tardy_tree=tt)
  tardy_models.append(tm)
  assert len(tm.bodies) == 1
  assert tm.bodies[0].joint.degrees_of_freedom == 6
  exercise_linear_velocity_manipulations(tardy_model=tm)
  #
  expected_e_kin_1 = [
    1.00009768395,
    1.00002522865,
    1.00000107257,
    1.0,
    1.0,
    1.0,
    0.0]
  rnd = random.Random(0)
  for i_fv,fixed_vertices in enumerate([[0], [1], [2],
                                        [0,1], [0,2], [1,2],
                                        [0,1,2]]):
    tt = tardy_tree.construct(
      sites=sites,
      edge_list=edge_list,
      fixed_vertex_lists=[fixed_vertices])
    assert tt.cluster_manager.clusters == [[0,1,2]]
    tm = construct_tardy_model(
      labels=labels, sites=sites, masses=masses, tardy_tree=tt)
    tardy_models.append(tm)
    assert len(tm.bodies) == 1
    assert tm.bodies[0].joint.degrees_of_freedom \
        == [3,1,0][len(fixed_vertices)-1]
    tm.assign_random_velocities(e_kin_target=1, random_gauss=rnd.gauss)
    if (len(fixed_vertices) != 3):
      assert approx_equal(tm.e_kin(), 1, eps=1e-10)
    else:
      assert approx_equal(tm.e_kin(), 0, eps=1e-10)
    tm.dynamics_step(delta_t=0.01)
    assert approx_equal(tm.e_kin(), expected_e_kin_1[i_fv], eps=1e-10)
    exercise_linear_velocity_manipulations(tardy_model=tm)
  #
  sites[2] = matrix.col([2,0,0])
  assert approx_equal((sites[0]-sites[1]).cos_angle(sites[2]-sites[1]), -1)
  for fixed_vertices in [[0,1], [0,2], [1,2]]:
    tt = tardy_tree.construct(
      sites=sites,
      edge_list=edge_list,
      fixed_vertex_lists=[fixed_vertices])
    assert tt.cluster_manager.clusters == [[0,1,2]]
    tm = construct_tardy_model(
      labels=labels, sites=sites, masses=masses, tardy_tree=tt)
    tardy_models.append(tm)
    assert len(tm.bodies) == 1
    assert tm.bodies[0].joint.degrees_of_freedom == 0
    exercise_linear_velocity_manipulations(tardy_model=tm)
  return tardy_models

def exercise_tardy_model(out, n_dynamics_steps, delta_t, tardy_model):
  tardy_model.check_d_e_pot_d_q()
  e_pots = flex.double([tardy_model.e_pot()])
  e_kins = flex.double([tardy_model.e_kin()])
  for i_step in xrange(n_dynamics_steps):
    tardy_model.dynamics_step(delta_t=delta_t)
    e_pots.append(tardy_model.e_pot())
    e_kins.append(tardy_model.e_kin())
  e_tots = e_pots + e_kins
  tardy_model.check_d_e_pot_d_q()
  print >> out, "degrees of freedom:", tardy_model.degrees_of_freedom
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

def exercise_dynamics_quick(
      out, tardy_model, n_dynamics_steps, delta_t=0.0001):
  relative_range = exercise_tardy_model(
    out=out,
    n_dynamics_steps=n_dynamics_steps,
    delta_t=delta_t,
    tardy_model=tardy_model)
  if (out is not sys.stdout):
    if (len(tardy_model.tardy_tree.cluster_manager.loop_edges) == 0):
      assert relative_range < 1.e-5
    else:
      assert relative_range < 2.e-4
  print >> out

def exercise_minimization_quick(out, tardy_model, max_iterations=5):
  print >> out, "Minimization:"
  print >> out, "  start e_pot:", tardy_model.e_pot()
  e_pot_start = tardy_model.e_pot()
  tardy_model.minimization(max_iterations=max_iterations)
  print >> out, "  final e_pot:", tardy_model.e_pot()
  e_pot_final = tardy_model.e_pot()
  if (out is not sys.stdout):
    assert e_pot_final < e_pot_start * 0.7
  print >> out

def construct_tardy_model(
      labels,
      sites,
      masses,
      tardy_tree):
  cm = tardy_tree.cluster_manager
  return tardy_model(
    labels=labels,
    sites=sites,
    masses=masses,
    tardy_tree=tardy_tree,
    potential_obj=potential_object(
      sites=sites,
      wells=sites,
      restraint_edges=cm.loop_edges+cm.loop_edge_bendings))

def exercise_with_tardy_model(out, tardy_model, n_dynamics_steps):
  tardy_model.tardy_tree.show_summary(out=out, vertex_labels=None)
  #
  ri = tardy_model.root_indices()
  assert len(ri) == 1
  nosiet = tardy_model.number_of_sites_in_each_tree()
  assert len(nosiet) == 1
  assert nosiet[0][0] == ri[0]
  assert nosiet[0][1] == len(tardy_model.sites)
  somiet = tardy_model.number_of_sites_in_each_tree()
  assert len(somiet) == 1
  assert somiet[0][0] == ri[0]
  assert somiet[0][1] == len(tardy_model.sites)
  #
  exercise_qd_e_kin_scales(tardy_model=tardy_model)
  exercise_random_velocities(tardy_model=tardy_model)
  exercise_linear_velocity_manipulations(tardy_model=tardy_model)
  tardy_model.assign_random_velocities(e_kin_target=1)
  assert approx_equal(tardy_model.e_kin(), 1)
  exercise_dynamics_quick(
    out=out, tardy_model=tardy_model, n_dynamics_steps=n_dynamics_steps)
  exercise_minimization_quick(out=out, tardy_model=tardy_model)

n_test_models = len(test_cases_tardy_pdb.test_cases)

def get_test_model_by_index(i, fixed_vertex_lists=[]):
  tc = test_cases_tardy_pdb.test_cases[i]
  tt = tc.tardy_tree_construct(fixed_vertex_lists=fixed_vertex_lists)
  return construct_tardy_model(
    labels=tc.labels,
    sites=tc.sites,
    masses=[1.0]*len(tc.sites),
    tardy_tree=tt)

test_case_5_fixed_vertices_expected_dof = [
  ([0,16,17], 2),
  ([16,17], 5),
  ([12], 8)]

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
  exercise_accumulate_in_each_tree()
  exercise_near_singular_hinges()
  exercise_fixed_vertices_special_cases()
  #
  if (1):
    for i in xrange(n_test_models):
      print >> out, "test model index:", i
      tardy_model = get_test_model_by_index(i=i)
      exercise_with_tardy_model(
        out=out, tardy_model=tardy_model, n_dynamics_steps=n_dynamics_steps)
      if (i == 0):
        assert tardy_model.degrees_of_freedom == 3
        fixed_vertices = [0]
        print >> out, "test model index:", i, \
          "fixed_vertices:", fixed_vertices
        tardy_model = get_test_model_by_index(
          i=i, fixed_vertex_lists=[fixed_vertices])
        assert tardy_model.degrees_of_freedom == 0
      elif (i == 5):
        assert tardy_model.degrees_of_freedom == 11
        for fixed_vertices,expected_dof in \
              test_case_5_fixed_vertices_expected_dof:
          print >> out, "test model index:", i, \
            "fixed_vertices:", fixed_vertices
          tardy_model = get_test_model_by_index(
            i=i, fixed_vertex_lists=[fixed_vertices])
          assert tardy_model.degrees_of_freedom == expected_dof
          exercise_with_tardy_model(
            out=out,
            tardy_model=tardy_model,
            n_dynamics_steps=n_dynamics_steps)
  #
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
