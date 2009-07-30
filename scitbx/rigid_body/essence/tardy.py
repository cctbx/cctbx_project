from __future__ import division
from scitbx.rigid_body.essence import featherstone
from scitbx.rigid_body.essence import body_lib
from scitbx.rigid_body.essence import joint_lib
import scitbx.lbfgs
from scitbx.array_family import flex
from scitbx import matrix
import random
import math

def construct_bodies(
      sites,
      masses,
      cluster_manager,
      near_singular_hinges_angular_tolerance_deg=5):
  assert len(masses) == len(sites)
  result = []
  abs_cos_limit = abs(math.cos(math.radians(
    near_singular_hinges_angular_tolerance_deg)))
  fvgci = cluster_manager.fixed_vertices_given_cluster_index_dict()
  for ic,cluster in enumerate(cluster_manager.clusters):
    body_sites = [matrix.col(sites[i]) for i in cluster]
    body_masses = [masses[i] for i in cluster]
    he = cluster_manager.hinge_edges[ic]
    fixed_vertices = fvgci.get(ic)
    if (fixed_vertices is not None):
      if (   len(fixed_vertices) > 2
          or len(fixed_vertices) == len(cluster)):
        body = body_lib.zero_dof(sites=body_sites, masses=body_masses)
      elif (len(fixed_vertices) == 1):
        body = body_lib.spherical(
          sites=body_sites,
          masses=body_masses,
          pivot=sites[fixed_vertices[0]])
      elif (len(fixed_vertices) == 2):
        normal_sites = [matrix.col(sites[i]) for i in fixed_vertices]
        pivot = normal_sites[1]
        axis = pivot - normal_sites[0]
        for site in body_sites:
          abs_cos = abs(axis.cos_angle(site - pivot, value_if_undefined=1))
          if (abs_cos < abs_cos_limit):
            body = body_lib.revolute(
              sites=body_sites,
              masses=body_masses,
              pivot=pivot,
              normal=axis.normalize())
            break
        else:
          body = body_lib.zero_dof(sites=body_sites, masses=body_masses)
      else:
        raise AssertionError
      body.parent = -1
    elif (he[0] == -1):
      if (len(body_sites) == 1):
        body = body_lib.translational(
          sites=body_sites, masses=body_masses)
      else:
        body = body_lib.six_dof(sites=body_sites, masses=body_masses)
      body.parent = -1
    else:
      normal_sites = [matrix.col(sites[i]) for i in he]
      body = body_lib.revolute(
        sites=body_sites,
        masses=body_masses,
        pivot=normal_sites[1],
        normal=(normal_sites[1]-normal_sites[0]).normalize())
      body.parent = cluster_manager.cluster_indices[he[1]]
    body.i_seqs = cluster
    result.append(body)
  return result

class model(object):

  def __init__(O,
        labels,
        sites,
        masses,
        tardy_tree,
        potential_obj,
        near_singular_hinges_angular_tolerance_deg=5):
    O.labels = labels
    O.sites = sites
    O.masses = masses
    O.tardy_tree = tardy_tree
    O.potential_obj = potential_obj
    O.near_singular_hinges_angular_tolerance_deg = \
      near_singular_hinges_angular_tolerance_deg
    O.bodies = construct_bodies(
      sites=sites,
      masses=masses,
      cluster_manager=tardy_tree.cluster_manager,
      near_singular_hinges_angular_tolerance_deg=
        near_singular_hinges_angular_tolerance_deg)
    O.number_of_trees = 0
    O.degrees_of_freedom = 0
    O.packed_q_size = 0
    O.packed_qd_size = 0
    for body in O.bodies:
      if (body.parent == -1): O.number_of_trees += 1
      O.degrees_of_freedom += body.joint.degrees_of_freedom
      O.packed_q_size += body.joint.q_size
    O.packed_qd_size = O.degrees_of_freedom
    O.flag_positions_as_changed()

  def flag_positions_as_changed(O):
    O.__featherstone_system_model = None
    O.__aja_array = None
    O.__jar_array = None
    O.__sites_moved = None
    O.__e_pot = None
    O.__d_e_pot_d_sites = None
    O.__f_ext_array = None
    O.flag_velocities_as_changed()

  def flag_velocities_as_changed(O):
    O.__qdd_array = None
    O.__e_kin = None

  def root_indices(O):
    result = []
    for ib,body in enumerate(O.bodies):
      if (body.parent == -1):
        result.append(ib)
    return result

  def pack_q(O):
    result = flex.double()
    result.reserve(O.packed_q_size)
    for body in O.bodies:
      result.extend(flex.double(body.joint.get_q()))
    assert result.size() == O.packed_q_size
    return result

  def unpack_q(O, packed_q):
    assert packed_q.size() == O.packed_q_size
    i = 0
    for body in O.bodies:
      n = body.joint.q_size
      body.joint = body.joint.new_q(q=packed_q[i:i+n])
      i += n
    assert i == O.packed_q_size
    O.flag_positions_as_changed()

  def pack_qd(O):
    result = flex.double()
    result.reserve(O.packed_qd_size)
    for body in O.bodies:
      result.extend(flex.double(body.qd))
    assert result.size() == O.packed_qd_size
    return result

  def unpack_qd(O, packed_qd):
    assert packed_qd.size() == O.packed_qd_size
    i = 0
    for body in O.bodies:
      n = body.joint.degrees_of_freedom
      body.qd = matrix.col(packed_qd[i:i+n])
      i += n
    assert i == O.packed_qd_size
    O.flag_velocities_as_changed()

  def _accumulate_in_each_tree(O, attr):
    result = []
    accu = [0] * len(O.bodies)
    for ib in xrange(len(O.bodies)-1,-1,-1):
      body = O.bodies[ib]
      accu[ib] += getattr(body, attr)
      if (body.parent == -1):
        result.append((ib, accu[ib]))
      else:
        accu[body.parent] += accu[ib]
    return result

  def number_of_sites_in_each_tree(O):
    return O._accumulate_in_each_tree(attr="number_of_sites")

  def sum_of_masses_in_each_tree(O):
    return O._accumulate_in_each_tree(attr="sum_of_masses")

  def mean_linear_velocity(O, number_of_sites_in_each_tree):
    if (number_of_sites_in_each_tree is None):
      number_of_sites_in_each_tree = O.number_of_sites_in_each_tree()
    sum_v = matrix.col((0,0,0))
    sum_n = 0
    for ib,n in number_of_sites_in_each_tree:
      body = O.bodies[ib]
      v = body.joint.get_linear_velocity(qd=body.qd)
      if (v is None): continue
      sum_v += v * n
      sum_n += n
    if (sum_n == 0):
      return None
    return sum_v / sum_n

  def subtract_from_linear_velocities(O, number_of_sites_in_each_tree, value):
    if (number_of_sites_in_each_tree is None):
      number_of_sites_in_each_tree = O.number_of_sites_in_each_tree()
    for ib,n in number_of_sites_in_each_tree:
      body = O.bodies[ib]
      v = body.joint.get_linear_velocity(qd=body.qd)
      if (v is None): continue
      body.qd = body.joint.new_linear_velocity(qd=body.qd, value=v-value)

  def featherstone_system_model(O):
    if (O.__featherstone_system_model is None):
      O.__featherstone_system_model = featherstone.system_model(
        bodies=O.bodies)
    return O.__featherstone_system_model

  def aja_array(O):
    if (O.__aja_array is None):
      O.__aja_array = []
      for body in O.bodies:
        aja = body.alignment.cb_b0 * body.joint.cb_sp * body.alignment.cb_0b
        if (body.parent != -1):
          aja = O.__aja_array[body.parent] * aja
        O.__aja_array.append(aja)
    return O.__aja_array

  def jar_array(O):
    if (O.__jar_array is None):
      O_aja = O.aja_array()
      O.__jar_array = []
      for body in O.bodies:
        jar = body.joint.cb_ps.r * body.alignment.cb_0b.r
        if (body.parent != -1):
          jar *= O_aja[body.parent].r.transpose()
        O.__jar_array.append(jar)
    return O.__jar_array

  def sites_moved(O):
    if (O.__sites_moved is None):
      O_aja = O.aja_array()
      O.__sites_moved = [None] * len(O.sites)
      n_done = 0
      clusters = O.tardy_tree.cluster_manager.clusters
      for ib in xrange(len(O.bodies)):
        aja = O_aja[ib]
        for i_seq in clusters[ib]:
          assert O.__sites_moved[i_seq] is None
          O.__sites_moved[i_seq] = aja * O.sites[i_seq]
          n_done += 1
      assert n_done == len(O.sites)
    return O.__sites_moved

  def e_pot(O):
    if (O.__e_pot is None):
      if (O.potential_obj is None):
        O.__e_pot = 0
      else:
        O.__e_pot = O.potential_obj.e_pot(
          sites_moved=O.sites_moved())
    return O.__e_pot

  def d_e_pot_d_sites(O):
    if (O.__d_e_pot_d_sites is None):
      if (O.potential_obj is None):
        O.__d_e_pot_d_sites = [matrix.col((0,0,0))] * len(O.sites)
      else:
        O.__d_e_pot_d_sites = O.potential_obj.d_e_pot_d_sites(
          sites_moved=O.sites_moved())
    return O.__d_e_pot_d_sites

  def f_ext_array(O):
    if (O.__f_ext_array is None):
      O_jar_array = O.jar_array()
      O_d_e_pot_d_sites = O.d_e_pot_d_sites()
      O.__f_ext_array = []
      clusters = O.tardy_tree.cluster_manager.clusters
      for ib,body in enumerate(O.bodies):
        cb_0b = body.alignment.cb_0b
        jar = O_jar_array[ib]
        f = matrix.col((0,0,0))
        nc = matrix.col((0,0,0))
        for i_seq in clusters[ib]:
          s = O.sites[i_seq]
          force_bf = -(jar * O_d_e_pot_d_sites[i_seq])
          f += force_bf
          nc += (cb_0b * s).cross(force_bf)
        O.__f_ext_array.append(matrix.col((nc, f)).resolve_partitions())
    return O.__f_ext_array

  def d_e_pot_d_q(O):
    return O.featherstone_system_model().d_e_pot_d_q(
      f_ext_array=O.f_ext_array())

  def d_e_pot_d_q_packed(O):
    result = flex.double()
    for v in O.featherstone_system_model().d_e_pot_d_q(
               f_ext_array=O.f_ext_array()):
      result.extend(flex.double(v))
    return result

  def qdd_array(O):
    if (O.__qdd_array is None):
      O.__qdd_array = O.featherstone_system_model().forward_dynamics_ab(
        tau_array=None, f_ext_array=O.f_ext_array())
    return O.__qdd_array

  def e_kin(O):
    if (O.__e_kin is None):
      O.__e_kin = O.featherstone_system_model().e_kin()
    return O.__e_kin

  def e_tot(O):
    return O.e_kin() + O.e_pot()

  def reset_e_kin(O, e_kin_target, e_kin_epsilon=1e-12):
    assert e_kin_target >= 0
    assert e_kin_epsilon > 0
    O_e_kin = O.e_kin()
    if (O_e_kin >= e_kin_epsilon):
      factor = (e_kin_target / O_e_kin)**0.5
      for body in O.bodies:
        body.qd *= factor
    O.flag_velocities_as_changed()

  def assign_zero_velocities(O):
    for body in O.bodies:
      body.qd = body.joint.qd_zero
    O.flag_velocities_as_changed()

  def assign_random_velocities(O,
        e_kin_target=None,
        e_kin_epsilon=1e-12,
        random_gauss=None):
    if (e_kin_target is None):
      work_e_kin_target = 1
    elif (e_kin_target == 0):
      O.assign_zero_velocities()
      return
    else:
      assert e_kin_target >= 0
      work_e_kin_target = e_kin_target
    qd_e_kin_scales = flex.double(
      O.featherstone_system_model().qd_e_kin_scales(
        e_kin_epsilon=e_kin_epsilon))
    if (O.degrees_of_freedom != 0):
      qd_e_kin_scales *= (work_e_kin_target / O.degrees_of_freedom)**0.5
    if (random_gauss is None):
      random_gauss = random.gauss
    i_qd = 0
    for body in O.bodies:
      qd_new = []
      for qd in body.joint.qd_zero:
        qd_new.append(qd + random_gauss(mu=0, sigma=qd_e_kin_scales[i_qd]))
        i_qd += 1
      body.qd = matrix.col(qd_new)
    assert i_qd == O.degrees_of_freedom
    O.flag_velocities_as_changed()
    if (e_kin_target is not None):
      O.reset_e_kin(e_kin_target=e_kin_target, e_kin_epsilon=e_kin_epsilon)
    return qd_e_kin_scales

  def dynamics_step(O, delta_t):
    O_qdd_array = O.qdd_array()
    for body in O.bodies:
      body.joint = body.joint.time_step_position(
        qd=body.qd, delta_t=delta_t)
    for body,qdd in zip(O.bodies, O_qdd_array):
      body.qd = body.joint.time_step_velocity(
        qd=body.qd, qdd=qdd, delta_t=delta_t)
    O.flag_positions_as_changed()

  def minimization(O, max_iterations=None, callback_after_step=None):
    return refinery(
      tardy_model=O,
      max_iterations=max_iterations,
      callback_after_step=callback_after_step)

class refinery(object):

  def __init__(O, tardy_model, max_iterations=None, callback_after_step=None):
    O.tardy_model = tardy_model
    O.__callback_after_step = callback_after_step
    O.x = tardy_model.pack_q()
    O.function_evaluations_total = 0
    O.lbfgs_steps_total = 0
    O.lbfgs_restarts = 0
    while True:
      lbfgs_steps_prev = O.lbfgs_steps_total
      scitbx.lbfgs.run(
        target_evaluator=O,
        termination_params=scitbx.lbfgs.termination_parameters(
          max_iterations=max_iterations),
        exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
          ignore_line_search_failed_step_at_lower_bound=True))
      if (O.lbfgs_steps_total == lbfgs_steps_prev):
        break
      if (max_iterations is not None):
        if (O.lbfgs_steps_total >= max_iterations):
          break
        max_iterations -= O.lbfgs_steps_total
      O.lbfgs_restarts += 1
    del O.x

  def callback_after_step(O, minimizer):
    O.lbfgs_steps_total += 1
    if (O.__callback_after_step is not None):
      O.__callback_after_step(minimizer=minimizer)

  def compute_functional_and_gradients(O):
    O.function_evaluations_total += 1
    O.tardy_model.unpack_q(packed_q=O.x)
    f = O.tardy_model.e_pot()
    g = O.tardy_model.d_e_pot_d_q_packed()
    return f, g
