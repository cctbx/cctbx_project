from __future__ import division
from scitbx.rigid_body.essence import featherstone
from scitbx.rigid_body.essence import body_lib
import scitbx.lbfgs
from scitbx.array_family import flex
from scitbx import matrix
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
  body_lib.set_cb_tree(bodies=result)
  return result

class model(featherstone.system_model):

  def __init__(O,
        labels,
        sites,
        masses,
        tardy_tree,
        potential_obj,
        near_singular_hinges_angular_tolerance_deg=5):
    super(model, O).__init__(bodies=construct_bodies(
      sites=sites,
      masses=masses,
      cluster_manager=tardy_tree.cluster_manager,
      near_singular_hinges_angular_tolerance_deg=
        near_singular_hinges_angular_tolerance_deg))
    O.labels = labels
    O.sites = sites
    O.masses = masses
    O.tardy_tree = tardy_tree
    O.potential_obj = potential_obj
    O.near_singular_hinges_angular_tolerance_deg = \
      near_singular_hinges_angular_tolerance_deg
    O.flag_positions_as_changed()

  def flag_positions_as_changed(O):
    O.__sites_moved = None
    O.__e_pot = None
    O.__d_e_pot_d_sites = None
    O.__f_ext_array = None
    super(model, O).flag_positions_as_changed()

  def flag_velocities_as_changed(O):
    O.__qdd_array = None
    super(model, O).flag_velocities_as_changed()

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
    """
Gradients of potential energy (defined via f_ext_array) w.r.t. positional
coordinates q.
    """
    return [body.joint.tau_as_d_e_pot_d_q(tau=tau)
      for body,tau in zip(O.bodies,
                          O.f_ext_as_tau(f_ext_array=O.f_ext_array()))]

  def d_e_pot_d_q_packed(O):
    result = flex.double()
    result.reserve(O.q_packed_size)
    for v in O.d_e_pot_d_q():
      result.extend(flex.double(v))
    assert result.size() == O.q_packed_size
    return result

  def qdd_array(O):
    if (O.__qdd_array is None):
      O.__qdd_array = O.forward_dynamics_ab(
        tau_array=None, f_ext_array=O.f_ext_array(), grav_accn=None)
    return O.__qdd_array

  def e_tot(O):
    return O.e_kin() + O.e_pot()

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
    O.tardy_model.unpack_q(q_packed=O.x)
    f = O.tardy_model.e_pot()
    g = O.tardy_model.d_e_pot_d_q_packed()
    return f, g
