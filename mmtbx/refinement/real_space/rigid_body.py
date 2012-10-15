from __future__ import division
from cctbx import maptbx
from cctbx.array_family import flex
import scitbx.rigid_body
import scitbx.graph.tardy_tree
import scitbx.lbfgs
from scitbx import matrix

def real_space_rigid_body_gradients_simple(
      unit_cell,
      density_map,
      sites_cart_0,
      center_of_mass,
      q,
      unit_quaternion_delta=0.01,
      translation_delta=0.3):
  result = flex.double()
  q_delta = q.deep_copy()
  def get(i, delta):
    fs = []
    for signed_delta in [delta, -delta]:
      q_delta[i] = q[i] + signed_delta
      aja = matrix.rt(scitbx.rigid_body.joint_lib_six_dof_aja_simplified(
        center_of_mass=center_of_mass,
        q=q_delta))
      sites_cart_delta = aja * sites_cart_0
      rs_f = maptbx.real_space_target_simple(
        unit_cell=unit_cell,
        density_map=density_map,
        sites_cart=sites_cart_delta,
        selection=flex.bool(sites_cart_delta.size(),True))
      fs.append(rs_f)
    result.append((fs[0]-fs[1])/(2*delta))
  for i in xrange(4): get(i=i, delta=unit_quaternion_delta)
  for i in xrange(3): get(i=i+4, delta=translation_delta)
  return result

class refine(object):

  def __init__(self,
        residue,
        density_map,
        geometry_restraints_manager,
        real_space_target_weight,
        real_space_gradients_delta,
        lbfgs_termination_params,
        unit_cell,
        cctbx_geometry_restraints_flags=None,
        states_collector=None):
    self.states_collector = states_collector
    self.cctbx_geometry_restraints_flags = cctbx_geometry_restraints_flags
    self.residue = residue
    self.density_map = density_map
    self.geometry_restraints_manager = geometry_restraints_manager
    self.real_space_gradients_delta = real_space_gradients_delta
    self.real_space_target_weight = real_space_target_weight
    #
    self.unit_cell = unit_cell
    self.sites_cart_residue_0 = residue.atoms().extract_xyz()
    self.residue_center_of_mass = self.sites_cart_residue_0.mean()
    residue_tardy_tree = scitbx.graph.tardy_tree.construct(
      n_vertices=self.sites_cart_residue_0.size(),
      edge_list="all_in_one_rigid_body") \
        .build_tree() \
        .fix_near_singular_hinges(sites=None)
    self.residue_tardy_model = scitbx.rigid_body.tardy_model(
      labels=None,
      sites=self.sites_cart_residue_0,
      masses=flex.double(self.sites_cart_residue_0.size(), 1),
      tardy_tree=residue_tardy_tree,
      potential_obj=self)
    self.x = self.residue_tardy_model.pack_q()
    assert self.x.size() == 7 # other cases not implemented
    #
    self.number_of_function_evaluations = -1
    self.f_start, self.g_start = self.compute_functional_and_gradients()
    self.rs_f_start = self.rs_f
    lbfgs_exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True,
      ignore_line_search_failed_step_at_upper_bound = True,
      ignore_line_search_failed_maxfev              = True)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      exception_handling_params = lbfgs_exception_handling_params)
    self.f_final, self.g_final = self.compute_functional_and_gradients()
    self.rs_f_final = self.rs_f
    del self.rs_f
    del self.x
    del self.residue_center_of_mass
    del self.sites_cart_residue_0
    del self.unit_cell

  def compute_functional_and_gradients(self):
    if (self.number_of_function_evaluations == 0):
      self.number_of_function_evaluations += 1
      return self.f_start, self.g_start
    self.number_of_function_evaluations += 1
    self.residue_tardy_model.unpack_q(q_packed=self.x)
    self.sites_cart_residue = self.residue_tardy_model.sites_moved()
    if(self.states_collector is not None):
      self.states_collector.add(sites_cart = self.sites_cart_residue)
    rs_f = maptbx.real_space_target_simple(
      unit_cell=self.unit_cell,
      density_map=self.density_map,
      sites_cart=self.sites_cart_residue,
      selection=flex.bool(self.sites_cart_residue.size(),True))
    rs_g = real_space_rigid_body_gradients_simple(
      unit_cell=self.unit_cell,
      density_map=self.density_map,
      sites_cart_0=self.sites_cart_residue_0,
      center_of_mass=self.residue_center_of_mass,
      q=self.x)
    self.rs_f = rs_f
    rs_f *= -self.real_space_target_weight
    rs_g *= -self.real_space_target_weight
    if (self.geometry_restraints_manager is None):
      f = rs_f
      g = rs_g
    else:
      gr_e = self.geometry_restraints_manager.energies_sites(
        sites_cart=self.sites_cart_residue,
        flags = self.cctbx_geometry_restraints_flags,
        compute_gradients=True)
      self.__d_e_pot_d_sites = gr_e.gradients
      f = rs_f + gr_e.target
      g = rs_g + self.residue_tardy_model.d_e_pot_d_q_packed()
    return f, g.as_double()

  def d_e_pot_d_sites(self, sites_moved):
    result = self.__d_e_pot_d_sites
    del self.__d_e_pot_d_sites
    return result
