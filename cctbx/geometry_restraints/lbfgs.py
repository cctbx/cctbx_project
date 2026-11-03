from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx.array_family import flex
import scitbx.lbfgs

class empty: pass

class lbfgs(object):

  def __init__(self,
      sites_cart,
      geometry_restraints_manager,
      riding_h_manager=None,
      correct_special_position_tolerance=1.0,
      geometry_restraints_flags=None,
      lbfgs_termination_params=None,
      lbfgs_core_params=None,
      lbfgs_exception_handling_params=None,
      disable_asu_cache=False,
      sites_cart_selection=None,
      site_labels=None,
      states_collector=None):
    if (lbfgs_termination_params is None):
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations=1000)
    self.riding_h_manager = riding_h_manager
    self.correct_special_position_tolerance=correct_special_position_tolerance
    self.site_labels = site_labels
    self.states_collector = states_collector
    self.tmp = empty()
    self.rmsd_bonds, self.rmsd_angles = None, None
    if sites_cart_selection:
      self.sites_cart_selection = flex.bool(sites_cart_selection)
    else:
      self.sites_cart_selection = None
    # The refinable coordinates: the selected atoms, less the H atoms when
    # they ride on their parents.
    self.tmp.refinable = flex.bool(sites_cart.size(), True)
    if(self.sites_cart_selection is not None):
      self.tmp.refinable = self.tmp.refinable & self.sites_cart_selection
    if(self.riding_h_manager is not None):
      self.tmp.refinable = \
        self.tmp.refinable & self.riding_h_manager.not_hd_selection
      self.tmp.refinable_not_hd = \
        self.tmp.refinable.select(self.riding_h_manager.not_hd_selection)
    self.x = flex.double(self.tmp.refinable.count(True)*3, 0)
    self.tmp.geometry_restraints_manager = geometry_restraints_manager
    self.tmp.geometry_restraints_flags = geometry_restraints_flags
    self.tmp.disable_asu_cache = disable_asu_cache
    self.tmp.sites_cart = sites_cart
    self.tmp.sites_shifted = sites_cart
    self.first_target_result = None
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params,
      exception_handling_params=lbfgs_exception_handling_params)
    self.apply_shifts()
    self.compute_target(compute_gradients=True)
    self.final_target_result = self.tmp.target_result
    sites_cart.clear()
    sites_cart.extend(self.tmp.sites_shifted)
    del self.tmp
    del self.x
    self.first_target_value = self.first_target_result.target
    self.final_target_value = self.final_target_result.target

  def apply_shifts(self):
    shifted = self.tmp.sites_cart.select(self.tmp.refinable) + \
      flex.vec3_double(self.x)
    self.tmp.sites_shifted = self.tmp.sites_cart.deep_copy()
    self.tmp.sites_shifted.set_selected(self.tmp.refinable, shifted)
    if(self.riding_h_manager is not None):
      self.riding_h_manager.idealize_riding_h_positions(
        sites_cart = self.tmp.sites_shifted)
    if(self.states_collector is not None):
      self.states_collector.add(sites_cart = self.tmp.sites_shifted)
    if (self.tmp.geometry_restraints_manager.crystal_symmetry is not None):
      crystal_symmetry = self.tmp.geometry_restraints_manager.crystal_symmetry
      site_symmetry_table \
        = self.tmp.geometry_restraints_manager.site_symmetry_table
      assert site_symmetry_table is not None
      for i_seq in site_symmetry_table.special_position_indices():
        self.tmp.sites_shifted[i_seq] = crystal.correct_special_position(
          crystal_symmetry=crystal_symmetry,
          tolerance=self.correct_special_position_tolerance,
          special_op=site_symmetry_table.get(i_seq).special_op(),
          site_cart=self.tmp.sites_shifted[i_seq])

  def compute_target(self, compute_gradients):
    self.tmp.target_result = \
      self.tmp.geometry_restraints_manager.energies_sites(
        sites_cart=self.tmp.sites_shifted,
        flags=self.tmp.geometry_restraints_flags,
        compute_gradients=compute_gradients,
        disable_asu_cache=self.tmp.disable_asu_cache,
        site_labels=self.site_labels)

  def compute_functional_and_gradients(self):
    if (self.first_target_result is None):
      assert self.x.all_eq(0)
    else:
      self.apply_shifts()
    self.compute_target(compute_gradients=True)
    self.f = self.tmp.target_result.target
    if (self.first_target_result is None):
      self.first_target_result = self.tmp.target_result
    gradients = self.tmp.target_result.gradients
    if(self.riding_h_manager is None):
      gradients = gradients.select(self.tmp.refinable)
    else:
      # The gradients on the riding H atoms are transferred to their parents;
      # the result holds the non-H atoms only.
      gradients = self.riding_h_manager.gradients_reduced_cpp(
        gradients    = gradients,
        sites_cart   = self.tmp.sites_shifted,
        hd_selection = self.riding_h_manager.hd_selection)
      gradients = gradients.select(self.tmp.refinable_not_hd)
    self.g = gradients.as_double()
    return self.f, self.g
