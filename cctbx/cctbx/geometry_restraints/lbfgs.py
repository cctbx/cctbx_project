from cctbx import crystal
from cctbx.array_family import flex
import scitbx.lbfgs

class empty: pass

class lbfgs:

  def __init__(self,
      sites_cart,
      geometry_restraints_manager,
      geometry_restraints_flags=None,
      lbfgs_termination_params=None,
      lbfgs_core_params=None,
      lbfgs_exception_handling_params=None,
      disable_asu_cache=False):
    if (lbfgs_termination_params is None):
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations=1000)
    self.n = sites_cart.size()*3
    self.x = flex.double(self.n, 0)
    self.tmp = empty()
    self.tmp.geometry_restraints_manager = geometry_restraints_manager
    self.tmp.geometry_restraints_flags = geometry_restraints_flags
    self.tmp.disable_asu_cache = disable_asu_cache
    self.tmp.sites_cart = sites_cart
    self.tmp.sites_shifted = sites_cart
    self.tmp.lock_pair_proxies = False
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
    self.first_target_value = self.first_target_result.target()
    self.final_target_value = self.final_target_result.target()

  def apply_shifts(self):
    self.tmp.sites_shifted = self.tmp.sites_cart + flex.vec3_double(self.x)
    if (self.tmp.geometry_restraints_manager.crystal_symmetry is not None):
      unit_cell \
        = self.tmp.geometry_restraints_manager.crystal_symmetry.unit_cell()
      site_symmetry_table \
        = self.tmp.geometry_restraints_manager.site_symmetry_table
      assert site_symmetry_table is not None
      for i_seq in site_symmetry_table.special_position_indices():
        self.tmp.sites_shifted[i_seq] = crystal.correct_special_position(
          unit_cell=unit_cell,
          special_op=site_symmetry_table.get(i_seq).special_op(),
          site_cart=self.tmp.sites_shifted[i_seq])

  def compute_target(self, compute_gradients):
    self.tmp.target_result = self.tmp.geometry_restraints_manager.energies(
      sites_cart=self.tmp.sites_shifted,
      flags=self.tmp.geometry_restraints_flags,
      compute_gradients=compute_gradients,
      disable_asu_cache=self.tmp.disable_asu_cache,
      lock_pair_proxies=self.tmp.lock_pair_proxies)
    self.tmp.lock_pair_proxies = True

  def callback_after_step(self, minimizer):
    self.tmp.lock_pair_proxies = False

  def __call__(self):
    if (self.first_target_result is None):
      assert self.x.all_eq(0)
    else:
      self.apply_shifts()
    self.compute_target(compute_gradients=True)
    self.f = self.tmp.target_result.target()
    if (self.first_target_result is None):
      self.first_target_result = self.tmp.target_result
    self.g = self.tmp.target_result.gradients.as_double()
    return self.x, self.f, self.g
