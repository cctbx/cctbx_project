from cctbx.array_family import flex
import scitbx.lbfgs

class empty: pass

class lbfgs:

  def __init__(self,
      sites_cart,
      restraints_manager,
      lbfgs_termination_params=None,
      lbfgs_core_params=None,
      disable_asu_cache=00000):
    if (lbfgs_termination_params is None):
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations=1000)
    self.n = sites_cart.size()*3
    self.x = flex.double(self.n, 0)
    self.tmp = empty()
    self.tmp.restraints_manager = restraints_manager
    self.tmp.disable_asu_cache = disable_asu_cache
    self.tmp.sites_cart = sites_cart
    self.tmp.sites_shifted = sites_cart
    self.tmp.lock_pair_proxies = 00000
    self.first_target_result = None
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    self.apply_shifts()
    self.compute_target(compute_gradients=0001)
    self.final_target_result = self.tmp.target_result
    sites_cart.clear()
    sites_cart.extend(self.tmp.sites_shifted)
    del self.tmp
    self.first_target_value = self.first_target_result.target()
    self.final_target_value = self.final_target_result.target()

  def apply_shifts(self):
    self.tmp.sites_shifted = self.tmp.sites_cart + flex.vec3_double(self.x)
    if (self.tmp.restraints_manager.crystal_symmetry is not None):
      unit_cell = self.tmp.restraints_manager.crystal_symmetry.unit_cell()
      site_symmetry_table = self.tmp.restraints_manager.site_symmetry_table
      assert site_symmetry_table is not None
      for i_seq in site_symmetry_table.special_position_indices():
        site_frac = unit_cell.fractionalize(self.tmp.sites_shifted[i_seq])
        site_special_frac = site_symmetry_table.get(i_seq).special_op() \
                          * site_frac
        distance_moved = unit_cell.distance(site_special_frac, site_frac)
        if (distance_moved > 1.e-2):
          raise AssertionError("Corrupt gradient calculations.")
        # it is essential to reset the site here because
        # otherwise rounding error accumulate over many cycles
        self.tmp.sites_shifted[i_seq] = unit_cell.orthogonalize(
          site_special_frac)

  def compute_target(self, compute_gradients):
    self.tmp.target_result = self.tmp.restraints_manager.energies(
      sites_cart=self.tmp.sites_shifted,
      compute_gradients=compute_gradients,
      disable_asu_cache=self.tmp.disable_asu_cache,
      lock_pair_proxies=self.tmp.lock_pair_proxies)
    self.tmp.lock_pair_proxies = 0001

  def callback_after_step(self, minimizer):
    self.tmp.lock_pair_proxies = 00000

  def __call__(self):
    if (self.first_target_result is None):
      assert self.x.all_eq(0)
    else:
      self.apply_shifts()
    self.compute_target(compute_gradients=0001)
    self.f = self.tmp.target_result.target()
    if (self.first_target_result is None):
      self.first_target_result = self.tmp.target_result
    self.g = self.tmp.target_result.gradients.as_double()
    return self.x, self.f, self.g
