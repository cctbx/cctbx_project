from cctbx import restraints
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.python_utils.misc import adopt_init_args
import scitbx.lbfgs
import sys

class energies:

  def __init__(self, sites_cart,
                     bond_proxies=None,
                     repulsion_proxies=None,
                     compute_gradients=0001,
                     disable_asu_cache=00000):
    if (compute_gradients):
      self.gradients = flex.vec3_double(sites_cart.size(), [0,0,0])
    else:
      self.gradients = None
    if (bond_proxies is None):
      self.n_bond_proxies = 0
      self.bond_residual_sum = 0
    else:
      self.n_bond_proxies = bond_proxies.n_total()
      self.bond_residual_sum = restraints.bond_residual_sum(
        sites_cart=sites_cart,
        sorted_asu_proxies=bond_proxies,
        gradient_array=self.gradients,
        disable_cache=disable_asu_cache)
    if (repulsion_proxies is None):
      self.n_repulsion_proxies = 0
      self.repulsion_residual_sum = 0
    else:
      self.n_repulsion_proxies = repulsion_proxies.n_total()
      self.repulsion_residual_sum = restraints.repulsion_residual_sum(
        sites_cart=sites_cart,
        sorted_asu_proxies=repulsion_proxies,
        gradient_array=self.gradients,
        function=restraints.repulsion_function(),
        disable_cache=disable_asu_cache)

  def target(self):
    return(self.bond_residual_sum
         + self.repulsion_residual_sum)

  def gradient_norm(self):
    if (self.gradients is not None):
      return flex.sum_sq(self.gradients.as_double())

  def show(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "target: %.6g" % self.target()
    if (self.gradients is not None):
      print >> f, "  norm of gradients: %.6g" % self.gradient_norm()
    print >> f, "  bond_residual_sum (n=%d): %.6g" % (
      self.n_bond_proxies, self.bond_residual_sum)
    print >> f, "  repulsion_residual_sum (n=%d): %.6g" % (
      self.n_repulsion_proxies, self.repulsion_residual_sum)

class lbfgs:

  def __init__(self, structure,
                     repulsion_params,
                     repulsion_types,
                     bond_params_table,
                     shell_sym_tables,
                     nonbonded_distance_cutoff,
                     nonbonded_buffer,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=None):
    adopt_init_args(self, locals())
    if (lbfgs_termination_params is None):
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations=1000)
    self._sites_cart = structure.sites_cart()
    self._sites_shifted = self._sites_cart
    self.n = self._sites_cart.size()*3
    self.x = flex.double(self.n, 0)
    self.first_target_result = None
    self._sites_cart_proxy_calculation = None
    self.activate_repulsion = 00000
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    structure.set_sites_cart(sites_cart=self._sites_shifted)
    structure.random_shift_sites(max_shift_cart=0.2)
    structure.apply_symmetry_sites()
    self._sites_cart = structure.sites_cart()
    self._sites_shifted = self._sites_cart
    self.x = flex.double(self.n, 0)
    self._sites_cart_proxy_calculation = None
    self.activate_repulsion = 0001
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    self.apply_shifts()
    self.compute_target(compute_gradients=0001)
    structure.set_sites_cart(sites_cart=self._sites_shifted)
    self.final_target_result = self._target_result
    del self.x
    del self.f
    del self.g
    if (self._sites_cart_proxy_calculation is not None):
      del self._bond_proxies
      del self._repulsion_proxies
    del self._sites_cart_proxy_calculation
    del self.activate_repulsion
    del self.minimizer
    del self._target_result
    del self._sites_shifted
    del self._sites_cart

  def apply_shifts(self):
    self._sites_shifted = self._sites_cart + flex.vec3_double(self.x)
    unit_cell = self.structure.unit_cell()
    site_symmetry_table = self.structure.site_symmetry_table()
    for i_seq in site_symmetry_table.special_position_indices():
      site_frac = unit_cell.fractionalize(self._sites_shifted[i_seq])
      site_special_frac = site_symmetry_table.get(i_seq).special_op() \
                        * site_frac
      distance_moved = unit_cell.distance(site_special_frac, site_frac)
      if (distance_moved > 1.e-2):
        print "WARNING: LARGE distance_moved: i_seq+1=%d, %.6g" % (
          i_seq+1, distance_moved)
      # it is essential to reset the site here because over
      # many cycles rounding error accumulate otherwise
      self._sites_shifted[i_seq] = unit_cell.orthogonalize(site_special_frac)

  def compute_target(self, compute_gradients):
    if (   self._sites_cart_proxy_calculation is None
        or self._sites_cart_proxy_calculation.max_distance(self._sites_shifted)
            >= self.nonbonded_buffer):
      shell_distance_cutoffs = [flex.max(crystal.get_distances(
        pair_sym_table=shell_sym_table,
        orthogonalization_matrix
          =self.structure.unit_cell().orthogonalization_matrix(),
        sites_frac=self.structure.unit_cell().fractionalization_matrix()
                  *self._sites_shifted)) * (1+1.e-6)
          for shell_sym_table in self.shell_sym_tables]
      asu_mappings = crystal.symmetry.asu_mappings(self.structure,
        buffer_thickness=max(
          max(shell_distance_cutoffs),
          self.nonbonded_distance_cutoff+self.nonbonded_buffer))
      asu_mappings.process_sites_cart(
        original_sites=self._sites_shifted,
        site_symmetry_table=self.structure.site_symmetry_table())
      shell_asu_tables = [crystal.pair_asu_table(
        asu_mappings=asu_mappings).add_pair_sym_table(
          sym_table=shell_sym_table)
            for shell_sym_table in self.shell_sym_tables]
      pair_proxies = restraints.pair_proxies(
        repulsion_params=self.repulsion_params,
        repulsion_types=self.repulsion_types,
        bond_params_table=self.bond_params_table,
        shell_asu_tables=shell_asu_tables,
        bonded_distance_cutoff=max(shell_distance_cutoffs),
        nonbonded_distance_cutoff=self.nonbonded_distance_cutoff,
        nonbonded_buffer=self.nonbonded_buffer)
      self._bond_proxies = pair_proxies.bond_proxies
      if (not self.activate_repulsion):
        self._repulsion_proxies = None
      else:
        self._repulsion_proxies = pair_proxies.repulsion_proxies
      self._sites_cart_proxy_calculation = self._sites_shifted
    self._target_result = energies(
      sites_cart=self._sites_shifted,
      bond_proxies=self._bond_proxies,
      repulsion_proxies=self._repulsion_proxies,
      compute_gradients=compute_gradients)

  def __call__(self):
    if (self.first_target_result is None):
      assert self.x.all_eq(0)
    else:
      self.apply_shifts()
    self.compute_target(compute_gradients=0001)
    if (self.first_target_result is None):
      self.first_target_result = self._target_result
    self.f = self._target_result.target()
    self.g = self._target_result.gradients.as_double()
    return self.x, self.f, self.g

  def callback_after_step(self, minimizer):
    if (0):
      print "minimizer f,iter,nfun:", self.f,minimizer.iter(),minimizer.nfun()
    if (0):
      self._target_result.show()
