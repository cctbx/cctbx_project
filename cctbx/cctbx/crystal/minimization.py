from cctbx import restraints
from cctbx.array_family import flex
from scitbx.python_utils.misc import adopt_init_args
from scitbx.python_utils.misc import time_log
import scitbx.lbfgs
import sys

class energies:

  def __init__(self, sites_cart,
                     asu_mappings,
                     bond_asu_proxies=None,
                     repulsion_asu_proxies=None,
                     compute_gradients=0001):
    if (compute_gradients):
      self.gradients = flex.vec3_double(sites_cart.size(), [0,0,0])
    else:
      self.gradients = None
    if (bond_asu_proxies is None):
      self.n_bond_asu_proxies = 0
      self.bond_asu_residual_sum = 0
    else:
      assert asu_mappings is not None
      self.n_bond_asu_proxies = len(bond_asu_proxies)
      self.bond_asu_residual_sum = restraints.bond_residual_sum(
        sites_cart=sites_cart,
        asu_mappings=asu_mappings,
        proxies=bond_asu_proxies,
        gradient_array=self.gradients,
        disable_cache=0001)
    if (repulsion_asu_proxies is None):
      self.n_repulsion_asu_proxies = 0
      self.repulsion_residual_sum = 0
    else:
      self.n_repulsion_asu_proxies = repulsion_asu_proxies.size()
      self.repulsion_residual_sum = restraints.repulsion_residual_sum(
        sites_cart=sites_cart,
        asu_mappings=asu_mappings,
        proxies=repulsion_asu_proxies,
        gradient_array=self.gradients,
        function=restraints.repulsion_function(),
        disable_cache=0001)

  def target(self):
    return(self.bond_asu_residual_sum
         + self.repulsion_residual_sum)

  def gradient_norm(self):
    if (self.gradients is not None):
      return flex.sum_sq(self.gradients.as_double())

  def show(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "target: %.6g" % self.target()
    if (self.gradients is not None):
      print >> f, "  norm of gradients: %.6g" % self.gradient_norm()
    print >> f, "  bond_asu_residual_sum (n=%d): %.6g" % (
      self.n_bond_asu_proxies, self.bond_asu_residual_sum)
    print >> f, "  repulsion_residual_sum (n=%d): %.6g" % (
      self.n_repulsion_asu_proxies, self.repulsion_residual_sum)

class lbfgs:

  def __init__(self, sites_cart,
                     asu_mappings=None,
                     bond_asu_proxies=None,
                     repulsion_asu_proxies=None,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=None):
    adopt_init_args(self, locals())
    self.time_total = time_log("lbfgs total").start()
    self.time_compute_target = time_log("lbfgs compute_target")
    if (lbfgs_termination_params is None):
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations=1000)
    self.n = sites_cart.size()*3
    self.x = flex.double(self.n, 0)
    self._sites_shifted = self.sites_cart
    self.first_target_value = None
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    self.apply_shifts()
    self.compute_target(compute_gradients=0001)
    sites_cart.clear()
    sites_cart.extend(self._sites_shifted)
    del self._sites_shifted
    self.final_target_value = self.target_result.target()
    self.time_total.stop()

  def apply_shifts(self):
    self._sites_shifted = self.sites_cart + flex.vec3_double(self.x)
    if (1):
      unit_cell = self.asu_mappings.unit_cell()
      site_symmetry_table = self.asu_mappings.site_symmetry_table()
      for i_seq in site_symmetry_table.special_position_indices():
        site_frac = unit_cell.fractionalize(self._sites_shifted[i_seq])
        site_special_frac = site_symmetry_table.get(i_seq).special_op() \
                          * site_frac
        distance_moved = unit_cell.distance(site_special_frac, site_frac)
        if (distance_moved > 1.e-2):
          print "LARGE distance_moved: i_seq+1=%d, %.6g" % (
            i_seq+1, distance_moved)

  def compute_target(self, compute_gradients):
    self.time_compute_target.start()
    self.target_result = energies(
      sites_cart=self._sites_shifted,
      asu_mappings=self.asu_mappings,
      bond_asu_proxies=self.bond_asu_proxies,
      repulsion_asu_proxies=self.repulsion_asu_proxies,
      compute_gradients=compute_gradients)
    self.time_compute_target.stop()

  def __call__(self):
    if (self.first_target_value is None):
      assert self.x.all_eq(0)
    else:
      self.apply_shifts()
    self.compute_target(compute_gradients=0001)
    self.f = self.target_result.target()
    if (self.first_target_value is None):
      self.first_target_value = self.f
    self.g = self.target_result.gradients.as_double()
    return self.x, self.f, self.g

  def callback_after_step(self, minimizer):
    if (0):
      print "minimizer f,iter,nfun:", self.f,minimizer.iter(),minimizer.nfun()
    if (0):
      self.target_result.show()

  def show_times(self):
    print self.time_total.report()
    print self.time_compute_target.report()
    print "time unaccounted for: %.2f" % (
      self.time_total.accumulation
      - self.time_compute_target.accumulation)
