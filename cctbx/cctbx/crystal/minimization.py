from cctbx import restraints
from cctbx.array_family import flex
from scitbx.python_utils.misc import adopt_init_args
from scitbx.python_utils.misc import time_log
import scitbx.lbfgs
import sys

class energies:

  def __init__(self, sites_cart,
                     asu_mappings=None,
                     bond_proxies=None,
                     bond_sym_proxies=None,
                     bond_sorted_proxies=None,
                     angle_proxies=None,
                     dihedral_proxies=None,
                     chirality_proxies=None,
                     planarity_proxies=None,
                     repulsion_proxies=None,
                     compute_gradients=0001):
    if (compute_gradients):
      self.gradients = flex.vec3_double(sites_cart.size(), [0,0,0])
    else:
      self.gradients = None
    if (bond_proxies is None):
      self.n_bond_proxies = 0
      self.bond_residual_sum = 0
    else:
      self.n_bond_proxies = len(bond_proxies)
      self.bond_residual_sum = restraints.bond_residual_sum(
        sites_cart=sites_cart,
        proxies=bond_proxies,
        gradient_array=self.gradients)
    if (bond_sym_proxies is None):
      self.n_bond_sym_proxies = 0
      self.bond_sym_residual_sum = 0
    else:
      assert asu_mappings is not None
      self.n_bond_sym_proxies = len(bond_sym_proxies)
      self.bond_sym_residual_sum=restraints.bond_residual_sum(
        sites_cart=sites_cart,
        asu_mappings=asu_mappings,
        proxies=bond_sym_proxies,
        gradient_array=self.gradients)
    if (bond_sorted_proxies is None):
      self.n_bond_sorted_proxies = 0
      self.bond_sorted_residual_sum = 0
    else:
      assert asu_mappings is not None
      self.n_bond_sorted_proxies = bond_sorted_proxies.n_total()
      self.bond_sorted_residual_sum=restraints.bond_residual_sum(
        sites_cart=sites_cart,
        sorted_proxies=bond_sorted_proxies,
        gradient_array=self.gradients)
    if (angle_proxies is None):
      self.n_angle_proxies = 0
      self.angle_residual_sum = 0
    else:
      self.n_angle_proxies = len(angle_proxies)
      self.angle_residual_sum = restraints.angle_residual_sum(
        sites_cart=sites_cart,
        proxies=angle_proxies,
        gradient_array=self.gradients)
    if (dihedral_proxies is None):
      self.n_dihedral_proxies = 0
      self.dihedral_residual_sum = 0
    else:
      self.n_dihedral_proxies = len(dihedral_proxies)
      self.dihedral_residual_sum = restraints.dihedral_residual_sum(
          sites_cart=sites_cart,
          proxies=dihedral_proxies,
          gradient_array=self.gradients)
    if (chirality_proxies is None):
      self.n_chirality_proxies = 0
      self.chirality_residual_sum = 0
    else:
      self.n_chirality_proxies = len(chirality_proxies)
      self.chirality_residual_sum = restraints.chirality_residual_sum(
          sites_cart=sites_cart,
          proxies=chirality_proxies,
          gradient_array=self.gradients)
    if (planarity_proxies is None):
      self.n_planarity_proxies = 0
      self.planarity_residual_sum = 0
    else:
      self.n_planarity_proxies = len(planarity_proxies)
      self.planarity_residual_sum = restraints.planarity_residual_sum(
          sites_cart=sites_cart,
          proxies=planarity_proxies,
          gradient_array=self.gradients)
    if (repulsion_proxies is None):
      self.n_repulsion_proxies = 0
      self.repulsion_residual_sum = 0
    else:
      assert asu_mappings is not None
      self.n_repulsion_proxies = repulsion_proxies.n_total()
      self.repulsion_residual_sum = restraints.repulsion_residual_sum(
        sites_cart=sites_cart,
        sorted_proxies=repulsion_proxies,
        gradient_array=self.gradients,
        function=restraints.repulsion_function(),
        disable_cache=00000)

  def target(self):
    return(self.bond_residual_sum
         + self.bond_sym_residual_sum
         + self.bond_sorted_residual_sum
         + self.angle_residual_sum
         + self.dihedral_residual_sum
         + self.chirality_residual_sum
         + self.planarity_residual_sum
         + self.repulsion_residual_sum)

  def gradient_norm(self):
    if (self.gradients is not None):
      return flex.sum_sq(self.gradients.as_double())

  def show(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "target:", self.target()
    print >> f, "  bond_residual_sum (n=%d):" % self.n_bond_proxies,\
      self.bond_residual_sum
    print >> f, "  bond_sym_residual_sum (n=%d):" % self.n_bond_sym_proxies,\
      self.bond_sym_residual_sum
    print >> f, "  bond_sorted_residual_sum (n=%d):" % \
      self.n_bond_sorted_proxies, self.bond_sorted_residual_sum
    print >> f, "  angle_residual_sum (n=%d):" % self.n_angle_proxies,\
      self.angle_residual_sum
    print >> f, "  dihedral_residual_sum (n=%d):" % self.n_dihedral_proxies,\
      self.dihedral_residual_sum
    print >> f, "  chirality_residual_sum (n=%d):" % self.n_chirality_proxies,\
      self.chirality_residual_sum
    print >> f, "  planarity_residual_sum (n=%d):" % self.n_planarity_proxies,\
      self.planarity_residual_sum
    print >> f, "  repulsion_residual_sum (n=%d):" % self.n_repulsion_proxies,\
      self.repulsion_residual_sum
    if (self.gradients is not None):
      print >> f, "  norm of gradients:", self.gradient_norm()

class lbfgs:

  def __init__(self, sites_cart,
                     site_symmetries=None,
                     asu_mappings=None,
                     bond_proxies=None,
                     bond_sym_proxies=None,
                     bond_sorted_proxies=None,
                     angle_proxies=None,
                     dihedral_proxies=None,
                     chirality_proxies=None,
                     planarity_proxies=None,
                     repulsion_proxies=None,
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
    if (self.site_symmetries is not None):
      sites_frac = self.site_symmetries[0].unit_cell() \
        .fractionalization_matrix() * self._sites_shifted
      sites_special = flex.vec3_double()
      for site_frac,site_symmetry in zip(sites_frac, self.site_symmetries):
        sites_special.append(site_symmetry.special_op()*site_frac)
        distance_moved = self.site_symmetries[0].unit_cell().distance(
          sites_special[-1], site_frac)
        if (distance_moved > 1.e-6):
          print "LARGE distance_moved:", distance_moved
      self._sites_shifted = self.site_symmetries[0].unit_cell() \
        .orthogonalization_matrix() * sites_special

  def compute_target(self, compute_gradients):
    self.time_compute_target.start()
    self.target_result = energies(
      sites_cart=self._sites_shifted,
      asu_mappings=self.asu_mappings,
      bond_proxies=self.bond_proxies,
      bond_sym_proxies=self.bond_sym_proxies,
      bond_sorted_proxies=self.bond_sorted_proxies,
      angle_proxies=self.angle_proxies,
      dihedral_proxies=self.dihedral_proxies,
      chirality_proxies=self.chirality_proxies,
      planarity_proxies=self.planarity_proxies,
      repulsion_proxies=self.repulsion_proxies,
      compute_gradients=compute_gradients)
    if (0):
      print "call:"
      self.target_result.show()
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
