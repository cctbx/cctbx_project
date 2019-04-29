from __future__ import absolute_import, division, print_function
from cctbx import adp_restraints
from cctbx.array_family import flex
from libtbx import adopt_init_args
import math
import sys

class energies(object):

  def __init__(self, u_cart,
                     u_iso=None,
                     use_u_aniso=None,
                     sites_cart=None,
                     adp_similarity_proxies=None,
                     rigid_bond_proxies=None,
                     isotropic_adp_proxies=None,
                     compute_gradients=True,
                     gradients_aniso_cart=None,
                     gradients_iso=None,
                     disable_asu_cache=False,
                     normalization=False):
    adopt_init_args(self, locals())
    self.number_of_restraints = 0
    self.residual_sum = 0
    self.normalization_factor = None
    if (adp_similarity_proxies is not None):
      assert u_iso is not None and use_u_aniso is not None
    if (rigid_bond_proxies is not None): assert sites_cart is not None
    if (sites_cart is not None): assert sites_cart.size() == u_cart.size()
    if (u_iso is not None): assert u_iso.size() == u_cart.size()
    if (use_u_aniso is not None): assert use_u_aniso.size() == u_cart.size()
    if (compute_gradients):
      if (self.gradients_aniso_cart is None):
        self.gradients_aniso_cart = flex.sym_mat3_double(
          sites_cart.size(), [0,0,0,0,0,0])
      else:
        assert self.gradients_aniso_cart.size() == sites_cart.size()
      if (u_iso is not None and self.gradients_iso is None):
        self.gradients_iso = flex.double(sites_cart.size(), 0)
      elif (u_iso is not None):
        assert self.gradients_iso.size() == sites_cart.size()
    if (adp_similarity_proxies is None):
      self.n_adp_similarity_proxies = None
      self.adp_similarity_residual_sum = 0
    else:
      self.n_adp_similarity_proxies = len(adp_similarity_proxies)
      self.adp_similarity_residual_sum = adp_restraints.adp_similarity_residual_sum(
        u_cart=u_cart,
        u_iso=u_iso,
        use_u_aniso=use_u_aniso,
        proxies=adp_similarity_proxies,
        gradients_aniso_cart=self.gradients_aniso_cart,
        gradients_iso=self.gradients_iso)
      self.number_of_restraints += 6 * self.n_adp_similarity_proxies
      self.residual_sum += self.adp_similarity_residual_sum
    if (rigid_bond_proxies is None):
      self.n_rigid_bond_proxies = None
      self.rigid_bond_residual_sum = 0
    else:
      self.n_rigid_bond_proxies = len(rigid_bond_proxies)
      self.rigid_bond_residual_sum = adp_restraints.rigid_bond_residual_sum(
        sites_cart=sites_cart,
        u_cart=u_cart,
        proxies=rigid_bond_proxies,
        gradients_aniso_cart=self.gradients_aniso_cart)
      self.number_of_restraints += self.n_rigid_bond_proxies
      self.residual_sum += self.rigid_bond_residual_sum
    if (isotropic_adp_proxies is None):
      self.n_isotropic_adp_proxies = None
      self.isotropic_adp_residual_sum = 0
    else:
      self.n_isotropic_adp_proxies = len(isotropic_adp_proxies)
      self.isotropic_adp_residual_sum = adp_restraints.isotropic_adp_residual_sum(
        u_cart=u_cart,
        proxies=isotropic_adp_proxies,
        gradients_aniso_cart=self.gradients_aniso_cart)
      self.number_of_restraints += self.n_isotropic_adp_proxies
      self.residual_sum += self.isotropic_adp_residual_sum
    self.finalize_target_and_gradients()

  def adp_similarity_deviation(self):
    if (self.n_adp_similarity_proxies is not None):
      adp_similarity_deltas_rms = adp_restraints.adp_similarity_deltas_rms(
        u_cart=self.u_cart,
        u_iso=self.u_iso,
        use_u_aniso=self.use_u_aniso,
        proxies=self.adp_similarity_proxies)
      a_sq = adp_similarity_deltas_rms * adp_similarity_deltas_rms
      a_ave = math.sqrt(flex.mean_default(a_sq, 0))
      a_max = math.sqrt(flex.max_default(a_sq, 0))
      a_min = math.sqrt(flex.min_default(a_sq, 0))
      return a_min, a_max, a_ave

  def rigid_bond_deviation(self):
    if (self.n_rigid_bond_proxies is not None):
      rigid_bond_deltas = adp_restraints.rigid_bond_deltas(
        sites_cart=self.sites_cart,
        u_cart=self.u_cart,
        proxies=self.rigid_bond_proxies)
      r_sq = rigid_bond_deltas * rigid_bond_deltas
      r_ave = math.sqrt(flex.mean_default(r_sq, 0))
      r_max = math.sqrt(flex.max_default(r_sq, 0))
      r_min = math.sqrt(flex.min_default(r_sq, 0))
      return r_min, r_max, r_ave

  def isotropic_adp_deviation(self):
    if (self.n_isotropic_adp_proxies is not None):
      isotropic_adp_deltas_rms = adp_restraints.isotropic_adp_deltas_rms(
        u_cart=self.u_cart,
        proxies=self.isotropic_adp_proxies)
      i_sq = isotropic_adp_deltas_rms * isotropic_adp_deltas_rms
      i_ave = math.sqrt(flex.mean_default(i_sq, 0))
      i_max = math.sqrt(flex.max_default(i_sq, 0))
      i_min = math.sqrt(flex.min_default(i_sq, 0))
      return i_min, i_max, i_ave

  def show(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print(prefix+"target: %.6g" % self.target, file=f)
    if (self.n_adp_similarity_proxies is not None):
      print(prefix+"  adp_similarity_residual_sum (n=%d): %.6g" % (
        self.n_adp_similarity_proxies, self.adp_similarity_residual_sum), file=f)
    if (self.n_rigid_bond_proxies is not None):
      print(prefix+"  rigid_bond_residual_sum (n=%d): %.6g" % (
        self.n_rigid_bond_proxies, self.rigid_bond_residual_sum), file=f)
    if (self.n_isotropic_adp_proxies is not None):
      print(prefix+"  isotropic_adp_residual_sum (n=%d): %.6g" % (
        self.n_isotropic_adp_proxies, self.isotropic_adp_residual_sum), file=f)

  def finalize_target_and_gradients(self):
    self.target = self.residual_sum
    if (self.normalization):
      self.normalization_factor = 1.0 / max(1, self.number_of_restraints)
      self.target *= self.normalization_factor
      if (self.gradients_aniso_cart is not None):
        self.gradients_aniso_cart *= self.normalization_factor
      if (self.gradients_iso is not None):
        self.gradients_iso *= self.normalization_factor
