from __future__ import division
from cctbx import geometry_restraints, adp_restraints
from cctbx.array_family import flex
from cctbx.xray import parameter_map
from smtbx.refinement import restraints
from smtbx.refinement.restraints.adp_restraints import\
     adp_similarity_restraints, isotropic_adp_restraints, rigid_bond_restraints
from smtbx.refinement.restraints.tests import trial_structure
from libtbx.test_utils import approx_equal
import libtbx

from scitbx import matrix

geom = geometry_restraints
adp = adp_restraints

rows_per_restraint = {
  geom.bond_similarity_proxy: 6,
  adp.adp_similarity_proxy: 6,
  adp.isotropic_adp_proxy: 6,
  }

class restraints_test_case:

  def __init__(self):
    xs = trial_structure()
    self.xray_structure = xs
    self.param_map = parameter_map(self.xray_structure.scatterers())
    assert self.proxies.size() > 0

  def run(self):
    self.exercise_ls_restraints()

  def exercise_ls_restraints(self):
    xs = self.xray_structure.deep_copy_scatterers()
    linearised_eqns = self.manager.build_linearised_eqns(xs)
    design_matrix = linearised_eqns.design_matrix.as_dense_matrix()
    fd_design = flex.double()
    for proxy in self.proxies:
      grads = self.fd_grads(proxy)
      for i, grad in enumerate(grads):
        fd_design.extend(grad)
    assert approx_equal(design_matrix, fd_design, 1e-5)
    assert approx_equal(
      linearised_eqns.n_restraints(),
      rows_per_restraint.get(self.proxies[0].__class__, 1) * self.proxies.size())


class geometry_restraints_test_case(restraints_test_case):

  def exercise_ls_restraints(self):
    xs = self.xray_structure.deep_copy_scatterers()
    match = restraints_test_case.exercise_ls_restraints(self)

  def fd_grads(self, proxy):
    grads = flex.double(self.param_map.n_parameters)
    eps = 1e-8
    uc = self.xray_structure.unit_cell()
    sites_cart = self.xray_structure.sites_cart().deep_copy()
    for i in xrange(self.param_map.n_scatterers):
      grad_site_cart = [0,0,0]
      for j in range(3):
        h = [0,0,0]
        h[j] = eps
        h = matrix.col(h)
        sites_cart[i] = matrix.col(sites_cart[i])+h
        r = self.restraint_t(uc, sites_cart, proxy)
        d1 = r.delta
        sites_cart[i] = matrix.col(sites_cart[i])-2*h
        r = self.restraint_t(uc, sites_cart, proxy)
        d2 = r.delta
        d_delta = (d1-d2)/(2*eps)
        grad_site_cart[j] = d_delta
      grad_site_frac = uc.fractionalize_gradient(grad_site_cart)
      for j in range(3):
        grads[self.param_map[i].site+j] = grad_site_frac[j]
    return [grads]

class bond_restraint_test_case(geometry_restraints_test_case):
  manager = restraints.manager(
    bond_proxies = geometry_restraints.shared_bond_simple_proxy([
      geom.bond_simple_proxy((0,30), 1.42, 1),
      geom.bond_simple_proxy((1,21), 1.42, 1)
    ]))
  proxies = manager.bond_proxies
  restraint_t = geom.bond

class angle_restraint_test_case(geometry_restraints_test_case):
  manager = restraints.manager(
    angle_proxies = geometry_restraints.shared_angle_proxy([
      geom.angle_proxy((30, 0, 19), 115, 1),
      geom.angle_proxy((21, 1, 2), 110, 1)
    ]))
  proxies = manager.angle_proxies
  restraint_t = geom.angle

class dihedral_restraint_test_case(geometry_restraints_test_case):
  manager = restraints.manager(
    dihedral_proxies = geometry_restraints.shared_dihedral_proxy([
      geom.dihedral_proxy((21, 19, 24, 26), 180, 1),
      geom.dihedral_proxy((5, 26, 28, 7), 60, 1)
    ]))
  proxies = manager.dihedral_proxies
  restraint_t = geom.dihedral

class adp_restraints_test_case(restraints_test_case):

  def __init__(self):
    restraints_test_case.__init__(self)

  def fd_grads(self, proxy):
    n_restraints = rows_per_restraint.get(proxy.__class__, 1)
    grads = [flex.double(self.param_map.n_parameters) for i in range(n_restraints)]
    eps = 1e-8
    uc = self.xray_structure.unit_cell()
    xs = self.xray_structure
    u_cart = xs.scatterers().extract_u_cart(uc).deep_copy()
    u_iso = xs.scatterers().extract_u_iso().deep_copy()
    for n in xrange(n_restraints):
      for i in xrange(self.param_map.n_scatterers):
        use_u_aniso = self.param_map[i].u_aniso != -1
        use_u_iso = self.param_map[i].u_iso != -1
        for j in range(6):
          if use_u_aniso:
            h = [0,0,0,0,0,0]
            h[j] = eps
            h = matrix.sym(sym_mat3=h)
            u_cart[i]=list((matrix.sym(sym_mat3=u_cart[i]) + h).as_sym_mat3())
            r = self.restraint(proxy, u_cart=u_cart)
            if isinstance(r, adp.rigid_bond):
              d1 = r.delta_z()
            else:
              d1 = r.deltas()[n]
            u_cart[i]=list((matrix.sym(sym_mat3=u_cart[i]) - 2*h).as_sym_mat3())
            r = self.restraint(proxy, u_cart=u_cart)
            if isinstance(r, adp.rigid_bond):
              d2 = r.delta_z()
            else:
              d2 = r.deltas()[n]
          elif use_u_iso:
            u_iso[i] += eps
            r = self.restraint(proxy, u_iso=u_iso)
            if isinstance(r, adp.rigid_bond):
              d1 = r.delta_z()
            else:
              d1 = r.deltas()[n]
            u_iso[i] -= 2*eps
            r = self.restraint(proxy, u_iso=u_iso)
            if isinstance(r, adp.rigid_bond):
              d2 = r.delta_z()
            else:
              d2 = r.deltas()[n]
          d_delta = (d1-d2)/(2*eps)
          if not isinstance(r, adp.rigid_bond) and j > 2:
            d_delta *= 2 # off diagonals count twice
          if use_u_aniso:
            grads[n][self.param_map[i].u_aniso+j] = d_delta
          elif use_u_iso:
            grads[n][self.param_map[i].u_iso] = d_delta
            break
    return grads

class isotropic_adp_test_case(adp_restraints_test_case):
  proxies = isotropic_adp_restraints(
    xray_structure=trial_structure()).proxies
  # no need to test all of them every time
  proxies = adp.shared_isotropic_adp_proxy(
    flex.select(proxies, flags=flex.random_bool(proxies.size(), 0.5)))
  manager = restraints.manager(isotropic_adp_proxies=proxies)

  def restraint(self, proxy, u_iso=None, u_cart=None):
    if u_cart is None:
      u_cart=self.xray_structure.scatterers().extract_u_cart(
        self.xray_structure.unit_cell())
    return adp.isotropic_adp(u_cart, proxy)

class adp_similarity_test_case(adp_restraints_test_case):
  proxies = adp_similarity_restraints(
    xray_structure=trial_structure()).proxies
  # no need to test all of them every time
  proxies = adp.shared_adp_similarity_proxy(
    flex.select(proxies, flags=flex.random_bool(proxies.size(), 0.5)))
  manager = restraints.manager(adp_similarity_proxies=proxies)

  def restraint(self, proxy, u_iso=None, u_cart=None):
    if u_cart is None:
      u_cart=self.xray_structure.scatterers().extract_u_cart(
        self.xray_structure.unit_cell())
    if u_iso is None:
      u_iso=self.xray_structure.scatterers().extract_u_iso()
    use_u_aniso=self.xray_structure.use_u_aniso()
    return adp.adp_similarity(u_cart, u_iso, use_u_aniso, proxy)

class rigid_bond_test_case(adp_restraints_test_case):
  proxies = rigid_bond_restraints(
    xray_structure=trial_structure()).proxies
  # no need to test all of them every time
  proxies = adp.shared_rigid_bond_proxy(
    flex.select(proxies, flags=flex.random_bool(proxies.size(), 0.3)))
  manager = restraints.manager(rigid_bond_proxies=proxies)

  def restraint(self, proxy, u_iso=None, u_cart=None):
    if u_cart is None:
      u_cart = self.xray_structure.scatterers().extract_u_cart(
        self.xray_structure.unit_cell())
    sites_cart = self.xray_structure.sites_cart()
    return adp.rigid_bond(sites_cart, u_cart, proxy)

def exercise_ls_restraints():
  bond_restraint_test_case().run()
  angle_restraint_test_case().run()
  dihedral_restraint_test_case().run()
  isotropic_adp_test_case().run()
  adp_similarity_test_case().run()
  rigid_bond_test_case().run()

def run():
  libtbx.utils.show_times_at_exit()
  exercise_ls_restraints()

if __name__ == '__main__':
  run()
