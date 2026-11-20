from __future__ import absolute_import, division, print_function
from scitbx.lstbx import normal_eqns_solving
from cctbx import geometry_restraints, adp_restraints, sgtbx, adptbx
from cctbx.array_family import flex
from cctbx.xray import parameter_map
from smtbx.refinement import restraints
from smtbx.refinement.restraints.adp_restraints import\
     adp_similarity_restraints, isotropic_adp_restraints,\
     rigid_bond_restraints, fixed_u_eq_adp_restraints,\
     adp_u_eq_similarity_restraints, adp_volume_similarity_restraints

from cctbx.adp_restraints import adp_restraint_params
import smtbx.utils
import smtbx.development
from smtbx.refinement import constraints, least_squares
from libtbx.test_utils import approx_equal
from libtbx.utils import wall_clock_time
import libtbx
from scitbx import matrix
from six.moves import range

geom = geometry_restraints
adp = adp_restraints

rows_per_restraint = {
  geom.bond_similarity_proxy: 6,
  adp.adp_similarity_proxy: 6,
  adp.isotropic_adp_proxy: 6,
  }

def set_gradient_flags(xray_structure):
  """ Refine everything, which is what these test cases are about.

  It matters to the proxy builders and not only to the refinement: the ADP
  restraint builders skip a scatterer whose ADPs are not being refined, so a
  structure whose flags are still at their defaults yields no ADP proxies at
  all. The test cases build their proxies in the class body, before any
  instance exists, so the structure they build them from has to be flagged
  here rather than in __init__.
  """
  for sc in xray_structure.scatterers():
    sc.flags.set_grad_site(True)
    if sc.flags.use_u_aniso(): sc.flags.set_grad_u_aniso(True)
    if sc.flags.use_u_iso(): sc.flags.set_grad_u_iso(True)
  return xray_structure


def refinable_sucrose():
  return set_gradient_flags(smtbx.development.sucrose())


class restraints_test_case:

  def __init__(self):
    self.xray_structure = set_gradient_flags(smtbx.development.sucrose())
    self.tolerance = 1e-4

    self.param_map = parameter_map(self.xray_structure.scatterers())
    assert self.proxies.size() > 0

  def run(self):
    self.exercise_ls_restraints()

  def exercise_ls_restraints(self):
    xs = self.xray_structure.deep_copy_scatterers()
    linearised_eqns = self.manager.build_linearised_eqns(xs, xs.parameter_map())
    design_matrix = linearised_eqns.design_matrix.as_dense_matrix()
    fd_design = flex.double()
    for proxy in self.proxies:
      grads = self.fd_grads(proxy)
      for i, grad in enumerate(grads):
        fd_design.extend(grad)
    assert approx_equal(design_matrix, fd_design, self.tolerance)

class geometry_restraints_test_case(restraints_test_case):

  def exercise_ls_restraints(self):
    match = restraints_test_case.exercise_ls_restraints(self)

  def fd_grads(self, proxy):
    grads = flex.double(self.param_map.n_parameters)
    eps = 1e-8
    uc = self.xray_structure.unit_cell()
    sites_cart = self.xray_structure.sites_cart().deep_copy()
    for i in range(self.param_map.n_scatterers):
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

class chirality_restraint_test_case(geometry_restraints_test_case):
  manager = restraints.manager(
    chirality_proxies = geometry_restraints.shared_chirality_proxy([
      geom.chirality_proxy((0, 19, 30, 21), volume_ideal=2.5,
                           both_signs=False, weight=1),
      geom.chirality_proxy((1, 2, 21, 22), volume_ideal=0.0,
                           both_signs=False, weight=1)
    ]))
  proxies = manager.chirality_proxies
  restraint_t = geom.chirality

def exercise_coincident_bond_similarity():
  """A SADI one of whose pairs has collapsed to a point.

  bond_similarity divides by each bond length to get the direction to restrain
  along. Two atoms of a pair at the same point made that 0/0, and six NaN went
  into the design matrix - which then spread through the normal matrix, so the
  Cholesky failure named an unrelated atom. Same shape as the coplanar FLAT
  case, in a much commoner restraint.

  Unlike the dihedral, nothing diverges on the way in: the numerator carries a
  factor of the bond vector, so it cancels the length and the entries stay at
  25 down to 1e-12 A. Only the exact coincidence is a problem, so the check is
  simply that it produces no NaN and still restrains the pair that is fine.
  """
  from cctbx import crystal, xray
  cs = crystal.symmetry(unit_cell=(50, 50, 50, 90, 90, 90),
                        space_group_symbol="P1")

  def row(sep):
    xs = xray.structure(crystal_symmetry=cs)
    for i, s in enumerate([(0., 0., 0.), (1.5, 0., 0.),
                           (5., 0., 0.), (5. + sep, 0., 0.)]):
      sc = xray.scatterer(label="C%d" % i, scattering_type="C",
                          site=[c / 50. for c in s])
      sc.flags.set_grad_site(True)
      xs.add_scatterer(sc)
    proxy = geometry_restraints.bond_similarity_proxy(
      i_seqs=[(0, 1), (2, 3)], weights=(1.0, 1.0))
    mgr = restraints.manager(
      bond_similarity_proxies=
        geometry_restraints.shared_bond_similarity_proxy([proxy]))
    eqns = mgr.build_linearised_eqns(xs, xs.parameter_map())
    return list(eqns.design_matrix.as_dense_matrix())

  for sep in (1.5, 1e-6, 1e-12, 0.):
    r = row(sep)
    assert [v for v in r if v != v] == [], "NaN at separation %g" % sep
    assert [v for v in r if abs(v) == float("inf")] == [], (
      "inf at separation %g" % sep)
    # the intact pair still has to be restrained, or the guard has thrown the
    # whole restraint away rather than the one direction it cannot define
    assert max(abs(v) for v in r) > 1e-6, (
      "separation %g left nothing of the restraint" % sep)


def exercise_degenerate_dihedral():
  """A torsion about an axis a terminal atom is sitting on.

  The row of a dihedral restraint scales as 1/(perpendicular distance of the
  terminal atom from the central bond axis), so it does not converge to zero
  at the degeneracy the way the comment in dihedral.h used to claim - it
  diverges. Measured on the linearised row, 5.7e3 for an ordinary geometry
  against 5.7e11 with a terminal atom 1e-8 A off the axis. Nothing is NaN, so
  nothing complains; the normal matrix is simply singular and the Cholesky
  failure names an unrelated atom.

  Two things have to hold. Where the geometry is sound the row must still be
  the true derivative, which is checked against finite differences of delta -
  a guard that quietly changed the mathematics would be worse than the
  divergence. Where it is degenerate the row must be refused outright rather
  than returned enormous.
  """
  from cctbx import crystal, xray
  cs = crystal.symmetry(unit_cell=(50, 50, 50, 90, 90, 90),
                        space_group_symbol="P1")

  def row_for(offset):
    xs = xray.structure(crystal_symmetry=cs)
    sites = [(0., offset, 0.), (0., 0., 0.), (1., 0., 0.), (1., offset, 0.5)]
    for i, s in enumerate(sites):
      sc = xray.scatterer(label="C%d" % i, scattering_type="C",
                          site=[c / 50. for c in s])
      sc.flags.set_grad_site(True)
      xs.add_scatterer(sc)
    proxy = geometry_restraints.dihedral_proxy(
      (0, 1, 2, 3), angle_ideal=0., weight=100)
    mgr = restraints.manager(
      dihedral_proxies=geometry_restraints.shared_dihedral_proxy([proxy]))
    eqns = mgr.build_linearised_eqns(xs, xs.parameter_map())
    return xs, proxy, list(eqns.design_matrix.as_dense_matrix())

  # sound geometry: the row is the derivative of delta, checked numerically
  xs, proxy, row = row_for(1.0)
  assert [v for v in row if v != v] == [], "NaN in a well conditioned row"
  assert max(abs(v) for v in row) > 1e-6, "a sound dihedral gave a zero row"

  eps = 1e-6
  sites_cart = xs.sites_cart().deep_copy()
  fd = []
  for i in range(len(sites_cart)):
    for j in range(3):
      deltas = []
      for sign in (1, -1):
        sc = sites_cart.deep_copy()
        s = list(sc[i]); s[j] += sign * eps; sc[i] = s
        deltas.append(geometry_restraints.dihedral(
          sites=[sc[k] for k in proxy.i_seqs], angle_ideal=0.,
          weight=100, periodicity=1).delta)
      fd.append((deltas[0] - deltas[1]) / (2 * eps))
  # the design matrix is in fractional coordinates, the finite differences in
  # Cartesian, so compare after the same change of basis the linearisation does
  orth = cs.unit_cell().orthogonalization_matrix()
  fd_frac = []
  for i in range(len(sites_cart)):
    g = fd[3 * i:3 * i + 3]
    fd_frac.extend([orth[0] * g[0], orth[4] * g[1], orth[8] * g[2]])
  assert approx_equal(row, fd_frac, 1e-3), (row, fd_frac)

  # degenerate geometry: refused, not enormous
  for offset in (1e-4, 1e-6, 1e-8, 0.):
    _, _, row = row_for(offset)
    assert [v for v in row if v != v] == [], "NaN at offset %g" % offset
    assert max(abs(v) for v in row) == 0, (
      "offset %g gave a row of %g, which the normal matrix cannot carry"
      % (offset, max(abs(v) for v in row)))


def exercise_coplanar_chirality():
  """A chirality restraint whose sites are exactly coplanar.

  FLAT reaches the refinement as chirality restraints of zero ideal volume,
  and atoms on a mirror are coplanar by symmetry, so the volume and every
  gradient are exactly zero. Recovering the design matrix row by dividing that
  out was 0/0, and the NaN reached the normal matrix as a Cholesky failure
  blaming an unrelated parameter.
  """
  from cctbx import crystal, xray
  from scitbx import matrix
  cs = crystal.symmetry(unit_cell=(8.32, 6.2744, 20.6559, 90, 90, 90),
                        space_group_symbol="P n m a")
  xs = xray.structure(crystal_symmetry=cs)
  for label, site in [("O1",  (0.511044, 0.75, 0.348495)),
                      ("N1",  (0.531067, 0.75, 0.480322)),
                      ("H1a", (0.521941, 0.75, 0.393406)),
                      ("C8",  (0.355582, 0.75, 0.333011))]:
    sc = xray.scatterer(label=label, site=site)
    sc.flags.set_grad_site(True)
    xs.add_scatterer(sc)

  proxy = geom.chirality_proxy((1, 3, 2, 0), volume_ideal=0.0,
                               both_signs=False, weight=100)
  uc = xs.unit_cell()
  r = geom.chirality(uc, xs.sites_cart(), proxy)
  assert r.volume_model == 0 and r.delta == 0, (r.volume_model, r.delta)

  mgr = restraints.manager(
    chirality_proxies=geometry_restraints.shared_chirality_proxy([proxy]))
  eqns = mgr.build_linearised_eqns(xs, xs.parameter_map())
  row = list(eqns.design_matrix.as_dense_matrix())
  assert [v for v in row if v != v] == [], "NaN in the design matrix row"
  assert max(abs(v) for v in row) > 1e-6, "row is all zero"

  # against finite differences of delta, which stay well defined at the
  # degeneracy even though the analytic scaling did not
  eps = 1e-8
  sites_cart = xs.sites_cart().deep_copy()
  fd = flex.double(xs.parameter_map().n_parameters)
  pm = xs.parameter_map()
  for i in range(pm.n_scatterers):
    g = [0, 0, 0]
    for j in range(3):
      h = [0, 0, 0]
      h[j] = eps
      sites_cart[i] = matrix.col(sites_cart[i]) + matrix.col(h)
      d1 = geom.chirality(uc, sites_cart, proxy).delta
      sites_cart[i] = matrix.col(sites_cart[i]) - 2*matrix.col(h)
      d2 = geom.chirality(uc, sites_cart, proxy).delta
      sites_cart[i] = matrix.col(sites_cart[i]) + matrix.col(h)
      g[j] = (d1 - d2)/(2*eps)
    gf = uc.fractionalize_gradient(g)
    for j in range(3):
      fd[pm[i].site + j] = gf[j]
  assert approx_equal(row, list(fd), 1e-4), (row, list(fd))

class adp_restraints_test_case(restraints_test_case):

  def __init__(self):
    restraints_test_case.__init__(self)

  def fd_grads(self, proxy):
    dynamic_restraint_proxy_classes = (
      adp.adp_u_eq_similarity_proxy,
      adp.adp_volume_similarity_proxy,
    )
    if isinstance(proxy, (dynamic_restraint_proxy_classes)):
      n_restraints = len(proxy.i_seqs)
    else:
      n_restraints = rows_per_restraint.get(proxy.__class__, 1)
    grads = [flex.double(self.param_map.n_parameters) for i in range(n_restraints)]
    eps = 1e-8
    uc = self.xray_structure.unit_cell()
    xs = self.xray_structure
    u_cart = xs.scatterers().extract_u_cart(uc).deep_copy()
    u_star = xs.scatterers().extract_u_star().deep_copy()
    u_iso = xs.scatterers().extract_u_iso().deep_copy()
    single_delta_classes = (
      adp.fixed_u_eq_adp,
    )
    for n in range(n_restraints):
      for i in range(self.param_map.n_scatterers):
        use_u_aniso = self.param_map[i].u_aniso > -1
        use_u_iso = self.param_map[i].u_iso > -1
        for j in range(6):
          if use_u_aniso:
            h = [0,0,0,0,0,0]
            h[j] = eps
            h = matrix.sym(sym_mat3=h)
            u_star[i]=list((matrix.sym(sym_mat3=u_star[i]) + h).as_sym_mat3())
            r = self.restraint(proxy, u_cart=flex.sym_mat3_double([
              adptbx.u_star_as_u_cart(uc, u) for u in u_star]))
            if isinstance(r, adp.rigid_bond):
              d1 = r.delta_z()
            elif isinstance(r, single_delta_classes):
              d1 = r.delta()
            else:
              d1 = r.deltas()[n]
            u_star[i]=list((matrix.sym(sym_mat3=u_star[i]) - 2*h).as_sym_mat3())
            r = self.restraint(proxy, u_cart=flex.sym_mat3_double([
              adptbx.u_star_as_u_cart(uc, u) for u in u_star]))
            if isinstance(r, adp.rigid_bond):
              d2 = r.delta_z()
            elif isinstance(r, single_delta_classes):
              d2 = r.delta()
            else:
              d2 = r.deltas()[n]
          elif use_u_iso:
            u_iso[i] += eps
            r = self.restraint(proxy, u_iso=u_iso)
            if isinstance(r, adp.rigid_bond):
              d1 = r.delta_z()
            elif isinstance(r, single_delta_classes):
              d1 = r.delta()
            else:
              d1 = r.deltas()[n]
            u_iso[i] -= 2*eps
            r = self.restraint(proxy, u_iso=u_iso)
            if isinstance(r, adp.rigid_bond):
              d2 = r.delta_z()
            elif isinstance(r, single_delta_classes):
              d2 = r.delta()
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
    xray_structure=refinable_sucrose()).proxies
  # no need to test all of them every time
  proxies = adp.shared_isotropic_adp_proxy(
    flex.select(proxies, flags=flex.random_bool(proxies.size(), 0.5)))
  manager = restraints.manager(isotropic_adp_proxies=proxies)

  def restraint(self, proxy, u_iso=None, u_cart=None):
    if u_cart is None:
      u_cart=self.xray_structure.scatterers().extract_u_cart(
        self.xray_structure.unit_cell())
    return adp.isotropic_adp(
      adp_restraint_params(u_cart=u_cart),
      proxy)

class fixed_u_eq_adp_test_case(adp_restraints_test_case):
  proxies = fixed_u_eq_adp_restraints(
    xray_structure=refinable_sucrose(),
    u_eq_ideal=0.025).proxies
  # no need to test all of them every time
  proxies = adp.shared_fixed_u_eq_adp_proxy(
    flex.select(proxies, flags=flex.random_bool(proxies.size(), 0.5)))
  manager = restraints.manager(fixed_u_eq_adp_proxies=proxies)

  def restraint(self, proxy, u_iso=None, u_cart=None):
    if u_cart is None:
      u_cart=self.xray_structure.scatterers().extract_u_cart(
        self.xray_structure.unit_cell())
    if u_iso is None:
      u_iso=self.xray_structure.scatterers().extract_u_iso()
    use_u_aniso=self.xray_structure.use_u_aniso()
    return adp.fixed_u_eq_adp(
      adp_restraint_params(u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso),
      proxy)

class adp_similarity_test_case(adp_restraints_test_case):
  proxies = adp_similarity_restraints(
    xray_structure=refinable_sucrose()).proxies
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
    return adp.adp_similarity(
      adp_restraint_params(u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso),
      proxy)

class adp_u_eq_similarity_test_case(adp_restraints_test_case):
  proxies = adp_u_eq_similarity_restraints(
    xray_structure=refinable_sucrose()).proxies
  # no need to test all of them every time
  #proxies = adp.shared_adp_u_eq_similarity_proxy(
    #flex.select(proxies, flags=flex.random_bool(proxies.size(), 0.5)))
  manager = restraints.manager(adp_u_eq_similarity_proxies=proxies)

  def restraint(self, proxy, u_iso=None, u_cart=None):
    if u_cart is None:
      u_cart=self.xray_structure.scatterers().extract_u_cart(
        self.xray_structure.unit_cell())
    if u_iso is None:
      u_iso=self.xray_structure.scatterers().extract_u_iso()
    use_u_aniso=self.xray_structure.use_u_aniso()
    return adp.adp_u_eq_similarity(
      adp_restraint_params(u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso),
      proxy)

class adp_volume_similarity_test_case(adp_restraints_test_case):
  proxies = adp_volume_similarity_restraints(
    xray_structure=refinable_sucrose()).proxies
  manager = restraints.manager(adp_volume_similarity_proxies=proxies)
  def __init__(self):
    adp_restraints_test_case.__init__(self)
    # eigen values and eigen vectors are dependent after all...
    # may need to make smaller
    self.tolerance = 0.3
  def restraint(self, proxy, u_iso=None, u_cart=None):
    if u_cart is None:
      u_cart=self.xray_structure.scatterers().extract_u_cart(
        self.xray_structure.unit_cell())
    if u_iso is None:
      u_iso=self.xray_structure.scatterers().extract_u_iso()
    use_u_aniso=self.xray_structure.use_u_aniso()
    return adp.adp_volume_similarity(
      adp_restraint_params(u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso),
      proxy)

class rigid_bond_test_case(adp_restraints_test_case):
  proxies = rigid_bond_restraints(
    xray_structure=refinable_sucrose()).proxies
  # no need to test all of them every time
  proxies = adp.shared_rigid_bond_proxy(
    flex.select(proxies, flags=flex.random_bool(proxies.size(), 0.3)))
  manager = restraints.manager(rigid_bond_proxies=proxies)

  def restraint(self, proxy, u_iso=None, u_cart=None):
    if u_cart is None:
      u_cart = self.xray_structure.scatterers().extract_u_cart(
        self.xray_structure.unit_cell())
    sites_cart = self.xray_structure.sites_cart()
    return adp.rigid_bond(
      adp_restraint_params(sites_cart=sites_cart, u_cart=u_cart),
      proxy)

def exercise_restrained_refinement(options):
  import random
  random.seed(1)
  flex.set_random_seed(1)
  xs0 = smtbx.development.random_xray_structure(
    sgtbx.space_group_info('P1'),
    n_scatterers=options.n_scatterers,
    elements="random")
  for sc in xs0.scatterers():
    sc.flags.set_grad_site(True)
  sc0 = xs0.scatterers()
  uc = xs0.unit_cell()

  mi = xs0.build_miller_set(anomalous_flag=False, d_min=options.resolution)
  fo_sq = mi.structure_factors_from_scatterers(
    xs0, algorithm="direct").f_calc().norm()
  fo_sq = fo_sq.customized_copy(sigmas=flex.double(fo_sq.size(), 1))

  i, j, k, l = random.sample(range(options.n_scatterers), 4)
  bond_proxies = geometry_restraints.shared_bond_simple_proxy()
  w = 1e9
  d_ij = uc.distance(sc0[i].site, sc0[j].site)*0.8
  bond_proxies.append(geom.bond_simple_proxy(
    i_seqs=(i, j),
    distance_ideal=d_ij,
    weight=w))
  d_jk = uc.distance(sc0[j].site, sc0[k].site)*0.85
  bond_proxies.append(geom.bond_simple_proxy(
    i_seqs=(j, k),
    distance_ideal=d_jk,
    weight=w))
  d_ki = min(uc.distance(sc0[k].site, sc0[i].site)*0.9, (d_ij + d_jk)*0.8)
  bond_proxies.append(geom.bond_simple_proxy(
    i_seqs=(k, i),
    distance_ideal=d_ki,
    weight=w))
  d_jl = uc.distance(sc0[j].site, sc0[l].site)*0.9
  bond_proxies.append(geom.bond_simple_proxy(
    i_seqs=(j, l),
    distance_ideal=d_jl,
    weight=w))
  d_lk = min(uc.distance(sc0[l].site, sc0[k].site)*0.8, 0.75*(d_jk + d_jl))
  bond_proxies.append(geom.bond_simple_proxy(
    i_seqs=(l, k),
    distance_ideal=d_lk,
    weight=w))
  restraints_manager = restraints.manager(bond_proxies=bond_proxies)

  xs1 = xs0.deep_copy_scatterers()
  xs1.shake_sites_in_place(rms_difference=0.1)

  def ls_problem():
    xs = xs1.deep_copy_scatterers()
    reparametrisation = constraints.reparametrisation(
      structure=xs,
      constraints=[],
      connectivity_table=smtbx.utils.connectivity_table(xs),
      temperature=20)
    return least_squares.crystallographic_ls(
      fo_sq.as_xray_observations(),
      reparametrisation=reparametrisation,
      restraints_manager=restraints_manager)

  gradient_threshold, step_threshold = 1e-6, 1e-6
  eps = 5e-3

  ls = ls_problem()
  t = wall_clock_time()
  cycles = normal_eqns_solving.naive_iterations(
    ls,
    gradient_threshold=gradient_threshold,
    step_threshold=step_threshold,
    track_all=True)
  if options.verbose:
    print("%i %s steps in %.6f s" % (cycles.n_iterations, cycles, t.elapsed()))
  sc = ls.xray_structure.scatterers()
  for p in bond_proxies:
    d = uc.distance(*[ sc[i_pair].site for i_pair in p.i_seqs ])
    assert approx_equal(d, p.distance_ideal, eps)

  ls = ls_problem()
  t = wall_clock_time()
  cycles = normal_eqns_solving.levenberg_marquardt_iterations(
    ls,
    gradient_threshold=gradient_threshold,
    step_threshold=step_threshold,
    tau=1e-3,
    track_all=True)
  if options.verbose:
    print("%i %s steps in %.6f s" % (cycles.n_iterations, cycles, t.elapsed()))
  sc = ls.xray_structure.scatterers()
  sc = ls.xray_structure.scatterers()
  for p in bond_proxies:
    d = uc.distance(*[ sc[i].site for i in p.i_seqs ])
    assert approx_equal(d, p.distance_ideal, eps)

def exercise_add_equation():
  linearised_eqns = restraints.linearised_eqns_of_restraint(10, 10)
  delta = 0.5
  grads = flex.double((0,0,1,0,0,2,0,0,-1, 0))
  w = 10
  linearised_eqns.add_equation(delta, grads, w)
  assert linearised_eqns.n_restraints() == 1
  linearised_eqns.add_equation(delta, grads, w)
  linearised_eqns.add_equation(delta, grads, w)
  assert linearised_eqns.n_restraints() == 3
  from scitbx import sparse
  assert approx_equal(
    linearised_eqns.design_matrix.as_dense_matrix(),
    sparse.matrix(rows=10, columns=10,
                  elements_by_columns=[ { 0: 0, 1: 0, 2: 0 },
                                        { 0: 0, 1: 0, 2: 0 },
                                        { 0: 1, 1: 1, 2: 1 },
                                        { 0: 0, 1: 0, 2: 0 },
                                        { 0: 0, 1: 0, 2: 0 },
                                        { 0: 2, 1: 2, 2: 2 },
                                        { 0: 0, 1: 0, 2: 0 },
                                        { 0: 0, 1: 0, 2: 0 },
                                        { 0: -1, 1: -1, 2: -1 },
                                        { 0: 0, 1: 0, 2: 0 }, ]).as_dense_matrix())


def exercise_ls_restraints(options):
  exercise_add_equation()
  exercise_restrained_refinement(options)
  bond_restraint_test_case().run()
  angle_restraint_test_case().run()
  dihedral_restraint_test_case().run()
  chirality_restraint_test_case().run()
  exercise_coplanar_chirality()
  exercise_degenerate_dihedral()
  exercise_coincident_bond_similarity()

  isotropic_adp_test_case().run()
  adp_similarity_test_case().run()
  rigid_bond_test_case().run()
  fixed_u_eq_adp_test_case().run()
  adp_u_eq_similarity_test_case().run()
  adp_volume_similarity_test_case().run()

def run():
  libtbx.utils.show_times_at_exit()
  import sys
  from libtbx.option_parser import option_parser
  command_line = (option_parser()
    .option(None, "--verbose",
            action="store_true")
    .option(None, "--scatterers",
            dest='n_scatterers',
            type="int",
            default=5)
    .option(None, "--resolution",
            type="float",
            default=0.2)
  ).process(args=sys.argv[1:])
  exercise_ls_restraints(command_line.options)

if __name__ == '__main__':
  run()
