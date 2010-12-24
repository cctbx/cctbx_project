from __future__ import division
from scitbx.lstbx import normal_eqns_solving
from cctbx import geometry_restraints, adp_restraints, sgtbx
from cctbx.array_family import flex
from cctbx.xray import parameter_map
from smtbx.refinement import restraints
from smtbx.refinement.restraints.adp_restraints import\
     adp_similarity_restraints, isotropic_adp_restraints, rigid_bond_restraints
import smtbx.utils
import smtbx.development
from smtbx.refinement import constraints, least_squares
from libtbx.test_utils import approx_equal
from libtbx.utils import wall_clock_time
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
    self.xray_structure = smtbx.development.sucrose()
    for sc in self.xray_structure.scatterers():
      sc.flags.set_grad_site(True)
      if sc.flags.use_u_aniso(): sc.flags.set_grad_u_aniso(True)
      if sc.flags.use_u_iso(): sc.flags.set_grad_u_iso(True)

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
    assert approx_equal(design_matrix, fd_design, 1e-4)
    assert approx_equal(
      linearised_eqns.n_restraints(),
      rows_per_restraint.get(self.proxies[0].__class__, 1) * self.proxies.size())

class geometry_restraints_test_case(restraints_test_case):

  def exercise_ls_restraints(self):
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
    xray_structure=smtbx.development.sucrose()).proxies
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
    xray_structure=smtbx.development.sucrose()).proxies
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
    xray_structure=smtbx.development.sucrose()).proxies
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

  i, j, k, l = random.sample(xrange(options.n_scatterers), 4)
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
    distance_ideal=d_jl,
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
      fo_sq=fo_sq,
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
    print "%i %s steps in %.6f s" % (cycles.n_iterations, cycles, t.elapsed())
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
    print "%i %s steps in %.6f s" % (cycles.n_iterations, cycles, t.elapsed())
  sc = ls.xray_structure.scatterers()
  sc = ls.xray_structure.scatterers()
  for p in bond_proxies:
    d = uc.distance(*[ sc[i].site for i in p.i_seqs ])
    assert approx_equal(d, p.distance_ideal, eps)

def exercise_ls_restraints(options):
  exercise_restrained_refinement(options)
  bond_restraint_test_case().run()
  angle_restraint_test_case().run()
  dihedral_restraint_test_case().run()
  isotropic_adp_test_case().run()
  adp_similarity_test_case().run()
  rigid_bond_test_case().run()

def run():
  libtbx.utils.show_times_at_exit()
  import sys
  from libtbx.option_parser import option_parser
  command_line = (option_parser()
    .option(None, "--verbose",
            action="store_true")
    .option(None, "--scatterers",
            dest='n_scatterers',
            type="int")
    .option(None, "--resolution",
            type="float")
  ).process(args=sys.argv[1:])
  exercise_ls_restraints(command_line.options)

if __name__ == '__main__':
  run()
