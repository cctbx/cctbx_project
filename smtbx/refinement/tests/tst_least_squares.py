from __future__ import division
from scitbx.linalg import eigensystem, svd
from scitbx import matrix
from cctbx import sgtbx, crystal, xray, miller, adptbx, uctbx
from cctbx import euclidean_model_matching as emma
from cctbx.array_family import flex
from smtbx.refinement import least_squares
from smtbx.refinement import constraints
import smtbx.utils
from libtbx.test_utils import approx_equal
import libtbx.utils
from stdlib import math
import sys
import random


class refinement_test(object):

  ls_cycle_repeats = 1

  def run(self):
    if self.ls_cycle_repeats == 1:
      self.exercise_ls_cycles()
    else:
      print "%s in %s" % (self.purpose, self.hall)
      for n in xrange(self.ls_cycle_repeats):
        self.exercise_ls_cycles()
        print '.',
        sys.stdout.flush()
      print


class site_refinement_test(refinement_test):

  debug = 1
  purpose = "site refinement"

  def run(self):
    self.exercise_floating_origin_restraints()
    refinement_test.run(self)

  def __init__(self):
    sgi = sgtbx.space_group_info("Hall: %s" % self.hall)
    cs = sgi.any_compatible_crystal_symmetry(volume=1000)
    xs = xray.structure(crystal.special_position_settings(cs))
    for sc in self.scatterers():
      sc.flags.set_use_u_iso(False).set_use_u_aniso(False)\
              .set_grad_site(True)
      xs.add_scatterer(sc)
    self.reference_xs = xs.as_emma_model()
    self.xray_structure = xs

    mi = cs.build_miller_set(d_min=0.5, anomalous_flag=False)
    ma = mi.structure_factors_from_scatterers(xs, algorithm="direct").f_calc()
    self.fo_sq = ma.norm().customized_copy(
      sigmas=flex.double(ma.size(), 1.))

  def exercise_floating_origin_restraints(self):
    n = self.n_independent_params
    eps_zero_rhs = 1e-6
    connectivity_table = smtbx.utils.connectivity_table(self.xray_structure)
    reparametrisation = constraints.reparametrisation(
      structure=self.xray_structure,
      geometrical_constraints=[],
      connectivity_table=connectivity_table)
    normal_eqns = least_squares.normal_equations(
      self.fo_sq,
      reparametrisation,
      weighting_scheme=least_squares.unit_weighting(),
      floating_origin_restraint_relative_weight=0)
    normal_eqns.build_up()
    assert normal_eqns.reduced.right_hand_side.all_approx_equal(
      0, eps_zero_rhs),\
           list(normal_eqns.reduced.right_hand_side)
    unrestrained_normal_matrix = normal_eqns.reduced.normal_matrix_packed_u
    assert len(unrestrained_normal_matrix) == n*(n+1)//2
    ev = eigensystem.real_symmetric(
      unrestrained_normal_matrix.matrix_packed_u_as_symmetric())
    unrestrained_eigenval = ev.values()
    unrestrained_eigenvec = ev.vectors()

    normal_eqns = least_squares.normal_equations(
      self.fo_sq,
      reparametrisation,
      weighting_scheme=least_squares.unit_weighting(),
    )
    normal_eqns.build_up()
    assert normal_eqns.reduced.right_hand_side.all_approx_equal(0, eps_zero_rhs),\
            list(normal_eqns.reduced.right_hand_side)
    restrained_normal_matrix = normal_eqns.reduced.normal_matrix_packed_u
    assert len(restrained_normal_matrix) == n*(n+1)//2
    ev = eigensystem.real_symmetric(
      restrained_normal_matrix.matrix_packed_u_as_symmetric())
    restrained_eigenval = ev.values()
    restrained_eigenvec = ev.vectors()

    # The eigendecomposition of the normal matrix
    # for the unrestrained problem is:
    #    A = sum_{0 <= i < n-p-1} lambda_i v_i^T v_i
    # where the eigenvalues lambda_i are sorted in decreasing order
    # and p is the dimension of the continous origin shift space.
    # In particular A v_i = 0, n-p <= i < n.
    # In the restrained case, it becomes:
    #    A' = A + sum_{n-p <= i < n} mu v_i^T v_i

    p = len(self.continuous_origin_shift_basis)
    assert approx_equal(restrained_eigenval[p:], unrestrained_eigenval[:-p],
                        eps=1e-12)
    assert unrestrained_eigenval[-p]/unrestrained_eigenval[-p-1] < 1e-12

    if p > 1:
      # eigenvectors are stored by rows
      unrestrained_null_space = unrestrained_eigenvec.matrix_copy_block(
        i_row=n-p, i_column=0,
        n_rows=p, n_columns=n)
      assert self.rank(unrestrained_null_space) == p

      restrained_space = restrained_eigenvec.matrix_copy_block(
        i_row=0, i_column=0,
        n_rows=p, n_columns=n)
      assert self.rank(restrained_space) == p

      singular = flex.double(
        self.continuous_origin_shift_basis)
      assert self.rank(singular) == p

      rank_finder = flex.double(n*3*p)
      rank_finder.resize(flex.grid(3*p, n))
      rank_finder.matrix_paste_block_in_place(unrestrained_null_space,
                                              i_row=0, i_column=0)
      rank_finder.matrix_paste_block_in_place(restrained_space,
                                              i_row=p, i_column=0)
      rank_finder.matrix_paste_block_in_place(singular,
                                              i_row=2*p, i_column=0)
      assert self.rank(rank_finder) == p
    else:
      # this branch handles the case p=1
      # it's necessary to work around a bug in the svd module
      # ( nx1 matrices crashes the code )
      assert approx_equal(
        restrained_eigenvec[0:n].angle(
          unrestrained_eigenvec[-n:]) % math.pi, 0)
      assert approx_equal(
        unrestrained_eigenvec[-n:].angle(
          flex.double(self.continuous_origin_shift_basis[0])) % math.pi, 0)

    # Do the floating origin restraints prevent the structure from floating?
    xs = self.xray_structure.deep_copy_scatterers()
    normal_eqns = least_squares.normal_equations(
      self.fo_sq,
      reparametrisation,
      weighting_scheme=least_squares.unit_weighting(),
    )
    barycentre_0 = xs.sites_frac().mean()
    while 1:
      xs.shake_sites_in_place(rms_difference=0.15)
      xs.apply_symmetry_sites()
      barycentre_1 = xs.sites_frac().mean()
      delta = matrix.col(barycentre_1) - matrix.col(barycentre_0)
      moved_far_enough = 0
      for singular in self.continuous_origin_shift_basis:
        e = matrix.col(singular[:3])
        if not approx_equal(delta.dot(e), 0, eps=0.01, out=None):
          moved_far_enough += 1
      if moved_far_enough: break

    # one refinement cycle
    normal_eqns.build_up()
    normal_eqns.solve()
    shifts = normal_eqns.shifts

    # That's what floating origin restraints are for!
    # Note that in the presence of special position, that's different
    # from the barycentre not moving along the continuous shift directions.
    # TODO: typeset notes about that subtlety.
    for singular in self.continuous_origin_shift_basis:
      assert approx_equal(shifts.dot(flex.double(singular)), 0, eps=1e-12)

  def rank(cls, a):
    """ row rank of a """
    rank_revealing = svd.real(a.deep_copy(),
                              accumulate_u=False, accumulate_v=False)
    return rank_revealing.numerical_rank(rank_revealing.sigma[0]*1e-9)
  rank = classmethod(rank)

  def exercise_ls_cycles(self):
    xs = self.xray_structure.deep_copy_scatterers()
    connectivity_table = smtbx.utils.connectivity_table(xs)
    reparametrisation = constraints.reparametrisation(
      structure=xs,
      geometrical_constraints=[],
      connectivity_table=connectivity_table)
    normal_eqns = least_squares.normal_equations(
      self.fo_sq, reparametrisation,
      weighting_scheme=least_squares.unit_weighting())
    emma_ref = xs.as_emma_model()
    xs.shake_sites_in_place(rms_difference=0.1)

    objectives = []
    scales = []
    fo_sq_max = flex.max(self.fo_sq.data())
    for i in xrange(5):
      normal_eqns.build_up()
      objectives.append(normal_eqns.objective)
      scales.append(normal_eqns.scale_factor)
      gradient_relative_norm = normal_eqns.gradient.norm()/fo_sq_max
      normal_eqns.solve_and_apply_shifts()
      shifts = normal_eqns.shifts

    assert approx_equal(normal_eqns.scale_factor, 1, eps=1e-5)
    assert approx_equal(normal_eqns.objective, 0)
    # skip next-to-last one to allow for no progress and rounding error
    assert objectives[0] >= objectives[1] >= objectives[3], objectives
    assert approx_equal(gradient_relative_norm, 0, eps=1e-9)

    match = emma.model_matches(emma_ref, xs.as_emma_model()).refined_matches[0]
    assert match.rt.r == matrix.identity(3)
    for pair in match.pairs:
      assert approx_equal(match.calculate_shortest_dist(pair), 0, eps=1e-4)


class adp_refinement_test(refinement_test):

  random_u_cart_scale = 0.2
  purpose = "ADP refinement"

  def __init__(self):
    sgi = sgtbx.space_group_info("Hall: %s" % self.hall)
    cs = sgi.any_compatible_crystal_symmetry(volume=1000)
    xs = xray.structure(crystal.special_position_settings(cs))
    for i, sc in enumerate(self.scatterers()):
      sc.flags.set_use_u_iso(False).set_use_u_aniso(True)\
              .set_grad_u_aniso(True)
      xs.add_scatterer(sc)
      site_symm = xs.site_symmetry_table().get(i)
      u_cart = adptbx.random_u_cart(u_scale=self.random_u_cart_scale)
      u_star = adptbx.u_cart_as_u_star(cs.unit_cell(), u_cart)
      xs.scatterers()[-1].u_star = site_symm.average_u_star(u_star)
    self.xray_structure = xs

    mi = cs.build_miller_set(d_min=0.5, anomalous_flag=False)
    ma = mi.structure_factors_from_scatterers(xs, algorithm="direct").f_calc()
    self.fo_sq = ma.norm().customized_copy(
      sigmas=flex.double(ma.size(), 1.))

  def exercise_ls_cycles(self):
    xs = self.xray_structure.deep_copy_scatterers()
    xs.shake_adp() # it must happen before the reparamtrisation is constructed
                   # because the ADP values are read then and only then.
    connectivity_table = smtbx.utils.connectivity_table(xs)
    reparametrisation = constraints.reparametrisation(
      structure=xs,
      geometrical_constraints=[],
      connectivity_table=connectivity_table)
    normal_eqns = least_squares.normal_equations(
      xs, self.fo_sq, reparametrisation,
      weighting_scheme=least_squares.unit_weighting())

    objectives = []
    gradient_relative_norms = []
    scales = []
    fo_sq_max = flex.max(self.fo_sq.data())
    for i in xrange(10):
      normal_eqns.build_up()
      objectives.append(normal_eqns.objective)
      scales.append(normal_eqns.scale_factor)
      gradient_relative_norms.append(normal_eqns.gradient.norm()/fo_sq_max)
      normal_eqns.solve_and_apply_shifts()
      shifts = normal_eqns.shifts

    assert approx_equal(normal_eqns.scale_factor, 1, eps=1e-4)
    assert approx_equal(normal_eqns.objective, 0)
    # skip next-to-last one to allow for no progress and rounding error
    n = len(objectives)
    assert objectives[0] > objectives[n-1], objectives
    assert approx_equal(gradient_relative_norms[-1], 0, eps=1e-6)

    for sc0, sc1 in zip(self.xray_structure.scatterers(), xs.scatterers()):
      assert approx_equal(sc0.u_star, sc1.u_star)


class p1_test(object):

  hall = "P 1"
  n_independent_params = 9
  continuous_origin_shift_basis = [ (1,0,0)*3,
                                    (0,1,0)*3,
                                    (0,0,1)*3 ]

  def scatterers(self):
    yield xray.scatterer("C1", (0.1, 0.2, 0.3))
    yield xray.scatterer("C2", (0.4, 0.7, 0.8))
    yield xray.scatterer("C3", (-0.1, -0.8, 0.6))


class p2_test(object):

  hall = "P 2x"
  n_independent_params = 7
  continuous_origin_shift_basis = [ (1,0,0, 1, 1,0,0) ]

  def scatterers(self):
    yield xray.scatterer("C1", (0.1, 0.2, 0.3))
    yield xray.scatterer("C2", (-0.3, 0, 0)) # on 2-axis
    yield xray.scatterer("C3", (0.4, 0.1, -0.1)) # near 2-axis

class pm_test(object):

  hall = "P -2x"
  n_independent_params = 8
  continuous_origin_shift_basis = [ (0,1,0, 1,0, 0,1,0),
                                    (0,0,1, 0,1, 0,0,1) ]

  def scatterers(self):
    yield xray.scatterer("C1", (0.1, 0.2, 0.3)) # near mirror plance
    yield xray.scatterer("C2", (0, -0.3, 0.4)) # on mirror plane
    yield xray.scatterer("C3", (0.7, 0.1, -0.1))


class site_refinement_in_p1_test(p1_test, site_refinement_test): pass

class site_refinement_in_p2_test(p2_test, site_refinement_test): pass

class site_refinement_in_pm_test(pm_test, site_refinement_test): pass


class adp_refinement_in_p1_test(p1_test, adp_refinement_test): pass

class adp_refinement_in_p2_test(p2_test, adp_refinement_test): pass

class adp_refinement_in_pm_test(pm_test, adp_refinement_test): pass


def exercise_normal_equations():
  adp_refinement_in_p1_test().run()
  adp_refinement_in_pm_test().run()
  adp_refinement_in_p2_test().run()

  site_refinement_in_p1_test().run()
  site_refinement_in_pm_test().run()
  site_refinement_in_p2_test().run()

class special_positions_test(object):

  delta_site   = 0.1 # angstrom
  delta_u_star = 0.1 # %

  def __init__(self, n_runs, **kwds):
    libtbx.adopt_optional_init_args(self, kwds)
    self.n_runs = n_runs
    self.crystal_symmetry = crystal.symmetry(
      unit_cell=uctbx.unit_cell((5.1534, 5.1534, 8.6522, 90, 90, 120)),
      space_group_symbol='Hall: P 6c')
    self.structure = xray.structure(
      self.crystal_symmetry.special_position_settings(),
      flex.xray_scatterer((
        xray.scatterer('K1',
                        site=(0, 0, -0.00195),
                        u=self.u_cif_as_u_star((0.02427, 0.02427, 0.02379,
                                                0.01214, 0.00000, 0.00000))),
        xray.scatterer('S1',
                       site=(1/3, 2/3, 0.204215),
                       u=self.u_cif_as_u_star((0.01423, 0.01423, 0.01496,
                                               0.00712, 0.00000, 0.00000 ))),
        xray.scatterer('Li1',
                       site=(1/3, 2/3, 0.815681),
                       u=self.u_cif_as_u_star((0.02132, 0.02132, 0.02256,
                                               0.01066, 0.00000, 0.00000 ))),
        xray.scatterer('O1',
                       site=(1/3, 2/3, 0.035931),
                       u=self.u_cif_as_u_star((0.06532, 0.06532, 0.01669,
                                               0.03266, 0.00000, 0.00000 ))),
        xray.scatterer('O2',
                       site=(0.343810, 0.941658, 0.258405),
                       u=self.u_cif_as_u_star((0.02639,  0.02079, 0.05284,
                                               0.01194, -0.00053,-0.01180 )))
      )))
    mi = self.crystal_symmetry.build_miller_set(anomalous_flag=False,
                                                d_min=0.5)
    fo_sq = mi.structure_factors_from_scatterers(
      self.structure, algorithm="direct").f_calc().norm()
    self.fo_sq = fo_sq.customized_copy(sigmas=flex.double(fo_sq.size(), 1))

  def u_cif_as_u_star(self, u_cif):
    return adptbx.u_cif_as_u_star(self.crystal_symmetry.unit_cell(), u_cif)

  def shake_point_group_3(self, sc):
    _, _, c, _, _, _ = self.crystal_symmetry.unit_cell().parameters()

    x, y, z = sc.site
    z += random.uniform(-self.delta_site, self.delta_site)/c
    sc.site = (x, y, z)

    u11, _, u33, _, _, _ = sc.u_star
    u11 *= 1 + random.uniform(-self.delta_u_star, self.delta_u_star)
    u33 *= 1 + random.uniform(-self.delta_u_star, self.delta_u_star)
    sc.u_star = (u11, u11, u33, u11/2, 0, 0)

  def run(self):
    if self.n_runs > 1:
      print 'small inorganic refinement'
      for i in xrange(self.n_runs):
        print '.',
        self.exercise()
      print
    else:
      self.exercise()

  def exercise(self):
    xs0 = self.structure
    xs = xs0.deep_copy_scatterers()
    k1, s1, li1, o1, o2 = xs.scatterers()
    self.shake_point_group_3(k1)
    self.shake_point_group_3(s1)
    self.shake_point_group_3(li1)
    self.shake_point_group_3(o1)
    o2.site = tuple(
      [ x*(1 + random.uniform(-self.delta_site, self.delta_site))
        for x in o2.site])
    o2.u_star = tuple(
      [ u*(1 + random.uniform(-self.delta_u_star, self.delta_u_star))
        for u in o2.u_star])

    for sc in xs.scatterers():
      sc.flags.set_use_u_iso(False).set_use_u_aniso(True)
      sc.flags.set_grad_site(True).set_grad_u_aniso(True)
      connectivity_table = smtbx.utils.connectivity_table(xs)
    reparametrisation = constraints.reparametrisation(
      structure=xs,
      geometrical_constraints=[],
      connectivity_table=connectivity_table)
    normal_eqns = least_squares.normal_equations(
      self.fo_sq, reparametrisation,
      weighting_scheme=least_squares.unit_weighting())

    objectives = []
    scales = []
    fo_sq_max = flex.max(self.fo_sq.data())
    shifts = []
    for i in xrange(10):
      normal_eqns.build_up()
      objectives.append(normal_eqns.objective)
      scales.append(normal_eqns.scale_factor)
      gradient_relative_norm = normal_eqns.gradient.norm()/fo_sq_max
      a = normal_eqns.reduced.normal_matrix_packed_u.deep_copy()
      normal_eqns.solve_and_apply_shifts()
      shifts.append(normal_eqns.shifts)

    ## Test whether refinement brought back the shaked structure to its
    ## original state
    match = emma.model_matches(xs0.as_emma_model(),
                               xs.as_emma_model()).refined_matches[0]
    assert match.rt.r == matrix.identity(3)
    assert not match.singles1 and not match.singles2
    assert match.rms < 1e-6

    delta_u_carts= (   xs.scatterers().extract_u_cart(xs.unit_cell())
                    - xs0.scatterers().extract_u_cart(xs.unit_cell())).norms()
    assert flex.abs(delta_u_carts) < 1e-6

    assert approx_equal(normal_eqns.scale_factor, 1, eps=1e-4)

    ## Test covariance matrix
    cov = normal_eqns.covariance_matrix(normalised_by_goof=False)\
        .matrix_packed_u_as_symmetric()
    m, n = cov.accessor().focus()
    # x,y for point group 3 sites are fixed: no variance or correlation
    for i in (0, 9, 18, 27,):
      assert cov.matrix_copy_block(i, 0, 2, n) == 0

    # u_star coefficients u13 and u23 for point group 3 sites are fixed
    # to 0: again no variance or correlation with any other param
    for i in (7, 16, 25, 34,):
      assert cov.matrix_copy_block(i, 0, 2, n).as_1d()\
             .all_approx_equal(0., 1e-20)

    # u_star coefficients u11, u22 and u12 for point group 3 sites
    # are totally correlated, with variances in ratios 1:1:1/2
    for i in (3, 12, 21, 30,):
      assert cov[i, i] != 0
      assert approx_equal(cov[i, i], cov[i+1, i+1], eps=1e-15)
      assert approx_equal(cov[i, i+1]/cov[i, i], 1, eps=1e-12)
      assert approx_equal(cov[i, i+3]/cov[i, i], 0.5, eps=1e-12)


def run():
  libtbx.utils.show_times_at_exit()
  import sys
  from libtbx.option_parser import option_parser
  command_line = (option_parser()
    .option(None, "--fix_random_seeds",
            action="store_true",
            default=False)
    .option(None, "--runs",
            type='int',
            default=1)
  ).process(args=sys.argv[1:])
  if command_line.options.fix_random_seeds:
    flex.set_random_seed(1)
    random.seed(2)
  n_runs = command_line.options.runs
  if n_runs > 1: refinement_test.ls_cycle_repeats = n_runs
  exercise_normal_equations()
  special_positions_test(n_runs).run()

if __name__ == '__main__':
  run()
