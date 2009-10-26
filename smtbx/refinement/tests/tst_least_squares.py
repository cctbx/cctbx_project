from scitbx.linalg import eigensystem, svd
from scitbx import matrix
from cctbx import sgtbx, crystal, xray, miller
from cctbx import euclidean_model_matching as emma
from cctbx.array_family import flex
from smtbx.refinement import least_squares
from smtbx import refinement
from libtbx.test_utils import approx_equal
import libtbx.utils
from stdlib import math
import sys

class test_case(object):

  debug = 1
  ls_cycle_repeats = 1

  def run(self):
    self.exercise_floating_origin_restraints()
    if self.ls_cycle_repeats == 1:
      self.exercise_ls_cycles()
    else:
      print self.hall
      for n in xrange(self.ls_cycle_repeats):
        self.exercise_ls_cycles()
        print '.',
        sys.stdout.flush()
      print

  def __init__(self):
    sgi = sgtbx.space_group_info("Hall: %s" % self.hall)
    cs = sgi.any_compatible_crystal_symmetry(volume=1000)
    xs = xray.structure(crystal.special_position_settings(cs))
    for sc in self.scatterers():
      xs.add_scatterer(sc)
    self.reference_xs = xs.as_emma_model()
    self.xray_structure = xs

    mi = cs.build_miller_set(d_min=0.5, anomalous_flag=False)
    ma = mi.structure_factors_from_scatterers(xs, algorithm="direct").f_calc()
    self.fo_sq = ma.norm().customized_copy(
      sigmas=flex.double(ma.size(), 1.))

  def exercise_floating_origin_restraints(self):
    n = self.n_independent_params
    if sys.version.find('GCC 4.0.0 20050519 (Red Hat 4.0.0-8)') == -1:
      eps_zero_rhs = 1e-12
    else:
      eps_zero_rhs = 5e-7    
    normal_eqns = least_squares.normal_equations(
      self.xray_structure,
      self.fo_sq,
      weighting_scheme=least_squares.unit_weighting(),
      floating_origin_restraint_relative_weight=0)
    normal_eqns.build_up()
    assert normal_eqns.reduced.right_hand_side.all_approx_equal(0, eps_zero_rhs),\
            list(normal_eqns.reduced.right_hand_side)
    unrestrained_normal_matrix = normal_eqns.reduced.normal_matrix_packed_u
    assert len(unrestrained_normal_matrix) == n*(n+1)//2
    ev = eigensystem.real_symmetric(
      unrestrained_normal_matrix.matrix_packed_u_as_symmetric())
    unrestrained_eigenval = ev.values()
    unrestrained_eigenvec = ev.vectors()

    normal_eqns = least_squares.normal_equations(
      self.xray_structure,
      self.fo_sq,
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
      xs,
      self.fo_sq,
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
    shifts = normal_eqns.solve()

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
    normal_eqns = least_squares.normal_equations(
      xs, self.fo_sq, weighting_scheme=least_squares.unit_weighting())
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
      shifts = normal_eqns.solve_and_apply_shifts()

    assert approx_equal(normal_eqns.scale_factor, 1, eps=1e-5)
    assert approx_equal(normal_eqns.objective, 0)
    # skip next-to-last one to allow for no progress and rounding error
    assert objectives[0] >= objectives[1] >= objectives[3], objectives
    assert approx_equal(gradient_relative_norm, 0, eps=1e-9)

    match = emma.model_matches(emma_ref, xs.as_emma_model()).refined_matches[0]
    assert match.rt.r == matrix.identity(3)
    for pair in match.pairs:
      assert approx_equal(match.calculate_shortest_dist(pair), 0, eps=1e-4)


class p1_test_case(test_case):

  hall = "P 1"
  n_independent_params = 9
  continuous_origin_shift_basis = [ (1,0,0)*3,
                                    (0,1,0)*3,
                                    (0,0,1)*3 ]

  def scatterers(self):
    sc = xray.scatterer("C1", (0.1, 0.2, 0.3))
    sc.flags.set_grad_site(True)
    yield sc
    sc = xray.scatterer("C2", (0.4, 0.7, 0.8))
    sc.flags.set_grad_site(True)
    yield sc
    sc = xray.scatterer("C3", (-0.1, -0.8, 0.6))
    sc.flags.set_grad_site(True)
    yield sc


class p2_test_case(test_case):

  hall = "P 2x"
  n_independent_params = 7
  continuous_origin_shift_basis = [ (1,0,0, 1, 1,0,0) ]

  def scatterers(self):
    sc = xray.scatterer("C1", (0.1, 0.2, 0.3))
    sc.flags.set_grad_site(True)
    yield sc
    sc = xray.scatterer("C2", (-0.3, 0, 0)) # on 2-axis
    sc.flags.set_grad_site(True)
    yield sc
    sc = xray.scatterer("C3", (0.4, 0.1, -0.1)) # near 2-axis
    sc.flags.set_grad_site(True)
    yield sc

class pm_test_case(test_case):

  hall = "P -2x"
  n_independent_params = 8
  continuous_origin_shift_basis = [ (0,1,0, 1,0, 0,1,0),
                                    (0,0,1, 0,1, 0,0,1) ]

  def scatterers(self):
    sc = xray.scatterer("C1", (0.1, 0.2, 0.3)) # near mirror plance
    sc.flags.set_grad_site(True)
    yield sc
    sc = xray.scatterer("C2", (0, -0.3, 0.4)) # on mirror plane
    sc.flags.set_grad_site(True)
    yield sc
    sc = xray.scatterer("C3", (0.7, 0.1, -0.1))
    sc.flags.set_grad_site(True)
    yield sc


class r3_test_case(test_case):

  hall = "R 3 (-y+z, x+z, -x+y+z)"
  n_independent_params = 7
  continuous_origin_shift_basis = [(1,1,1, 1, 1,1,1)]

  def scatterers(self):
    sc = xray.scatterer("C1", (0.1, 0.2, 0.3))
    sc.flags.set_grad_site(True)
    yield sc
    sc = xray.scatterer("C2", (0.4, 0.4, 0.4)) # on 3-axis
    sc.flags.set_grad_site(True)
    yield sc
    sc = xray.scatterer("C3", (0.3, 0.4, 0.4)) # near 3-axis
    sc.flags.set_grad_site(True)
    yield sc


def exercise_normal_equations():
  p1_test_case().run()
  pm_test_case().run()
  p2_test_case().run()
  r3_test_case().run()


def run():
  libtbx.utils.show_times_at_exit()
  import sys
  if sys.argv[1:]:
    test_case.ls_cycle_repeats = int(sys.argv[1])
  exercise_normal_equations()

if __name__ == '__main__':
  run()
