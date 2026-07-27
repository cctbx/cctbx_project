from __future__ import absolute_import, division, print_function

from cctbx import crystal, sgtbx, xray
from cctbx import euclidean_model_matching as emma
from cctbx.array_family import flex
from scitbx import matrix
from scitbx.lstbx import normal_eqns_solving
from smtbx.refinement import constraints, least_squares
from smtbx.refinement.restraints import origin_fixing_restraints
import smtbx.utils
from libtbx.test_utils import approx_equal
import math

# Small enough that the whole file runs in seconds, big enough that the
# parameters differ in scale the way a real refinement's do: the curvature with
# respect to a fractional coordinate and to a U_iso are orders of magnitude
# apart, which is what the preconditioning in scipy_iterations is there for.
unit_cell = (9.5, 10.5, 11.5, 88., 95., 91.)
space_group_symbol = "P1"
scatterer_sites = ((0.10, 0.20, 0.30), (0.40, 0.70, 0.80),
                   (-0.10, -0.80, 0.60), (0.25, 0.45, 0.15))
u_isos = (0.02, 0.035, 0.025, 0.03)


def structure(site_delta=0., u_delta=0.):
  """ The test structure, optionally moved off the point it is defined at.

  Everything about it is fixed, cell included, and the displacement is a
  function of the parameter index rather than random: a test which fails only
  sometimes is worse than no test.
  """
  cs = crystal.symmetry(unit_cell=unit_cell,
                        space_group_symbol=space_group_symbol)
  xs = xray.structure(crystal.special_position_settings(cs))
  for i, site in enumerate(scatterer_sites):
    site = tuple(x + site_delta*math.sin(3.1*(3*i + j))
                 for j, x in enumerate(site))
    sc = xray.scatterer("C%i" % (i + 1), site,
                        u=u_isos[i] + u_delta*(-1)**i)
    sc.flags.set_use_u_iso(True).set_use_u_aniso(False)
    sc.flags.set_grad_site(True).set_grad_u_iso(True)
    xs.add_scatterer(sc)
  return xs


def reference_structure():
  return structure()


def shaken_structure():
  return structure(site_delta=0.03, u_delta=0.004)


def observations_from(xs, d_min=0.8, noise=0.03):
  """ The intensities xs would give rise to, with a little noise on them.

  The noise matters. A fit which can be made perfect is atypically forgiving of
  the one approximation in the gradient: the SHELX weights depend on Fc and the
  gradient treats them as constants, but the error that introduces is
  multiplied by the residual, which a noiseless fit drives to zero. Real data
  leave a residual behind and then the error tells.

  The perturbation is a fixed function of the reflection index rather than
  random, so that a failure here is always the same failure.
  """
  mi = xs.crystal_symmetry().build_miller_set(d_min=d_min,
                                              anomalous_flag=False)
  ma = mi.structure_factors_from_scatterers(xs, algorithm="direct").f_calc()
  fo_sq = ma.norm()
  sigmas = noise*(fo_sq.data() + flex.mean(fo_sq.data()))
  perturbation = flex.double(
    [math.sin(7.3*i) for i in range(fo_sq.size())])
  fo_sq = fo_sq.customized_copy(data=fo_sq.data() + sigmas*perturbation,
                                sigmas=sigmas)
  return fo_sq.as_xray_observations()


def least_squares_for(xs, observations, weighting_scheme):
  reparametrisation = constraints.reparametrisation(
    structure=xs, constraints=[],
    connectivity_table=smtbx.utils.connectivity_table(xs))
  return least_squares.crystallographic_ls(
    observations, reparametrisation,
    weighting_scheme=weighting_scheme,
    origin_fixing_restraints_type=
    origin_fixing_restraints.atomic_number_weighting)


def exercise_parameter_vector():
  """ parameter_vector() shall report where the parameters really are.

  apply_shifts() does not always move the parameters by the shifts it is given:
  validate() constrains some of them, a U_iso driven negative being reset to
  1e-4. A minimiser which places the parameters itself works out its next shift
  from where it believes the previous one left them, so believing the request
  rather than the outcome puts every subsequent shift out by the difference --
  and the minimiser is never told, so it goes on refining a structure which is
  not the one it is being given objectives for.
  """
  observations = observations_from(reference_structure())
  ls = least_squares_for(shaken_structure(), observations,
                         least_squares.unit_weighting())
  ls.build_up()
  before = ls.parameter_vector()
  assert before.size() == ls.reparametrisation.n_independents, before.size()

  # a shift small enough that nothing is constrained is applied in full
  small = flex.double(before.size(), 0)
  small[0] = 1e-3
  ls.apply_shifts(small)
  after = ls.parameter_vector()
  assert approx_equal(after, before + small, eps=1e-12), \
    (list(after), list(before + small))
  ls.apply_shifts(-small)

  # a shift which drives a U_iso negative is not, and that is visible
  u_iso_indices = [i for i, annotation
                   in enumerate(ls.reparametrisation.component_annotations)
                   if annotation.endswith('.uiso')]
  assert u_iso_indices, ls.reparametrisation.component_annotations
  i = u_iso_indices[0]
  before = ls.parameter_vector()
  clamping = flex.double(before.size(), 0)
  clamping[i] = -10*before[i]
  ls.apply_shifts(clamping)
  after = ls.parameter_vector()
  assert after[i] > 0, after[i]
  assert after[i] != before[i] + clamping[i], (after[i], before[i])
  assert approx_equal(after[i], 1e-4, eps=1e-12), after[i]


def exercise_gradient_consistency():
  """ The gradient of the crystallographic objective shall be its gradient.

  This is the premise the whole module rests on, and worth checking on a real
  structure rather than only on the test problems of scitbx: the objective has
  the overall scale factor optimised away, carries the restraints, and is put
  together by a good deal of C++.

  Weights which do not depend on Fc -- unit, sigma, and the SHELX scheme with
  a = b = 0, which is just 1/sigma^2 -- give a gradient which is exact. Where a
  is not zero the weights depend on Fc and the gradient treats them as
  constants, so it cannot be exact; how far off it is, is measured rather than
  assumed, since it decides how much a line search can be asked to do.

  Every evaluation is made on a freshly built L.S. object, so that all of them
  derive the scale factor the weights need in the same way. Reusing one object
  would fold in the path dependence which exercise_weighting_path_dependence
  measures separately.

  Only a few parameters are differenced: each costs two passes over the data,
  and they are all of a kind. The error is measured against the norm of the
  whole gradient rather than component by component, a component which happens
  to be near zero saying nothing about how well a line search will fare.
  """
  observations = observations_from(reference_structure())
  xs = shaken_structure()
  n_differenced = 6
  h = 1e-6

  for name, weighting_scheme, tolerance in (
    ('unit', least_squares.unit_weighting(), 1e-6),
    ('sigma', least_squares.sigma_weighting(), 1e-6),
    ('shelx a=0', least_squares.mainstream_shelx_weighting(a=0, b=0), 1e-6),
    ('shelx a=0.1',
     least_squares.mainstream_shelx_weighting(a=0.1, b=0), 0.1)):
    ls = least_squares_for(xs.deep_copy_scatterers(), observations,
                           weighting_scheme)
    ls.build_up()
    gradient = -ls.opposite_of_gradient()

    def objective_at(shifts):
      fresh = least_squares_for(xs.deep_copy_scatterers(), observations,
                                weighting_scheme)
      fresh.apply_shifts(shifts)
      fresh.build_up()
      return fresh.objective()

    worst = 0
    for i in range(n_differenced):
      shifts = flex.double(gradient.size(), 0)
      shifts[i] = h
      forward = objective_at(shifts)
      shifts[i] = -h
      backward = objective_at(shifts)
      finite_difference = (forward - backward)/(2*h)
      worst = max(worst, abs(finite_difference - gradient[i]))
    worst /= gradient.norm()
    print("\t%-12s weighting: gradient differs from finite differences by "
          "%.2e of its norm" % (name, worst))
    assert worst < tolerance, (name, worst)


def exercise_weighting_path_dependence():
  """ How much does the objective depend on how its point was arrived at?

  build_up() hands the weighting scheme the scale factor optimised at the point
  evaluated *before* it. Where the weights depend on that scale factor, the
  objective is therefore not quite a function of the parameters alone. Repeated
  builds at one point are a fixed-point iteration which settles quickly, so it
  is a transient rather than an ambiguity, but a line search probing a distant
  point and coming back sees it. Gauss-Newton never does, walking towards the
  minimum and never back.

  The unit and sigma schemes are free of it, their weights depending on neither
  Fc nor the scale factor, and are asserted to be reproducible to within the
  round-off of taking the parameters away and bringing them back.
  """
  observations = observations_from(reference_structure())
  detour = 0.05

  for name, weighting_scheme, tolerance in (
    ('unit', least_squares.unit_weighting(), 1e-11),
    ('sigma', least_squares.sigma_weighting(), 1e-11),
    ('shelx', least_squares.mainstream_shelx_weighting(a=0.1, b=0), 0.1)):
    ls = least_squares_for(shaken_structure(), observations, weighting_scheme)
    ls.build_up()
    ls.build_up() # let the scale factor settle before measuring anything
    direct = ls.objective()

    excursion = flex.double(ls.opposite_of_gradient().size(), 0)
    excursion.set_selected(flex.size_t(range(4)),
                           flex.double([detour]*4))
    ls.apply_shifts(excursion)
    ls.build_up()
    ls.apply_shifts(-excursion)
    ls.build_up()
    roundabout = ls.objective()
    ls.build_up()
    settled = ls.objective()

    difference = abs(roundabout - direct)/direct
    recovered = abs(settled - direct)/direct
    print("\t%-5s weighting: after a detour of %.2f the objective at the same "
          "point is out by %.2e relative, %.2e after one more build"
          % (name, detour, difference, recovered))
    assert difference < tolerance, (name, difference)
    # whatever the detour did, one further build very nearly undoes it: the
    # objective is a function of the parameters in the limit of repeated builds
    assert recovered <= max(difference, tolerance/10), (name, recovered)


def exercise_site_and_adp_refinement():
  """ The scipy minimisers shall refine what Gauss-Newton refines.

  Gauss-Newton sets the standard, and it is the right standard to use rather
  than the structure the data were generated from: the data carry noise, so the
  refined structure is not the generating one, and what is being tested is the
  minimiser and not the crystallography.

  Both weighting schemes are exercised, and the difference between them is the
  point of doing so. With weights that do not depend on Fc the gradient is
  exact and every method reaches Gauss-Newton's objective to the last digit.
  With a = 0.1 it is out by some 5%, and every method then stops a little
  short -- by 0.007% for L-BFGS-B up to 1.7% for trust-ncg -- which is what the
  looser tolerance allows for, and no fault of the minimisers.
  """
  from smtbx.refinement import minimisers
  observations = observations_from(reference_structure())

  for name, weighting_scheme, tolerance in (
    ('shelx a=0', least_squares.mainstream_shelx_weighting(a=0), 1e-6),
    ('shelx a=0.1',
     least_squares.mainstream_shelx_weighting(a=0.1), 3e-2)):
    reference_xs = shaken_structure()
    reference_ls = least_squares_for(reference_xs, observations,
                                     weighting_scheme)
    normal_eqns_solving.naive_iterations(
      reference_ls, gradient_threshold=1e-12, step_threshold=1e-7)
    best = reference_ls.objective()
    # the scale factor is not refined but optimised away at every evaluation,
    # so it is an outcome of the minimisation like any other and is held to
    # the same standard: Gauss-Newton's value, not the nominal 1
    best_scale = reference_ls.scale_factor()
    emma_reference = reference_xs.as_emma_model()

    for method in ('trust-ncg', 'L-BFGS-B', 'CG', 'BFGS'):
      xs = shaken_structure()
      ls = least_squares_for(xs, observations, weighting_scheme)
      cycles = minimisers.crystallographic_scipy_iterations(
        ls, method=method, n_max_iterations=200, max_evaluations=5000,
        g_tolerance=1e-10, f_tolerance=1e-14, x_tolerance=1e-10)
      achieved = ls.objective()
      assert achieved <= best*(1 + tolerance), (name, method, achieved, best)
      assert approx_equal(ls.scale_factor(), best_scale, eps=2*tolerance), \
        (name, method, ls.scale_factor(), best_scale)

      match = emma.model_matches(
        emma_reference, xs.as_emma_model()).refined_matches[0]
      assert match.rt.r == matrix.identity(3), (name, method)
      for pair in match.pairs:
        distance = match.calculate_shortest_dist(pair)
        assert approx_equal(distance, 0, eps=1e-2), (name, method, distance)
      for refined, sc in zip(reference_xs.scatterers(), xs.scatterers()):
        assert approx_equal(sc.u_iso, refined.u_iso, eps=1e-3), \
          (name, method, sc.label, sc.u_iso, refined.u_iso)

      print("\t%-12s %-10s %3d iterations, %4d evaluations, objective "
            "%.6e (Gauss-Newton %.6e)"
            % (name, method, cycles.n_iterations, cycles.n_evaluations,
               achieved, best))


def exercise_standard_uncertainties():
  """ Standard uncertainties shall be available once a minimisation is over.

  Whoever refines a structure wants s.u. afterwards, and they come from the
  normal matrix. A minimiser which does not build normal equations to take its
  steps still has to leave them describing where it stopped.
  """
  from smtbx.refinement import minimisers
  observations = observations_from(reference_structure())
  weighting_scheme = least_squares.mainstream_shelx_weighting(a=0)

  xs = shaken_structure()
  ls = least_squares_for(xs, observations, weighting_scheme)
  minimisers.crystallographic_scipy_iterations(
    ls, method='trust-ncg', n_max_iterations=200, max_evaluations=5000)

  covariance = ls.covariance_matrix_and_annotations()
  variances = covariance.diagonal()
  assert variances.size() == ls.reparametrisation.n_independents, \
    (variances.size(), ls.reparametrisation.n_independents)
  for i, variance in enumerate(variances):
    assert variance > 0, (i, variance, covariance.annotations[i])
    assert not math.isnan(variance), (i, covariance.annotations[i])
  print("\tstandard uncertainties span %.2e to %.2e"
        % (math.sqrt(flex.min(variances)), math.sqrt(flex.max(variances))))


def run():
  try:
    import scipy.optimize # noqa: F401
  except ImportError:
    print('Skipping tst_scipy_minimisers: scipy is not available')
    return
  exercise_parameter_vector()
  exercise_gradient_consistency()
  exercise_weighting_path_dependence()
  exercise_site_and_adp_refinement()
  exercise_standard_uncertainties()
  print('OK')


if __name__ == '__main__':
  run()
