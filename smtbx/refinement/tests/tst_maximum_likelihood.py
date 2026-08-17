from __future__ import absolute_import, division, print_function

""" Maximum-likelihood refinement: the pieces that can be checked on their own.

The amplitude target itself is not reimplemented here - cctbx/xray/targets/mlf.h
already has it, and cctbx.xray.mlf_target_and_gradients exposes it. What needs
checking is the step that lets that target drive smtbx's existing Gauss-Newton
accumulator, which only ever receives (y_calc, grad y_calc, y_obs, weight) and
knows nothing about likelihoods.

The trick is to hand it an *effective* observation and weight:

    w_eff  = 2a^2/(eb)  acentric,  a^2/(eb)  centric
    yo_eff = (fo/a)*I1/I0(2a*fo*|Fc|/(eb))   acentric
             (fo/a)*tanh(a*fo*|Fc|/(eb))     centric

chosen so that w_eff*(yo_eff - |Fc|) == -dT/d|Fc| exactly. The accumulator then
produces the true ML gradient while forming a positive-definite matrix, because
w_eff is the positive part of the curvature.

That identity is what this file checks, against cctbx's own implementation
rather than against a second copy of the same reasoning - the two conventions
it is easy to get wrong being that the target class stores conj(dT/dFc)/n_work,
and that alpha and beta are pre-multiplied by k and k^2 inside it.
"""

import math
import random

from cctbx import xray
from cctbx.array_family import flex
from scitbx.math import bessel_i1_over_i0
from libtbx.test_utils import approx_equal


def d_target_d_modulus_fcalc(fo, fc, a, b, k, e, centric):
  """ dT/d|Fc| for one reflection, taken from cctbx's own MLF target.

  The class returns conj(dT/dFc)/n_work; with a single work reflection n_work
  is 1, and dT/dFc = (dT/d|Fc|)*conj(Fc)/|Fc|, so undoing the conjugation and
  the phase leaves the real derivative wrt the modulus.
  """
  t = xray.mlf_target_and_gradients(
    f_obs=flex.double([fo]),
    r_free_flags=flex.bool([False]),
    f_calc=flex.complex_double([fc]),
    alpha=flex.double([a]),
    beta=flex.double([b]),
    scale_factor=k,
    epsilons=flex.double([e]),
    centric_flags=flex.bool([centric]),
    compute_gradients=True)
  g = t.gradients_work()[0]
  return (g.conjugate()*fc/abs(fc)).real


def effective_observation(fo, modulus_fc, a, b, k, e, centric):
  """ The (weight, observation) pair handed to the accumulator. """
  a *= k
  b *= k*k
  eb = e*b
  if centric:
    w = a*a/eb
    yo = (fo/a)*math.tanh(a*fo*modulus_fc/eb)
  else:
    w = 2.*a*a/eb
    yo = (fo/a)*bessel_i1_over_i0(2.*a*fo*modulus_fc/eb)
  return w, yo


def exercise_effective_observation_reproduces_the_gradient():
  random.seed(0)
  worst = 0.
  for centric in (False, True):
    for _ in range(400):
      fo = random.uniform(0.5, 60.)
      modulus_fc = random.uniform(0.5, 60.)
      phase = random.uniform(0., 2.*math.pi)
      fc = complex(modulus_fc*math.cos(phase), modulus_fc*math.sin(phase))
      a = random.uniform(0.2, 1.)
      b = random.uniform(0.5, 40.)
      k = random.uniform(0.5, 2.)
      e = random.choice([1., 2., 3., 4.])
      d = d_target_d_modulus_fcalc(fo, fc, a, b, k, e, centric)
      w, yo = effective_observation(fo, modulus_fc, a, b, k, e, centric)
      worst = max(worst, abs(w*(yo - modulus_fc) + d)/max(abs(d), 1e-12))
  assert worst < 1e-9, worst
  print("\tthe effective observation reproduces dT/d|Fc|")


def exercise_the_weight_is_positive():
  """ The accumulator cannot take a negative weight.

  rank_n_update and store_design_row both assume w >= 0, and linear_ls::solve
  factorises by Cholesky, so a negative weight is not merely inaccurate - it
  breaks the solve. w_eff is the positive part of the curvature by
  construction; this is the guard that it stays that way.
  """
  for centric in (False, True):
    for a in (0.05, 0.5, 1.):
      for b in (1e-2, 1., 100.):
        w, _ = effective_observation(12., 9., a, b, 1., 2., centric)
        assert w > 0, (centric, a, b, w)
  print("\tthe effective weight is positive")


def effective_observation_in_intensity_space(fo_sq, observable, a, b, k, e,
                                             centric):
  """ The same reduction expressed against F^2, which is what the loop has.

  The build loop computes u = |Fc|^2 and du/dp, not |Fc|. Because
  du/dp = 2|Fc| d|Fc|/dp, the |Fc|-space equations can be formed from the F^2
  gradients by rescaling the weight and shifting the observation, which is what
  smtbx/refinement/ml_target.h does and what this checks.
  """
  fo = math.sqrt(fo_sq)
  modulus_fc = math.sqrt(observable)
  w_amp, yo_amp = effective_observation(fo, modulus_fc, a, b, k, e, centric)
  w_u = w_amp/(4.*observable)
  yo_u = 2.*modulus_fc*yo_amp - observable
  return w_u, yo_u


def exercise_intensity_space_form_is_the_same_equations():
  """ w_u (yo_u - u) must equal -dT/du, with dT/du = (dT/d|Fc|)/(2|Fc|).

  If that holds, the normal equations built from the F^2 observable and its
  gradients are bit-for-bit the |Fc|-space maximum-likelihood ones, and the
  observable never has to change.
  """
  random.seed(1)
  worst = 0.
  for centric in (False, True):
    for _ in range(400):
      fo = random.uniform(0.5, 60.)
      modulus_fc = random.uniform(0.5, 60.)
      phase = random.uniform(0., 2.*math.pi)
      fc = complex(modulus_fc*math.cos(phase), modulus_fc*math.sin(phase))
      a = random.uniform(0.2, 1.)
      b = random.uniform(0.5, 40.)
      k = random.uniform(0.5, 2.)
      e = random.choice([1., 2., 3., 4.])
      d_dmod = d_target_d_modulus_fcalc(fo, fc, a, b, k, e, centric)
      d_du = d_dmod/(2.*modulus_fc)
      u = modulus_fc*modulus_fc
      w_u, yo_u = effective_observation_in_intensity_space(
        fo*fo, u, a, b, k, e, centric)
      worst = max(worst, abs(w_u*(yo_u - u) + d_du)/max(abs(d_du), 1e-12))
  assert worst < 1e-9, worst
  print("\tthe F^2 form gives the same equations, so the observable can stay")


def exercise_cxx_agrees_with_the_python_reduction():
  """ The C++ in ml_target.h against the arithmetic written out above. """
  from smtbx.refinement import least_squares as ls_ext
  ext = ls_ext.ext
  random.seed(2)
  worst_w = worst_yo = 0.
  for centric in (False, True):
    for _ in range(200):
      fo = random.uniform(0.5, 60.)
      modulus_fc = random.uniform(0.5, 60.)
      a = random.uniform(0.2, 1.)
      b = random.uniform(0.5, 40.)
      k = random.uniform(0.5, 2.)
      e = random.choice([1., 2., 3., 4.])
      u = modulus_fc*modulus_fc
      # sigma_fo_sq only sets the scale the near-zero guard measures against,
      # and every |Fc| here is well clear of it; a value of zero would make the
      # guard fall back on fo_sq alone, which is what it did before the guard
      # was scaled by the measurement.
      ok, yo_c, w_c = ext.ml_effective_observation(
        fo_sq=fo*fo, sigma_fo_sq=fo, observable=u, alpha=a, beta=b,
        scale_factor=k, epsilon=e, centric=centric)
      assert ok
      w_p, yo_p = effective_observation_in_intensity_space(
        fo*fo, u, a, b, k, e, centric)
      worst_w = max(worst_w, abs(w_c - w_p)/max(abs(w_p), 1e-12))
      worst_yo = max(worst_yo, abs(yo_c - yo_p)/max(abs(yo_p), 1e-12))
  assert worst_w < 1e-12, worst_w
  assert worst_yo < 1e-10, worst_yo
  print("\tthe C++ reduction matches the arithmetic")


def exercise_mli_approaches_mlf_as_the_error_vanishes():
  """ Intensity likelihood is the Rice model convolved with the measurement
  error, and that convolution is the whole of the difference.

  A plain change of variable to Io = Fo^2 would multiply the density by
  1/(2Fo), which contains no Fc and therefore changes no derivative - such a
  target refines exactly as MLF does, which is why cctbx's mli.h, a copy of
  mlf.h, is a placeholder rather than a nearly-finished target.

  So the test that matters is that the convolved gradient tends to the
  amplitude one as sigma tends to zero, and departs from it when the error is
  real.
  """
  from smtbx.refinement import least_squares as ls_ext
  ext = ls_ext.ext
  a, b, k, e = 0.85, 6., 1., 1.
  fo, modulus_fc = 9., 7.
  u = modulus_fc*modulus_fc
  exact_mlf = d_target_d_modulus_fcalc(
    fo, complex(modulus_fc, 0.), a, b, k, e, False)
  previous = None
  for sigma in (8., 4., 2., 1., 0.5, 0.25):
    ok, d = ext.mli_d_target_d_modulus_fc(
      fo_sq=fo*fo, sigma_fo_sq=sigma, observable=u, alpha=a, beta=b,
      scale_factor=k, epsilon=e, centric=False)
    assert ok
    gap = abs(d - exact_mlf)
    if previous is not None:
      # halving sigma must bring it closer, and it does so quadratically
      assert gap < previous, (sigma, gap, previous)
    previous = gap
  assert previous < 1e-3, previous
  # and with a real error the two genuinely differ
  ok, d_big = ext.mli_d_target_d_modulus_fc(
    fo_sq=fo*fo, sigma_fo_sq=8., observable=u, alpha=a, beta=b,
    scale_factor=k, epsilon=e, centric=False)
  assert ok and abs(d_big - exact_mlf) > 1e-2, (d_big, exact_mlf)
  print("\tthe intensity target converges on the amplitude one as sigma -> 0")


def exercise_mli_handles_a_negative_intensity():
  """ The reason for an intensity target: a measured Io below zero is data.

  The quadrature integrates over the true intensity, which cannot be negative;
  a negative measurement simply shifts where the error distribution sits. It
  must produce a finite gradient rather than a domain error.
  """
  from smtbx.refinement import least_squares as ls_ext
  ext = ls_ext.ext
  ok, d = ext.mli_d_target_d_modulus_fc(
    fo_sq=-3.5, sigma_fo_sq=4., observable=25., alpha=0.8, beta=9.,
    scale_factor=1., epsilon=1., centric=False)
  assert ok, "a negative measured intensity was refused"
  assert d == d and abs(d) < 1e6, d      # not NaN, not divergent
  print("\ta negative measured intensity is usable, not a domain error")


def exercise_a_refinement_builds_end_to_end():
  """ The driver assembles and produces usable normal equations.

  Not a claim about refinement quality - only that alpha/beta estimation, the
  fixed-scale accumulator and the loop fit together, that both targets give
  finite non-trivial equations, and that neither of them is silently the
  least-squares ones.
  """
  from cctbx.development import random_structure
  from cctbx.array_family import flex
  from scitbx.lstbx import normal_eqns
  from smtbx.refinement import least_squares, constraints, sigma_a
  import smtbx.utils

  if not sigma_a.is_available():
    print("\tskipped: the alpha/beta estimator is not available here")
    return

  xs = random_structure.xray_structure(
    space_group_symbol="P 21 21 21",
    elements=["C"]*14 + ["O"]*3, u_iso=0.045, random_u_iso=True)
  fo_sq = xs.structure_factors(d_min=1.0).f_calc().norm()
  fo_sq = fo_sq.customized_copy(
    sigmas=flex.double([max(0.5, 0.04*v) for v in fo_sq.data()]))
  fo_sq.set_observation_type_xray_intensity()
  free = fo_sq.generate_r_free_flags(fraction=0.1, max_free=None).data()

  def build(engine, target):
    xs2 = xs.deep_copy_scatterers()
    random.seed(7)
    for sc in xs2.scatterers():
      sc.site = tuple(c + random.uniform(-0.004, 0.004) for c in sc.site)
      sc.flags.set_grad_site(True)
      sc.flags.set_use_u_iso(True)
      sc.flags.set_grad_u_iso(True)
    rep = constraints.reparametrisation(
      structure=xs2, constraints=[],
      connectivity_table=smtbx.utils.connectivity_table(xs2))
    klass = least_squares.crystallographic_ls_class(
      non_linear_ls_with_separable_scale_factor=engine,
      n_parameters=rep.n_independents, may_parallelise=False)
    kwds = dict(weighting_scheme=least_squares.unit_weighting())
    if target is not None:
      kwds.update(ml_target=target, ml_free_flags=free,
                  ml_free_reflections_per_bin=40)
    ls = klass(fo_sq.as_xray_observations(), rep, **kwds)
    ls.build_up()
    return list(ls.step_equations().right_hand_side())

  b_ls = build(normal_eqns.non_linear_ls_with_separable_scale_factor_BLAS_2,
               None)
  assert len(b_ls) and max(abs(v) for v in b_ls) > 0
  for target in ("mlf", "mli"):
    b = build(normal_eqns.non_linear_ls_with_fixed_scale_factor_BLAS_2, target)
    assert len(b) == len(b_ls)
    assert all(v == v and abs(v) < 1e30 for v in b), target
    assert max(abs(v) for v in b) > 0, target
    departure = (max(abs(x - y) for x, y in zip(b, b_ls))
                 /max(max(abs(v) for v in b_ls), 1e-12))
    assert departure > 1e-6, (target, "identical to least squares")
  print("\tboth targets build end to end and differ from least squares")


def exercise_blas_3_gives_the_same_equations():
  """ The level-3 accumulator must form the same equations as level 2

  The two reach them differently: level 2 applies a rank-1 update per
  reflection, level 3 buffers rows of sqrt(W)J and folds them in with a syrk,
  summing the same terms in a different order and through a different library.
  The maximum-likelihood path is the more exacting test of that, its weights
  spanning a far wider range than a least-squares weighting scheme's because
  the F^2-space weight carries a factor 1/(4|Fc|^2).

  fast_linalg initialises itself on import, so a test which does not import it
  gets the level-2 accumulator whatever it asks for.
  """
  from cctbx.development import random_structure
  from cctbx.array_family import flex
  from scitbx.lstbx import normal_eqns
  from smtbx.refinement import least_squares, constraints, sigma_a
  import smtbx.utils

  if not sigma_a.is_available():
    print("\tskipped: the alpha/beta estimator is not available here")
    return
  try:
    from fast_linalg import env
  except ImportError:
    print("\tskipped: fast_linalg is not available here")
    return
  if not env.initialised:
    print("\tskipped: no OpenBLAS behind fast_linalg")
    return

  xs = random_structure.xray_structure(
    space_group_symbol="P 21 21 21",
    elements=["C"]*14 + ["O"]*3, u_iso=0.045, random_u_iso=True)
  fo_sq = xs.structure_factors(d_min=1.0).f_calc().norm()
  fo_sq = fo_sq.customized_copy(
    sigmas=flex.double([max(0.5, 0.04*v) for v in fo_sq.data()]))
  fo_sq.set_observation_type_xray_intensity()
  free = fo_sq.generate_r_free_flags(fraction=0.1, max_free=None).data()

  def build(engine, target, parallelise):
    xs2 = xs.deep_copy_scatterers()
    random.seed(11)
    for sc in xs2.scatterers():
      sc.site = tuple(c + random.uniform(-0.004, 0.004) for c in sc.site)
      sc.flags.set_grad_site(True)
      sc.flags.set_use_u_iso(True)
      sc.flags.set_grad_u_iso(True)
    rep = constraints.reparametrisation(
      structure=xs2, constraints=[],
      connectivity_table=smtbx.utils.connectivity_table(xs2))
    klass = least_squares.crystallographic_ls_class(
      non_linear_ls_with_separable_scale_factor=engine,
      n_parameters=rep.n_independents, may_parallelise=parallelise)
    kwds = dict(weighting_scheme=least_squares.unit_weighting(),
                may_parallelise=parallelise)
    if target is not None:
      kwds.update(ml_target=target, ml_free_flags=free,
                  ml_free_reflections_per_bin=40)
    ls = klass(fo_sq.as_xray_observations(), rep, **kwds)
    ls.build_up()
    eqs = ls.step_equations()
    return (list(eqs.right_hand_side()),
            list(eqs.normal_matrix_packed_u()))

  for target in ("mlf", "mli"):
    b2, a2 = build(normal_eqns.non_linear_ls_with_fixed_scale_factor_BLAS_2,
                   target, False)
    b3, a3 = build(normal_eqns.non_linear_ls_with_fixed_scale_factor_BLAS_3,
                   target, False)
    assert len(b2) == len(b3) and len(a2) == len(a3), target
    scale_b = max(max(abs(v) for v in b2), 1e-12)
    scale_a = max(max(abs(v) for v in a2), 1e-12)
    db = max(abs(x - y) for x, y in zip(b2, b3))/scale_b
    da = max(abs(x - y) for x, y in zip(a2, a3))/scale_a
    # a reordered sum in double precision, not a bit-identical one
    assert db < 1e-10, (target, "right hand side", db)
    assert da < 1e-10, (target, "normal matrix", da)

    # and threading must not change them either
    b3t, a3t = build(normal_eqns.non_linear_ls_with_fixed_scale_factor_BLAS_3,
                     target, True)
    db = max(abs(x - y) for x, y in zip(b3, b3t))/scale_b
    da = max(abs(x - y) for x, y in zip(a3, a3t))/scale_a
    assert db < 1e-10, (target, "threaded right hand side", db)
    assert da < 1e-10, (target, "threaded normal matrix", da)
  print("\tthe BLAS 3 accumulator gives the same equations, threaded or not")


def exercise_the_accumulator_family_is_chosen_and_enforced():
  """ ml_target must select a fixed-scale accumulator, and refuse the other

  A separable-scale accumulator under a likelihood recovers a scale which is
  not free, alpha and beta already containing it, and converges somewhere
  plausible and wrong rather than failing: hence both halves of the check.
  """
  from scitbx.lstbx import normal_eqns
  from smtbx.refinement import least_squares

  for target in ("mlf", "mli"):
    klass = least_squares.crystallographic_ls_class(
      n_parameters=100000, may_parallelise=False, ml_target=target)
    chosen = klass.non_linear_ls_engine.__name__
    assert "fixed_scale" in chosen, (target, chosen)
  # and without it, the family is unchanged
  klass = least_squares.crystallographic_ls_class(
    n_parameters=100000, may_parallelise=False)
  assert "separable" in klass.non_linear_ls_engine.__name__

  # an explicitly passed engine is still honoured, and then the guard is what
  # stands between the caller and a quietly wrong refinement
  klass = least_squares.crystallographic_ls_class(
    non_linear_ls_with_separable_scale_factor=
      normal_eqns.non_linear_ls_with_separable_scale_factor_BLAS_2,
    n_parameters=100, may_parallelise=False)
  instance = klass.__new__(klass)
  instance.ml_target = "mlf"
  try:
    instance.check_maximum_likelihood_accumulator()
  except RuntimeError as e:
    assert "fixed-scale" in str(e), str(e)
  else:
    raise AssertionError("a separable-scale accumulator was accepted")
  print("\tml_target picks the fixed-scale accumulator and refuses the other")


def run():
  exercise_effective_observation_reproduces_the_gradient()
  exercise_the_weight_is_positive()
  exercise_intensity_space_form_is_the_same_equations()
  exercise_cxx_agrees_with_the_python_reduction()
  exercise_mli_approaches_mlf_as_the_error_vanishes()
  exercise_mli_handles_a_negative_intensity()
  exercise_a_refinement_builds_end_to_end()
  exercise_blas_3_gives_the_same_equations()
  exercise_the_accumulator_family_is_chosen_and_enforced()
  print("OK")


if __name__ == '__main__':
  run()
