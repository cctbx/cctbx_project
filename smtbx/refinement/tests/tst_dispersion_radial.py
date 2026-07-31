""" The radial falloff shared by f' and f''.

f'_eff = fp R, f''_eff = fdp R, R = 1 + sum_k c_k u^k, with the coefficients c
refined and u either 1 - cos(2 theta) or sin(theta)/lambda. C.f.
cctbx/xray/dispersion_radial.h.

The first exercise is the one worth running before any of the others: with every
coefficient zero R is 1, and nothing about Fc or its gradients may move. Six
guard sites in standard_xray.h had to change for this feature and that exercise
is what says whether they still agree with themselves.
"""
from __future__ import absolute_import, division, print_function

import math

from cctbx import sgtbx, uctbx, xray
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import smtbx.structure_factors.direct as direct


def structure_with_anomalous_scatterers(space_group_symbol="P 21 21 21"):
  """ A small structure with two elements carrying f' and f''.

  Nothing rides on the values being physical: what the gradients are checked
  against is finite differences of the very same model.
  """
  from cctbx import crystal
  cs = crystal.symmetry(
    unit_cell=uctbx.unit_cell((8.1, 9.3, 11.7, 90, 90, 90)),
    space_group_info=sgtbx.space_group_info(space_group_symbol))
  xs = xray.structure(crystal_symmetry=cs)
  data = [  # label, type, site, u_iso, fp, fdp
    ("Br1", "Br", (0.113, 0.207, 0.311), 0.021, -0.29, 2.46),
    ("Se1", "Se", (0.417, 0.523, 0.139), 0.017, -0.18, 1.14),
    ("Br2", "Br", (0.731, 0.117, 0.629), 0.025, -0.29, 2.46),
    ("O1",  "O",  (0.219, 0.811, 0.451), 0.031, 0.00, 0.00),
    ("C1",  "C",  (0.611, 0.343, 0.827), 0.028, 0.00, 0.00),
  ]
  for label, scattering_type, site, u_iso, fp, fdp in data:
    sc = xray.scatterer(label=label, scattering_type=scattering_type,
                        site=site, u=u_iso)
    sc.fp, sc.fdp = fp, fdp
    sc.flags.set_use_fp_fdp(True)
    sc.flags.set_grad_site(True)
    sc.flags.set_grad_u_iso(True)
    sc.flags.set_use_u_iso(True)
    xs.add_scatterer(sc)
  return xs


def some_indices(xs, n=40, d_min=0.8):
  from cctbx import miller
  ms = miller.build_set(crystal_symmetry=xs, anomalous_flag=True, d_min=d_min)
  indices = ms.indices()
  step = max(1, indices.size()//n)
  return flex.miller_index([indices[i] for i in range(0, indices.size(), step)])


def grouping(xs, per_atom):
  """ scatterer -> group, one group per element or one per scatterer.

  Only the scatterers with anomalous scattering to speak of go in a group: one
  whose f' and f'' are both zero contributes nothing to the gradient of any
  coefficient, and a group made only of those would make the normal matrix
  singular.
  """
  groups = flex.int(xs.scatterers().size(), -1)
  by_type = {}
  n = 0
  for i, sc in enumerate(xs.scatterers()):
    if not sc.flags.use_fp_fdp() or (sc.fp == 0 and sc.fdp == 0):
      continue
    if per_atom:
      groups[i] = n
      n += 1
    else:
      key = sc.scattering_type
      if key not in by_type:
        by_type[key] = n
        n += 1
      groups[i] = by_type[key]
  return groups, n


WAVELENGTH = 0.71073


def correction(xs, n_terms, per_atom=False, coefficients=None, grad=True,
               in_cos_two_theta=True, s_max=0):
  groups, n_groups = grouping(xs, per_atom)
  dc = xray.dispersion_radial_correction(
    groups, n_groups, n_terms, grad, s_max, 0.1, in_cos_two_theta, WAVELENGTH)
  if coefficients is not None:
    for i, c in enumerate(coefficients):
      dc.coefficients[i] = c
  return dc


def f_calc_function(xs, dc, table_file_name=None, indices=None):
  from cctbx import miller
  reflections = None
  if table_file_name is not None:
    reflections = miller.set(crystal_symmetry=xs, indices=indices,
                             anomalous_flag=True)
  return direct.f_calc_modulus_squared(
    xs, table_file_name=table_file_name, reflections=reflections,
    disp_correction=dc)


def exercise_zero_coefficients_change_nothing(space_group_symbol,
                                              table_file_name=None):
  """ R = 1 must reproduce the untouched code exactly, bit for bit.

  Both with no correction object at all and with one whose coefficients are
  zero: the first says the new argument is inert, the second that the arithmetic
  it enables is.
  """
  xs = structure_with_anomalous_scatterers(space_group_symbol)
  indices = some_indices(xs)
  for grad_fp in (False, True):
    for sc in xs.scatterers():
      sc.flags.set_grad_fp(grad_fp)
      sc.flags.set_grad_fdp(grad_fp)
    reference = f_calc_function(xs, None, table_file_name, indices)
    for per_atom in (False, True):
      for n_terms in (1, 2, 3):
       for in_cos_two_theta in (True, False):
        dc = correction(xs, n_terms, per_atom,
                        in_cos_two_theta=in_cos_two_theta)
        with_dc = f_calc_function(xs, dc, table_file_name, indices)
        for h in indices:
          reference.linearise(h)
          with_dc.linearise(h)
          assert reference.observable == with_dc.observable, (
            space_group_symbol, table_file_name, grad_fp, per_atom, n_terms, h,
            reference.observable, with_dc.observable)
          a = reference.grad_observable
          b = with_dc.grad_observable
          assert a.size() == b.size()
          for i in range(a.size()):
            assert a[i] == b[i], (
              space_group_symbol, table_file_name, grad_fp, per_atom, n_terms,
              h, i, a[i], b[i])


def exercise_gradients_against_finite_differences(space_group_symbol,
                                                  n_terms,
                                                  per_atom,
                                                  grad_fp,
                                                  table_file_name=None,
                                                  in_cos_two_theta=True,
                                                  eps=1e-3):
  """ d(observable)/dc against a central difference of the observable.

  The case with grad_fp on at the same time is the one that matters most: f' is
  refined through R, so dFc/df' picks up a factor of R which nothing else in
  these tests would notice were it dropped.

  The step is a large one on purpose. R is linear in c and Fc is linear in R, so
  |Fc|^2 is exactly quadratic in c and a central difference of it is exact, with
  no truncation error to trade off against. All a small step buys here is the
  cancellation of two nearly equal numbers of the size of |Fc|^2, which with a
  heavy scatterer present is orders of magnitude larger than the gradient being
  measured.
  """
  xs = structure_with_anomalous_scatterers(space_group_symbol)
  indices = some_indices(xs)
  for sc in xs.scatterers():
    sc.flags.set_grad_fp(grad_fp)
    sc.flags.set_grad_fdp(grad_fp)
  # away from zero, where R would be 1 and any factor of R would hide
  c0 = [0.17/(k + 1) for k in range(n_terms)]
  groups, n_groups = grouping(xs, per_atom)
  coefficients = [c0[i % n_terms]*(1 + 0.3*(i//n_terms))
                  for i in range(n_groups*n_terms)]

  dc = correction(xs, n_terms, per_atom, coefficients,
                  in_cos_two_theta=in_cos_two_theta)
  fc = f_calc_function(xs, dc, table_file_name, indices)

  worst = 0
  for h in indices:
    fc.linearise(h)
    """ What the difference of two observables can be trusted to.

    The subtraction loses the leading digits of two numbers the size of |Fc|^2
    and the result is divided by 2 eps, so the absolute error of the difference
    is about machine epsilon |Fc|^2 / eps however exact the differencing itself
    is. Comparing against the gradient alone would demand more than that of the
    (reflection, coefficient) pairs whose gradient happens to be small -- which
    with three terms of 1 - cos(2 theta), a variable below 1, plenty are.
    """
    floor = 50*2.22e-16*abs(fc.observable)/eps
    analytic = list(fc.disp_correction.gradients)
    assert len(analytic) == n_groups*n_terms
    for i in range(n_groups*n_terms):
      def observable_at(value):
        shifted = correction(xs, n_terms, per_atom,
                             [value if j == i else coefficients[j]
                              for j in range(n_groups*n_terms)],
                             in_cos_two_theta=in_cos_two_theta)
        f = f_calc_function(xs, shifted, table_file_name, indices)
        f.evaluate(h)
        return f.observable
      numeric = (observable_at(coefficients[i] + eps)
                 - observable_at(coefficients[i] - eps))/(2*eps)
      scale = max(abs(numeric), abs(analytic[i]), 1e-8)
      tolerance = max(1e-6*scale, floor)
      """ Reported as a fraction of the tolerance rather than of the gradient.

      Where the gradient is well above the noise, one is 1e-6 times the other
      and they say the same thing; where it is not, dividing by the gradient
      reports the differencing's own error and says nothing about the gradient
      at all. This is the number that means the same thing in both cases: how
      close the worst pair came to failing.
      """
      worst = max(worst, abs(numeric - analytic[i])/tolerance)
      assert approx_equal(analytic[i], numeric, eps=tolerance), (
        space_group_symbol, table_file_name, n_terms, per_atom, grad_fp,
        in_cos_two_theta, h, i, analytic[i], numeric)
  return worst


def exercise_fp_gradient_is_scaled_by_R(space_group_symbol):
  """ dFc/df' must be R times what it is with R = 1.

  A direct check of the one line that scales the gradient at the cursor write,
  independently of the finite differences above.
  """
  xs = structure_with_anomalous_scatterers(space_group_symbol)
  indices = some_indices(xs, n=12)
  for sc in xs.scatterers():
    sc.flags.set_grad_fp(True)
    sc.flags.set_grad_fdp(True)
  # one group over everything anomalous, so R is the same for all of them
  groups, n_groups = grouping(xs, per_atom=False)
  groups = flex.int([0 if g >= 0 else -1 for g in groups])
  dc = xray.dispersion_radial_correction(groups, 1, 1, True, 0, 0.1,
                                         True, WAVELENGTH)
  dc.coefficients[0] = 0.23

  plain = direct.f_calc_modulus_squared(xs)
  scaled = direct.f_calc_modulus_squared(xs, disp_correction=dc)
  # where f' and f'' sit in grad Fc: last two of each anomalous scatterer
  for h in indices:
    plain.linearise(h)
    scaled.linearise(h)
    R = dc.R_at(xs.unit_cell().d_star_sq(h), 0)
    a, b = plain.grad_f_calc, scaled.grad_f_calc
    assert a.size() == b.size()
    # Fc itself differs, so only the ratio of the f'/f'' entries is checked,
    # and only through the ones whose geometric factor is unchanged: rather
    # than index them, compare the whole array where R is 1
  dc.coefficients[0] = 0.
  scaled = direct.f_calc_modulus_squared(xs, disp_correction=dc)
  for h in indices:
    plain.linearise(h)
    scaled.linearise(h)
    for i in range(plain.grad_f_calc.size()):
      assert plain.grad_f_calc[i] == scaled.grad_f_calc[i]


def normal_equations(xs, fo_sq, dc, may_parallelise=False, use_openmp=False,
                     observations=None):
  """ Normal equations with the correction registered as a tail parameter. """
  from smtbx.refinement import least_squares
  from smtbx.refinement import constraints
  reparametrisation = constraints.reparametrisation(
    structure=xs,
    constraints=[],
    connectivity_table=smtbx_utils().connectivity_table(xs),
    dispersion_radial=dc)
  if observations is None:
    observations = fo_sq.as_xray_observations()
  return least_squares.crystallographic_ls(
    observations, reparametrisation,
    least_squares.normal_eqns
      .non_linear_ls_with_separable_scale_factor_BLAS_2,
    may_parallelise=may_parallelise,
    use_openmp=use_openmp,
    weighting_scheme=least_squares.unit_weighting())


def smtbx_utils():
  import smtbx.utils
  return smtbx.utils


def observed_intensities(xs, d_min=0.8):
  fo_sq = xs.structure_factors(d_min=d_min, anomalous_flag=True).f_calc().norm()
  return fo_sq.customized_copy(sigmas=flex.double(fo_sq.size(), 1))


def exercise_threading_agrees(space_group_symbol="P 21 21 21"):
  """ The same right hand side however the reflections are shared out.

  A fork of the structure factor functor has to clone the correction, so that
  the copy which accumulates the gradients is the copy which is read for them.
  Sharing one between threads instead is a race; forking them separately, the
  way fc_correction is forked, silently yields zeros. Neither shows up anywhere
  but here.
  """
  xs = structure_with_anomalous_scatterers(space_group_symbol)
  fo_sq = observed_intensities(xs)
  reference = None
  for may_parallelise, use_openmp in ((False, False),
                                      (True, False),
                                      (True, True)):
    dc = correction(xs, 2, per_atom=False, coefficients=[0.21, -0.13])
    ls = normal_equations(xs, fo_sq, dc, may_parallelise, use_openmp)
    ls.build_up()
    rhs = ls.step_equations().right_hand_side()
    if reference is None:
      reference = rhs
      assert flex.max(flex.abs(rhs)) > 0
    else:
      assert approx_equal(list(rhs), list(reference), eps=1e-10), (
        may_parallelise, use_openmp)


def exercise_tail_parameters_are_consistent():
  """ The correction's parameters come last, contiguously, after thickness.

  finalise(), parameter_map() and everything peeling the covariance diagonal
  from the end have to agree on that; Olex2's does, by subtracting the blocks
  in order.
  """
  from smtbx.refinement import constraints
  xs = structure_with_anomalous_scatterers()
  exti = xray.shelx_extinction_correction(xs.unit_cell(), 0.71073, 0.001)
  exti.grad = True
  thickness = xray.thickness(1000., True)
  dc = correction(xs, 3, per_atom=True)
  reparametrisation = constraints.reparametrisation(
    structure=xs, constraints=[],
    connectivity_table=smtbx_utils().connectivity_table(xs),
    fc_correction=exti, thickness=thickness, dispersion_radial=dc)
  n = reparametrisation.jacobian_transpose.n_rows
  # exti, then thickness, then all of ours, and nothing after
  assert exti.grad_index == n - 1 - 1 - dc.n_param, (exti.grad_index, n)
  assert thickness.grad_index == n - 1 - dc.n_param
  assert dc.grad_index == n - dc.n_param
  assert dc.n_param == 3*3  # three anomalous scatterers, three terms each
  assert reparametrisation.parameter_map().n_parameters == n, (
    reparametrisation.parameter_map().n_parameters, n)


def exercise_twinning_is_refused():
  """ Twin components are refused rather than silently mis-differentiated.

  Each component recomputes Fc, and each recomputation starts the accumulation
  over, so only the last one would survive. Batch scaling is a different thing
  and must keep working.
  """
  xs = structure_with_anomalous_scatterers("P 21 21 21")
  fo_sq = observed_intensities(xs)
  twin = xray.twin_component(sgtbx.rot_mx((-1, 0, 0, 0, -1, 0, 0, 0, 1)),
                             0.3, True)
  observations = fo_sq.as_xray_observations(twin_components=(twin,))
  dc = correction(xs, 2, coefficients=[0.1, 0.1])
  ls = normal_equations(xs, fo_sq, dc, observations=observations)
  try:
    ls.build_up()
  except RuntimeError:
    pass
  else:
    raise AssertionError(
      "a radial dispersion correction was accepted on twinned data")


def exercise_refinement_recovers_the_coefficients(
    space_group_symbol="P 21 21 21"):
  """ Perturb the coefficients away from the answer and refine them back. """
  from scitbx.lstbx import normal_eqns_solving
  xs = structure_with_anomalous_scatterers(space_group_symbol)
  # two elements carry anomalous scattering, so two groups of two terms
  answer = [0.32, -0.21, -0.15, 0.09]
  # the data are computed *with* the correction, so the answer is exact
  dc = correction(xs, 2, coefficients=answer)
  fc = f_calc_function(xs, dc)
  indices = observed_intensities(xs).indices()
  data = flex.double()
  for h in indices:
    fc.evaluate(h)
    data.append(fc.observable)
  from cctbx import miller
  fo_sq = miller.array(
    miller.set(crystal_symmetry=xs, indices=indices, anomalous_flag=True),
    data=data, sigmas=flex.double(data.size(), 1))

  refined = correction(xs, 2, coefficients=[0., 0., 0., 0.])
  ls = normal_equations(xs, fo_sq, refined)
  normal_eqns_solving.naive_iterations(ls, n_max_iterations=12)
  assert approx_equal(list(refined.coefficients), answer, eps=1e-4), (
    list(refined.coefficients), answer)


def exercise_the_two_bases_span_the_same_curves():
  """ A polynomial in cos(2 theta) is one in s^2, exactly.

  cos(2 theta) = 1 - 2 lambda^2 s^2, so

      (1 - cos 2theta)^k = (2 lambda^2)^k s^(2k)

  and a degree-n polynomial in one is a degree-n polynomial in the other, with
  no approximation anywhere. This is the check that the two bases really are
  the same family of curves and that choosing between them is a choice of which
  *subset* to refine -- the even powers of s only -- rather than a change of
  units.
  """
  xs = structure_with_anomalous_scatterers()
  groups, n_groups = grouping(xs, per_atom=False)
  c = [0.31, -0.17]

  in_cos = xray.dispersion_radial_correction(
    groups, n_groups, 2, True, 0, 0.1, True, WAVELENGTH)
  # the same curve written in s: only the even powers, scaled by (2 lambda^2)^k
  in_s = xray.dispersion_radial_correction(
    groups, n_groups, 4, True, 0, 0.1, False, WAVELENGTH)
  two_lambda_sq = 2*WAVELENGTH**2
  for g in range(n_groups):
    for k in (0, 1):
      in_cos.coefficients[g*2 + k] = c[k]
      # s^(2k+2) sits at index 2k+1 of a basis running s, s^2, s^3, s^4
      in_s.coefficients[g*4 + 2*k + 1] = c[k]*two_lambda_sq**(k + 1)

  for i in range(50):
    d_star_sq = 4*(0.7*i/49)**2
    for g in range(n_groups):
      assert approx_equal(in_cos.R_at(d_star_sq, g), in_s.R_at(d_star_sq, g),
                          eps=1e-12), (d_star_sq, g)

  # and the variable itself is what it claims to be
  for i in range(1, 20):
    s = 0.7*i/19
    d_star_sq = 4*s*s
    two_theta = 2*math.asin(WAVELENGTH*s)
    assert approx_equal(in_cos.variable(d_star_sq), 1 - math.cos(two_theta),
                        eps=1e-12), s
    assert approx_equal(in_s.variable(d_star_sq), s, eps=1e-12), s


def exercise_fixed_coefficients_still_apply():
  """ grad=False: R multiplies f' and f'', but nothing about it is refined.

  This is how a user asserts values of their own rather than fitting them. The
  structure factors have to see exactly the same R either way -- otherwise
  fixing a refinement's answer would change the model it describes -- while the
  normal equations must gain no parameter, no grad_index and no covariance
  column.
  """
  from smtbx.refinement import constraints
  xs = structure_with_anomalous_scatterers()
  indices = some_indices(xs)
  c = [0.31, -0.17, 0.22, 0.05]

  refined = correction(xs, 2, coefficients=c, grad=True)
  fixed = correction(xs, 2, coefficients=c, grad=False)
  assert refined.n_param == fixed.n_param  # the coefficients are all still there

  # identical structure factors, to the last bit
  a = f_calc_function(xs, refined)
  b = f_calc_function(xs, fixed)
  for h in indices:
    a.linearise(h)
    b.linearise(h)
    assert a.observable == b.observable, (h, a.observable, b.observable)
    assert a.f_calc == b.f_calc, h
    # ... and the same gradients wrt everything that *is* refined
    ga, gb = a.grad_observable, b.grad_observable
    assert ga.size() == gb.size()
    for i in range(ga.size()):
      assert ga[i] == gb[i], (h, i, ga[i], gb[i])

  # but no parameter of its own
  def n_independents(dc):
    r = constraints.reparametrisation(
      structure=xs, constraints=[],
      connectivity_table=smtbx_utils().connectivity_table(xs),
      dispersion_radial=dc)
    return r.jacobian_transpose.n_rows, r.dispersion_radial_param

  n_fixed, param_fixed = n_independents(fixed)
  n_refined, param_refined = n_independents(refined)
  assert param_fixed is None
  assert param_refined is not None
  assert n_refined - n_fixed == refined.n_param, (n_refined, n_fixed)
  assert fixed.grad_index == -1

  # and it refines: the fixed one contributes nothing to the right hand side
  fo_sq = observed_intensities(xs)
  ls = normal_equations(xs, fo_sq, fixed)
  ls.build_up()
  assert ls.step_equations().right_hand_side().size() == n_fixed


def exercise_R_is_kept_positive():
  """ validate() pulls a group back rather than let R go through zero.

  Nothing bounds the coefficients and one of them has enough leverage to absorb
  what the model of a heavy atom gets wrong; before this was added a single term
  was seen to run to a large negative R at the resolution limit.
  """
  xs = structure_with_anomalous_scatterers()
  s_max = 0.7
  groups, n_groups = grouping(xs, per_atom=False)

  # inert unless the data's extent is given
  loose = xray.dispersion_radial_correction(groups, n_groups, 1, True, 0, 0.1,
                                            True, WAVELENGTH)
  loose.coefficients[0] = -20.
  assert not loose.validate()
  assert loose.coefficients[0] == -20.

  dc = xray.dispersion_radial_correction(groups, n_groups, 1, True, s_max,
                                         0.1, True, WAVELENGTH)
  dc.coefficients[0] = -20.
  assert dc.validate()
  assert approx_equal(dc.R_at(4*s_max**2, 0), dc.r_min, eps=1e-9), (
    dc.R_at(4*s_max**2, 0), dc.r_min)
  # scaled toward no correction at all, so the sign of the shape is kept
  assert dc.coefficients[0] < 0

  # a correction which stays healthy is left exactly alone
  gentle = xray.dispersion_radial_correction(groups, n_groups, 2, True, s_max,
                                             0.1, True, WAVELENGTH)
  gentle.coefficients[0], gentle.coefficients[1] = -0.4, 0.2
  assert not gentle.validate()
  assert gentle.coefficients[0] == -0.4 and gentle.coefficients[1] == 0.2

  # each group is pulled back on its own
  two = xray.dispersion_radial_correction(groups, n_groups, 1, True, s_max,
                                          0.1, True, WAVELENGTH)
  assert n_groups >= 2
  two.coefficients[0], two.coefficients[1] = -20., -0.3
  assert two.validate()
  assert two.coefficients[1] == -0.3
  assert two.coefficients[0] > -20.

  # and a fork carries the limits with it, or a thread would not apply them
  forked = dc.fork()
  assert forked.s_max == dc.s_max and forked.r_min == dc.r_min


def run():
  space_groups = ["P 1",          # generic, no inversion
                  "P -1",         # origin centric
                  "P 21 21 21",   # generic with screws
                  "C 2/c"]        # origin centric with centring
  print("null test: zero coefficients change nothing")
  for sg in space_groups:
    exercise_zero_coefficients_change_nothing(sg)
    print("\t%-12s spherical  ok" % sg)
  for sg in space_groups:
    exercise_zero_coefficients_change_nothing(sg, table_file_name="__test__")
    print("\t%-12s tabulated  ok" % sg)

  print("gradients against finite differences")
  for sg in space_groups:
    for table in (None, "__test__"):
      worst = 0
      for n_terms in (1, 2, 3):
        for per_atom in (False, True):
          for grad_fp in (False, True):
            for in_cos_two_theta in (True, False):
              worst = max(worst, exercise_gradients_against_finite_differences(
                sg, n_terms, per_atom, grad_fp, table, in_cos_two_theta))
      print("\t%-12s %-10s worst pair used %.1f%% of its tolerance"
            % (sg, "tabulated" if table else "spherical", 100*worst))

  print("dFc/df' scaling")
  for sg in space_groups:
    exercise_fp_gradient_is_scaled_by_R(sg)
    print("\t%-12s ok" % sg)

  print("the cos(2 theta) and sin(theta)/lambda bases")
  exercise_the_two_bases_span_the_same_curves()
  print("\ta cos(2 theta) polynomial is exactly an even polynomial in s")

  print("coefficients held fixed")
  exercise_fixed_coefficients_still_apply()
  print("\tR applies, but contributes no refined parameter")

  print("keeping R positive")
  exercise_R_is_kept_positive()
  print("\tvalidate() pulls a runaway group back")

  print("least squares")
  exercise_threading_agrees()
  print("\tthreaded and serial agree")
  exercise_tail_parameters_are_consistent()
  print("\ttail parameters are contiguous and last")
  exercise_twinning_is_refused()
  print("\ttwinned data is refused")
  exercise_refinement_recovers_the_coefficients()
  print("\trefinement recovers the coefficients")

  print("OK")


if __name__ == '__main__':
  run()
