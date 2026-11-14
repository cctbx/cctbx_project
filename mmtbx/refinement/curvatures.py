"""
Second derivatives (curvatures) of reciprocal-space refinement targets.

Layer 1: per-reflection derivatives with respect to the real and imaginary
parts of the calculated structure factor, in the (da, db) / (daa, dbb, dab)
convention of cctbx/xray/targets/common_results.h, so that the results feed
cctbx.xray.structure.grads_and_curvs_target_simple.
"""
from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex

def rice_radial_derivatives(f_obs, f_calc_abs, a, sigma, centric_flags):
  """
  First and second derivatives, with respect to |Fcalc|, of the Rice
  (acentric) and folded-normal (centric) negative log-likelihood

    acentric: a^2 |Fc|^2 / sigma - ln I0(2 a |Fo| |Fc| / sigma) + const
    centric:  a^2 |Fc|^2 / (2 sigma) - ln cosh(a |Fo| |Fc| / sigma) + const

  which is the per-reflection form of the mlf target
  (cctbx/xray/targets/mlf.h, target_one_h) with a = k alpha and
  sigma = epsilon k^2 beta.

  Parameters
  ----------
  f_obs : flex.double
    Observed amplitudes.
  f_calc_abs : flex.double
    Calculated amplitudes |Fcalc|.
  a : flex.double
    Multiplier of Fcalc in the mean of the distribution.
  sigma : flex.double
    Variance parameter of the distribution.
  centric_flags : flex.bool
    True for centric reflections.

  Returns
  -------
  t1, t2 : flex.double
    d T / d|Fc| and d^2 T / d|Fc|^2 per reflection. Both are zero where
    a <= 0, sigma <= 0 or |Fc| <= 0.
  """
  import scitbx.math
  n = f_obs.size()
  assert f_calc_abs.size() == n
  assert a.size() == n
  assert sigma.size() == n
  assert centric_flags.size() == n
  ok = (a > 0) & (sigma > 0) & (f_calc_abs > 0)
  safe_sigma = sigma.deep_copy()
  safe_sigma.set_selected(~ok, 1.0)
  fo_over_sigma = f_obs / safe_sigma
  # acentric
  x = 2 * a * fo_over_sigma * f_calc_abs
  m = scitbx.math.bessel_i1_over_i0(x)
  m_prime = scitbx.math.bessel_i1_over_i0_derivative(x)
  t1_acentric = 2 * a * a * f_calc_abs / safe_sigma \
              - 2 * a * fo_over_sigma * m
  t2_acentric = 2 * a * a / safe_sigma \
              - flex.pow2(2 * a * fo_over_sigma) * m_prime
  # centric
  x = a * fo_over_sigma * f_calc_abs
  th = flex.tanh(x)
  t1_centric = a * a * f_calc_abs / safe_sigma - a * fo_over_sigma * th
  t2_centric = a * a / safe_sigma \
             - flex.pow2(a * fo_over_sigma) * (1 - th * th)
  t1 = t1_acentric.deep_copy()
  t1.set_selected(centric_flags, t1_centric.select(centric_flags))
  t2 = t2_acentric.deep_copy()
  t2.set_selected(centric_flags, t2_centric.select(centric_flags))
  t1.set_selected(~ok, 0.0)
  t2.set_selected(~ok, 0.0)
  return t1, t2

def gradients_from_radial(f_calc, t1):
  """
  Gradient (da, db) of a target T(|Fc|) with respect to (Re Fc, Im Fc),
  as the complex number da + i db, from the radial derivative t1 = dT/d|Fc|.
  Zero where |Fc| = 0.
  """
  fc_abs = flex.abs(f_calc)
  ok = fc_abs > 0
  safe = fc_abs.deep_copy()
  safe.set_selected(~ok, 1.0)
  factor = t1 / safe
  factor.set_selected(~ok, 0.0)
  return f_calc * factor

def hessians_from_radial(f_calc, t1, t2):
  """
  Second derivatives (daa, dbb, dab) of a target T(|Fc|) with respect to
  (Re Fc, Im Fc) from the radial derivatives t1 = dT/d|Fc| and
  t2 = d^2 T/d|Fc|^2:

    daa = t2 A^2/|Fc|^2 + t1 B^2/|Fc|^3
    dbb = t2 B^2/|Fc|^2 + t1 A^2/|Fc|^3
    dab = (t2 - t1/|Fc|) A B/|Fc|^2

  Zero where |Fc| = 0.

  Returns
  -------
  flex.vec3_double
    One (daa, dbb, dab) triple per reflection.
  """
  fc_abs = flex.abs(f_calc)
  ok = fc_abs > 0
  safe = fc_abs.deep_copy()
  safe.set_selected(~ok, 1.0)
  aa = flex.pow2(flex.real(f_calc) / safe)
  bb = flex.pow2(flex.imag(f_calc) / safe)
  ab = flex.real(f_calc) * flex.imag(f_calc) / (safe * safe)
  t1_over_fc = t1 / safe
  daa = t2 * aa + t1_over_fc * bb
  dbb = t2 * bb + t1_over_fc * aa
  dab = (t2 - t1_over_fc) * ab
  daa.set_selected(~ok, 0.0)
  dbb.set_selected(~ok, 0.0)
  dab.set_selected(~ok, 0.0)
  return flex.vec3_double(daa, dbb, dab)

def mlf_hessians_work(f_obs, r_free_flags, f_calc, alpha, beta, scale_factor,
                      epsilons, centric_flags):
  """
  Per-reflection second derivatives of the mlf target for the working set,
  normalised like the gradients of cctbx.xray.mlf_target_and_gradients
  (divided by the number of working reflections). The arguments are those
  of mlf_target_and_gradients. Contributions are zero where
  mlf.h's d_target_one_h_over_fc returns zero (alpha <= 0, beta <= 1e-3 or
  |Fc| = 0).

  Returns
  -------
  flex.vec3_double
    (daa, dbb, dab) for each working reflection, in reflection order.
  """
  k = scale_factor
  if k <= 0:
    k = 1.0
  a = alpha * k
  sigma = epsilons * beta * (k * k)
  ok = (alpha > 0) & (beta > 1.e-3)
  a.set_selected(~ok, 0.0)
  fc_abs = flex.abs(f_calc)
  t1, t2 = rice_radial_derivatives(
    f_obs=f_obs, f_calc_abs=fc_abs, a=a, sigma=sigma,
    centric_flags=centric_flags)
  hessians = hessians_from_radial(f_calc=f_calc, t1=t1, t2=t2)
  if r_free_flags.size() == 0:
    work = flex.bool(f_obs.size(), True)
  else:
    work = ~r_free_flags
  n_work = work.count(True)
  assert n_work > 0
  return hessians.select(work) * (1.0 / n_work)
