"""
Second derivatives (curvatures) of reciprocal-space refinement targets.

Layer 1: per-reflection derivatives with respect to the real and imaginary
parts of the calculated structure factor, in the (da, db) / (daa, dbb, dab)
convention of cctbx/xray/targets/common_results.h, so that the results feed
cctbx.xray.structure.grads_and_curvs_target_simple.

Layer 2: per-atom Gauss-Newton (expected) curvatures of the data target
with respect to atomic parameters, from the per-reflection weights of
gauss_newton_weights and the structure-factor derivatives of
grads_and_curvs_target_simple.
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

def gauss_newton_weights(fmodel):
  """
  Per-reflection weight of the Gauss-Newton curvature of fmodel's data
  target: the expected curvature with respect to an atomic parameter p is

    sum_h  weights_h |d F_calc,h / d p|^2

  with weights_h = w_h hbar_h k_h^2, where w_h is the target's reflection
  weight, hbar_h the expected radial curvature (half the per-reflection
  Fisher information for acentric reflections, the full one for centric)
  and k_h the overall scale of F_model. Zero on the test set. Implemented
  for the targets ls_wunit_k1 and ml, which drops reflections with
  alpha <= 0 or beta <= 1e-3 as cctbx/xray/targets/mlf.h does.

  Parameters
  ----------
  fmodel : mmtbx.f_model.manager

  Returns
  -------
  flex.double
    One weight per reflection, in the order of fmodel.f_model().
  """
  from libtbx.utils import Sorry
  n = fmodel.f_model().data().size()
  work = ~fmodel.r_free_flags().data()
  n_work = work.count(True)
  assert n_work > 0
  if fmodel.target_name == "ls_wunit_k1":
    fo = fmodel.f_obs_work().data()
    k = fmodel.scale_k1_w()
    w = 1.0 / flex.sum(fo * fo)
    hbar = flex.double(n, w * k * k)
    hbar.set_selected(fmodel.f_obs().centric_flags().data(), 2 * w * k * k)
  elif fmodel.target_name == "ml":
    alpha, beta = fmodel.alpha_beta()
    k = fmodel.scale_ml_wrapper()
    if k is None or k <= 0:
      k = 1.0
    a = alpha.data() * k
    sigma = fmodel.epsilons * beta.data() * (k * k)
    ok = (alpha.data() > 0) & (beta.data() > 1.e-3) & (sigma > 0)
    safe_sigma = sigma.deep_copy()
    safe_sigma.set_selected(~ok, 1.0)
    hbar = a * a / safe_sigma * (1.0 / n_work)
    hbar.set_selected(~ok, 0.0)
  else:
    raise Sorry("Gauss-Newton weights are not available for target %s" %
      fmodel.target_name)
  k_eff = fmodel.k_isotropic() * fmodel.k_anisotropic()
  weights = hbar * k_eff * k_eff
  weights.set_selected(~work, 0.0)
  return weights

def gauss_newton_diagonal_sites_u_iso(xray_structure, miller_indices,
                                      weights):
  """
  Diagonal of the Gauss-Newton curvature sum_h weights_h |dF_h/dp|^2 for
  the Cartesian coordinates and the u_iso of every scatterer.

  The 3x3 block of each atom is accumulated in fractional coordinates by
  cctbx.xray.structure.grads_and_curvs_target_simple and transformed to
  Cartesian coordinates, so the result is exact for any cell. Scatterers
  must be isotropic and on general positions.

  Parameters
  ----------
  xray_structure : cctbx.xray.structure
  miller_indices : flex.miller_index
  weights : flex.double
    One weight per reflection, see gauss_newton_weights.

  Returns
  -------
  sites, u_iso : flex.double
    The 3 n curvatures of the Cartesian coordinates (x, y, z of atom 0,
    then atom 1, ...) and the n curvatures of u_iso.
  """
  n = xray_structure.scatterers().size()
  assert weights.size() == miller_indices.size()
  assert xray_structure.site_symmetry_table().special_position_indices() \
    .size() == 0
  for sc in xray_structure.scatterers():
    assert sc.flags.use_u_iso() and not sc.flags.use_u_aniso()
  zero = flex.double(weights.size(), 0)
  gact = xray_structure.grads_and_curvs_target_simple(
    miller_indices=miller_indices,
    da_db=flex.complex_double(zero, zero),
    daa_dbb_dab=flex.vec3_double(weights, weights, zero))
  n_params = 7  # site (3), u_iso, occupancy, fp, fdp
  curvs = gact.curvs
  assert curvs.size() == n_params * n
  curvs.reshape(flex.grid(n, n_params))
  u_iso = curvs.matrix_copy_column(3)
  # Cartesian diagonal d_i = sum_ab F_ai F_bi H_ab from the packed
  # fractional block (H_11, H_12, H_13, H_22, H_23, H_33), where x_frac =
  # F x_cart
  blocks = gact.site_curvs
  assert blocks.size() == 6 * n
  blocks.reshape(flex.grid(n, 6))
  hf = [blocks.matrix_copy_column(k) for k in range(6)]
  f = xray_structure.unit_cell().fractionalization_matrix()
  pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
  columns = []
  for i in range(3):
    d = flex.double(n, 0)
    for k, (a, b) in enumerate(pairs):
      c = f[3 * a + i] * f[3 * b + i]
      if a != b: c *= 2
      d += hf[k] * c
    columns.append(d)
  sites = flex.vec3_double(columns[0], columns[1], columns[2]).as_double()
  return sites, u_iso
