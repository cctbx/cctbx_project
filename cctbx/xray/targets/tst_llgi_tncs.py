from __future__ import absolute_import, division, print_function
""" Translational NCS in the LLGI functors (llgi.h tncs_frame): the model
term has mean a*ec and error variance (1 - a^2)*u on the E scale
normalised by RESN, and the data are normalised by t, the tNCS modulation
of the expected intensity. """
from cctbx.xray import ext
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from scitbx.math import bessel_ln_of_i0
from six.moves import range
import math

def rice_dobs_eeff(e_obs_sq, sig_e_obs_sq, centric_flags):
  """ Moment-matched (Dobs, Eeff), or (0.01, posterior rms) with no Rice
  solution, as mmtbx.refinement.llgi_hybrid.rice_dobs_eeff. """
  rm = ext.llgi_rice_moments(e_obs_sq=e_obs_sq, sig_e_obs_sq=sig_e_obs_sq,
    centric_flags=centric_flags)
  dobs = flex.double(); eeff = flex.double()
  for i in range(e_obs_sq.size()):
    if rm.valid[i]:
      dobs.append(math.sqrt(rm.dsqr[i])); eeff.append(rm.eeff[i])
    else:
      dobs.append(0.01); eeff.append(math.sqrt(max(rm.mu2[i], 0.)))
  return dobs, eeff, ~rm.valid

def make_case(n, centric, seed=0):
  """ Intensities and a model with tNCS: per reflection a (sigmaA), u
  (error-variance modulation), resn, fc, I, sig(I). """
  def frac(i, k): return ((seed + 7 * i + k) * 0.6180339887) % 1.0
  case = dict(n=n, centric=flex.bool(n, centric))
  case["a"] = flex.double([0.35 + 0.5 * frac(i, 1) for i in range(n)])
  case["u"] = flex.double([0.15 + 1.7 * frac(i, 2) for i in range(n)])
  case["resn"] = flex.double([0.8 + 0.6 * frac(i, 3) for i in range(n)])
  case["f_calc"] = flex.complex_double([
    complex((0.2 + 2.0 * frac(i, 4)) * math.cos(2 * math.pi * frac(i, 5)),
            (0.2 + 2.0 * frac(i, 4)) * math.sin(2 * math.pi * frac(i, 5)))
    for i in range(n)])
  case["i_obs"] = flex.double([(0.02 + 3.0 * frac(i, 6)) * case["resn"][i] ** 2
                               for i in range(n)])
  case["sig_i"] = flex.double([(0.05 + 0.8 * frac(i, 7)) * case["resn"][i] ** 2
                               for i in range(n)])
  return case

def t_consistent(case, scale=1.0, offset=0.0):
  """ a normalisation t > (1 - a^2) u: the expected-intensity modulation
  for a model term that is unmodulated, scaled and offset to give another
  valid t """
  return flex.double([scale * (case["a"][i] ** 2 + (1 - case["a"][i] ** 2) *
                      case["u"][i]) + offset for i in range(case["n"])])

def normalised_data(case, t, rice_kappa):
  """ data normalised by RESN^2 t: hybrid, Dobs, Feff = Eeff RESN sqrt(t) """
  esn = case["resn"] * case["resn"] * t
  e_obs_sq = case["i_obs"] / esn
  sig = case["sig_i"] / esn
  dobs, eeff, force_exact = rice_dobs_eeff(e_obs_sq, sig, case["centric"])
  hybrid = ext.llgi_hybrid(e_obs_sq=e_obs_sq, sig_e_obs_sq=sig,
    force_exact=force_exact, centric_flags=case["centric"],
    rice_kappa=rice_kappa)
  f_eff = eeff * case["resn"] * flex.sqrt(t)
  return dict(hybrid=hybrid, dobs=dobs, eeff=eeff, f_eff=f_eff,
              e_obs_sq=e_obs_sq, sig=sig)

def f_scale_target(case, t, data, hybrid=True, f_calc=None, sigmaa=None,
                   scatfrac=None, k=1.0, compute_gradients=False, tncs=True):
  n = case["n"]
  return ext.llgi_target_and_gradients(
    f_eff=data["f_eff"], r_free_flags=flex.bool(n, False),
    f_calc=case["f_calc"] if f_calc is None else f_calc, dobs=data["dobs"],
    sigmaa=case["a"] if sigmaa is None else sigmaa,
    scatfrac=flex.double(n, 1.) if scatfrac is None else scatfrac,
    scale_factor=k, teps=flex.double(n, 1.), resn=case["resn"],
    centric_flags=case["centric"], compute_gradients=compute_gradients,
    hybrid=data["hybrid"] if hybrid else None,
    tncs_u=case["u"] if tncs else flex.double(), tncs_t=t if tncs else flex.double())

def absolute_log_likelihood(case, target_per_reflection, data):
  """ log p(I) from a target computed in the t-normalised frame (exact
  path): LLGI + M(0,b) - ln(sqrt(2 pi) sig_I) [- ln(2 pi)/2 centric] """
  out = flex.double()
  null = data["hybrid"].null_log_z
  for i in range(case["n"]):
    v = -target_per_reflection[i] + null[i] - math.log(
      math.sqrt(2 * math.pi) * case["sig_i"][i])
    if case["centric"][i]: v -= 0.5 * math.log(2 * math.pi)
    out.append(v)
  return out

def reference_log_likelihood(case):
  """ log p(I | mean a*|fc|, variance resn^2 (1 - a^2) u) by the exact
  kernel with the per-reflection rescaling c = V + nu^2 (independent of
  tncs_frame) """
  n = case["n"]
  nu = flex.abs(case["f_calc"]) * case["a"]
  V = case["resn"] * case["resn"] * (1 - case["a"] * case["a"]) * case["u"]
  c = V + nu * nu
  a = flex.sqrt(nu * nu / c)
  ec = flex.double([1. if a[i] > 0 else 0. for i in range(n)])
  e_obs_sq = case["i_obs"] / c; sig = case["sig_i"] / c
  r = ext.llgi_exact_evaluate(e_obs_sq=e_obs_sq, sig_e_obs_sq=sig,
    e_calc=ec, sigmaa=a, centric_flags=case["centric"])
  h = ext.llgi_hybrid(e_obs_sq=e_obs_sq, sig_e_obs_sq=sig,
    force_exact=flex.bool(n, False), centric_flags=case["centric"], rice_kappa=0.)
  return absolute_log_likelihood(case, -r.ll, dict(hybrid=h))

def assert_gradient(ana, fd, label):
  """ finite-difference agreement to 1e-5 relative: the exact kernel's
  quadrature window follows the parameters, which limits the precision of
  a difference quotient """
  assert abs(ana - fd) <= 1e-5 * max(1., abs(fd)), (label, ana, fd)

def exercise_identity_without_modulation(centric):
  """ u = t = 1 reproduces the functors without tNCS """
  case = make_case(20, centric, seed=1)
  n = case["n"]
  one = flex.double(n, 1.)
  case["u"] = one
  sigmaa = case["a"] * 0.9; scatfrac = flex.double([0.9 + 0.005 * i for i in range(n)])
  assert flex.max(sigmaa) < 0.999   # the hybrid's exact path needs sigmaA < 0.999
  for rice_kappa in (0.0, 0.1, None):
    data = normalised_data(case, one, rice_kappa if rice_kappa is not None else 0.)
    use_h = rice_kappa is not None
    r0 = f_scale_target(case, one, data, hybrid=use_h, sigmaa=sigmaa,
      scatfrac=scatfrac, k=1.1, compute_gradients=True, tncs=False)
    r1 = f_scale_target(case, one, data, hybrid=use_h, sigmaa=sigmaa,
      scatfrac=scatfrac, k=1.1, compute_gradients=True, tncs=True)
    assert approx_equal(r0.target_per_reflection(), r1.target_per_reflection(), eps=1e-12)
    assert approx_equal(r0.gradients_work(), r1.gradients_work(), eps=1e-12)
    assert r0.n_exact() == r1.n_exact()
    sel = flex.bool([i % 3 != 0 for i in range(n)])
    args = dict(f_eff=data["f_eff"], selection=sel, f_calc=case["f_calc"],
      dobs=data["dobs"], sigmaa=sigmaa, scatfrac=scatfrac, scale_factor=1.1,
      teps=one, resn=case["resn"], centric_flags=case["centric"],
      hybrid=data["hybrid"] if use_h else None)
    s0 = ext.llgi_sigmaa_scatfrac_target_and_gradients(**args)
    s1 = ext.llgi_sigmaa_scatfrac_target_and_gradients(tncs_u=one, tncs_t=one, **args)
    assert approx_equal(s0.target(), s1.target(), eps=1e-12)
    assert approx_equal(s0.d_target_by_dsigmaa(), s1.d_target_by_dsigmaa(), eps=1e-12)
    assert approx_equal(s0.d_target_by_dscatfrac(), s1.d_target_by_dscatfrac(), eps=1e-12)
    assert s0.n_exact() == s1.n_exact()
    e_args = dict(e_eff=data["eeff"], selection=sel,
      e_model=flex.abs(case["f_calc"]), dobs=data["dobs"], sigmaa=sigmaa,
      centric_flags=case["centric"], hybrid=data["hybrid"] if use_h else None)
    e0 = ext.llgi_e_sigmaa_target_and_gradients(**e_args)
    e1 = ext.llgi_e_sigmaa_target_and_gradients(tncs_u=one, tncs_t=one, **e_args)
    assert approx_equal(e0.target(), e1.target(), eps=1e-12)
    assert approx_equal(e0.d_target_by_dsigmaa(), e1.d_target_by_dsigmaa(), eps=1e-12)
    assert e0.n_exact() == e1.n_exact()

def exercise_exact_path_is_absolute_likelihood(centric):
  """ exact path: the absolute log p(I) does not depend on the
  normalisation t and agrees with an independent evaluation """
  case = make_case(40, centric, seed=2)
  ref = reference_log_likelihood(case)
  for t in (t_consistent(case), t_consistent(case, scale=1.4, offset=0.3)):
    data = normalised_data(case, t, rice_kappa=0.)
    r = f_scale_target(case, t, data)
    absolute = absolute_log_likelihood(case, r.target_per_reflection(), data)
    assert r.n_exact() == case["n"]
    assert approx_equal(absolute, ref, eps=1e-8)

def ln_rice_over_wilson(f, mu, w, sigma_sq, centric):
  """ ln[Rice(f; mu, w) / Wilson(f; sigma_sq)], amplitude densities """
  if not centric:
    return (math.log(sigma_sq / w) - (f * f + mu * mu) / w + f * f / sigma_sq
            + bessel_ln_of_i0(2 * f * mu / w))
  x = f * mu / w
  return (0.5 * math.log(sigma_sq / w) - (f * f + mu * mu) / (2 * w)
          + f * f / (2 * sigma_sq) + x + math.log((1 + math.exp(-2 * x)) / 2))

def exercise_rice_path(centric):
  """ Rice path: Feff ~ Rice(Dobs a |fc|, resn^2 (t (1 - Dobs^2) +
  Dobs^2 (1 - a^2) u)) against the Wilson null with variance resn^2 t """
  case = make_case(40, centric, seed=3)
  t = t_consistent(case, scale=1.2, offset=0.1)
  data = normalised_data(case, t, rice_kappa=0.)
  r = f_scale_target(case, t, data, hybrid=False)
  for i in range(case["n"]):
    rn2 = case["resn"][i] ** 2; d = data["dobs"][i]; a = case["a"][i]
    w = rn2 * (t[i] * (1 - d * d) + d * d * (1 - a * a) * case["u"][i])
    gain = ln_rice_over_wilson(data["f_eff"][i], d * a * abs(case["f_calc"][i]),
      w, rn2 * t[i], centric)
    assert approx_equal(-r.target_per_reflection()[i], gain, eps=1e-6)

def exercise_fcalc_gradients(centric):
  case = make_case(16, centric, seed=4)
  t = t_consistent(case, scale=1.1)
  for rice_kappa in (0., None):
    data = normalised_data(case, t, 0. if rice_kappa is None else rice_kappa)
    use_h = rice_kappa is not None
    r = f_scale_target(case, t, data, hybrid=use_h, compute_gradients=True)
    g = r.gradients_work()
    eps = 1e-6
    for i in range(case["n"]):
      for part, delta in (("re", complex(eps, 0)), ("im", complex(0, eps))):
        vals = []
        for sign in (1, -1):
          fc = case["f_calc"].deep_copy(); fc[i] = fc[i] + sign * delta
          vals.append(f_scale_target(case, t, data, hybrid=use_h, f_calc=fc).target_work())
        fd = (vals[0] - vals[1]) / (2 * eps)
        ana = g[i].real if part == "re" else g[i].imag
        assert_gradient(ana, fd, (i, part))

def exercise_sigmaa_scatfrac_gradients(centric):
  case = make_case(16, centric, seed=5)
  n = case["n"]
  t = t_consistent(case, scale=1.15)
  sigmaa = case["a"] * 0.8
  scatfrac = flex.double([0.6 + 0.03 * i for i in range(n)])
  k = 1.05
  assert flex.max(sigmaa) < 0.999   # the hybrid's exact path needs sigmaA < 0.999
  sel = flex.bool([i % 4 != 1 for i in range(n)])
  for rice_kappa in (0., None):
    data = normalised_data(case, t, 0. if rice_kappa is None else rice_kappa)
    h = None if rice_kappa is None else data["hybrid"]
    def run(sa, sf):
      return ext.llgi_sigmaa_scatfrac_target_and_gradients(
        f_eff=data["f_eff"], selection=sel, f_calc=case["f_calc"],
        dobs=data["dobs"], sigmaa=sa, scatfrac=sf, scale_factor=k,
        teps=flex.double(n, 1.), resn=case["resn"],
        centric_flags=case["centric"], hybrid=h, tncs_u=case["u"], tncs_t=t)
    r = run(sigmaa, scatfrac)
    eps = 1e-6
    for i in range(n):
      for which in ("sigmaa", "scatfrac"):
        vals = []
        for sign in (1, -1):
          sa = sigmaa.deep_copy(); sf = scatfrac.deep_copy()
          if which == "sigmaa": sa[i] += sign * eps
          else: sf[i] += sign * eps
          vals.append(run(sa, sf).target())
        fd = (vals[0] - vals[1]) / (2 * eps)
        ana = (r.d_target_by_dsigmaa() if which == "sigmaa"
               else r.d_target_by_dscatfrac())[i]
        assert_gradient(ana, fd, (i, which))
    # E scale: e_model = |fc|/resn, eeff = Feff/(resn sqrt(t))
    e_model = flex.abs(case["f_calc"]) / case["resn"]
    e_eff = data["f_eff"] / (case["resn"] * flex.sqrt(t))
    def run_e(sa):
      return ext.llgi_e_sigmaa_target_and_gradients(e_eff=e_eff, selection=sel,
        e_model=e_model, dobs=data["dobs"], sigmaa=sa,
        centric_flags=case["centric"], hybrid=h, tncs_u=case["u"], tncs_t=t)
    re_ = run_e(sigmaa)
    for i in range(n):
      vals = []
      for sign in (1, -1):
        sa = sigmaa.deep_copy(); sa[i] += sign * eps
        vals.append(run_e(sa).target())
      fd = (vals[0] - vals[1]) / (2 * eps)
      assert_gradient(re_.d_target_by_dsigmaa()[i], fd, (i, "e-scale sigmaa"))
    # the E-scale and F-scale functors agree when k = scatfrac = 1
    rf = ext.llgi_sigmaa_scatfrac_target_and_gradients(
      f_eff=data["f_eff"], selection=sel, f_calc=case["f_calc"],
      dobs=data["dobs"], sigmaa=sigmaa, scatfrac=flex.double(n, 1.),
      scale_factor=1., teps=flex.double(n, 1.), resn=case["resn"],
      centric_flags=case["centric"], hybrid=h, tncs_u=case["u"], tncs_t=t)
    assert approx_equal(rf.target(), re_.target(), eps=1e-10)
    assert approx_equal(rf.d_target_by_dsigmaa(), re_.d_target_by_dsigmaa(), eps=1e-10)

def exercise_invalid_frame_contributes_nothing():
  """ t <= (1 - a^2) u: no contribution (like the negative-variance guard) """
  case = make_case(6, False, seed=6)
  n = case["n"]
  t = flex.double([(1 - case["a"][i] ** 2) * case["u"][i] * (0.5 if i % 2 else 1.5)
                   for i in range(n)])
  data = normalised_data(case, flex.double(n, 1.), rice_kappa=0.)
  r = f_scale_target(case, t, data, compute_gradients=True)
  for i in range(n):
    if i % 2:
      assert r.target_per_reflection()[i] == 0
      assert r.gradients_work()[i] == 0
    else:
      assert r.target_per_reflection()[i] != 0

def exercise_frame_arrays():
  """ llgi_tncs_frame (tncs_frame for arrays): a' and s, 0 where the frame
  is invalid (a outside (0, 1), u or t not positive, a'^2 not positive) """
  a = flex.double([0.5, 0.5, 0.5, 1.0, 0.0, 0.9, 0.5])
  u = flex.double([1.0, 2.0, 0.0, 1.0, 1.0, 1.5, 1.0])
  t = flex.double([1.0, 2.5, 1.0, 1.0, 1.0, 0.0, 0.5])
  r = ext.llgi_tncs_frame(a=a, u=u, t=t)
  assert list(r.valid) == [True, True, False, False, False, False, False]
  assert approx_equal(r.a_p[0], 0.5) and approx_equal(r.s[0], 1.)
  a_p = math.sqrt(1 - 0.75 * 2.0 / 2.5)
  assert approx_equal(r.a_p[1], a_p)
  assert approx_equal(r.s[1], 0.5 / (math.sqrt(2.5) * a_p))
  for i in range(2, 7):
    assert r.a_p[i] == 0 and r.s[i] == 0

def exercise():
  for centric in (False, True):
    exercise_identity_without_modulation(centric)
    exercise_exact_path_is_absolute_likelihood(centric)
    exercise_rice_path(centric)
    exercise_fcalc_gradients(centric)
    exercise_sigmaa_scatfrac_gradients(centric)
  exercise_invalid_frame_contributes_nothing()
  exercise_frame_arrays()
  print("OK")

if __name__ == "__main__":
  exercise()
