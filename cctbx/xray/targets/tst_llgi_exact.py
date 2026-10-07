from __future__ import absolute_import, division, print_function
from cctbx.xray import ext
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import math

# Reference values computed independently with mpmath (30 digits, adaptive
# quadrature over J = |E|^2 with dense breakpoints) for
# (E_obs^2, sigma(E_obs^2), E_calc, sigmaA, centric) -> exact LLGI.
# The first eight are the hybrid-LLGI handoff fixtures (sec. 13), whose
# own LLGI_exact values agree to their quoted 1e-8 except
# acentric_moderate sigmaA=0.9 (1.36870903 there, from a 32-node fixed
# Gauss-Laguerre rule).
reference = [
  ( 5.0,  2.0, 1.5, 0.6, False,  0.3799826485),
  ( 5.0,  2.0, 2.0, 0.9, False,  1.36870907),
  (46.5, 30.0, 1.5, 0.6, False,  0.0229181622),
  (46.5, 30.0, 2.0, 0.9, False,  0.1196025348),
  ( 5.0,  2.0, 1.5, 0.6, True,   0.339885682),
  ( 5.0,  2.0, 1.0, 0.3, True,   0.001630498345),
  ( 1.0,  1.0, 1.5, 0.6, True,  -0.07630331786),
  ( 1.0,  1.0, 1.0, 0.3, True,   8.9555402e-5),
  # Strongly Bessel-tilted cases.
  (400.0, 3.0, 4.0, 0.98, False, -4405.31226933),
  (400.0, 3.0, 8.0, 0.98, True,  -1417.13144187),
  (150.0, 10.0, 4.0, 0.95, False,   7.88554544938),
]

def evaluate(cases):
  return ext.llgi_exact_evaluate(
    e_obs_sq=flex.double([c[0] for c in cases]),
    sig_e_obs_sq=flex.double([c[1] for c in cases]),
    e_calc=flex.double([c[2] for c in cases]),
    sigmaa=flex.double([c[3] for c in cases]),
    centric_flags=flex.bool([c[4] for c in cases]))

def exercise_reference_values():
  r = evaluate(reference)
  for c, ll in zip(reference, r.ll):
    assert approx_equal(ll, c[5], eps=2e-7 * max(1, abs(c[5]))), (c, ll)

def exercise_sigmaa_zero_limit():
  # LLGI vanishes identically as sigmaA -> 0 (no information from model)
  cases = [(e, s, ec, 1e-9, c) for e, s in [(5, 2), (0.5, 3), (-1, 1)]
           for ec in [0.5, 2.5] for c in [False, True]]
  for ll in evaluate(cases).ll:
    assert abs(ll) < 1e-12, ll

def exercise_derivatives():
  cases = []
  for e, s in [(5, 2), (46.5, 30), (1, 1), (20, 0.5), (-2, 1), (42.2, 5.3)]:
    for ec in [0.3, 1.5, 4.0]:
      for a in [0.3, 0.75, 0.9]:
        for c in [False, True]:
          cases.append((e, s, ec, a, c))
  r = evaluate(cases)
  h = 1e-5
  def shifted(dec=0, da=0):
    return evaluate([(e, s, ec + dec, a + da, c) for e, s, ec, a, c in cases])
  fd_ec = (shifted(dec=h).ll - shifted(dec=-h).ll) / (2 * h)
  fd_a = (shifted(da=h).ll - shifted(da=-h).ll) / (2 * h)
  for i in range(len(cases)):
    scale = max(1, abs(fd_a[i]))
    assert abs(r.d_ll_d_ec[i] - fd_ec[i]) < 1e-6 * max(1, abs(fd_ec[i]))
    assert abs(r.d_ll_d_a[i] - fd_a[i]) < 1e-6 * scale, (cases[i],)
    # posterior mean E identity: dLLGI/dEc = (2 b a / S)(<E> - a Ec)
    e, s, ec, a, c = cases[i]
    b = 0.5 if c else 1.0
    assert approx_equal(r.d_ll_d_ec[i],
      2 * b * a / (1 - a * a) * (r.e_expected[i] - a * ec), eps=1e-9)

def exercise_rice_moments():
  # phasertng math::rice_from_intensity values for the handoff fixtures
  # (agree with 50-digit references except where noted there).
  e = flex.double([5, 46.5, 42.209152, 5, 1, 5])
  s = flex.double([2, 30, 5.3, 2, 1, 0.3])
  c = flex.bool([False, False, False, True, True, False])
  r = ext.llgi_rice_moments(e_obs_sq=e, sig_e_obs_sq=s, centric_flags=c)
  expected_mu4 = [6.01832086767407, 2.210249120908386, 228.302664458899,
                  9.47658702987293, 0.822616135807296]
  expected_d2 = [0.440760138337296, 0.0001089276402577, 0.0051725709101173,
                 0.647492483117873, 0.816480797938071]
  for i in range(5):
    assert r.valid[i]
    assert approx_equal(r.mu4[i], expected_mu4[i], eps=1e-7*expected_mu4[i])
    # D^2 is ill-conditioned near the validity boundary (case 2)
    assert approx_equal(r.dsqr[i], expected_d2[i], eps=2e-4*expected_d2[i])
  # Strong, precise reflection: D^2 close to 1, Eeff close to sqrt(E^2)
  assert r.valid[5]
  assert r.dsqr[5] > 0.99 and abs(r.eeff[5] - math.sqrt(5)) < 0.02
  # A reflection with no Rice solution (very large E^2 relative to sigma
  # in the low-information regime)
  r = ext.llgi_rice_moments(
    e_obs_sq=flex.double([75.6]), sig_e_obs_sq=flex.double([44.2]),
    centric_flags=flex.bool([False]))
  assert not r.valid[0]

def hybrid_inputs(n=40, seed=3):
  import random
  rnd = random.Random(seed)
  d = dict(
    e_obs_sq=flex.double([rnd.expovariate(1) + rnd.gauss(0, 0.5)
                          for i in range(n)]),
    sig_e_obs_sq=flex.double([0.1 + 2 * rnd.random() for i in range(n)]),
    e_eff=flex.double([0.2 + 2 * rnd.random() for i in range(n)]),
    e_model=flex.double([0.1 + 2 * rnd.random() for i in range(n)]),
    dobs=flex.double([0.5 + 0.5 * rnd.random() for i in range(n)]),
    sigmaa=flex.double([0.4 + 0.5 * rnd.random() for i in range(n)]),
    centric_flags=flex.bool([rnd.random() < 0.3 for i in range(n)]),
    selection=flex.bool(n, True))
  return d

def exercise_hybrid_e_scale():
  d = hybrid_inputs()
  n = d["e_eff"].size()
  args = dict((k, d[k]) for k in
    ["e_eff", "selection", "e_model", "dobs", "sigmaa", "centric_flags"])
  # rice_kappa: 0 = exact everywhere, 1e9 = never (Rice everywhere)
  def hybrid(threshold, force=None):
    return ext.llgi_hybrid(
      e_obs_sq=d["e_obs_sq"], sig_e_obs_sq=d["sig_e_obs_sq"],
      force_exact=force if force is not None else flex.bool(n, False),
      centric_flags=d["centric_flags"], rice_kappa=threshold)
  # A hybrid that never triggers reproduces the Rice target exactly
  cls = ext.llgi_e_sigmaa_target_and_gradients
  r0 = cls(**args)
  r1 = cls(hybrid=hybrid(1e9), **args)
  assert r1.n_exact() == 0
  assert r0.target() == r1.target()
  # All-exact: target is the mean exact LLGI, gradients match FD
  h = hybrid(0.0)
  rs = ext.llgi_e_sigmaa_target_and_gradients(hybrid=h, **args)
  assert rs.n_exact() == n
  ex = ext.llgi_exact_evaluate(d["e_obs_sq"], d["sig_e_obs_sq"],
    d["e_model"], d["sigmaa"], d["centric_flags"])
  assert approx_equal(rs.target(), -flex.mean(ex.ll), eps=1e-12)
  delta = 1e-6
  for i in [0, 7, 19]:
    r = cls(hybrid=h, **args)
    t = []
    for s in [delta, -delta]:
      a2 = dict(args)
      a2["sigmaa"] = args["sigmaa"].deep_copy()
      a2["sigmaa"][i] += s
      t.append(cls(hybrid=h, **a2).target())
    fd = (t[0] - t[1]) / (2 * delta)
    assert approx_equal(r.d_target_by_dsigmaa()[i], fd, eps=1e-6), i
  # Mixed: force_exact on some reflections only
  force = flex.bool([i % 5 == 0 for i in range(n)])
  rm = ext.llgi_e_sigmaa_target_and_gradients(
    hybrid=hybrid(1e9, force), **args)
  assert rm.n_exact() == force.count(True)

def exercise_hybrid_f_scale():
  d = hybrid_inputs()
  n = d["e_eff"].size()
  import cmath, random
  rnd = random.Random(5)
  resn = flex.double([1.5 + rnd.random() for i in range(n)])
  phases = [2 * math.pi * rnd.random() for i in range(n)]
  f_calc = flex.complex_double([d["e_model"][i] * resn[i] * 0.9
    * cmath.exp(1j * phases[i]) for i in range(n)])
  scatfrac = flex.double([0.8 + 0.3 * rnd.random() for i in range(n)])
  args = dict(f_eff=d["e_eff"] * resn, f_calc=f_calc, dobs=d["dobs"],
    sigmaa=d["sigmaa"] * 0.9, scatfrac=scatfrac, scale_factor=1.1,
    teps=flex.double(n, 1), resn=resn, centric_flags=d["centric_flags"])
  h = ext.llgi_hybrid(e_obs_sq=d["e_obs_sq"],
    sig_e_obs_sq=d["sig_e_obs_sq"], force_exact=flex.bool(n, False),
    centric_flags=d["centric_flags"], rice_kappa=0.0)
  rff = flex.bool(n, False)
  r0 = ext.llgi_target_and_gradients(r_free_flags=rff,
    compute_gradients=True, **args)
  r1 = ext.llgi_target_and_gradients(r_free_flags=rff,
    compute_gradients=True, hybrid=ext.llgi_hybrid(e_obs_sq=d["e_obs_sq"],
      sig_e_obs_sq=d["sig_e_obs_sq"], force_exact=flex.bool(n, False),
      centric_flags=d["centric_flags"], rice_kappa=1e9), **args)
  assert r1.n_exact() == 0 and r0.target_work() == r1.target_work()
  r = ext.llgi_target_and_gradients(r_free_flags=rff,
    compute_gradients=True, hybrid=h, **args)
  assert r.n_exact() == n
  # gradient w.r.t. f_calc (real and imaginary parts), work-set mean
  delta = 1e-6
  g = r.gradients_work()
  for i in [0, 11, 23]:
    for part in [1, 1j]:
      t = []
      for s in [delta, -delta]:
        fc = f_calc.deep_copy()
        fc[i] += s * part
        a2 = dict(args); a2["f_calc"] = fc
        t.append(ext.llgi_target_and_gradients(r_free_flags=rff,
          compute_gradients=False, hybrid=h, **a2).target_work())
      fd = (t[0] - t[1]) / (2 * delta)
      analytic = g[i].real if part == 1 else g[i].imag
      assert approx_equal(analytic, fd, eps=1e-6), (i, part, analytic, fd)
  # sigmaa / scatfrac gradients
  sel = flex.bool(n, True)
  rs = ext.llgi_sigmaa_scatfrac_target_and_gradients(
    selection=sel, hybrid=h, **args)
  assert rs.n_exact() == n
  for name, getter in [("sigmaa", "d_target_by_dsigmaa"),
                       ("scatfrac", "d_target_by_dscatfrac")]:
    for i in [2, 17]:
      t = []
      for s in [delta, -delta]:
        a2 = dict(args); a2[name] = args[name].deep_copy(); a2[name][i] += s
        t.append(ext.llgi_sigmaa_scatfrac_target_and_gradients(
          selection=sel, hybrid=h, **a2).target())
      fd = (t[0] - t[1]) / (2 * delta)
      assert approx_equal(getattr(rs, getter)()[i], fd, eps=1e-6), (name, i)

def exercise_french_wilson_inverse():
  """ (I, sigI, <I>) -> exact French-Wilson (F, SIGF) -> back. """
  import random
  rnd = random.Random(11)
  cases = []
  for k in range(30):
    c = rnd.random() < 0.4
    h = rnd.uniform(-3, 30)
    sig_i = math.exp(rnd.uniform(math.log(5), math.log(5000)))
    mean_i = sig_i * math.exp(rnd.uniform(0, 4))
    cases.append((c, h, sig_i, mean_i))
  # forward French-Wilson moments (exact posterior), by mpmath quadrature
  import mpmath as mp
  f, sigf, mean_i, cen, i_true, s_true = [], [], [], [], [], []
  for c, h, s, mi in cases:
    p = -0.5 if c else 0
    w = lambda u: (u**p) * mp.e**(-(u - h)**2 / 2)
    pts = [0, max(h, 0) + 0.5, max(h, 0) + 12, mp.inf]
    z = mp.quad(w, pts)
    m1 = mp.quad(lambda u: w(u) * mp.sqrt(u), pts) / z
    m2 = mp.quad(lambda u: w(u) * u, pts) / z
    f.append(float(m1) * math.sqrt(s))
    sigf.append(float(mp.sqrt(m2 - m1 * m1)) * math.sqrt(s))
    mean_i.append(mi)
    cen.append(c)
    i_true.append(s * (h + (0.5 if c else 1.0) * s / mi))
    s_true.append(s)
  inv = ext.llgi_french_wilson_inverse(f=flex.double(f),
    sigf=flex.double(sigf), mean_intensity=flex.double(mean_i),
    centric_flags=flex.bool(cen))
  for i in range(len(f)):
    assert inv.valid[i] and not inv.prior_dominated[i]
    assert abs(inv.sig_i_obs[i] / s_true[i] - 1) < 1e-5, i
    assert abs(inv.i_obs[i] - i_true[i]) < 1e-5 * s_true[i] * max(1, abs(
      cases[i][1])), i
  # SIGF/F at or beyond the large-sigma limit: prior-dominated
  inv = ext.llgi_french_wilson_inverse(f=flex.double([1.0, 1.0]),
    sigf=flex.double([0.6, 0.8]), mean_intensity=flex.double([1.0, 1.0]),
    centric_flags=flex.bool([False, True]))
  assert inv.prior_dominated[0] and inv.prior_dominated[1]

def exercise_rice_moments_large_sigma():
  """ Near the edge of the Rice family (D -> 0, large Eeff: weak
  information, typically large sigma) the solution is very sensitive to
  <E^2>, and <E^4> = m <E^2> + k sig^2 is a difference of terms of order
  sig^2. References from 60-digit integration. """
  cases = [
    # e_obs_sq, sig, centric, mu2, D, Eeff
    (18.33,  9.73,   False, 1.2027432978797986, 0.01484041782, 30.35728178),
    (41.31,  26.34,  False, 1.0598740303535263, 0.01403904474, 17.4580336),
    (796.17, 563.10, False, 1.0025108963288873, 0.003631826334, 13.83333973),
    (173.26, 138.73, True,  1.0180059498475505, 0.01683120269, 8.034940596),
    (75.67,  57.19,  True,  1.0464130606123748, 0.01859118521, 11.63118874)]
  r = ext.llgi_rice_moments(
    e_obs_sq=flex.double([c[0] for c in cases]),
    sig_e_obs_sq=flex.double([c[1] for c in cases]),
    centric_flags=flex.bool([c[2] for c in cases]))
  for i, c in enumerate(cases):
    assert abs(r.mu2[i] / c[3] - 1) < 1e-10, (i, r.mu2[i], c[3])
    assert r.valid[i], i
    assert abs(math.sqrt(r.dsqr[i]) / c[4] - 1) < 0.01, (i, r.dsqr[i])
    assert abs(r.eeff[i] / c[5] - 1) < 0.01, (i, r.eeff[i])

def exercise_small_sigma_fixtures():
  """ Hybrid LLGI handoff, revision 2, sec. 13: small-sigma fixtures. The
  exact LLGI and gradient, D^2, and the hybrid rule: exact where
  1 - D^2 > rice_kappa (1 - sigmaA^2)^2 (rice_kappa = 0.1 for both acentric
  and centric). With sigmaA = 0.9 the limit is 1 - D^2 > 0.0036: the two
  strong reflections stay on Rice, the two weak ones (Rice errors 0.24 and
  0.044 at ec = 3) are exact. """
  fixtures = [
    # e_obs_sq, sig, centric, D2, sigmaA, ec, llgi_exact, grad_exact, exact
    (1.0,      0.05, False, 0.9987460790, 0.9, 2.0, -3.07807616, -7.78898929, False),
    (0.063337, 0.05, False, 0.9866387871, 0.9, 3.0, -30.03548112, -22.21102762, True),
    (1.0,      0.05, True,  0.9993724380, 0.9, 2.0, -1.04799042, -3.78147179, False),
    (0.105036, 0.05, True,  0.9924916055, 0.9, 3.0, -14.60839593, -11.05683607, True)]
  e2 = flex.double([f[0] for f in fixtures])
  sig = flex.double([f[1] for f in fixtures])
  cen = flex.bool([f[2] for f in fixtures])
  sa = flex.double([f[4] for f in fixtures])
  r = ext.llgi_exact_evaluate(e_obs_sq=e2, sig_e_obs_sq=sig,
    e_calc=flex.double([f[5] for f in fixtures]), sigmaa=sa, centric_flags=cen)
  for i, f in enumerate(fixtures):
    assert approx_equal(r.ll[i], f[6], eps=1e-6), (i, r.ll[i], f[6])
    assert approx_equal(r.d_ll_d_ec[i], f[7], eps=1e-6), (i, r.d_ll_d_ec[i])
  def hybrid(rice_kappa):
    return ext.llgi_hybrid(e_obs_sq=e2, sig_e_obs_sq=sig,
      force_exact=flex.bool(len(fixtures), False), centric_flags=cen,
      rice_kappa=rice_kappa)
  h = hybrid(0.1)
  assert approx_equal(h.dsqr, [f[3] for f in fixtures], eps=1e-8)
  assert list(h.exact_selection(sa)) == [f[8] for f in fixtures]
  # lower sigmaA: more model error, so the weak reflections go back to Rice
  assert h.exact_selection(flex.double(4, 0.5)).count(True) == 0
  # rice_kappa = 0 (sigmaA fits): exact wherever possible
  assert hybrid(0).exact_selection(sa).count(True) == len(fixtures)
  assert h.with_rice_kappa(0).exact_selection(sa).count(True) == len(fixtures)
  sel = h.select(flex.bool([False, True, True, False]))
  assert list(sel.exact_selection(sa.select(flex.size_t([1, 2])))) \
    == [True, False]

def run():
  exercise_reference_values()
  exercise_sigmaa_zero_limit()
  exercise_derivatives()
  exercise_rice_moments()
  exercise_hybrid_e_scale()
  exercise_hybrid_f_scale()
  exercise_french_wilson_inverse()
  exercise_rice_moments_large_sigma()
  exercise_small_sigma_fixtures()
  print("OK")

if (__name__ == "__main__"):
  run()
