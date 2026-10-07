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

# (centric, h = I/sigI, sigI, <I>, F, SIGF): exact French-Wilson posterior
# moments F = <sqrt(J)>, SIGF = sd(sqrt(J)), computed with mpmath (30
# digits, adaptive quadrature over u = J/sigI with weight
# u^p exp(-(u - h)^2/2), p = -1/2 centric, 0 acentric; then scaled by
# sqrt(sigI)). <I> is chosen so that the prior term sigI/<I> is O(1).
french_wilson_reference = [
  (False, 15.4724887407, 2962.11391564, 19077.3821844, 213.970176194, 6.93093844489),
  (False, 16.383699352, 17.9040398534, 138.748221383, 17.1190042807, 0.52354424971),
  (False, 23.1682367932, 9.57944554571, 32.2405431498, 14.8941372725, 0.321772680024),
  (True, 23.7182696341, 601.569005297, 711.27746253, 119.369586483, 2.52231709079),
  (False, 28.8370067771, 457.865166965, 5371.30171853, 114.889021708, 1.99339251616),
  (True, -2.50497568066, 192.359572332, 244.098251375, 4.59344735172, 3.35966686543),
  (True, 4.98411945095, 6.15485423483, 39.3690085536, 5.44748136753, 0.581619222368),
  (False, 24.8000952411, 180.443959637, 2336.91529412, 66.8820289838, 1.34966058514),
  (False, 18.8608345524, 117.74965657, 358.243318068, 47.1093637548, 1.25085482032),
  (False, 29.8578241747, 1658.12265381, 28132.4725791, 222.472554497, 3.72788856501),
  (True, 4.57897475306, 36.8205139565, 48.7618932879, 12.7260389089, 1.50148956433),
  (False, 10.2131935623, 1732.6900042, 8131.36855568, 132.866549674, 6.54051787761),
  (False, 24.961222519, 5.01885692834, 11.6123810774, 11.1904632229, 0.224360057399),
  (False, 12.5095801085, 4365.62902154, 21401.4748304, 233.504697648, 9.36707772586),
  (True, 17.7720121037, 1082.6861479, 3185.30318056, 138.548123081, 3.9143272794),
  (True, 7.97532564029, 3901.20393964, 80918.8884617, 175.318100878, 11.2331806565),
  (True, 5.13080231347, 10.04867803, 12.768927726, 7.06951648873, 0.730101028086),
  (False, 2.86337823055, 238.152655901, 1425.97545652, 25.6974419272, 4.81054177886),
  (True, 21.1525091186, 12.3558109706, 162.225439742, 16.1529161629, 0.382950041795),
  (True, 10.8849353691, 21.7553975263, 64.01015852, 15.3390273724, 0.712669073932),
  (False, 23.5125796017, 40.8700768134, 1407.90579572, 30.9923488595, 0.659732201545),
  (True, 10.0110630234, 1828.52418721, 23827.7725155, 134.781905727, 6.82337952499),
  (True, 29.6469560178, 21.812232089, 61.2878212639, 25.4187619849, 0.429333870214),
  (False, 7.85552904308, 38.7207973754, 51.9338643209, 17.4046481138, 1.11828994),
  (True, 16.2302483395, 26.7922239521, 296.855969116, 20.823030919, 0.64473190882),
  (True, 11.9558674559, 3770.28362445, 26102.9507536, 211.749207835, 8.93902713882),
  (False, 25.5953470375, 17.6788150207, 32.750175045, 21.2678747886, 0.415821648935),
  (False, 23.9874643322, 28.0198433932, 59.8664121403, 25.9197322888, 0.540806862716),
  (False, 28.0333615711, 19.4419306899, 869.54972366, 23.3420028463, 0.416624282876),
  (False, 16.9166290879, 91.9055368348, 139.228991917, 39.4128084823, 1.16722118375),
]

def exercise_french_wilson_inverse():
  """ (I, sigI, <I>) -> exact French-Wilson (F, SIGF) -> back. """
  f, sigf, mean_i, cen, i_true, s_true = [], [], [], [], [], []
  for c, h, s, mi, f_fw, sigf_fw in french_wilson_reference:
    f.append(f_fw)
    sigf.append(sigf_fw)
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
      french_wilson_reference[i][1])), i
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
