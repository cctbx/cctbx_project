from __future__ import absolute_import, division, print_function
""" mmtbx.refinement.llgi_normalisation: the BEST curve, resolution bins,
the SigmaN model with its restraints, and the joint data-only fit of SigmaN
and the tNCS coherence. """
import math
import time
import numpy
from cctbx import adptbx, crystal, miller
from cctbx.array_family import flex
from cctbx.xray import ext as xray_ext
from libtbx.test_utils import approx_equal
from six.moves import cStringIO as StringIO
import mmtbx.refinement.llgi_normalisation as ln
from mmtbx.refinement.llgi_tncs import U_MIN

def make_set(symbol, cell, d_min):
  cs = crystal.symmetry(unit_cell=cell, space_group_symbol=symbol)
  return miller.build_set(cs, anomalous_flag=False, d_min=d_min)

def exercise_best_curve_matches_phasertng():
  assert len(ln.BEST_TABLE) == 300
  k = 8.563714792513845e-06
  assert ln.BEST_K == k and ln.BEST_SSQR0 == 0.009
  assert ln.BEST_STEP == 0.004091973
  t0, t1, t299 = 119420.094, 101659.641, 6939.107
  assert (ln.BEST_TABLE[0], ln.BEST_TABLE[1], ln.BEST_TABLE[299]) == (t0, t1, t299)
  ssqr = flex.double([0.009, 0.009 + 0.004091973, 0.009 + 299 * 0.004091973,
    0.009 + 0.5 * 0.004091973, 0.0, 0.009 + 400 * 0.004091973])
  v = ln.best_curve_values(ssqr)
  expected = [k * t0, k * t1, k * t299, k * (t0 + t1) / 2, k * t0, k * t299]
  for a, b in zip(v, expected):
    assert abs(a - b) <= 1e-12 * abs(b), (a, b)
  # linear interpolation between table entries, as best_curve_at_ssqr::get
  s = 0.009 + 37.3 * 0.004091973
  e = k * (ln.BEST_TABLE[37] + 0.3 * (ln.BEST_TABLE[38] - ln.BEST_TABLE[37]))
  assert abs(ln.best_curve_values(flex.double([s]))[0] - e) < 1e-12 * e

def exercise_bin_restraint_sigmas():
  s = ln.bin_restraint_sigmas(flex.double([10, 6, 3, 2]))
  assert approx_equal(list(s), [1.15, 0.366, 0.177, 0.158], eps=1e-12)

def exercise_bins():
  ms = make_set("P 1", (50, 60, 70, 90, 90, 90), 2.0)
  dss = ms.d_star_sq().data()
  n = dss.size()
  bins = ln.resolution_bins(dss)
  assert bins.n_bins == max(6, min(n // 500, 50)) == 50, bins.n_bins
  assert bins.d_mid.size() == bins.n_bins
  idx = numpy.array(list(bins.index))
  counts = numpy.bincount(idx, minlength=bins.n_bins)
  assert counts.min() >= 10
  order = numpy.argsort(numpy.array(list(dss)))
  assert (numpy.diff(idx[order]) >= 0).all()
  # d_mid decreases with bin number, and lies within the data range
  dm = list(bins.d_mid)
  assert all(a > b for a, b in zip(dm[:-1], dm[1:]))
  assert dm[-1] > 2.0 and dm[0] < ms.d_max_min()[0]
  # equal width in d*^2: the mid d*^2 steps are constant
  mid = 1 / numpy.array(dm) ** 2
  assert approx_equal(numpy.diff(mid), numpy.diff(mid)[0] * numpy.ones(len(mid) - 1),
    eps=1e-5)
  # few reflections: bins are merged until each has at least 10
  sub = dss.select(flex.size_t(range(0, 40)))
  b2 = ln.resolution_bins(sub, n_min=6)
  c2 = numpy.bincount(numpy.array(list(b2.index)), minlength=b2.n_bins)
  assert b2.n_bins < 6 and c2.min() >= 10

def random_x(model, seed=1):
  rs = numpy.random.RandomState(seed)
  return rs.normal(size=model.n_parameters()) * 0.05

def exercise_beta_constraints():
  def model_for(symbol, cell):
    ms = make_set(symbol, cell, 3.0)
    return ms, ln.sigman_model(ms, ln.resolution_bins(ms.d_star_sq().data()))
  ms, m = model_for("P 43 21 2", (60, 60, 80, 90, 90, 90))
  rs = numpy.random.RandomState(0)
  x = m.start(flex.double(ms.size(), 1.0))
  x[1 + m.n_bins:] = rs.normal(size=m.n_parameters() - 1 - m.n_bins)
  b = m.beta(x)
  assert abs(b[0] - b[1]) < 1e-14 and abs(b[0]) > 1e-3
  assert b[3] == b[4] == b[5] == 0
  assert m.n_parameters() == 1 + m.n_bins + 2
  ms, m = model_for("C 2", (60, 40, 50, 90, 100, 90))
  x = numpy.zeros(m.n_parameters())
  x[1 + m.n_bins:] = rs.normal(size=m.n_parameters() - 1 - m.n_bins)
  b = m.beta(x)
  assert b[3] == 0 and b[5] == 0 and abs(b[4]) > 1e-3
  assert m.n_parameters() == 1 + m.n_bins + 4
  ms, m = model_for("P 61 2 2", (50, 50, 90, 90, 90, 120))
  assert m.n_parameters() == 1 + m.n_bins + 2
  x = numpy.zeros(m.n_parameters())
  x[1 + m.n_bins:] = rs.normal(size=2)
  b = m.beta(x)
  # phasertng's convention: b12 multiplies hk, so b11 = b22 = b12
  assert abs(b[0] - b[1]) < 1e-14 and abs(b[3] - b[0]) < 1e-14
  assert b[4] == b[5] == 0
  # with resolution dependence the independent parameters are doubled
  m2 = ln.sigman_model(ms, ln.resolution_bins(ms.d_star_sq().data()),
    resolution_dependent=True)
  assert m2.n_parameters() == 1 + m.n_bins + 4

def exercise_derivatives_fd():
  ms = make_set("P 1", (50, 60, 70, 80, 95, 100), 3.5)
  bins = ln.resolution_bins(ms.d_star_sq().data())
  for symbol, cell in (("P 1", (50, 60, 70, 80, 95, 100)),
                       ("P 61 2 2", (50, 50, 90, 90, 90, 120))):
    ms = make_set(symbol, cell, 3.5)
    bins = ln.resolution_bins(ms.d_star_sq().data())
    m = ln.sigman_model(ms, bins, resolution_dependent=True)
    n = m.n_parameters()
    x = random_x(m)
    x[:1 + bins.n_bins] *= 5
    g = m.d_ln_sigman_dx(x)
    assert g.shape == (n, ms.size())
    f0 = numpy.array(list(m.ln_sigman(x)))
    h = 1e-4
    for i in range(n):
      xp = x.copy(); xm = x.copy()
      xp[i] += h; xm[i] -= h
      fd = (numpy.array(list(m.ln_sigman(xp))) -
            numpy.array(list(m.ln_sigman(xm)))) / (2 * h)
      assert abs(fd - g[i]).max() < 1e-6 * max(1, abs(g[i]).max()), (symbol, i)
    # the model is linear in x apart from the BEST constant
    lnbest = numpy.log(numpy.array(list(
      ln.best_curve_values(ms.d_star_sq().data()))))
    assert abs(f0 - (g.T.dot(x) + lnbest)).max() < 1e-9
    v, gr, hs = m.restraint(x)
    assert gr.shape == (n,) and hs.shape == (n, n)
    for i in range(n):
      xp = x.copy(); xm = x.copy()
      xp[i] += h; xm[i] -= h
      fd = (m.restraint(xp)[0] - m.restraint(xm)[0]) / (2 * h)
      assert abs(fd - gr[i]) < 1e-6 * max(1, abs(gr[i])), (symbol, i, fd, gr[i])
      fdh = (m.restraint(xp)[1] - m.restraint(xm)[1]) / (2 * h)
      assert abs(fdh - hs[i]).max() < 1e-6 * max(1, abs(hs[i]).max()), (symbol, i)
    assert abs(hs - hs.T).max() < 1e-12 * abs(hs).max()
    assert numpy.linalg.eigvalsh(hs).min() > -1e-9 * abs(hs).max()

def exercise_restraint_values():
  ms = make_set("P 21 21 21", (50, 60, 70, 90, 90, 90), 2.5)
  bins = ln.resolution_bins(ms.d_star_sq().data())
  m = ln.sigman_model(ms, bins)
  x = numpy.zeros(m.n_parameters())
  assert m.restraint(x)[0] == 0
  # bin term only: 0.5 (ln b / sigma)^2
  x[1 + 2] = 0.3
  sig = list(ln.bin_restraint_sigmas(bins.d_mid))
  assert approx_equal(m.restraint(x)[0], 0.5 * (0.3 / sig[2]) ** 2, eps=1e-12)
  # sphericity: sigma = 10 (d_min/2.5)^2 on a* a* etc.; beta11 only, with
  # the isotropic part removed (orthogonal cell: coefficient (4/3) A^2)
  x = numpy.zeros(m.n_parameters())
  x[1 + bins.n_bins] = 0.01
  a, b, c = [1 / 50., 1 / 60., 1 / 70.]
  A, B, C = 50., 60., 70.
  sigma = 10 * (ms.d_min() / 2.5) ** 2
  sb = [sigma * a * a / 4, sigma * b * b / 4, sigma * c * c / 4]
  iso = (4 / 3.) * 0.01 * A * A
  p = [0.01 - iso * a * a / 4, -iso * b * b / 4, -iso * c * c / 4]
  expect = 0.5 * sum((p[i] / sb[i]) ** 2 for i in range(3))
  assert approx_equal(m.restraint(x)[0], expect, eps=1e-10), (m.restraint(x)[0], expect)

def exercise_isotropic_beta():
  for symbol, cell, constr in (
      ("P 1", (50, 60, 70, 80, 95, 100), False),
      ("P 61 2 2", (50, 50, 90, 90, 90, 120), True),
      ("P 21 21 21", (50, 60, 70, 90, 90, 90), True)):
    ms = make_set(symbol, cell, 3.0)
    bins = ln.resolution_bins(ms.d_star_sq().data())
    m = ln.sigman_model(ms, bins, resolution_dependent=True)
    uc = ms.unit_cell()
    rp = uc.reciprocal_parameters()
    cg, cb, ca = [math.cos(math.radians(v)) for v in (rp[5], rp[4], rp[3])]
    B = 17.0
    iso = [B / 4 * rp[0] ** 2, B / 4 * rp[1] ** 2, B / 4 * rp[2] ** 2,
           B / 2 * rp[0] * rp[1] * cg, B / 2 * rp[0] * rp[2] * cb,
           B / 2 * rp[1] * rp[2] * ca]
    # solve for the independent parameters reproducing the isotropic beta
    nb = bins.n_bins
    x = numpy.zeros(m.n_parameters())
    x[1 + nb:1 + nb + m.n_beta_ind()] = m.independent_from_beta(iso)
    assert approx_equal(list(m.beta(x)), iso, eps=1e-12)
    assert abs(numpy.array(list(m.anisotropy_term(x))) - 1).max() < 1e-12
    # beta1 isotropic as well
    x[1 + nb + m.n_beta_ind():] = m.independent_from_beta(iso)
    assert abs(numpy.array(list(m.anisotropy_term(x))) - 1).max() < 1e-12
    # but a real anisotropy is not removed
    x[1 + nb] += 0.5 * iso[0]
    assert abs(numpy.array(list(m.anisotropy_term(x))) - 1).max() > 1e-4

def exercise_start():
  ms = make_set("P 21 21 21", (50, 60, 70, 90, 90, 90), 2.5)
  bins = ln.resolution_bins(ms.d_star_sq().data())
  m = ln.sigman_model(ms, bins)
  # data exactly following BEST (times eps) with a scale: start recovers it
  best = ln.best_curve_values(ms.d_star_sq().data())
  eps = ms.epsilons().data().as_double()
  i_obs = 7.0 * best * eps
  x = m.start(i_obs)
  assert x.shape == (m.n_parameters(),)
  fit = numpy.array(list(m.ln_sigman(x)))
  assert abs(fit - numpy.log(7.0 * numpy.array(list(best)))).max() < 1e-6
  assert (x[1 + bins.n_bins:] == 0).all()

# ---------------------------------------------------------------------------
# The joint fit (fit, normalisation_params)

P212121 = ("P 21 21 21", (60, 70, 80, 90, 90, 90))
B_ANISO = (-5., 0., 5., 0., 0., 0.)   # B_cart, max - min eigenvalue 10 A^2
TNCS_T = [(0, 0, 0.47), (0.31, 0, 0)]
SCALE_TRUE = 50.

def _true_sigman(ms, beta_aniso_b):
  """ SigmaN_true = scale BEST exp(-h.beta.h) wiggle(bin), numpy """
  dss = ms.d_star_sq().data()
  best = ln.best_curve_values(dss).as_numpy_array()
  u_star = adptbx.u_cart_as_u_star(ms.unit_cell(),
    adptbx.b_as_u(tuple(beta_aniso_b)))
  hkl = ms.indices().as_vec3_double().as_numpy_array()
  h, k, l = hkl[:, 0], hkl[:, 1], hkl[:, 2]
  u = u_star
  q = (u[0] * h * h + u[1] * k * k + u[2] * l * l + 2 * u[3] * h * k +
    2 * u[4] * h * l + 2 * u[5] * k * l)
  aniso = numpy.exp(-2 * math.pi ** 2 * q)
  bins = ln.resolution_bins(dss)
  b = bins.index.as_numpy_array()
  wiggle = 1 + 0.05 * numpy.cos(2.3 * b)
  return SCALE_TRUE * best * aniso * wiggle

def _a_group(ms, t, amplitude=1.):
  """ consensus pair term of two copies related by the translation t """
  hkl = ms.indices().as_vec3_double().as_numpy_array()
  dss = ms.d_star_sq().data().as_numpy_array()
  return amplitude * numpy.cos(2 * math.pi * hkl.dot(numpy.array(t))) * \
    numpy.exp(-10 * dss / 4)

def _simulate(space_group, cell, d_min, beta_aniso_b, gamma, v,
    sigma_fraction, seed, n_groups=1):
  """ intensities following Wilson statistics with SigmaN_true (BEST, an
  anisotropic B, a fixed +-5% bin pattern), tNCS modulation
  max(1 + gamma exp(-2 pi^2 v^2 s^2) sum A_g, U_MIN) and Gaussian errors of
  sigma_fraction Sigma (a number, or one per reflection). Returns (i_obs,
  a_groups) """
  rs = numpy.random.RandomState(seed)
  ms = make_set(space_group, cell, d_min)
  dss = ms.d_star_sq().data().as_numpy_array()
  eps = ms.epsilons().data().as_numpy_array().astype(float)
  cen = ms.centric_flags().data().as_numpy_array()
  amplitude = 2 * 0.5 if n_groups == 1 else 0.5
  a = [_a_group(ms, TNCS_T[g], amplitude) for g in range(n_groups)]
  t = 1 + gamma * numpy.exp(-2 * math.pi ** 2 * v ** 2 * dss) * sum(a)
  t = numpy.maximum(t, U_MIN)
  sigma = eps * _true_sigman(ms, beta_aniso_b) * t
  n = ms.size()
  j = numpy.where(cen, sigma * rs.normal(size=n) ** 2,
    sigma * rs.exponential(size=n))
  sig = numpy.asarray(sigma_fraction) * sigma
  i = j + rs.normal(size=n) * sig
  i_obs = ms.array(data=flex.double(i), sigmas=flex.double(sig))
  i_obs.set_observation_type_xray_intensity()
  return i_obs, [flex.double(x) for x in a]

def _resn_true(i_obs):
  eps = i_obs.epsilons().data().as_numpy_array()
  return numpy.sqrt(eps * _true_sigman(i_obs, B_ANISO))

def _check_sigman(r, i_obs, tol_db=1., tol_resn=0.02):
  assert abs(r.delta_b_aniso - 10.) < tol_db, r.delta_b_aniso
  ratio = r.resn.data().as_numpy_array() / _resn_true(i_obs)
  assert numpy.median(abs(ratio - 1)) < tol_resn, numpy.median(abs(ratio - 1))
  assert math.isfinite(r.log_likelihood)

def exercise_fit_recovers_truth():
  i_obs, a = _simulate(P212121[0], P212121[1], 2.0, B_ANISO, 0.9, 0.4,
    0.1, seed=1)
  log = StringIO()
  r = ln.fit(i_obs, tncs_a_groups=a, log=log)
  assert abs(r.coherence.gamma - 0.9) < 0.02, r.coherence.gamma
  assert abs(r.coherence.v - 0.4) < 0.1, r.coherence.v
  _check_sigman(r, i_obs)
  assert r.gamma_groups is None
  assert r.warnings == [], r.warnings
  assert r.resn.indices().all_eq(i_obs.indices())
  assert r.teps.size() == r.anisotropy.size() == i_obs.size()
  # teps follows the fitted coherence
  t = 1 + r.coherence.at(i_obs.d_star_sq().data()) * a[0]
  t.set_selected(t < U_MIN, U_MIN)
  assert approx_equal(r.teps.data(), t, eps=1e-12)
  assert r.n_teps_floored == (1 + r.coherence.at(
    i_obs.d_star_sq().data()) * a[0] < U_MIN).count(True)
  assert r.selection.count(True) == i_obs.size() - r.n_outliers
  assert len(log.getvalue()) > 0

def exercise_no_tncs():
  i_obs, a = _simulate(P212121[0], P212121[1], 2.0, B_ANISO, 0., 0.,
    0.1, seed=2)
  r = ln.fit(i_obs, tncs_a_groups=None)
  assert r.coherence is None and r.gamma_groups is None
  assert (r.teps.data() == 1).all_eq(True)
  assert r.n_teps_floored == 0
  _check_sigman(r, i_obs)

def exercise_gamma_bound():
  i_obs, a = _simulate(P212121[0], P212121[1], 2.0, B_ANISO, 1.2, 0.4,
    0.1, seed=3)
  r = ln.fit(i_obs, tncs_a_groups=a)
  assert r.coherence.gamma == 1.05, r.coherence.gamma
  w = ("tNCS coherence gamma = 1.050 exceeds 1: the data are more strongly "
    "modulated than identical copies would give")
  assert w in r.warnings, r.warnings
  # a lower gamma_max is honoured
  p = ln.normalisation_params.extract()
  p.gamma_max = 0.8
  r = ln.fit(i_obs, tncs_a_groups=a, params=p)
  assert r.coherence.gamma == 0.8 and r.warnings == []

def exercise_outlier_rejected():
  i_obs, a = _simulate(P212121[0], P212121[1], 2.0, B_ANISO, 0.9, 0.4,
    0.1, seed=4)
  r0 = ln.fit(i_obs, tncs_a_groups=a)
  cen = i_obs.centric_flags().data()
  i = [k for k in range(i_obs.size()) if not cen[k] and r0.selection[k]][100]
  sigma = r0.resn.data()[i] ** 2 * r0.teps.data()[i]
  data = i_obs.data().deep_copy()
  data[i] = 40 * sigma
  r = ln.fit(i_obs.customized_copy(data=data), tncs_a_groups=a)
  assert not r.selection[i]
  assert r.n_outliers >= 1
  assert r.n_outlier_rounds >= 1

def exercise_weak_shells():
  ms = make_set(P212121[0], P212121[1], 2.0)
  bins = ln.resolution_bins(ms.d_star_sq().data())
  b = bins.index.as_numpy_array()
  top = b >= bins.n_bins - 2
  sf = numpy.where(top, 10., 0.1)
  i_obs, a = _simulate(P212121[0], P212121[1], 2.0, B_ANISO, 0.9, 0.4,
    sf, seed=5)
  r = ln.fit(i_obs, tncs_a_groups=a)
  assert math.isfinite(r.log_likelihood)
  for k in (bins.n_bins - 2, bins.n_bins - 1):
    c = math.exp(r.x[1 + k])
    assert 0.7 <= c <= 1.4, (k, c)

def exercise_zero_sigma_excluded():
  i_obs, a = _simulate(P212121[0], P212121[1], 2.0, B_ANISO, 0.9, 0.4,
    0.1, seed=6)
  sig = i_obs.sigmas().deep_copy()
  zero = [7, 70, 700, 7000, 17000]
  for k in zero: sig[k] = 0
  sig[zero[-1]] = -1.
  r = ln.fit(i_obs.customized_copy(sigmas=sig), tncs_a_groups=a)
  assert r.n_zero_sigma == 5
  for k in zero: assert not r.selection[k]
  # a given selection that includes them still excludes them
  sel = flex.bool(i_obs.size(), True)
  r = ln.fit(i_obs.customized_copy(sigmas=sig), tncs_a_groups=a,
    selection=sel)
  assert r.n_zero_sigma == 5
  for k in zero: assert not r.selection[k]

def exercise_epsilon():
  i_obs, a = _simulate("P 61 2 2", (50, 50, 90, 90, 90, 120), 2.5,
    (0, 0, 0, 0, 0, 0), 0., 0., 0.1, seed=7)
  eps = i_obs.epsilons().data().as_numpy_array()
  assert (eps > 1).sum() > 10
  r = ln.fit(i_obs)
  sigman = numpy.exp(r.model.ln_sigman(r.x).as_numpy_array())
  ratio = r.resn.data().as_numpy_array() ** 2 / sigman
  assert abs(ratio - eps).max() < 1e-10 * eps.max()
  # beta obeys the symmetry: b11 = b22 = b12, b13 = b23 = 0
  b = r.model.beta(r.x)
  assert abs(b[0] - b[1]) < 1e-14 and abs(b[0] - b[3]) < 1e-14
  assert b[4] == b[5] == 0

def exercise_fixed_selection():
  i_obs, a = _simulate(P212121[0], P212121[1], 2.0, B_ANISO, 0.9, 0.4,
    0.1, seed=8)
  data = i_obs.data().deep_copy()
  data[123] = 50 * data[123] + 1e4    # an outlier, so the selection matters
  i_obs = i_obs.customized_copy(data=data)
  r0 = ln.fit(i_obs, tncs_a_groups=a)
  assert r0.n_outlier_rounds >= 1 and not r0.selection[123]
  r = ln.fit(i_obs, tncs_a_groups=a, selection=r0.selection)
  assert r.selection.all_eq(r0.selection)
  assert r.n_outlier_rounds == 0
  assert r.n_would_change == 0
  assert approx_equal(r.log_likelihood, r0.log_likelihood, eps=1e-3)
  # with everything selected the outlier would change status
  r = ln.fit(i_obs, tncs_a_groups=a, selection=flex.bool(i_obs.size(), True))
  assert r.selection[123] and r.n_would_change >= 1 and r.n_outliers == 0

def exercise_hierarchical():
  i_obs, a = _simulate(P212121[0], P212121[1], 2.0, B_ANISO, 0.9, 0.4,
    0.1, seed=9, n_groups=2)
  p = ln.normalisation_params.extract()
  p.hierarchical_coherence = True
  r = ln.fit(i_obs, tncs_a_groups=a, params=p)
  assert len(r.gamma_groups) == 2
  assert abs(r.gamma_groups[0] - r.gamma_groups[1]) < 0.05, r.gamma_groups
  assert abs(r.coherence.gamma - 0.9) < 0.05, r.coherence.gamma
  # without it, the shared gamma
  r = ln.fit(i_obs, tncs_a_groups=a)
  assert r.gamma_groups is None and abs(r.coherence.gamma - 0.9) < 0.03

def exercise_resolution_dependent_anisotropy():
  i_obs, a = _simulate(P212121[0], P212121[1], 2.0, B_ANISO, 0.9, 0.4,
    0.1, seed=10)
  p = ln.normalisation_params.extract()
  p.resolution_dependent_anisotropy = "off"
  r_off = ln.fit(i_obs, tncs_a_groups=a, params=p)
  assert not r_off.resolution_dependent
  p.resolution_dependent_anisotropy = "on"
  r_on = ln.fit(i_obs, tncs_a_groups=a, params=p)
  assert r_on.resolution_dependent
  assert r_on.model.n_parameters() == r_off.model.n_parameters() + 3
  r_auto = ln.fit(i_obs, tncs_a_groups=a)
  # no true resolution dependence: the gain is below the threshold
  # (chi-squared 99.9% point for 2 parameters, halved)
  assert r_auto.resolution_dependent_threshold == 13.816 / 2
  assert not r_auto.resolution_dependent
  assert r_auto.resolution_dependent_gain < r_auto.resolution_dependent_threshold
  # a resolution-dependent anisotropy in the data is detected
  ms = i_obs
  x = ms.d_star_sq().data().as_numpy_array()
  x = x / x.max()
  hkl = ms.indices().as_vec3_double().as_numpy_array()
  uc = ms.unit_cell()
  u1 = adptbx.u_cart_as_u_star(uc, adptbx.b_as_u((8., -8., 0., 0., 0., 0.)))
  q1 = u1[0] * hkl[:, 0] ** 2 + u1[1] * hkl[:, 1] ** 2 + u1[2] * hkl[:, 2] ** 2
  f = numpy.exp(-2 * math.pi ** 2 * q1 * x)
  r_true = ln.fit(i_obs.customized_copy(
    data=i_obs.data() * flex.double(f), sigmas=i_obs.sigmas() *
    flex.double(f)), tncs_a_groups=a)
  assert r_true.resolution_dependent, r_true.resolution_dependent_gain

def exercise_derivatives_of_target():
  """ the target's gradient and exact Hessian by finite differences, with
  the resolution-dependent anisotropy, tNCS and hierarchical eta """
  i_obs, a = _simulate(P212121[0], P212121[1], 3.5, B_ANISO, 0.9, 0.4,
    0.3, seed=11, n_groups=2)
  data = ln._fit_data(i_obs, a, ln.normalisation_params.extract())
  model = ln.sigman_model(i_obs, data.bins, resolution_dependent=True)
  t = ln._target(data, model, flex.bool(i_obs.size(), True).as_numpy_array(),
    hierarchical=True, sigma_eta=0.1)
  rs = numpy.random.RandomState(3)
  p = t.start(i_obs)
  p[t.i_tncs] = 0.85
  p[t.i_tncs + 1] = 0.3
  p[t.i_tncs + 2:] = [0.05, -0.08]
  p[t.n_q - 2:t.n_q] = rs.normal(size=2) * 1e-4
  p[1 + data.bins.n_bins] = 1e-4   # beta11
  f, g, h = t.evaluate(p, hessian="exact")
  hf = t.evaluate(p, hessian="fisher")[2]
  assert abs(hf - hf.T).max() < 1e-9 * abs(hf).max()
  n = p.size
  assert g.shape == (n,) and h.shape == (n, n)
  for i in range(n):
    step = 1e-5 * max(abs(p[i]), 1e-3)
    pp = p.copy(); pm = p.copy()
    pp[i] += step; pm[i] -= step
    fp, gp = t.evaluate(pp, hessian=None)[:2]
    fm, gm = t.evaluate(pm, hessian=None)[:2]
    fd = (fp - fm) / (2 * step)
    assert abs(fd - g[i]) < 1e-5 * max(1., abs(g[i])), (i, fd, g[i])
    fdh = (gp - gm) / (2 * step)
    assert abs(fdh - h[i]).max() < 1e-4 * max(1., abs(h[i]).max()), (i,
      abs(fdh - h[i]).max(), abs(h[i]).max())

def exercise_newton_fallback():
  """ an indefinite Hessian gives no Newton step, and the Fisher form is
  used instead; the fit converges from the start, where that happens (with
  per-group coherence: gamma and eta_g are degenerate at eta = 0) """
  h = numpy.array([[2., 0.], [0., -1.]])
  g = numpy.array([1., 1.])
  free = numpy.array([True, True])
  assert ln._newton_step(h, g, free) is None
  s = ln._newton_step(h, g, free, ridge=True)
  assert s is not None and g.dot(s) < 0
  s = ln._newton_step(numpy.diag([2., 4.]), g, free)
  assert approx_equal(s, [-0.5, -0.25])
  s = ln._newton_step(numpy.diag([2., 4.]), g, numpy.array([True, False]))
  assert approx_equal(s, [-0.5, 0])
  i_obs, a = _simulate(P212121[0], P212121[1], 2.5, B_ANISO, 0.9, 0.4,
    0.1, seed=14, n_groups=2)
  data = ln._fit_data(i_obs, a, ln.normalisation_params.extract())
  model = ln.sigman_model(i_obs, data.bins)
  t = ln._target(data, model, data.valid, hierarchical=True)
  lower, upper = t.bounds(1.05)
  m = ln._minimise(t, t.start(i_obs), lower, upper, 50)
  assert m.n_fisher >= 1 and m.converged
  assert abs(m.p[t.i_tncs] - 0.9) < 0.05

def exercise_outlier_prefilter():
  """ the cheap bound used to skip tail probabilities gives the same
  outliers as evaluating them everywhere """
  rs = numpy.random.RandomState(12)
  n = 200000
  sigma = numpy.exp(rs.uniform(-3, 3, size=n))
  s = numpy.exp(rs.uniform(-9, 4, size=n))
  e = numpy.concatenate([rs.uniform(-30, 0, size=n // 4),
    rs.uniform(0, 60, size=n // 4), rs.exponential(size=n // 2) * 3])
  e[:n // 8] = -s[:n // 8] * rs.uniform(0, 8, size=n // 8)
  i = e * sigma
  cen = rs.uniform(size=n) < 0.5
  sig = s * sigma
  full = xray_ext.wilson_intensity_tail_probability(i_obs=flex.double(i),
    sig_i_obs=flex.double(sig), sigma=flex.double(sigma),
    centric_flags=flex.bool(cen)).as_numpy_array() < 1e-6
  fast, n_eval = ln._wilson_outliers(i, sig, sigma, cen, 1e-6)
  assert (fast == full).all(), (fast != full).sum()
  assert full.sum() > 1000 and n_eval < n

def exercise_fit_speed():
  i_obs, a = _simulate(P212121[0], P212121[1], 2.0, B_ANISO, 0.9, 0.4,
    0.1, seed=13)
  t0 = time.time()
  ln.fit(i_obs, tncs_a_groups=a)
  dt = time.time() - t0
  assert dt < 10, dt

def run():
  exercise_best_curve_matches_phasertng()
  exercise_bin_restraint_sigmas()
  exercise_bins()
  exercise_beta_constraints()
  exercise_derivatives_fd()
  exercise_restraint_values()
  exercise_isotropic_beta()
  exercise_start()
  exercise_outlier_prefilter()
  exercise_newton_fallback()
  exercise_derivatives_of_target()
  exercise_fit_recovers_truth()
  exercise_no_tncs()
  exercise_gamma_bound()
  exercise_outlier_rejected()
  exercise_weak_shells()
  exercise_zero_sigma_excluded()
  exercise_epsilon()
  exercise_fixed_selection()
  exercise_hierarchical()
  exercise_resolution_dependent_anisotropy()
  exercise_fit_speed()
  print("OK")

if __name__ == "__main__":
  run()
