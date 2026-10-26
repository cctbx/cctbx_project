from __future__ import absolute_import, division, print_function
""" mmtbx.refinement.llgi_normalisation: the BEST curve, resolution bins and
the SigmaN model with its restraints. """
import math
import numpy
from cctbx import crystal, miller
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import mmtbx.refinement.llgi_normalisation as ln

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

def run():
  exercise_best_curve_matches_phasertng()
  exercise_bin_restraint_sigmas()
  exercise_bins()
  exercise_beta_constraints()
  exercise_derivatives_fd()
  exercise_restraint_values()
  exercise_isotropic_beta()
  exercise_start()
  print("OK")

if __name__ == "__main__":
  run()
