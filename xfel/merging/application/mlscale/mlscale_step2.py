from __future__ import absolute_import, division, print_function

"""
mlscale_step2 -- joint maximum-likelihood scaling and merging.

WHAT IT DOES

Given observations I_jh = g_j * theta_h * p + noise, estimate every theta_h
(the full intensity) and the global shape parameters that set the partiality
distribution.  The per-frame scale factors g_j are INTEGRATED OUT, not fitted.

WHY THE SCALE FACTORS ARE INTEGRATED OUT

There is one g_j per frame and the number of observations per frame never
grows -- around 19 for small-unit-cell serial data.  If g_j is maximised and
plugged in (a "profile" likelihood), the noise it absorbs is counted as
evidence, and theta stays biased however many frames are added.  This is the
Neyman-Scott incidental-parameter problem, and in simulation at 15
reflections per frame it makes the profile estimator WORSE than a plain
Monte Carlo merge.  Integrating g_j against a lognormal prior removes that
bias and additionally shrinks badly-determined frame scales toward the
population instead of letting them run away.

THE LIKELIHOOD TABLE

    I = m*p + eps,   m = g*theta,   eps ~ N(0, sigma)
    x = I/m = p + eps/m

so the density of x depends only on (rho, sigma/m), and

    log f(I ; m, sigma, rho) = log f_x(x ; rho, sigma/m) - ln m

Three table axes therefore suffice, and per-observation sigma is handled
exactly rather than approximated by a shell average.

IDENTIFIABILITY -- READ THE SHAPE SCAN IN THE LOG

theta and rho are degenerate in SCALE: raising rho and lowering every theta
in a shell leaves the intensities alone.  They are not degenerate in SHAPE,
and the shape responds to rho once rho is around 1 or more.  So the shape
parameters are determined by the high-resolution shells and only weakly, if
at all, by the low-resolution ones.  The scan table reports the log
likelihood against the shape parameter; if that curve is flat, the data are
not constraining it and the result rests on the assumed physical form.
"""

import math
import numpy as np

from xfel.merging.application.worker import worker
from xfel.merging.application.merge.merger import merger
from xfel.merging.application.reflection_table_utils import \
    reflection_table_utils as rt_util
from dials.array_family import flex

SQ2PI = math.sqrt(2.0*math.pi)
C_WINDOW = 2.0          # prediction window: C RLP widths past the shell edge

# numpy renamed trapz -> trapezoid in 2.0; cctbx environments still ship 1.x
_trapz = getattr(np, 'trapezoid', None) or np.trapz


# --------------------------------------------------------------------------
# partiality model
# --------------------------------------------------------------------------
def p_of_s(s, rho):
  """Fraction of a Gaussian RLP caught by a top-hat Ewald shell."""
  from scipy.special import erf
  return 0.5*(erf(s + rho/2.0) - erf(s - rho/2.0))


def logp_density(rho, ny=1500):
  """Smooth density of ln p for s0 uniform over the prediction window."""
  S = rho/2.0 + C_WINDOW
  pmax, pmin = p_of_s(0.0, rho), p_of_s(S, rho)
  y = np.linspace(math.log(pmin) - 0.3, math.log(pmax) + 0.3, ny)
  dy = y[1] - y[0]
  sg = np.linspace(0.0, S, 60001)
  pg = p_of_s(sg, rho)
  yy = np.clip(y, math.log(pmin), math.log(pmax))
  F = 1.0 - np.interp(-np.exp(yy), -pg, sg)/S
  F[y <= math.log(pmin)] = 0.0
  F[y >= math.log(pmax)] = 1.0
  return y, np.maximum(np.gradient(F, dy), 0.0)


def mean_partiality(rho):
  """<p> for a single rho.  Does a numerical integration - do not call this
  per observation; use mean_partiality_array()."""
  y, f = logp_density(rho)
  return float(_trapz(np.exp(y)*f, y))


_PBAR_GRID = None


def mean_partiality_array(rho, n=160):
  """<p> for an array of rho, by interpolation on a cached log-spaced grid.

  Calling mean_partiality() once per observation would run a numerical
  integration 10^5 times.  <p>(rho) is smooth and monotonic, so a grid of a
  couple of hundred points is ample.
  """
  global _PBAR_GRID
  rho = np.clip(np.asarray(rho, float), 1e-3, 60.0)
  if _PBAR_GRID is None:
    g = np.exp(np.linspace(math.log(1e-3), math.log(60.0), n))
    _PBAR_GRID = (g, np.array([mean_partiality(r) for r in g]))
  g, v = _PBAR_GRID
  return np.interp(rho, g, v)


def partiality_moments(rho):
  """<p> and var(p) for a single rho."""
  y, f = logp_density(rho)
  pv = np.exp(y)
  m1 = float(_trapz(pv*f, y))
  m2 = float(_trapz(pv*pv*f, y))
  return m1, max(m2 - m1*m1, 1e-12)


class LikelihoodTable(object):
  """log f(I ; m, sigma, rho), tabulated in x = I/m and r = sigma/m.

  TWO THINGS THIS HAS TO GET RIGHT

  1. The x range must scale with r.  x = p + eps/m has tails of width r, so
     a fixed range covers the density only while r is small.

  2. There must be NO hard floor.  A floor is a one-sided constraint: it is
     escaped by making m larger, so the fit inflates theta.  Outside the
     tabulated region we fall back to a moment-matched Gaussian,

         f(I) ~ N(I ; m*<p>,  sqrt(sigma^2 + m^2 var(p)))

     which is exact in the limit that noise dominates -- precisely the
     regime that falls off the table -- and smooth everywhere else.  Weak
     observations then contribute what they should (very little) instead of
     pushing the fit.
  """

  N_SIGMA = 7.0

  def __init__(self, rhos, log_r, nx=400):
    self.rhos = np.asarray(rhos, float)
    self.log_r = np.asarray(log_r, float)
    r = np.exp(self.log_r)
    self.x_lo = -0.35 - self.N_SIGMA*r
    self.x_hi = 1.35 + self.N_SIGMA*r
    self.nx = nx
    mom = [partiality_moments(rho) for rho in self.rhos]
    self.pbar = np.array([m[0] for m in mom])
    self.pvar = np.array([m[1] for m in mom])
    self.tab = np.empty((len(self.rhos), len(self.log_r), nx))
    self.xgs = np.empty((len(self.log_r), nx))
    for jr in range(len(self.log_r)):
      self.xgs[jr] = np.linspace(self.x_lo[jr], self.x_hi[jr], nx)
    for ir, rho in enumerate(self.rhos):
      yb, fb = logp_density(rho)
      ys, fs = yb[::3], fb[::3]
      dy = ys[1] - ys[0]
      pv = np.exp(ys)
      for jr in range(len(self.log_r)):
        rr = r[jr]
        xg = self.xgs[jr]
        K = np.exp(-0.5*((xg[:, None] - pv[None, :])/rr)**2)/(rr*SQ2PI)
        self.tab[ir, jr] = np.log(np.maximum(K @ (fs*dy), 1e-300))
    self.dxs = (self.x_hi - self.x_lo)/(nx - 1.0)
    self.dlr = self.log_r[1] - self.log_r[0]
    self.nr = len(self.log_r)
    self.n_fallback = 0
    self.n_total = 0

  def gaussian_branch(self, I, sigma, m, ridx):
    """Moment-matched Gaussian density of I.  Already a density in I."""
    mu = m*self.pbar[ridx]
    S2 = sigma*sigma + (m*m)*self.pvar[ridx]
    return -0.5*(I - mu)**2/S2 - 0.5*np.log(S2) - math.log(SQ2PI)

  def __call__(self, I, sigma, lnm, ridx):
    m = np.exp(lnm)
    x = I/m
    lr = np.log(sigma) - lnm

    in_r = (lr >= self.log_r[0]) & (lr <= self.log_r[-1])
    jrf = np.clip((lr - self.log_r[0])/self.dlr, 0, self.nr - 1.001)
    j0 = jrf.astype(np.int32)
    fr = jrf - j0

    def at(j):
      ix = (x - self.x_lo[j])/self.dxs[j]
      ok = (ix >= 0) & (ix <= self.nx - 1.001)
      ixc = np.clip(ix, 0, self.nx - 1.001)
      i0 = ixc.astype(np.int32)
      fx = ixc - i0
      a = self.tab[ridx, j, i0]
      b = self.tab[ridx, j, i0 + 1]
      return a + fx*(b - a), ok

    v0, ok0 = at(j0)
    v1, ok1 = at(j0 + 1)
    use_table = in_r & ok0 & ok1
    self.n_fallback += int((~use_table).sum())
    self.n_total += use_table.size
    tabval = v0 + fr*(v1 - v0) - lnm            # -lnm: Jacobian x -> I
    return np.where(use_table, tabval,
                    self.gaussian_branch(I, sigma, m, ridx))

  def fallback_fraction(self):
    return self.n_fallback/max(self.n_total, 1)


# --------------------------------------------------------------------------
# worker
# --------------------------------------------------------------------------
class ml_merger(merger):
  """Marginal-likelihood scaling and merging.

  Inherits the MTZ output machinery from the standard merger.
  """

  def __repr__(self):
    return "Maximum-likelihood scaling and merging (marginal)"

  # ---------------------------------------------------------------- setup

  def validate(self):
    p = self.params.mlscale.fit
    assert p.n_nodes >= 5, "mlscale.fit.n_nodes should be at least 5"
    assert p.max_em_iter >= 1

  def run(self, experiments, reflections):
    self.logger.log_step_time("MLSCALE_FIT")
    comm = self.mpi_helper.comm
    rank = self.mpi_helper.rank

    need = ['mlscale.frame', 'mlscale.qlen', 'mlscale.wavelength',
            'mlscale.bandwidth', 'intensity.sum.value',
            'intensity.sum.variance', 'miller_index_asymmetric']
    if reflections is not None and reflections.size():
      for k in need:
        if k not in reflections:
          raise KeyError("mlscale_step2 needs column '%s'. Run mlscale_step1 "
                         "earlier in dispatch.step_list, and make sure the "
                         "mlscale columns are in input.persistent_refl_cols."
                         % k)

    # ---- gather everything to rank 0.  For small-unit-cell data this is a
    #      few tens of MB at most, and it removes all MPI logic from the fit.
    payload = self.local_payload(reflections)
    parts = comm.gather(payload, root=0)

    merged = rt_util.merged_reflection_table()
    if rank == 0:
      merged = self.fit_and_merge(parts)

    self.logger.log_step_time("MLSCALE_FIT", True)

    # rank 0 holds the whole merged table; the others contribute nothing
    self.gather_and_output_reflections(merged, 'all')
    return None, reflections

  def local_payload(self, reflections):
    if reflections is None or reflections.size() == 0:
      return None
    var = reflections['intensity.sum.variance'].as_numpy_array()
    ok = np.isfinite(var) & (var > 0)
    inten = reflections['intensity.sum.value'].as_numpy_array()
    ok &= np.isfinite(inten)
    hkl = np.array([tuple(h) for h in reflections['miller_index_asymmetric']],
                   dtype=np.int32)
    return dict(
        hkl=hkl[ok],
        I=inten[ok],
        sig=np.sqrt(var[ok]),
        frame=reflections['mlscale.frame'].as_numpy_array()[ok],
        q=reflections['mlscale.qlen'].as_numpy_array()[ok],
        wav=reflections['mlscale.wavelength'].as_numpy_array()[ok],
        bw=reflections['mlscale.bandwidth'].as_numpy_array()[ok])

  # ---------------------------------------------------------------- fit

  def rho_of(self, q, wav, bw, amp, power):
    """rho = shell thickness / RLP width, as a power law in q.

    Physically, T ~ bw*q^2*lambda/2 and w ~ q*eta (mosaicity-limited) or
    w ~ 1/D (size-limited), giving power 1 or 2 respectively.  Only the
    combination is identifiable, so it is parameterised directly:

        rho(q) = amp * (q/q_ref)^power

    with q_ref the high-resolution limit of the data.
    """
    return amp*(q/self.q_ref)**power

  def fit_and_merge(self, parts):
    from scipy.special import logsumexp
    from numpy.polynomial.hermite_e import hermegauss

    fp = self.params.mlscale.fit
    log = self.logger.main_log

    parts = [p for p in parts if p]
    hkl = np.concatenate([p['hkl'] for p in parts])
    I = np.concatenate([p['I'] for p in parts])
    sig = np.concatenate([p['sig'] for p in parts])
    frame = np.concatenate([p['frame'] for p in parts])
    q = np.concatenate([p['q'] for p in parts])
    wav = np.concatenate([p['wav'] for p in parts])
    bw = np.concatenate([p['bw'] for p in parts])

    # dense indices
    uh, ref = np.unique(hkl, axis=0, return_inverse=True)
    uf, fr = np.unique(frame, return_inverse=True)
    NREF, NFRAME, NOBS = len(uh), len(uf), len(I)
    mult = np.bincount(ref, minlength=NREF)
    self.q_ref = float(np.max(q))

    log("")
    log("=" * 74)
    log("mlscale_step2: maximum-likelihood scaling and merging")
    log("=" * 74)
    log("  observations            %d" % NOBS)
    log("  unique reflections      %d" % NREF)
    log("  frames                  %d" % NFRAME)
    log("  reflections per frame   %.1f" % (NOBS/max(NFRAME, 1)))
    log("  multiplicity  mean %.1f  median %d  min %d  max %d"
        % (mult.mean(), int(np.median(mult)), mult.min(), mult.max()))
    for t in (1, 2, 5, 10, 30):
      log("    reflections with multiplicity >= %-4d %6d  (%.1f%%)"
          % (t, int((mult >= t).sum()), 100.0*(mult >= t).mean()))
    log("  d range                 %.2f to %.2f A" % (1.0/q.min(), 1.0/q.max()))
    log("  median I/sigma          %.2f" % float(np.median(I/sig)))
    log("  fraction I < 0          %.3f" % float(np.mean(I < 0)))
    log("")

    # ---- conventional baseline, for comparison
    w = 1.0/sig**2
    t_conv = self.conventional(I, sig, ref, fr, NREF, NFRAME,
                               self.rho_of(q, wav, bw, fp.rho_amp_start,
                                           fp.rho_power))

    # ---- outer loop: scan the shape amplitude
    amps = np.array(fp.rho_amp_scan) if fp.rho_amp_scan else \
        fp.rho_amp_start*np.array([0.5, 0.7, 0.85, 1.0, 1.2, 1.5, 2.0])
    log("  shape scan:  rho(q) = amp * (q/%.3f)^%.2f" % (self.q_ref, fp.rho_power))
    log("  %8s %8s %8s %14s %10s %10s"
        % ("amp", "rho@lo", "rho@hi", "marginal logL", "sd(ln g)",
           "frac Gauss"))
    log("  (frac Gauss is the fraction handled by the moment-matched")
    log("   Gaussian fallback rather than the tabulated shape.  These are")
    log("   noise-dominated observations that carry little shape")
    log("   information; a large fraction is not an error, but it does mean")
    log("   the shape is being determined by a subset of the data.)")
    results = []
    for amp in amps:
      out = self.em_fit(I, sig, ref, fr, NREF, NFRAME,
                        self.rho_of(q, wav, bw, amp, fp.rho_power), t_conv,
                        n_iter=fp.scan_em_iter)
      results.append((amp, out))
      log("  %8.3f %8.3f %8.3f %14.1f %10.3f %10.4f"
          % (amp, amp*(q.min()/self.q_ref)**fp.rho_power, amp,
             out['logL'], out['sg'], out['floor_frac']))

    best_i = int(np.argmax([r[1]['logL'] for r in results]))
    best_amp = results[best_i][0]
    log("")
    log("  best amp = %.3f  (refitting with %d EM iterations)"
        % (best_amp, fp.max_em_iter))
    best = self.em_fit(I, sig, ref, fr, NREF, NFRAME,
                       self.rho_of(q, wav, bw, best_amp, fp.rho_power), t_conv,
                       n_iter=fp.max_em_iter)
    self.report_curvature(amps, [r[1]['logL'] for r in results], best_i, log)

    rho = self.rho_of(q, wav, bw, best_amp, fp.rho_power)
    # recompute the baseline at the SAME rho, or the ratio column mixes two
    # different partiality corrections
    t_conv = self.conventional(I, sig, ref, fr, NREF, NFRAME, rho)
    log("")
    log("  final EM: fraction using the Gaussian fallback: %.4f"
        % best['floor_frac'])
    if best['floor_frac'] > 0.5:
      log("  *** NOTE: most observations are noise-dominated and handled by")
      log("  *** the Gaussian branch. The partiality shape is then being")
      log("  *** determined by a minority of the data.")
    self.report_shells(q, rho, ref, mult, best, t_conv, log)
    self.report_frames(best, log)
    self.report_convergence(best, log)
    self.compare(best, t_conv, q, ref, NREF, log)

    # ---- build the merged table
    theta = np.exp(best['lnt'])
    sigma_t = theta*best['rel_err']
    keep = mult >= max(1, self.params.merging.minimum_multiplicity)
    log("")
    log("  output: %d reflections with multiplicity >= %d"
        % (int(keep.sum()), self.params.merging.minimum_multiplicity))

    idx = np.where(keep)[0]
    if idx.size == 0:
      log("  *** WARNING: no reflection met the multiplicity cutoff.")
      return rt_util.merged_reflection_table()

    # NOTE: rt_util.merged_reflection_table() comes back with its four columns
    # already declared and zero rows, so assigning a length-N column to it
    # trips the flex row-count assertion.  Build a fresh table instead: the
    # first assignment sets the row count, and the schema ends up identical.
    table = flex.reflection_table()
    table['miller_index'] = flex.miller_index(
        [(int(a), int(b), int(c)) for a, b, c in uh[idx]])
    table['intensity'] = flex.double(np.ascontiguousarray(theta[idx],
                                                          dtype=float))
    table['sigma'] = flex.double(np.ascontiguousarray(sigma_t[idx],
                                                      dtype=float))
    table['multiplicity'] = flex.int(np.ascontiguousarray(mult[idx],
                                                          dtype=np.int32))
    return table

  # ---------------------------------------------------------------- pieces

  def conventional(self, I, sig, ref, fr, NREF, NFRAME, rho):
    """Ratio-of-sums scaling plus inverse-variance weighted merge."""
    pbar = mean_partiality_array(rho)
    w = 1.0/sig**2
    t = (np.bincount(ref, w*I, NREF)/np.maximum(np.bincount(ref, w, NREF), 1e-30))
    t = t/np.maximum(np.bincount(ref, w*pbar, NREF) /
                     np.maximum(np.bincount(ref, w, NREF), 1e-30), 1e-6)
    for _ in range(30):
      num = np.bincount(fr, w*I*t[ref]*pbar, NFRAME)
      den = np.bincount(fr, w*(t[ref]*pbar)**2, NFRAME)
      g = np.maximum(num/np.maximum(den, 1e-30), 1e-8)
      g /= math.exp(float(np.mean(np.log(g))))
      ws = w*(g[fr]*pbar)**2
      t = (np.bincount(ref, ws*I/(g[fr]*pbar), NREF) /
           np.maximum(np.bincount(ref, ws, NREF), 1e-30))
    return np.maximum(t, 1e-8)

  def em_fit(self, I, sig, ref, fr, NREF, NFRAME, rho, t_start,
             n_iter=None):
    from scipy.special import logsumexp
    from numpy.polynomial.hermite_e import hermegauss
    fp = self.params.mlscale.fit

    rho = np.clip(rho, 1e-3, 60.0)
    edges = np.linspace(rho.min(), rho.max(), fp.n_rho_bins + 1)
    centres = 0.5*(edges[1:] + edges[:-1])
    ridx = np.clip(np.digitize(rho, edges) - 1, 0,
                   fp.n_rho_bins - 1).astype(np.int32)

    lnt = np.log(np.maximum(t_start, 1e-8))
    lr = np.log(sig) - (lnt[ref] + 0.0)
    log_r = np.linspace(np.percentile(lr, 0.1) - 3.0,
                        np.percentile(lr, 99.9) + 3.0, fp.n_noise_bins)
    lik = LikelihoodTable(centres, log_r)

    t_gh, w_gh = hermegauss(fp.n_nodes)
    w_gh = w_gh/w_gh.sum()

    sg = fp.sg_start
    trace = []
    for it in range(n_iter or fp.max_em_iter):
      nodes = sg*t_gh
      LL = np.empty((fp.n_nodes, NFRAME))
      for k in range(fp.n_nodes):
        LL[k] = np.bincount(fr, lik(I, sig, lnt[ref] + nodes[k], ridx), NFRAME)
      LL += np.log(w_gh)[:, None]
      norm = logsumexp(LL, axis=0)
      post = np.exp(LL - norm[None, :])

      # M-step.  A fixed fine grid here would cost n_nodes * n_steps table
      # lookups per iteration and dominate the run.  Instead use a coarse grid
      # whose range SHRINKS as the fit settles: early iterations need reach,
      # later ones need resolution, and no single fixed grid gives both
      # cheaply.  A grid search (unlike a Newton step on a bilinearly
      # interpolated table) cannot decrease the objective, so EM stays
      # monotone.
      cur = lnt.copy()
      rng_it = max(fp.step_floor, fp.step_range*(fp.step_decay**it))
      STEP = np.linspace(-rng_it, rng_it, fp.n_steps)
      sc = np.zeros((fp.n_steps, NREF))
      for k in range(fp.n_nodes):
        wk = post[k][fr]
        base = nodes[k] + cur[ref]
        for j, st in enumerate(STEP):
          sc[j] += np.bincount(ref, wk*lik(I, sig, base + st, ridx), NREF)
      kk = np.clip(np.argmax(sc, 0), 1, fp.n_steps - 2)
      ar = np.arange(NREF)
      y0, y1, y2 = sc[kk - 1, ar], sc[kk, ar], sc[kk + 1, ar]
      den = y0 - 2*y1 + y2
      sh = np.where(den < 0, -0.5*(y2 - y0)/np.where(den == 0, -1e-9, den), 0.0)
      lnt = cur + STEP[kk] + np.clip(sh, -1, 1)*(STEP[1] - STEP[0])
      curv = den/((STEP[1] - STEP[0])**2)          # negative at a maximum

      m1 = (post*nodes[:, None]).sum(0)
      shift = float(m1.mean())
      lnt += shift                                # pin the global scale
      m2 = (post*(nodes**2)[:, None]).sum(0) - 2*shift*m1 + shift**2
      sg = float(np.clip(math.sqrt(max(float(m2.mean()), 1e-6)),
                         0.02, 3.0))
      trace.append((float(norm.sum()), sg))
      if it > 3 and abs(trace[-1][0] - trace[-2][0]) < fp.tolerance:
        break

    # sigma on ln theta from the curvature of the marginal likelihood
    rel_err = 1.0/np.sqrt(np.maximum(-curv, 1e-8))

    return dict(lnt=lnt, sg=sg, logL=trace[-1][0], trace=trace,
                post=post, nodes=nodes, rel_err=np.clip(rel_err, 1e-4, 10.0),
                n_iter=len(trace), rho=rho, ridx=ridx,
                floor_frac=lik.fallback_fraction())

  # ---------------------------------------------------------------- reports

  def report_curvature(self, amps, lls, best_i, log):
    lls = np.asarray(lls)
    span = lls.max() - lls.min()
    log("  likelihood span across the scan: %.1f log units" % span)
    if span < 3.0:
      log("  *** WARNING: the scan is FLAT. The data do not constrain the")
      log("  *** shape amplitude; the result rests entirely on the assumed")
      log("  *** physical form. Treat the resolution dependence as an")
      log("  *** assumption, not a measurement.")
    elif best_i in (0, len(amps) - 1):
      log("  *** WARNING: the optimum is at the edge of the scan range.")
      log("  *** Widen mlscale.fit.rho_amp_scan.")
    else:
      # least-squares parabola through the three points nearest the peak;
      # the scan spacing is uneven, so an equal-spacing formula is wrong
      a = np.asarray(amps[best_i - 1:best_i + 2], float)
      y = np.asarray(lls[best_i - 1:best_i + 2], float)
      c2, c1, c0 = np.polyfit(a, y, 2)
      if c2 < 0:
        peak = -c1/(2*c2)
        sigma_amp = 1.0/math.sqrt(-2*c2)      # logL falls by 1/2 at 1 sigma
        log("  parabolic optimum amp = %.3f   1-sigma = %.3f" % (peak, sigma_amp))
      else:
        log("  the three points around the optimum are not concave; "
            "treat the maximum as unresolved")

  def report_shells(self, q, rho, ref, mult, best, t_conv, log):
    nb = self.params.statistics.n_bins
    qr = q  # per observation
    edges = np.linspace(qr.min(), qr.max(), nb + 1)
    idx = np.clip(np.digitize(qr, edges) - 1, 0, nb - 1)
    log("")
    log("  per-shell shape and outcome")
    log("  %-16s %8s %8s %8s %8s %10s %10s"
        % ("d (A)", "N obs", "rho", "p_max", "<p>", "ML/conv", "regime"))
    theta = np.exp(best['lnt'])
    for i in range(nb):
      sel = idx == i
      if sel.sum() < 5:
        continue
      r = float(np.median(rho[sel]))
      rr = np.unique(ref[sel])
      ratio = float(np.median(theta[rr]/np.maximum(t_conv[rr], 1e-12)))
      pmax = float(p_of_s(0.0, r))
      log("  %6.2f - %-7.2f %8d %8.3f %8.4f %8.4f %10.4f   %s"
          % (1.0/edges[i], 1.0/edges[i + 1], int(sel.sum()), r, pmax,
             float(mean_partiality_array([r])[0]), ratio,
             "ok" if pmax > 0.95 else "rho weakly determined"))

  def report_frames(self, best, log):
    lng = (best['post']*best['nodes'][:, None]).sum(0)
    spread = np.sqrt(np.maximum(
        (best['post']*(best['nodes']**2)[:, None]).sum(0) - lng**2, 0))
    log("")
    log("  frame scale factors")
    log("    fitted prior width sd(ln g)  %.4f" % best['sg'])
    log("    posterior mean ln g:  sd %.4f  range %.3f to %.3f"
        % (lng.std(), lng.min(), lng.max()))
    log("    per-frame posterior width:  median %.4f  p90 %.4f"
        % (float(np.median(spread)), float(np.percentile(spread, 90))))
    frac = float(np.mean(spread > 0.7*best['sg']))
    log("    frames whose scale is barely determined "
        "(posterior width > 0.7*prior): %.1f%%" % (100*frac))
    if frac > 0.5:
      log("    *** Most frames carry little scale information. The shrinkage")
      log("    *** is doing the work; check reflections per frame.")

  def report_convergence(self, best, log):
    log("")
    log("  EM convergence (%d iterations)" % best['n_iter'])
    tr = best['trace']
    show = list(range(min(3, len(tr)))) + \
        list(range(3, len(tr), max(1, len(tr)//6)))
    for i in sorted(set(show)):
      log("    iter %3d   logL %14.2f   sd(ln g) %.4f" % (i + 1, tr[i][0], tr[i][1]))
    if len(tr) > 2 and abs(tr[-1][0] - tr[-2][0]) > self.params.mlscale.fit.tolerance:
      log("    *** WARNING: not converged. Raise mlscale.fit.max_em_iter.")

  def compare(self, best, t_conv, q, ref, NREF, log):
    theta = np.exp(best['lnt'])
    ok = (t_conv > 0) & (theta > 0)
    lr = np.log(theta[ok]/t_conv[ok])
    lr = lr - np.median(lr)
    log("")
    log("  ML versus conventional merge")
    log("    ln(I_ML / I_conv):  sd %.4f   p10 %+.4f   p90 %+.4f"
        % (lr.std(), float(np.percentile(lr, 10)), float(np.percentile(lr, 90))))
    log("    correlation of ln I: %.4f"
        % float(np.corrcoef(np.log(theta[ok]), np.log(t_conv[ok]))[0, 1]))
    log("    median relative sigma from the ML curvature: %.4f"
        % float(np.median(best['rel_err'])))
    log("")
    log("    A systematic trend of ln(I_ML/I_conv) with resolution IS the")
    log("    partiality correction. It will look like a B-factor change in")
    log("    refinement, so compare ADPs between the two merges - that is")
    log("    where the evidence lives, not in R1 or in flat K values.")


if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(ml_merger)
