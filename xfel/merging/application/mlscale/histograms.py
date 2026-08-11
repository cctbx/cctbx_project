from __future__ import absolute_import, division, print_function

"""
Per-reflection observed intensity histograms, sampled over resolution and model
intensity.  Called from mlscale_step1 when mlscale.histograms.enable=True.

WHAT THIS IS FOR

The shape of the observed intensity distribution for a single Miller index
carries the partiality information.  Two regimes bracket the behaviour:

  monochromatic / thin slice   rho << 1   the Ewald shell cuts a thin slice
                                          through the RLP, no observation is
                                          ever fully integrated, the
                                          distribution piles up near zero
  pink beam / full integration rho >> 1   the shell contains the whole RLP,
                                          most observations are full, the
                                          distribution is unimodal near I_full

The crossover resolution decides whether a partiality correction matters over
the range being refined, and whether the absolute intensity is identifiable at
all (it is not, wherever p_max < 1).

RESOLUTION SAMPLING

Rows are chosen by TARGET d-SPACING, not by equal-count bins.  The number of
unique reflections in a shell grows as q^2, so a decile of a d-sorted list is
a thin sliver at the high-resolution limit: for d_min = 0.8 A, 90 percent of
unique reflections lie above 1.7 A and equal-count deciles 1-9 never reach low
resolution at all.  Targets equally spaced in 1/d sample the actual range.

WHY THE X AXIS IS I_obs / I_calc

Raw intensities make panels incomparable, since every reflection has its own
I_full.  Dividing by I_calc leaves a single global constant shared by all
panels, so the upper edge of each distribution is proportional to p_max for
that shell, and the variation of the edge with resolution reads out the
crossover directly.

This is the one place in the mlscale pipeline where I_calc is used, and it is
used only for display and selection, never for estimation.

READING THE ZERO PEAK -- IMPORTANT

A pile-up near zero has two possible causes and they must not be confused:

  partiality   the reflection was far off the Ewald sphere; a real, physical
               small value.  The signature is that the pile is WIDE compared
               with the measurement error.
  noise        the reflection is simply too weak to measure; the pile is the
               background-subtraction error distribution centred on zero.
               The signature is that the pile has the width of sigma.

The 'sig/q99' and 'I/sig' columns below separate them.  Where sig/q99 is
comparable to the width of the zero peak, the zero peak is noise and carries
no partiality information.
"""

import os
import numpy as np


def _model_lookup(i_model):
  """dict from asu miller index tuple -> (model intensity, d spacing)."""
  d = i_model.d_spacings().data()
  out = {}
  for idx, val, dd in zip(i_model.indices(), i_model.data(), d):
    out[tuple(idx)] = (float(val), float(dd))
  return out


def select_reflections(params, i_model, global_counts, logger):
  """Choose one Miller index per (resolution target, intensity percentile).

  Runs on rank 0 only.  Returns a list of dicts, or [] if nothing qualifies.
  """
  hp = params.mlscale.histograms
  lookup = _model_lookup(i_model)

  cand = []
  for hkl, count in global_counts.items():
    if count < hp.min_multiplicity:
      continue
    entry = lookup.get(hkl)
    if entry is None:
      continue
    icalc, dspacing = entry
    if not np.isfinite(icalc) or icalc <= 0 or dspacing <= 0:
      continue
    cand.append((hkl, icalc, dspacing, count))

  if not cand:
    logger.main_log("mlscale histograms: no reflection has multiplicity >= %d "
                    "and a positive model intensity; nothing to plot. Lower "
                    "mlscale.histograms.min_multiplicity." % hp.min_multiplicity)
    return []

  q = np.array([1.0 / t[2] for t in cand])
  order = np.argsort(q)
  cand = [cand[i] for i in order]
  q = q[order]
  n = len(cand)

  logger.main_log("mlscale histograms: %d reflections eligible "
                  "(multiplicity >= %d), d from %.2f to %.2f A"
                  % (n, hp.min_multiplicity, 1.0 / q[0], 1.0 / q[-1]))

  # ---- resolution targets, equally spaced in 1/d unless given explicitly
  if hp.target_d_spacings:
    targets = [float(d) for d in hp.target_d_spacings if d and d > 0]
  else:
    nrow = max(1, int(hp.n_resolution_rows))
    qt = np.linspace(q[0], q[-1], nrow + 2)[1:-1] if nrow > 1 \
        else np.array([0.5 * (q[0] + q[-1])])
    targets = [1.0 / v for v in qt]
  targets = sorted(targets, reverse=True)     # low resolution first

  # Pool by a WINDOW IN 1/d, not by a fixed count. An equal-count pool spans a
  # huge d range at low resolution (few reflections there) and a sliver at high
  # resolution, so count-based pools for different targets overlap almost
  # completely -- which is exactly what happened with the decile version.
  qspan = q[-1] - q[0]
  halfwin = 0.5 * hp.resolution_pool_fraction * qspan

  selected = []
  used = set()                                  # global: no HKL in two rows
  for target_d in targets:
    qt = 1.0 / target_d
    lo = int(np.searchsorted(q, qt - halfwin, side='left'))
    hi = int(np.searchsorted(q, qt + halfwin, side='right'))
    if hi - lo < len(hp.intensity_percentiles):
      # widen symmetrically until there is something to sample
      centre = int(np.searchsorted(q, qt))
      need = max(len(hp.intensity_percentiles), 8)
      lo = max(0, centre - need); hi = min(n, centre + need)
    block = [t for t in cand[lo:hi] if t[0] not in used]
    if not block:
      continue
    d_vals = [t[2] for t in block]
    block = sorted(block, key=lambda t: t[1])   # by model intensity
    m = len(block)
    seen = set()
    for pct in hp.intensity_percentiles:
      j = max(0, min(m - 1, int(round(pct / 100.0 * (m - 1)))))
      hkl, icalc, dspacing, count = block[j]
      if hkl in seen or hkl in used:
        continue                                # pool too small; skip dupes
      seen.add(hkl); used.add(hkl)
      selected.append(dict(hkl=hkl, icalc=icalc, d=dspacing, count=count,
                           target_d=target_d, percentile=pct,
                           pool_d_lo=min(d_vals), pool_d_hi=max(d_vals),
                           pool_size=m))
  return selected


def gather_observations(selected_hkls, reflections, mpi_helper):
  """Collect every observation of each selected HKL from all ranks at rank 0.

  Returns, on rank 0, a dict hkl -> (intensities, sigmas).
  """
  comm = mpi_helper.comm
  wanted = set(selected_hkls)

  loc_i = {h: [] for h in wanted}
  loc_s = {h: [] for h in wanted}
  if reflections is not None and reflections.size() > 0 and wanted:
    key = ('miller_index_asymmetric' if 'miller_index_asymmetric' in reflections
           else 'miller_index')
    inten = reflections['intensity.sum.value'].as_numpy_array()
    var = reflections['intensity.sum.variance'].as_numpy_array()
    for i, hkl in enumerate(reflections[key]):
      t = tuple(hkl)
      if t in wanted:
        loc_i[t].append(inten[i])
        loc_s[t].append(var[i])

  local = {h: (np.asarray(loc_i[h], dtype=float),
               np.asarray(loc_s[h], dtype=float))
           for h in wanted if loc_i[h]}
  allparts = comm.gather(local, root=0)

  if mpi_helper.rank != 0:
    return None
  merged = {}
  for h in wanted:
    ii, vv = [], []
    for part in allparts:
      if part and h in part:
        ii.append(part[h][0]); vv.append(part[h][1])
    if ii:
      merged[h] = (np.concatenate(ii), np.sqrt(np.maximum(np.concatenate(vv), 0)))
    else:
      merged[h] = (np.zeros(0), np.zeros(0))
  return merged


def make_plot(params, selected, observations, logger):
  """Draw the panel grid and log the quantitative summary.  Rank 0 only."""
  hp = params.mlscale.histograms
  outdir = params.output.output_dir
  prefix = params.output.prefix

  panels = []
  for s in selected:
    obs, sig = observations.get(s['hkl'], (np.zeros(0), np.zeros(0)))
    ok = np.isfinite(obs) & np.isfinite(sig)
    obs, sig = obs[ok], sig[ok]
    scale = s['icalc'] if (hp.normalize_by_model and s['icalc'] > 0) else 1.0
    panels.append(dict(s, x=obs / scale, sx=sig / scale))

  if not panels:
    logger.main_log("mlscale histograms: no observations gathered.")
    return

  logger.main_log("")
  logger.main_log("mlscale histograms: distribution shape vs resolution")
  logger.main_log("  x = I_obs / I_calc; the overall scale is one global "
                  "constant, so q99 tracks p_max.")
  logger.main_log("  sig/q99 is the measurement error on the same axis. If the "
                  "zero peak is no wider")
  logger.main_log("  than sig/q99, that peak is NOISE, not partiality, and "
                  "carries no shape information.")
  logger.main_log("")
  logger.main_log("  %-15s %4s %7s %6s %8s %8s %8s %8s %8s"
                  % ("hkl", "pct", "d (A)", "N", "median", "q99",
                     "sig/q99", "I/sig", "frac<0"))
  for p in panels:
    x, sx = p['x'], p['sx']
    if x.size == 0:
      logger.main_log("  %-15s %4d %7.2f %6d   (no observations)"
                      % (str(p['hkl']), p['percentile'], p['d'], 0))
      continue
    q99 = float(np.percentile(x, 99))
    msig = float(np.median(sx)) if sx.size else float('nan')
    isig = (float(np.median(x)) / msig) if msig > 0 else float('nan')
    logger.main_log("  %-15s %4d %7.2f %6d %8.3f %8.3f %8.3f %8.2f %8.3f"
                    % (str(p['hkl']).replace(' ', ''), p['percentile'], p['d'],
                       x.size, float(np.median(x)), q99,
                       (msig / q99 if q99 > 0 else float('nan')), isig,
                       float(np.mean(x < 0))))
  logger.main_log("")

  by_row = {}
  for p in panels:
    by_row.setdefault(p['target_d'], []).append(p)

  logger.main_log("  pooled by resolution row:")
  logger.main_log("  %-9s %14s %7s %8s %8s %8s"
                  % ("target d", "pool d range", "N", "med/q99", "sig/q99",
                     "frac<0"))
  for target_d in sorted(by_row, reverse=True):
    ps = by_row[target_d]
    xs = [p['x'] for p in ps if p['x'].size]
    ss = [p['sx'] for p in ps if p['sx'].size]
    if not xs:
      continue
    xs = np.concatenate(xs); ss = np.concatenate(ss)
    q99 = float(np.percentile(xs, 99))
    logger.main_log("  %-9.2f %6.2f - %-6.2f %7d %8.3f %8.3f %8.3f"
                    % (target_d, ps[0]['pool_d_hi'], ps[0]['pool_d_lo'],
                       xs.size, (float(np.median(xs)) / q99 if q99 > 0 else 0),
                       (float(np.median(ss)) / q99 if q99 > 0 else 0),
                       float(np.mean(xs < 0))))
  logger.main_log("")
  logger.main_log("  med/q99 near 0.1 -> thin slice, piled at zero; "
                  "0.6-0.8 -> approaching full integration.")
  logger.main_log("  But check sig/q99 first: if it is of order 1 the data are "
                  "noise-dominated and")
  logger.main_log("  the shape says nothing about partiality.")
  logger.main_log("")

  if hp.save_npz:
    npz_path = os.path.join(outdir, "%s_mlscale_histograms.npz" % prefix)
    payload = {}
    for i, p in enumerate(panels):
      payload['x_%02d' % i] = p['x']
      payload['sx_%02d' % i] = p['sx']
      payload['meta_%02d' % i] = np.array(
          [p['hkl'][0], p['hkl'][1], p['hkl'][2], p['d'], p['icalc'],
           p['target_d'], p['percentile']], dtype=float)
    np.savez_compressed(npz_path, **payload)
    logger.main_log("  raw panel data: %s" % os.path.abspath(npz_path))

  try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
  except ImportError:
    logger.main_log("  matplotlib unavailable; skipped the plot "
                    "(the .npz still holds the data)")
    return

  rows = sorted(by_row, reverse=True)
  pcts = list(hp.intensity_percentiles)
  nrow, ncol = len(rows), len(pcts)
  fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 2.2 * nrow),
                           squeeze=False)

  xlabel = (r'$I_{\rm obs}\ /\ I_{\rm calc}$' if hp.normalize_by_model
            else r'$I_{\rm obs}$')
  for r, target_d in enumerate(rows):
    row = {p['percentile']: p for p in by_row[target_d]}
    for c, pct in enumerate(pcts):
      ax = axes[r][c]
      p = row.get(pct)
      if p is None or p['x'].size == 0:
        ax.text(0.5, 0.5, 'no data', transform=ax.transAxes,
                ha='center', va='center', fontsize=8, color='0.5')
        ax.set_xticks([]); ax.set_yticks([])
        continue
      x, sx = p['x'], p['sx']
      hi = np.percentile(x, 99.5)
      lo = min(0.0, np.percentile(x, 0.5))
      ax.hist(x, bins=np.linspace(lo, hi * 1.1, hp.n_bins),
              color='#4878A8', lw=0)
      ax.axvline(0.0, color='0.15', lw=1.0)
      # measurement error, drawn on the same axis: if the zero peak is no
      # wider than this bar, that peak is noise rather than partiality
      if sx.size:
        msig = float(np.median(sx))
        ylim = ax.get_ylim()[1]
        ax.plot([-msig, msig], [ylim * 0.92] * 2, color='#C44E52', lw=1.8,
                solid_capstyle='butt')
        ax.text(0.02, 0.93, r'$\pm\sigma$', transform=ax.transAxes,
                fontsize=6.5, color='#C44E52', va='top')
      ax.tick_params(labelsize=7)
      ax.set_yticks([])
      ax.spines[['top', 'right', 'left']].set_visible(False)
      ax.set_title("%s  d=%.2f\u00c5  N=%d"
                   % (str(p['hkl']).replace(' ', ''), p['d'], x.size),
                   fontsize=7, pad=3)
      if r == 0:
        ax.text(0.5, 1.34, "I percentile %d" % pct, transform=ax.transAxes,
                ha='center', fontsize=8, color='0.3')
      if c == 0:
        ax.text(-0.20, 0.5, "d \u2248 %.1f \u00c5" % target_d,
                transform=ax.transAxes, rotation=90, va='center', ha='right',
                fontsize=8, color='0.3')
      if r == nrow - 1:
        ax.set_xlabel(xlabel, fontsize=8)

  fig.suptitle("Observed intensity distributions by resolution and model "
               "intensity", fontsize=10)
  fig.tight_layout(rect=[0.02, 0, 1, 0.95])
  png = os.path.join(outdir, "%s_mlscale_histograms.png" % prefix)
  fig.savefig(png, dpi=180)
  plt.close(fig)
  logger.main_log("  plot: %s" % os.path.abspath(png))
