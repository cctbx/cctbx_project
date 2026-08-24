from __future__ import absolute_import, division, print_function
import numpy as np

METRICS = ('cellparams', 'l1', 'l2', 'ncdist', 's6', 'dc7unsrt', 'v7')
DEFAULT_METRIC = 'ncdist'

# SAUC scaling weights for dc7unsrt: [1,1,1,.5,.5,.5,1/3] applied to (d*d)
_DC7_W = np.array([1., 1., 1., .5, .5, .5, 1. / 3.])


def features(metric, params):
  """
  Compute the feature vector(s) used by `distance` for a given metric.

  Parameters
  ----------
  metric : str
      One of `METRICS`.
  params : array-like, shape (..., 6)
      Reduced cell parameters (a, b, c, alpha, beta, gamma), with
      a, b, c in Angstrom and alpha, beta, gamma in DEGREES.

  Returns
  -------
  ndarray, shape (..., k)
      k = 6 for cellparams/l1/l2/ncdist/s6, k = 7 for dc7unsrt/v7.
      Works unchanged for a single cell (6,), a query's settings (S,6)
      and the whole database (N,6).

  `v7` loops over rows through `cctbx.uctbx` (~9 us/row) because it needs
  two real lattice reductions (Niggli + reciprocal minimum cell); every
  other metric is pure array arithmetic on `params`.
  """
  p = np.asarray(params, dtype=float)
  a, b, c = p[..., 0], p[..., 1], p[..., 2]

  if metric == 'cellparams':
    return p

  if metric == 'l1' or metric == 'l2':
    al, be, ga = np.radians(p[..., 3]), np.radians(p[..., 4]), np.radians(p[..., 5])
    return np.stack(
        [a, b, c, al * (b + c) / 2, be * (a + c) / 2, ga * (a + b) / 2],
        axis=-1)

  if metric in ('ncdist', 's6', 'dc7unsrt'):
    al, be, ga = np.radians(p[..., 3]), np.radians(p[..., 4]), np.radians(p[..., 5])
    g1, g2, g3 = a * a, b * b, c * c
    g4 = 2 * b * c * np.cos(al)
    g5 = 2 * a * c * np.cos(be)
    g6 = 2 * a * b * np.cos(ga)
    g = np.stack([g1, g2, g3, g4, g5, g6], axis=-1)

    if metric == 'ncdist':
      return g

    if metric == 's6':
      return np.stack([
          g4 / 2, g5 / 2, g6 / 2,
          -g1 - g5 / 2 - g6 / 2,
          -g2 - g4 / 2 - g6 / 2,
          -g3 - g4 / 2 - g5 / 2,
      ], axis=-1)

    # dc7unsrt
    body = g1 + g2 + g3 + np.minimum.reduce([
        g4 + g5 + g6, -g4 - g5 + g6, -g4 + g5 - g6, g4 - g5 - g6])
    return np.stack([
        g1, g2, g3,
        g2 + g3 - np.abs(g4),
        g1 + g3 - np.abs(g5),
        g1 + g2 - np.abs(g6),
        body,
    ], axis=-1)

  if metric == 'v7':
    from cctbx import uctbx
    flat = p.reshape(-1, 6)
    out = np.empty((flat.shape[0], 7))
    for i, row in enumerate(flat):
      uc = uctbx.unit_cell(list(row))
      d = sorted(uc.niggli_cell().parameters()[:3])
      r = sorted(1.0 / x for x in uc.reciprocal().minimum_cell().parameters()[:3])[::-1]
      out[i] = d + r + [uc.volume() ** (1. / 3.)]
    return out.reshape(p.shape[:-1] + (7,))

  raise ValueError(metric)


def distance(metric, f1, f2):
  """
  Distance between two feature arrays produced by `features`.

  Parameters
  ----------
  metric : str
      One of `METRICS`.
  f1, f2 : array-like
      Feature arrays (from `features`), broadcastable against each other,
      with the feature index on the last axis.

  Returns
  -------
  ndarray
      The distance with the last axis reduced -- a 0-d ndarray for two
      (k,) inputs, (N,) for (N,k) against (k,). Values carry SAUC's
      reported scaling, so they are directly comparable with the SAUC
      web tool.
  """
  f1 = np.asarray(f1, dtype=float)
  f2 = np.asarray(f2, dtype=float)
  d = f2 - f1

  if metric == 'cellparams':
    return np.abs(d).sum(-1)

  if metric == 'l1':
    return np.abs(d).sum(-1) / np.sqrt(6)

  if metric == 'l2':
    return np.sqrt((d * d).sum(-1))

  if metric == 'ncdist' or metric == 's6':
    return 0.1 * ((d * d).sum(-1)) ** 0.25

  if metric == 'dc7unsrt':
    return 0.1 * ((_DC7_W * d * d).sum(-1)) ** 0.25

  if metric == 'v7':
    return np.sqrt(6 / 7) * np.sqrt((d * d).sum(-1))

  raise ValueError(metric)


# Known consequence: with SAUC's 24-permutation (S6) and 6-permutation
# (DC7) minimisations dropped -- they are redundant once settings are
# matched -- `s6` is a linear reweighting of the same delta-G6 that
# `ncdist` uses, so the two metrics rank similarly. Adding the 24-perm
# minimum back later is a ~5-line running-min loop; not built here.
