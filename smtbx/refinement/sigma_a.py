""" alpha and beta (sigma_A) for maximum-likelihood refinement

An adapter onto the estimator in mmtbx/max_lik, which implements the
Lunin-Skovoroda algorithm in C++ with the binning, smoothing and interpolation
around it in maxlik.py. Returns alpha and beta per reflection, ordered as the
observations are.

mmtbx is imported lazily, within the functions that need it: smtbx does not
depend on mmtbx and must import without it, so a caller requesting maximum
likelihood where it is unavailable receives a single diagnostic rather than an
ImportError at module load. The smtbx build can supply the extension without
the mmtbx module being configured; see smtbx/max_lik/SConscript.

References:
 - Lunin & Skovoroda (1995) Acta Cryst. A51, 880-887.
 - Read (1986) Acta Cryst. A42, 140-149.
"""
from __future__ import absolute_import, division, print_function

from cctbx.array_family import flex


class missing_mmtbx(RuntimeError):
  pass


def _maxlik():
  try:
    from mmtbx.max_lik import maxlik
  except ImportError as e:
    raise missing_mmtbx(
      "maximum-likelihood refinement needs the alpha/beta estimator from "
      "mmtbx.max_lik, which is not available in this installation (%s). "
      "Least-squares refinement is unaffected." % e)
  return maxlik


def alpha_beta(f_obs, f_calc, r_free_flags,
               free_reflections_per_bin=140,
               interpolation=True,
               add_sigma_squared_to_beta=False):
  """ Per-reflection alpha and beta, as two flex.double aligned with f_obs.

  `f_obs` and `f_calc` are miller arrays over the same indices; `f_calc` may be
  complex and is reduced to amplitudes here. `r_free_flags` is a flex.bool with
  True for the test set.

  `add_sigma_squared_to_beta` folds the experimental variance into beta, the
  convention Refmac uses. Without it the maximum-likelihood weight
  2*alpha^2/(epsilon*beta) depends only on the resolution shell, so every
  reflection in a shell is weighted alike and the measured sigmas play no part
  at all. Leave it off for the intensity target, which models the experimental
  error explicitly by convolution and would otherwise count it twice.
  """
  maxlik = _maxlik()
  assert f_obs.indices().size() == f_calc.indices().size()
  assert r_free_flags.size() == f_obs.indices().size()
  amplitudes = abs(f_calc)
  epsilons = f_obs.epsilons().data().as_double()
  manager = maxlik.alpha_beta_est_manager(
    f_obs=f_obs,
    f_calc=amplitudes,
    free_reflections_per_bin=free_reflections_per_bin,
    flags=r_free_flags,
    interpolation=interpolation,
    epsilons=epsilons)
  alpha, beta = manager.alpha_beta()
  a = alpha.data().deep_copy()
  b = beta.data().deep_copy()
  if add_sigma_squared_to_beta and f_obs.sigmas() is not None:
    b += f_obs.sigmas()*f_obs.sigmas()
  return a, b


def centric_flags_and_epsilons(f_obs):
  """ The other two per-reflection quantities the targets need.

  Both come from the symmetry rather than from the data, and both are asked for
  by every likelihood target, so they are gathered in one place instead of at
  each call site.
  """
  return (f_obs.centric_flags().data(),
          f_obs.epsilons().data().as_double())


def deterministic_free_flags(f_obs, fraction=0.1):
  """ A test set derived from the data, reproducible without being stored

  The set is a pure function of the Miller indices, the unit cell and the space
  group, so the same data always yields the same free reflections and no flag
  column need be written. It follows that different reflections give a
  different set, as a test set carried over onto other data would not be one.

  The global flex random generator is not used: seeding it would affect every
  other consumer of randomness in the process, and the current seed cannot be
  read back to restore it. Hashing each index requires no state.

  f_obs is a unique set in the asymmetric unit, so symmetry-related reflections
  cannot be split between the work and test sets.
  """
  assert 0 < fraction < 1
  uc = f_obs.unit_cell().parameters()
  # a seed from the crystal rather than from the clock
  seed = hash((f_obs.space_group().type().number(),
               tuple(round(p, 4) for p in uc),
               f_obs.indices().size())) & 0x7fffffff
  flags = flex.bool(f_obs.indices().size(), False)
  for i, h in enumerate(f_obs.indices()):
    # a cheap integer mix; the constants are the usual odd multipliers
    x = (h[0]*73856093) ^ (h[1]*19349663) ^ (h[2]*83492791) ^ seed
    x &= 0xffffffff
    x = (x ^ (x >> 16))*0x45d9f3b & 0xffffffff
    x = (x ^ (x >> 16))*0x45d9f3b & 0xffffffff
    x = x ^ (x >> 16)
    flags[i] = ((x % 1000000)/1000000.0) < fraction
  return flags


def is_available():
  """ Whether alpha/beta can be estimated at all in this installation. """
  try:
    _maxlik()
  except missing_mmtbx:
    return False
  return True
