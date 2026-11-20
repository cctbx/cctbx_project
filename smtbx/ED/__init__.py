from __future__ import absolute_import, division, print_function

"""Helpers for dynamical electron diffraction that belong in Python.

The calculation itself is all in C++ (smtbx/ED/*.h). What lives here is the
part that has to make a judgement rather than a computation.
"""

import math

from cctbx.array_family import flex


def weighted_intensity_change(a, b, sigmas):
  """How far apart two sets of calculated intensities are, as a fit would see.

  The obvious comparison -- how much the intensities moved, relative to their
  own size -- is the wrong one, and misleadingly so: on a real dataset it reads
  1e-6 while the refinement it feeds is still 18% away from convergence. The
  intensities barely move in absolute terms, but they move systematically, and
  a least-squares fit weighs them by 1/sigma^2 against the observations rather
  than by how large they are.

  So this is the change in the shape a wR2 is built from,

      sqrt( sum w (a - b)^2 / sum w b^2 ),   w = 1/sigma^2

  which does track the fit: over three decades on TyrosineED it stays
  proportional to the error in the refined objective, to within a factor of
  about two.
  """
  num = den = 0.0
  for x, y, s in zip(a, b, sigmas):
    if s <= 0:
      continue
    w = 1.0/(s*s)
    num += w*(x - y)**2
    den += w*y*y
  return math.sqrt(num/den) if den > 0 else 0.0


def estimate_int_points(data, params, indices, sigmas, n_start=5, n_max=160,
                        n_sample=400, tol=1e-3, log=None):
  """How finely the rocking curve has to be sampled, for this crystal.

  The cost of a dynamical refinement is very nearly linear in the number of
  integration points, and the number needed is not a constant: the rocking
  curve oscillates at a rate set by the specimen thickness, so a thicker
  crystal needs a finer integration and any fixed default is wrong for
  somebody. Left at a default that is too coarse the error does not announce
  itself -- on TyrosineED the usual setting of 10 leaves the refined thickness
  about 20% out, with nothing in the output to say so.

  A closed form does not work either. Two-beam theory gives the rocking curve
  as sin^2(pi t s)/(pi s)^2, so sampling it should need about 2*span*t points;
  on the same data that is out by a factor of three, because many-beam coupling
  makes the real oscillation faster than the two-beam one. So the number is
  found by measurement: the density is doubled until doing so stops changing
  the answer.

  The density is doubled until doing so stops mattering, and then the answer is
  interpolated rather than taken from the ladder. Doubling brackets the
  requirement but overshoots it by up to a factor of two, and the cost of the
  refinement is very nearly linear in this number, so the overshoot is real
  time. The error falls as a power of the density, so two rungs give the
  exponent and the crossing can be solved for directly, at no extra cost.

  @param data an N_beam_shared_data, which computes the intensities
  @param params its refinement_params, borrowed and restored -- the integration
    step is the one thing this changes, and it puts it back
  @param indices the reflections to judge on
  @param sigmas their standard uncertainties, for the weights
  @param n_sample how many reflections to actually use, taken at an even
    stride so the choice is reproducible. The measure is statistical and a few
    hundred settle it, where using every reflection makes deciding cost about
    as much as the refinement it is meant to speed up. None uses all of them.
  @param tol the weighted change to converge to. 1e-3 corresponds to roughly
    1% in the refined objective, measured on TyrosineED.
  @param log a file-like object to report the ladder to, or None
  @return the number of integration points to use
  """
  span = params.int_span
  original_step = params.int_step

  if n_sample is not None and n_sample < len(indices):
    stride = len(indices)//n_sample
    sel = list(range(0, len(indices), stride))
    judged = flex.miller_index([indices[i] for i in sel])
    judged_sigmas = [sigmas[i] for i in sel]
  else:
    judged, judged_sigmas = indices, sigmas

  # The width cache sets how far each beam is integrated over; the step only
  # sets how many samples fall inside that, so the cache does not depend on the
  # density and is built once rather than at every rung. Rebuilding it each
  # time made deciding cost about as much as the refinement it informs.
  data.build_width_cache()

  def intensities(n):
    params.int_step = span/n
    return list(data.compute_dynI(judged))

  try:
    n = n_start
    previous = intensities(n)
    ladder = []
    while n*2 <= n_max:
      finer = intensities(n*2)
      change = weighted_intensity_change(previous, finer, judged_sigmas)
      ladder.append((n, change))
      if log is not None:
        log.write("  integration %3d -> %3d points: weighted change %.2e\n"
                  % (n, 2*n, change))
      if change < tol:
        break
      n *= 2
      previous = finer

    if ladder[-1][1] >= tol:      # never converged; the cap is the answer
      return n
    if len(ladder) < 2:           # converged immediately, nothing to fit
      return ladder[-1][0]

    # change(n) ~ n^-p over the last two rungs; solve change(n) = tol
    (n1, e1), (n2, e2) = ladder[-2], ladder[-1]
    if not (e1 > 0 and e2 > 0 and e1 > e2):
      return n2
    p = math.log(e1/e2)/math.log(n2/n1)
    wanted = n1*(e1/tol)**(1.0/p)
    return int(max(n_start, min(n_max, math.ceil(wanted))))
  finally:
    # the caller's setting is not ours to change; the choice is the return value
    params.int_step = original_step
