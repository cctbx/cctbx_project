# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

from __future__ import absolute_import, division, print_function

import math

from cctbx.array_family import flex
from scitbx.lbfgs.tst_curvatures import lbfgs_with_curvatures_mix_in
from six.moves import range


class find_scale(lbfgs_with_curvatures_mix_in):
  def __init__(self, scaler, params):

    """This function is largely redundant, because it duplicates what is
    done during mark1 scaling.

    @param scaler Database structure of scaling input
    @param params work_params
    """

    # Extract an ordered union of all Miller indices observed on all
    # frames, and the database structure of observations.
    self._millers = scaler.millers['merged_asu_hkl']
    self._observations = scaler._observations

    # XXX Could be more clever about this here, because this will
    # determine scale factors for rejected frames as well!  Better
    # named selected_frames?
    self._subset = scaler.frames['data_subset']

    self._data = self._observations.get_double('i')
    self._hkl = self._observations.get_int('hkl_id')
    self._sigmas = self._observations.get_double('sigi')
    self._frames = self._observations.get_int('frame_id')

    # XXX Useless assert?
    assert len(self._hkl) == len(self._data) \
      and  len(self._hkl) == len(self._sigmas)

    # Initialise all per-frame scale factors to one.
    n_frames = len(self._subset)
    self.x = flex.double(n_frames + len(self._millers))
    for i in range(n_frames):
      self.x[i] = 1

    # For each Miller index, the weighted (XXX) average intensity of
    # all the observations serves as an initial estimate of the merged
    # intensity.  This is all Monte Carlo scaling would do.
    assert len(self._millers) == len(scaler.summed_wt_I) \
      and  len(self._millers) == len(scaler.summed_weight)

    for i in range(len(self._millers)):
      if scaler.summed_weight[i] > 0:
        self.x[n_frames + i] = scaler.summed_wt_I[i] / scaler.summed_weight[i]

    # The weight of each observation is (1 / sigma)**2, where sigma is
    # the standard deviation of the observation as determined during
    # integration.  An observation is assigned a weight of zero if
    #
    #   The observation was made on a rejected frame
    #
    #   The integrated intensity of the observation is non-positive
    #
    #   The variance of the observation, s**2, as determined during
    #   integration, is non-positive
    #
    #   The d-spacing of the observation lies outside the
    #   user-supplied resolution limits
    #
    # XXX Check Bolotovsky et al.: use sigma**2 or sigma for the
    # weighting?
    self.w = flex.double(len(self._hkl))
    for i in range(len(self.w)):
      if not self._subset[self._frames[i]]:
        continue

      if not params.include_negatives and self._data[i] <= 0:
        continue

      # XXX Should compare against sqrt(eps) instead?  See also
      # scales_non_positive below.
      v = self._sigmas[i]**2
      if v <= 0:
        continue

      # Test d_min first, because it is more likely to have a lower
      # resolution limit than an upper resolution limit.  XXX Is this
      # ever enforced in practice, i.e. is this the first time the
      # limits are applied?
      d = scaler.params.target_unit_cell.d(self._millers[self._hkl[i]])
      if (params.d_min is not None and d < params.d_min) or \
         (params.d_max is not None and d > params.d_max):
        continue

      self.w[i] = 1 / v

    # Should be the last call in the application-specific minimizer
    # class.  This will call lbfgs's run() function and perform
    # optimization.
    super(find_scale, self).__init__() #max_iterations=2000


  def compute_functional_and_gradients(self):
    """The compute_functional_and_gradients() function

    @return Two-tuple of the value of the functional, and an
            <code>n</code>-long vector with the values of the
            gradients at the current position
    """

    #from libtbx.development.timers import Profiler
    from xfel import compute_functional_and_gradients

    n_frames = len(self._subset)

    #p = Profiler("compute_functional_and_gradients [C++]")
    (f, g) = compute_functional_and_gradients(
      self.x, self.w, n_frames, self._observations)
    #del p

    # XXX Only output this every 100 iterations or so.
    scales = self.x[0:len(self._subset)]
    stats = flex.mean_and_variance(scales)
    print("* f =% 10.4e, g =% f+/-%f" % (
      math.sqrt(f),
      stats.mean(),
      stats.unweighted_sample_standard_deviation()))

    # Warn if there are non_positive per-frame scaling factors.
    scales_non_positive = scales.select(scales <= 1e-6) # XXX Or just zero!
    if len(scales_non_positive) > 0:
      stats = flex.mean_and_variance(scales_non_positive)
      if len(scales_non_positive) > 1:
        sigma = stats.unweighted_sample_standard_deviation()
      else:
        sigma = 0
      print("Have %d non-positive per-frame scaling factors: " \
        "%f+/-%f [%f, %f]" % (
          len(scales_non_positive),
          stats.mean(),
          sigma,
          flex.min(scales_non_positive),
          flex.max(scales_non_positive)))

    return (f, g)


  def curvatures(self):
    from xfel import curvatures
    n_frames = len(self._subset)
    return curvatures(
      self.x, self.w, n_frames, self._observations)


  def get_scaling_results(self, results, scaler):
    from xfel import get_scaling_results_mark2

    return get_scaling_results_mark2(
      self.x, self.w, results, scaler.params.target_unit_cell)
