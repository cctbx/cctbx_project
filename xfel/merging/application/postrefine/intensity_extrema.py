from itertools import chain

import numpy as np


class IntensityExtremumOutsideThreshold(ValueError):
  pass


class IntensityExtrema(object):
  """Collect, and reject experiments based on, extrema of their intensities"""
  def __init__(self, comm, params):
    self.iqr_lim = params.postrefinement.intensity_extrema_iqr_dist_threshold
    self.comm = comm
    self.maxima_upper_lim = None
    self.minima_lower_lim = None

  def find_limits(self, experiments, reflections):
    maxima = []
    minima = []
    for iid, expt in enumerate(experiments):
      refl = reflections.select(reflections["id"] == iid)
      maxima.append(max(refl['intensity.sum.value']))
      minima.append(min(refl['intensity.sum.value']))
    maxima = np.array(list(chain.from_iterable(self.comm.allgather(maxima))))
    minima = np.array(list(chain.from_iterable(self.comm.allgather(minima))))
    maxima_quartiles = np.nanpercentile(maxima, [25, 50, 75])
    minima_quartiles = np.nanpercentile(minima, [25, 50, 75])
    maxima_iqr = maxima_quartiles[2] - maxima_quartiles[0]
    minima_iqr = maxima_quartiles[2] - maxima_quartiles[0]
    self.maxima_upper_lim = maxima_quartiles[1] + self.iqr_lim * maxima_iqr
    self.minima_lower_lim = minima_quartiles[1] - self.iqr_lim * minima_iqr

  def raise_if_outside_threshold(self):
    ...





