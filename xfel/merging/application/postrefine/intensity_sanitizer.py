from collections import Counter
from itertools import chain

import numpy as np

from dials.array_family import flex
from dxtbx.model import ExperimentList


class IntensityExtremumOutsideThreshold(ValueError):
  pass


class IntensitySanitizer(object):
  """Collect, and reject experiments based on, extrema of their intensities"""
  def __init__(self, comm, postrefinement_params):
    self.iqr_lim = postrefinement_params.intensity_extrema_iqr_dist_threshold
    self.comm = comm
    self.exception_counter = Counter()
    self.sanitized_expts = ExperimentList()
    self.sanitized_refls = flex.reflection_table()

  def sanitize(self, experiments, reflections):
    """Generate sanitized_expts, _refls & exception counter based on input"""
    maxima = []
    minima = []
    for iid, expt in enumerate(experiments):
      refl = reflections.select(reflections["id"] == iid)
      maxima.append(max(refl['intensity.sum.value']))
      minima.append(min(refl['intensity.sum.value']))
    all_maxima = np.array(list(chain.from_iterable(self.comm.allgather(maxima))))
    all_minima = np.array(list(chain.from_iterable(self.comm.allgather(minima))))
    maxima_quartiles = np.nanpercentile(all_maxima, [25, 50, 75])
    minima_quartiles = np.nanpercentile(all_minima, [25, 50, 75])
    maxima_iqr = maxima_quartiles[2] - maxima_quartiles[0]
    minima_iqr = minima_quartiles[2] - minima_quartiles[0]
    maxima_upper_lim = maxima_quartiles[1] + self.iqr_lim * maxima_iqr
    minima_lower_lim = minima_quartiles[1] - self.iqr_lim * minima_iqr
    for iid, expt in enumerate(experiments):
      if minima[iid] < minima_lower_lim or maxima[iid] > maxima_upper_lim:
        self.exception_counter[repr(IntensityExtremumOutsideThreshold())] += 1
        continue
      refl = reflections.select(reflections["id"] == iid)
      self.sanitized_expts.append(expt)
      self.sanitized_refls.extend(refl)
