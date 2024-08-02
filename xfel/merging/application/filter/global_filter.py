from __future__ import division

from collections import Counter
from enum import Enum
from itertools import chain

import numpy as np

from dials.array_family import flex
from dxtbx.model import ExperimentList
from xfel.merging.application.utils.data_counter import data_counter
from xfel.merging.application.worker import worker


def flat_array(iterable):
  """Return the contents of all input iterables flattened into a 1d np.array"""
  return np.array(list(chain.from_iterable(iterable)))


def uniques(*iterables):
  """Return a set of unique elements across all input iterables"""
  return set(sum([list(i) for i in iterables], []))


class FilterReasons(Enum):
  """Enumerator documenting all possible reasons for filtering expts/refls"""
  intensity_extremum_iqr_dist = "Intensity extremum outside IQR dist threshold"
  # add subsequent global filtering reasons here (step 1/3)

  @classmethod
  def max_reason_len(cls):
    return max(len(r.value) for r in cls)

  @classmethod
  def report_line(cls, reason, filtered_expts, filtered_refls):
    """Return a line for the report listing all input in nice format"""
    fmt = '- {:' + str(cls.max_reason_len() + 1) + '} {:6d} expts, {:9d} refls'
    return fmt.format(str(reason) + ':', filtered_expts, filtered_refls)


class GlobalFilter(worker):
  """Filter experiments & reflections based on their aggregated statistics"""

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    self.expt_filter_reasons = Counter()
    self.refl_filter_reasons = Counter()
    super(GlobalFilter, self).__init__(params=params, mpi_helper=mpi_helper,
                                       mpi_logger=mpi_logger)

  def __repr__(self):
    return """Filter expts & refls based on their aggregated statistics"""

  def filter_intensity_extrema(self, expts, refls):
    """Filter expts whose refls' intensity extrema don't fit the population"""
    iqr_lim = self.params.filter_global.intensity_extrema_iqr_dist_threshold
    maxima = []
    minima = []
    refl_list = []
    for iid, expt in enumerate(expts):
      refl = refls.select(refls["id"] == iid)
      maxima.append(max(refl['intensity.sum.value']))
      minima.append(min(refl['intensity.sum.value']))
      refl_list.append(refl)
    all_maxima = flat_array(self.mpi_helper.comm.allgather(maxima))
    all_minima = flat_array(self.mpi_helper.comm.allgather(minima))
    maxima_quartiles = np.nanpercentile(all_maxima, [25, 50, 75])
    minima_quartiles = np.nanpercentile(all_minima, [25, 50, 75])
    maxima_iqr = maxima_quartiles[2] - maxima_quartiles[0]
    minima_iqr = minima_quartiles[2] - minima_quartiles[0]
    maxima_upper_lim = maxima_quartiles[1] + iqr_lim * maxima_iqr
    minima_lower_lim = minima_quartiles[1] - iqr_lim * minima_iqr
    filtered_expts = ExperimentList()
    filtered_refls = flex.reflection_table()
    for expt, refl, minimum, maximum in zip(expts, refl_list, minima, maxima):
      if minimum < minima_lower_lim or maximum > maxima_upper_lim:
        filter_reason = FilterReasons.intensity_extremum_iqr_dist
        self.expt_filter_reasons[filter_reason] += 1
        self.refl_filter_reasons[filter_reason] += refl.size()
      else:
        filtered_expts.append(expt)
        filtered_refls.extend(refl)
    return filtered_expts, filtered_refls

  # implement subsequent global filtering algorithms here (step 2/3)

  def report_filter_reasons(self):
    self.logger.log('Experiments/reflections filtered on this rank due to:')
    for r in uniques(self.expt_filter_reasons, self.refl_filter_reasons):
      te = self.expt_filter_reasons[r]
      tr = self.refl_filter_reasons[r]
      self.logger.log(FilterReasons.report_line(r.value, te, tr))
    te = sum(self.expt_filter_reasons.values())
    tr = sum(self.refl_filter_reasons.values())
    self.logger.log(FilterReasons.report_line('TOTAL', te, tr))
    if self.mpi_helper.rank == 0:
      self.logger.main_log('Experiments/reflections filtered due to:')
      expt_filter_reasons = self.mpi_helper.count(self.expt_filter_reasons)
      refl_filter_reasons = self.mpi_helper.count(self.refl_filter_reasons)
      for r in uniques(expt_filter_reasons, refl_filter_reasons):
        te = expt_filter_reasons[r]
        tr = refl_filter_reasons[r]
        self.logger.main_log(FilterReasons.report_line(r.value, te, tr))
      te = sum(expt_filter_reasons.values())
      tr = sum(refl_filter_reasons.values())
      self.logger.main_log(FilterReasons.report_line('TOTAL', te, tr))

  def run(self, experiments, reflections):
    expts, refls = self.filter_intensity_extrema(experiments, reflections)
    # call subsequent global filtering algorithms here (step 3/3)
    self.report_filter_reasons()
    data_counter(self.params).count(expts, refls)
    return expts, refls
