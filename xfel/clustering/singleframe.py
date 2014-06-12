from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import math
from libtbx import easy_pickle
import logging
from cctbx.array_family import flex
import cPickle

class SingleFrame:
  """ Class that creates single-image agregate metrics/scoring that can then be
  used in downstream clustering or filtering procedures.
  """

  def __init__(self, path, filename, crystal_num=0):
    try:
      # Warn on error, but continue directory traversal.
      d = easy_pickle.load(path)
      self.is_polarization_corrected = False
      self.miller_array = d['observations'][crystal_num]
      self.path = path
      self.name = filename
      self.pg = d['pointgroup']
      self.uc = d['current_orientation'][crystal_num].unit_cell() \
        .niggli_cell() \
        .parameters()
      self.orientation = d['current_orientation'][crystal_num]
      self.total_i = d['observations'][crystal_num].sum()
      self.wavelength = d['wavelength']

      # Do polarization correction
      self.polarization_correction()
      self.minus_2B, self.G, self.log_i, \
          self.one_over_d_square, self.wilson_err = self.calc_wilson()
      if logging.Logger.root.level < logging.DEBUG:  # Extreme debug!
        self.plot_wilson()
      logging.debug("Extracted image {}".format(filename))
    except KeyError:
      logging.warning("Could not extract point group and unit cell from %s" % path)
    except (cPickle.UnpicklingError, ValueError, EOFError):
      logging.warning("Could not read %s. It may not be a pickle file." % path)

  def trim_res_limit(self, d_min=None, d_max=None):
    """ Remove all miller indicies outside the range of _d_min, _d_max.
    Changes the object in place.
    """
    if d_min is None:
      d_min = self.miller_array.d_min()
    if d_max is None:
      d_max = self.miller_array.d_max_min()[0]
    self.miller_array = self.miller_array.resolution_filter(d_max, d_min)

  def calc_wilson(self):
    """ Do a linear regression to fit G and B. Returns the coeficients minus_2B,
    G, the transformed data log_i, and one_over_d_sqare. Also returns fit_stats,
    which is a dictionairy.
    """
    log_i = flex.log(self.miller_array.sort().data())\
        .as_numpy_array()
    one_over_d_square = self.miller_array.sort().sin_theta_over_lambda_sq().data()\
        .as_numpy_array()

     # Discard negatives ToDo: one could get the mod of the negative values,
     # then plot them as negative in the linear fit.
    log_i, one_over_d_square = zip(*[i for i in zip(log_i, one_over_d_square)
                                  if i[0] >=0])

    minus_2B, G, r_val, _, std_err = linregress(one_over_d_square, log_i)

    # ignore p_val since this will be insanely small
    logging.debug("G: {}, -2B: {}, r: {}, std_err: {}".
      format(G, minus_2B, r_val, std_err))
    return minus_2B, G, log_i, one_over_d_square, {"R": r_val,
                                                   "Standard Error": std_err}

  def plot_wilson(self, width=30, ax=None):
    """ Makes a log(I) vs 1/d**2 plot, displaying the raw partial data, a
    rolling average of the data, and the Wilson model fit to the data.

    Params:
    - width: smoothing window size
    """

    if ax is None:
      fig = plt.figure()
      ax = fig.gca()
      direct_visualisation = True
    else:
      direct_visualisation = False

    smooth = self._moving_average(self.log_i, n=width)
    ax.plot(self.one_over_d_square[width - 1:], smooth,
          '--r', lw=3)
    ax.plot(self.one_over_d_square, self.log_i, 'bo', ms=2)
    ax.plot([0, -1 * self.G / self.minus_2B], [self.G, 0], 'y-', lw=2)
    plt.xlim(0, max(self.one_over_d_square))
    plt.xlabel("(sin(theta)/lambda)^2")
    plt.ylabel("ln(I)")
    plt.title("Single frame Wilson fit\n{}\nG: {}, B: {}, r: {}, std_err: {}".
              format(self.name, self.G, -1 * self.minus_2B / 2,
                     self.wilson_err['R'], self.wilson_err['Standard Error']))

    if direct_visualisation:
      plt.show()
    return ax

    """ Spline method removed because it will be v.slow
    from scipy.interpolate import UnivariateSpline as Spline
    from numpy import linspace
    xs = linspace(min(self.one_over_d_square), max(self.one_over_d_square), 100)
    spl = Spline(self.one_over_d_square, self.log_i, s=10000)
    ys = spl(xs)
    plt.plot(xs, ys, '--g', lw=3)
    """
    """ idiomatic CCTBX method removed because I want more fine-grained detail
     _d_star_p = 1.618034  # Golden ratio distribution for d-spacings
     binner = self.miller_array.setup_binner(n_bins=nbins)
     #logging.debug(str("{}".format(binner.show_summary())))
     bin_selections = [binner.selection(i) for i in binner.range_used()]
     means = [self.miller_array.select(sel).mean() for sel in bin_selections]
     log_means = [math.log(mil) if mil > 0 else 0 for mil in means]
     centers = binner.bin_centers(_d_star_p)
     d_centers = centers ** (-1 / _d_star_p)
     plt.plot(1/(d_centers**2), log_means)
     plt.show()
     """

  def polarization_correction(self):
    """ Perform basic polarization correction in place, and change the
    is_polarization_corrected flag to True.

    I_corrected = 2*I_uncorrected/(1 + cos(two_theta)**2)
    """
    two_theta = self.miller_array.two_theta(wavelength=self.wavelength).data()
    one_over_P = 2/(1 + (flex.cos(two_theta) ** 2))
    self.miller_array = self.miller_array.customized_copy(
      data=self.miller_array.data() * one_over_P)
    self.is_polarization_corrected = True

  @staticmethod
  def _moving_average(array, n=50):
    """ quick method for moving average, needed for smoothing plots. Implements
    a summer area table approach."""
    tmp = np.cumsum(array, dtype=float)
    tmp[n:] = tmp[n:] - tmp[:-n]
    return tmp[n - 1:] / n
