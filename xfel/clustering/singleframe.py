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

  def __init__(self, path, filename, crystal_num=0, remove_negative=False):
    try:
      # Warn on error, but continue directory traversal.
      d = easy_pickle.load(path)
      self.is_polarization_corrected = False
      # Miller arrays
      self.miller_array = d['observations'][crystal_num].sort()
      self.mapped_predictions = d['mapped_predictions'][0]
      # Image pickle info
      self.path = path
      self.name = filename
      # Unit cell info
      self.crystal_system = self.miller_array.crystal_symmetry()\
        .space_group().crystal_system()
      self.pg = d['pointgroup']
      self.uc = d['current_orientation'][crystal_num].unit_cell() \
        .niggli_cell() \
        .parameters()
      self.orientation = d['current_orientation'][crystal_num]
      # Agregate info
      self.total_i = d['observations'][crystal_num].sum()
      self.xbeam = d['xbeam']
      self.ybeam = d['ybeam']
      self.wavelength = d['wavelength']

      if remove_negative:
        self.filter_negative_intensities()

      # Do polarization correction
      self.polarization_correction()
      self.minus_2B, self.G, self.log_i, \
          self.sinsqtheta_over_lambda_sq, self.wilson_err = self.init_calc_wilson()

      if logging.Logger.root.level < logging.DEBUG:  # Extreme debug!
        self.plot_wilson()
      logging.debug("Extracted image {}".format(filename))
    except KeyError:
      logging.warning("Could not extract point group and unit cell from %s" % path)
    except (cPickle.UnpicklingError, ValueError, EOFError):
      logging.warning("Could not read %s. It may not be a pickle file." % path)

  def trim_res_limit(self, d_min=None, d_max=None):
    """
    Remove all miller indicies outside the range of _d_min, _d_max.
    Changes the object in place.

    :param d_min: min res of new miller array. Defaults to current value.
    :param d_max: max res of new miller array. Defaults to current value.
    """
    if d_min is None:
      d_min = self.miller_array.d_min()
    if d_max is None:
      d_max = self.miller_array.d_max_min()[0]
    self.miller_array = self.miller_array.resolution_filter(d_max, d_min).sort()

  def filter_negative_intensities(self):
    """
    Filters negative intensities from the Miller array. Acts in place.
    :return: acts in place.
    """
    i_I_positive = (self.miller_array.data() > 0)
    self.miller_array = self.miller_array.select(i_I_positive).sort()
    self.mapped_predictions = self.mapped_predictions.select(i_I_positive)


  def init_calc_wilson(self):
    """ Do a linear regression to fit G and B. Returns the coeficients minus_2B,
    G, the transformed data log_i, and one_over_d_sqare. Also returns fit_stats,
    which is a dictionairy.

    :return: minus_2B (gradient of fit),
             G (intercept of fit),
             log_i (dependent variable of fit),
             one_over_d_square (independent variable of fit).
    """
    log_i = flex.log(self.miller_array.sort().data())\
        .as_numpy_array()
    sinsqtheta_over_labmdasq = self.miller_array.sort()\
      .sin_theta_over_lambda_sq().data().as_numpy_array()

     # Discard negatives ToDo: one could get the mod of the negative values,
     # then plot them as negative in the linear fit.
    log_i, sinsqtheta_over_labmdasq = zip(*[i for i in zip(log_i, sinsqtheta_over_labmdasq)
                                  if i[0] >=0])

    minus_2B, G, r_val, _, std_err = linregress(sinsqtheta_over_labmdasq, log_i)

    # ignore p_val since this will be insanely small
    logging.debug("G: {}, -2B: {}, r: {}, std_err: {}".
      format(G, minus_2B, r_val, std_err))
    return minus_2B, G, log_i, sinsqtheta_over_labmdasq, {"R": r_val,
                                                   "Standard Error": std_err}

  def plot_wilson(self, width=30, ax=None):
    """ Makes a log(I) vs 1/d**2 plot, displaying the raw partial data, a
    rolling average of the data, and the Wilson model fit to the data.

    :param: width: smoothing window size
    :param: ax: optional axes object to ve used for plotting
    """

    if ax is None:
      fig = plt.figure()
      ax = fig.gca()
      direct_visualisation = True
    else:
      direct_visualisation = False

    smooth = self._moving_average(self.log_i, n=width)
    ax.plot(self.sinsqtheta_over_lambda_sq[width - 1:], smooth,
          '--r', lw=3)
    ax.plot(self.sinsqtheta_over_lambda_sq, self.log_i, 'bo', ms=2)
    ax.plot([0, -1 * self.G / self.minus_2B], [self.G, 0], 'y-', lw=2)
    plt.xlim(0, max(self.sinsqtheta_over_lambda_sq))
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

  def distance_from(self, other_uc):
    """
    Calculares to distance using NCDist from Andrews and Bernstein J. Appl.
    Cryst. 2014 between this frame and some other unit cell.
    :param other_uc: a 6-tuple of a, b, c, alpha, beta, gamma for some unit cell
    :return: the NCDist in A^2 to other_uc
    """
    from cctbx.uctbx.determine_unit_cell import NCDist
    self_g6 = self.make_g6(self.uc)
    other_g6 = self.make_g6(other_uc)
    return NCDist(self_g6, other_g6)

  @staticmethod
  def _moving_average(array, n=50):
    """ quick method for moving average, needed for smoothing plots. Implements
    a summer area table approach."""
    tmp = np.cumsum(array, dtype=float)
    tmp[n:] = tmp[n:] - tmp[:-n]
    return tmp[n - 1:] / n

  @staticmethod
  def make_g6(uc):
      """ Take a reduced Niggli Cell, and turn it into the G6 representation """
      a = uc[0] ** 2
      b = uc[1] ** 2
      c = uc[2] ** 2
      d = 2 * uc[1] * uc[2] * math.cos(uc[3])
      e = 2 * uc[0] * uc[2] * math.cos(uc[4])
      f = 2 * uc[0] * uc[1] * math.cos(uc[5])
      return [a, b, c, d, e, f]
