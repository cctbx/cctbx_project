""" Module for working with single images in a serial crystallography
dataset"""
from __future__ import absolute_import, division, print_function
from libtbx import easy_pickle
import numpy as np
import math
import logging
from cctbx.array_family import flex
from six.moves import cPickle as pickle
from .api import InputFrame
logger = logging.getLogger('sf')

class SingleFrame(InputFrame):
  """ Class that creates single-image agregate metrics/scoring that can then be
  used in downstream clustering or filtering procedures.
  """
  ANGSTROMS_TO_EV = 12398.425

  def __init__(self, path=None, filename=None, crystal_num=0,
               remove_negative=False, use_b=True, scale=True, dicti=None,
               pixel_size=None):
    """
    Constructor for SingleFrame object, using a cctbx.xfel integration pickle.

    :param path: path to integration pickle
    :param filename: the file name alone (used as a label)
    :param crystal_num: if multiple lattices present, the latice number.
    :param remove_negative: Boolean for removal of negative intensities
    :param use_b: if True, initialise scale and B, if false, use only mean-intensity scaling.
    :param dicti: optional. If a dictionairy is supplied here, will create object from that rather than attempting to read the file specified in path, filename.
    :param pixel_size: the size of pixels in mm. Defaults to a MAR detector with a warning at debug level of logging.
    :param scale: if False, will intialise scales to G=1, B=0.


    :return: a SingleFrame object, with the following Object attributes:


    Object attributes are:
        - `is_polarization_corrected`: Boolean flag indicatinf if polarization correction has been applied
        - `miller_array`: the cctbx.miller miller array of spot intensities.
        - `mapped_predictions`: the mapped_predictions locations
        - `path`: full path to the original file
        - `name`: file-name, used as an identifier
        - `crystal_system:
        - `pg`: point group of pickle
        - `uc`: Niggli unit cell as a tuple
        - `orientation`: cctbx crystal_orientation object
        - `total_i`: the total integrated intensity for this frame
        - `xbeam`: x-location of beam centre
        - `ybeam`: y-location of beam centre
        - `wavelength:
        - `spot_offset`: the mean offset between observed spots and predicted centroids. Only created if integration was performed using verbose_cv=True. Otherwise None.
        - `minus_2B`: the gradient of the ln(i) vs. sinsqtheta_over_lambda_sq plot
        - `G`: intercept of the of the ln(i) vs. sinsqtheta_over_lambda_sq plot
        - `log_i`: list of log_i intensities
        - `sinsqtheta_over_lambda_sq`: list of sinsqtheta_over_lambda_sq
        - `wilson_err`: standard error on the fit of ln(i) vs. sinsqtheta_over_lambda_sq
        - `miller_fullies`: a cctbx.miller array of fully recorded intensites.
    """
    if dicti is not None:
      d = dicti
    else:
      try:
        d = easy_pickle.load(path)
      except (pickle.UnpicklingError, ValueError, EOFError, IOError):
        d = {}
        logger.warning("Could not read %s. It may not be a pickle file." % path)
    if 'observations' not in d or len(d['observations'][crystal_num].data()) == 0:
      return
    try:
      if pixel_size:
        self.pixel_size = pixel_size
      else:
        logger.debug("No pixel size specified, defaulting to MAR (0.079346). "
                        "Bad times if this is not the correct detector!")
        self.pixel_size = 0.079346
      # Warn on error, but continue directory traversal.
      self.is_polarization_corrected = False
      # Miller arrays
      self.miller_array = d['observations'][crystal_num]
      self.mapped_predictions = d['mapped_predictions'][crystal_num]
      # Image pickle info
      self.path = path or d['path']
      self.name = filename
      # Unit cell info
      self.crystal_system = self.miller_array.crystal_symmetry()\
        .space_group().crystal_system()
      self.pg = d['pointgroup'].replace(' ', '')  # enforce consistency
      # XXX major bug here??? niggli cell not produced with knowledge of the centring symbol???
      self.uc = d['current_orientation'][crystal_num].unit_cell() \
        .niggli_cell() \
        .parameters()
      self.orientation = d['current_orientation'][crystal_num]
      # Agregate info
      self.total_i = d['observations'][crystal_num].sum()
      self.xbeam = d['xbeam']
      self.ybeam = d['ybeam']
      self.wavelength = d['wavelength']
      self.distance = d['distance']
      if 'correction_vectors' in d:
        all_corrections = []
        for spot in d['correction_vectors'][crystal_num]:
          dta = np.sqrt((spot['refinedcenter'][0] - spot['obscenter'][0]) ** 2
                      + (spot['refinedcenter'][1] - spot['obscenter'][1]) ** 2)
          all_corrections.append(dta)
        self.spot_offset = np.mean(all_corrections)
      else:
        self.spot_offset = None

      if remove_negative:
        self.filter_negative_intensities()

      # Do polarization correction
      self.polarization_correction()
      self.minus_2B, self.G, self.log_i, \
          self.sinsqtheta_over_lambda_sq, \
          self.wilson_err = self.init_calc_wilson(use_b)
      if not scale:
        self.minus_2B = 0
        self.G = 1
      if logger.root.level < logging.DEBUG:  # Extreme debug!
        self.plot_wilson()
      logger.debug("Extracted image {}".format(filename))
    except KeyError:
      logger.warning("Could not extract point group and unit cell from %s" % path)

    self.miller_fullies = None

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

  def n_reflections_by_sigi(self, sig_i_cuttoff):
    """
    Currently a placeholder that returns None.

    This method should return the number of reflection in the frame that have an
    I/sig(I) > sig_i_cuttoff
    """
    reflections_above_cuttoff = None
    return len(reflections_above_cuttoff)

  def init_calc_wilson(self, use_b_factor, i_corrections=None):
    """ If use_b_factor is
    :param i_corrections: allows flex array of correction factors (e.g. partialities) to be specified
    :param use_b_factor: if True, do a linear regression to fit G and B and returns the coeficients minus_2B, G, the transformed data log_i, and one_over_d_sqare. Also returns fit_stats, which is a dictionairy. If use_b_factor is False, then B is 0, and G is the mean intensity of the image. The r_value is then 0 (by definition), and the std_err is the standard error on the mean.

    :return minus_2B, G, log_i, on_over_d_square: `minus_2B`: gradient of fit; `G`: intercept of fit; `log_i`: dependent variable of fit; `one_over_d_square`: independent variable of fit.
    """
    if i_corrections:
      inten = (self.miller_array.sort().data() * i_corrections).as_numpy_array()
    else:
      inten = self.miller_array.sort().data().as_numpy_array()
    sinsqtheta_over_labmdasq = self.miller_array.sort()\
      .sin_theta_over_lambda_sq().data().as_numpy_array()

     # then plot them as negative in the linear fit.
    inten, sinsqtheta_over_labmdasq = zip(*[i for i
                                            in zip(inten,
                                                   sinsqtheta_over_labmdasq)
                                            if i[0] >= 0])

    if use_b_factor:
      from scipy.stats import linregress
      minus_2B, G, r_val, _, std_err = linregress(sinsqtheta_over_labmdasq,
                                                  np.log(inten))
    else:
      # If the model is a constant value, r_val = 0, and
      from scipy.stats import sem
      minus_2B, G, r_val, std_err = 0, np.mean(inten), 0, sem(inten)

    # ignore p_val since this will be insanely small
    logger.debug("G: {}, -2B: {}, r: {}, std_err: {}".
      format(G, minus_2B, r_val, std_err))
    return minus_2B, G, np.log(inten), sinsqtheta_over_labmdasq, {"R": r_val,
                                                   "Standard Error": std_err}

  def plot_wilson(self, width=30, ax=None):
    """ Makes a log(I) vs 1/d**2 plot, displaying the raw partial data, a
    rolling average of the data, and the Wilson model fit to the data.

    :param: width: smoothing window size
    :param: ax: optional axes object to ve used for plotting
    """

    import matplotlib.pyplot as plt
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
     #logger.debug(str("{}".format(binner.show_summary())))
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
    Calculates distance using NCDist from Andrews and Bernstein J. Appl.
    Cryst. 2014 between this frame and some other unit cell.
    :param:other_uc: a 6-tuple of a, b, c, alpha, beta, gamma for some unit cell
    :return: the NCDist in A^2 to other_uc
    """
    from cctbx.uctbx.determine_unit_cell import NCDist
    self_g6 = self.make_g6(self.uc)
    other_g6 = self.make_g6(other_uc)
    return NCDist(self_g6, other_g6)

  def to_panda(self):
    """ Returns the object attributes as a pandas series """
    import pandas as pd
    return pd.Series({'path': self.path,
                      'name': self.name,
                      'crystal_system': self.crystal_system,
                      'point group': self.pg,
                      'a': self.uc[0],
                      'b': self.uc[1],
                      'c': self.uc[2],
                      'alpha': self.uc[3],
                      'beta': self.uc[4],
                      'gamma': self.uc[5],
                      'total_i': self.total_i,
                      'wavelength': self.wavelength,
                      'spot_offset': self.spot_offset,
                      'minus_2B': self.minus_2B,
                      'G': self.G,
                      'willson_err': self.wilson_err})


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

class SingleDialsFrame(SingleFrame):
  def __init__(self, refl=None, expt=None, id=None, **kwargs):
    from xfel.command_line.frame_extractor import ConstructFrame
    frame = ConstructFrame(refl, expt).make_frame()
    SingleFrame.__init__(self, dicti=frame, path=str(id), **kwargs)
    self.experiment = expt
    self.reflections = refl

class SingleDialsFrameFromFiles(SingleFrame):
  def __init__(self, refls_path=None, expts_path=None, **kwargs):
    from xfel.command_line.frame_extractor import ConstructFrameFromFiles
    frame = ConstructFrameFromFiles(refls_path, expts_path).make_frame()
    SingleFrame.__init__(self, dicti=frame, path=" ".join((refls_path, expts_path)), **kwargs)

class CellOnlyFrame(SingleFrame):
  def __init__(self, crystal_symmetry, path=None, name=None, lattice_id=None):
    from six.moves import cStringIO as StringIO
    f = StringIO()
    self.crystal_symmetry = crystal_symmetry
    self.crystal_symmetry.show_summary(f=f)
    self.niggli_cell = self.crystal_symmetry.niggli_cell()
    self.niggli_cell.show_summary(f=f, prefix="   niggli-->")
    logger.info(f.getvalue())
    self.uc = self.niggli_cell.unit_cell().parameters()
    self.mm = self.niggli_cell.unit_cell().metrical_matrix()
    self.pg = "".join(self.crystal_symmetry.space_group().type().lookup_symbol().split())
    self.path = path
    self.name = name
    self.lattice_id = lattice_id

class SingleDialsFrameFromJson(SingleFrame):
  def __init__(self, expts_path=None, **kwargs):
    from dials.util.options import Importer, flatten_experiments
    importer = Importer([expts_path], read_experiments=True, read_reflections=False, check_format=False)
    if importer.unhandled:
      # in python 2: raise Exception("unable to process:"), importer.unhandled
      raise Exception("unable to process:")
    experiments_l = flatten_experiments(importer.experiments)
    assert len(experiments_l)==1, "Sorry, only supports one experiment per json at present."
    tcrystal = experiments_l[0].crystal
    from cctbx import crystal
    group = tcrystal.get_space_group()
    self.crystal_symmetry = crystal.symmetry(unit_cell=tcrystal.get_unit_cell(),
                                             space_group=group)
    self.crystal_symmetry.show_summary()
    self.niggli_cell = self.crystal_symmetry.niggli_cell()
    self.niggli_cell.show_summary(prefix="   niggli-->")
    self.uc = self.niggli_cell.unit_cell().parameters()
    self.mm = self.niggli_cell.unit_cell().metrical_matrix()
    self.pg = "".join(group.type().lookup_symbol().split())
    self.path = expts_path
