from __future__ import division
from libtbx import easy_pickle
import logging
import math


class SingleFrame:
  """ Class that creates single-image agregate metrics/scoring that can then be
  used in downstream clustering or filtering procedures.
  """


  def __init__(self, path, filename, crystal_num=0):
    try:
      # Warn on error, but continue directory traversal.
      d = easy_pickle.load(path)
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
      self.unmerged_b = self.unmerged_b()
      logging.debug("Extracted image {}".format(filename))
    except KeyError:
      logging.warning("Could not extract point group and unit cell from %s" % path)
    except:
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


  def unmerged_b(self, nbins=20):
    _d_star_p = 1.618034  # Golden ratio distribution for d-spacings
    binner = self.miller_array.setup_binner(n_bins=nbins)
    #logging.debug(str("{}".format(binner.show_summary())))
    bin_selections = [binner.selection(i) for i in binner.range_used()]
    means = [self.miller_array.select(sel).mean() for sel in bin_selections]
    log_means = [math.log(mil) if mil > 0 else 0 for mil in means]
    centers = binner.bin_centers(_d_star_p)
    d_centers = centers ** (-1 / _d_star_p)
    return (d_centers, log_means)
    if logging.Logger.root.level < logging.DEBUG:  # Extreme debug!
      import matplotlib.pyplot as plt
      plt.plot(d_centers, log_means)
      plt.gca().invert_xaxis()
      plt.show()


