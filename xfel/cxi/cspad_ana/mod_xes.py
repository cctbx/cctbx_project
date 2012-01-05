# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$
"""First- and second-order statistics for CS-PAD images

The mod_xes user analysis module performs the following sequence of
analysis:  dark correction (using average dark previously computed
by mod_average); removal of inactive pixels; common mode correction;
removal of pixels with high stddev ("hot pixels"); removal of
noise < 5 ADUs; selection of scan rows known to contain the
spectrum; summation-reduction so as to form a single spectrum.

XXX mod_xes must be run as a single process--guard against it!
"""

__version__ = "$Revision$"

import math
import os

from libtbx import easy_pickle
import scitbx.math
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import average_tbx


class mod_xes(average_tbx.average_mixin):
  """Class for generating first- and second-order statistics within
  the pyana framework

  XXX Maybe this module should be renamed to mod_stat12, mod_sstat or
  some such?
  """


  def __init__(self,
               address,
               pickle_dirname=".",
               pickle_basename="",
               roi=None,
               **kwds):
    """The mod_average class constructor stores the parameters passed
    from the pyana configuration file in instance variables.  All
    parameters, except @p address are optional, and hence need not be
    defined in pyana.cfg.

    @param address         Address string XXX Que?!
    @param pickle_dirname     Directory portion of output pickle file
                           XXX mean, mu?
    @param pickle_basename    Filename prefix of output pickle file
                           image XXX mean, mu?
    @param calib_dir       Directory with calibration information
    @param dark_path       Path to input dark image
    @param dark_stddev     Path to input dark standard deviation
    @param flags           inactive:  Eliminate the inactive pixels
                           noelastic: Eliminate elastic scattering
                           nohot:     Eliminate the hot pixels.
                           nonoise:   Eliminate nosiy pixels.
    @param n               The number of shots to process, or as many
                           as possible if undefined XXX Sort of
                           redundant with pyana
    """
    super(mod_xes, self).__init__(
      address=address,
      **kwds
    )
    self.pickle_dirname = cspad_tbx.getOptString(pickle_dirname)
    self.pickle_basename = cspad_tbx.getOptString(pickle_basename)
    self.roi = cspad_tbx.getOptROI(roi)

  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_xes, self).event(evt, env)
    if (evt.get("skip_event")):
      return

    if self.roi is not None:
      pixels = self.cspad_img[self.roi[2]:self.roi[3], self.roi[0]:self.roi[1]]
      dark_mask = self.dark_mask[self.roi[2]:self.roi[3], self.roi[0]:self.roi[1]]
      pixels = pixels.as_1d().select(dark_mask.as_1d())
    else:
      pixels = self.cspad_img.as_1d().select(self.dark_mask.as_1d())
    stats = scitbx.math.basic_statistics(pixels.as_double())
    s, ms = cspad_tbx.evt_time(evt)
    evt_time = s + ms/1000
    self.stats_logger.info("SKEWNESS %.3f %s" %(evt_time, stats.skew))
    self.stats_logger.info("KURTOSIS %.3f %s" %(evt_time, stats.kurtosis))

    #if self.nmemb % 1000 == 0 or math.log(self.nmemb, 2) % 1 == 0:
      #self.endjob(env)

  def endjob(self, env):
    """The endjob() function finalises the mean and standard deviation
    images and writes them to disk.

    @param env Environment object
    """

    super(mod_xes, self).endjob(env)
    if (self.nmemb > 0):
      if (self.pickle_dirname  is not None or
          self.pickle_basename is not None):
        if (not os.path.isdir(self.pickle_dirname)):
          os.makedirs(self.pickle_dirname)
        d = dict(
          sum_img = self.sum_img,
          sumsq_img = self.sumsq_img,
          nmemb = self.nmemb,
          sifoil = self.sifoil,
        )
        pickle_path = os.path.join(self.pickle_dirname,
                                   self.pickle_basename+str(env.subprocess())+".pickle")
        easy_pickle.dump(pickle_path, d)
        self.logger.info(
          "Pickle written to %s" % self.pickle_dirname)

    if (self.nfail == 0):
      self.logger.info(
        "%d images processed" % self.nmemb)
    else:
      self.logger.warn(
        "%d images processed, %d failed" % (self.nmemb, self.nfail))
