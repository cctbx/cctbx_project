# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$
"""First- and second-order statistics for CS-PAD images

The mod_average user analysis module computes the mean and the
standard deviation from the images in an XTC stream.  On successful
completion, the mean and standard deviation images are written to disk
as pickled dictionaries.
"""

__version__ = "$Revision$"

from xfel.cxi.cspad_ana import average_tbx
from xfel.cxi.cspad_ana import cspad_tbx


class mod_average(average_tbx.average_mixin):
  """Class for generating first- and second-order statistics within
  the pyana framework

  XXX Maybe this module should be renamed to mod_stat12, mod_sstat or
  some such?
  """


  def __init__(self,
               address,
               **kwds):
    """The mod_average class constructor stores the parameters passed
    from the pyana configuration file in instance variables.  All
    parameters, except @p address are optional, and hence need not be
    defined in pyana.cfg.

    @param address Address string XXX Que?!
    """

    super(mod_average, self).__init__(
      address=address,
      **kwds
    )


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_average, self).event(evt, env)
    if (evt.get("skip_event")):
      return

    self.logger.info("shot number %i" %self.nmemb)


  def endjob(self, env):
    """The endjob() function finalises the mean and standard deviation
    images and writes them to disk.  The distance and wavelength
    written to the standard deviation image is actually the average
    distance and wavelength, since standard deviations of those
    quantities do not make much sense in visualisation.

    @param env Environment object
    """

    super(mod_average, self).endjob(env)
    if (self.nmemb > 0):
      if (self.avg_dirname  is not None or
          self.avg_basename is not None):
        d = cspad_tbx.dpack(
          active_areas    = self.active_areas,
          beam_center_x   = cspad_tbx.pixel_size * self.beam_center[0],
          beam_center_y   = cspad_tbx.pixel_size * self.beam_center[1],
          data            = self.avg_img,
          distance        = self.avg_distance,
          wavelength      = self.avg_wavelength)
        p = cspad_tbx.dwritef(d, self.avg_dirname, self.avg_basename)
        self.logger.info(
          "Average written to %s" % p)

      if (self.stddev_dirname  is not None or
          self.stddev_basename is not None):
        d = cspad_tbx.dpack(
          active_areas    = self.active_areas,
          beam_center_x   = cspad_tbx.pixel_size * self.beam_center[0],
          beam_center_y   = cspad_tbx.pixel_size * self.beam_center[1],
          data            = self.stddev_img,
          distance        = self.avg_distance,
          wavelength      = self.avg_wavelength)
        p = cspad_tbx.dwritef(d, self.stddev_dirname, self.stddev_basename)
        self.logger.info(
          "Standard deviation written to %s" % p)

    if (self.nfail == 0):
      self.logger.info(
        "%d images processed" % self.nmemb)
    else:
      self.logger.warn(
        "%d images processed, %d failed" % (self.nmemb, self.nfail))
