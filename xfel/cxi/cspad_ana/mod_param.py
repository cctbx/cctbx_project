# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$
"""Simple analysis of per-shot parameters of a run

XXX mod_param must be run as a single process--guard against it!
"""

__version__ = "$Revision$"

import logging
import math

from xfel.cxi.cspad_ana import cspad_tbx


class mod_param(object):
  """Class for simple analysis of certain per-shot parameters of a
  run.  Currently attenuation, sample-detector distance, and
  wavelength are considered.  Since attenuation and distance is not
  expected to vary, their values are output only when changed.  For
  the wavelength, the average and standard deviation is computed.
  Note that incompleteness due to lacking image data is not detected,
  because it may slow down this module considerably.
  """


  def __init__(self):
    """The mod_param class constructor initialises the logger.
    """

    self.m_logger = logging.getLogger(__name__)
    self.m_logger.setLevel(logging.INFO)


  def __del__(self):
    logging.shutdown()


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    self.m_nfail = 0
    self.m_nmemb = 0

    self.m_no_detz       = 0
    self.m_no_sifoil     = 0
    self.m_no_wavelength = 0

    self.m_detz       = None
    self.m_detz_max   = None
    self.m_detz_min   = None
    self.m_sifoil     = None
    self.m_sifoil_max = None
    self.m_sifoil_min = None


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.
    The event() function does not log shots skipped due to
    incompleteness in order to keep the output streams clean.
    Instead, the number of skipped shots is reported by endjob().

    @param evt Event data object, a configure object
    @param env Environment object
    """

    if (evt.get("skip_event")):
      return

    detz = cspad_tbx.env_detz(env)
    if (detz is None):
      self.m_no_detz += 1

    sifoil = cspad_tbx.env_sifoil(env)
    if (sifoil is None):
      self.m_no_sifoil += 1

    timestamp = cspad_tbx.evt_timestamp(evt)
    if (timestamp is None):
      self.m_no_timestamp += 1

    wavelength = cspad_tbx.evt_wavelength(evt)
    if (wavelength is None):
      self.m_no_wavelength += 1

    if (detz       is None or
        sifoil     is None or
        timestamp  is None or
        wavelength is None):
      self.m_nfail += 1
      return

    if (detz != self.m_detz):
      if (self.m_detz is None):
        self.m_logger.info("%s: initial detz     % 8.4f" %
                           (timestamp, detz))
      else:
        self.m_logger.info("%s: detz     % 8.4f -> % 8.4f" %
                           (timestamp, self.m_detz, detz))

      self.m_detz = detz
      if (self.m_detz_max is None or detz > self.m_detz_max):
        self.m_detz_max = detz
      if (self.m_detz_min is None or detz < self.m_detz_min):
        self.m_detz_min = detz

    if (sifoil != self.m_sifoil):
      if (self.m_sifoil is None):
        self.m_logger.info("%s: initial Si-foil  % 8d" %
                           (timestamp, sifoil))
      else:
        self.m_logger.info("%s: Si-foil  % 8d -> % 8d" %
                           (timestamp, self.m_sifoil, sifoil))

      self.m_sifoil = sifoil
      if (self.m_sifoil_max is None or sifoil > self.m_sifoil_max):
        self.m_sifoil_max = sifoil
      if (self.m_sifoil_min is None or sifoil < self.m_sifoil_min):
        self.m_sifoil_min = sifoil

    # Accumulate the sum and the squared sum of the shifted the
    # wavelength.  The shift is taken as the first wavelength
    # encountered.  This may be more accurate than accumulating raw
    # values [Chan et al. (1983) Am. Stat. 37, 242-247].
    if (self.m_nmemb == 0):
      self.m_wavelength_shift  = wavelength
      self.m_wavelength_sum    = (wavelength - self.m_wavelength_shift)
      self.m_wavelength_sumsq  = (wavelength - self.m_wavelength_shift)**2
      self.m_nmemb             = 1
    else:
      self.m_wavelength_sum   += (wavelength - self.m_wavelength_shift)
      self.m_wavelength_sumsq += (wavelength - self.m_wavelength_shift)**2
      self.m_nmemb            += 1


  def endjob(self, env):
    """The endjob() function finalises the mean and standard deviation
    calculations, and reports on the total number of skipped shots.

    @param env Environment object
    """

    if (self.m_nmemb > 0):
      # The shift has to be added back to the average wavelength, but
      # not the standard deviation.
      avg_wavelength     = self.m_wavelength_sum / self.m_nmemb
      stddev_wavelength  = self.m_wavelength_sumsq \
          -                self.m_wavelength_sum * avg_wavelength
      avg_wavelength    += self.m_wavelength_shift

      if (stddev_wavelength < 0):
        stddev_wavelength = 0
      elif (self.m_nmemb > 1):
        stddev_wavelength = math.sqrt(stddev_wavelength / (self.m_nmemb - 1))
      else:
        stddev_wavelength = math.sqrt(stddev_wavelength / (self.m_nmemb - 0))

      self.m_logger.info("Det-z:      min % 12.6f, max   % 12.6f" %
                         (self.m_detz_min, self.m_detz_max))
      self.m_logger.info("Si-foil:    min % 12d, max   % 12d" %
                         (self.m_sifoil_min, self.m_sifoil_max))
      self.m_logger.info("Wavelength: mu  % 12.6f, sigma % 12.6f" %
                         (avg_wavelength, stddev_wavelength))

    self.m_logger.info("%5d images processed, %5d images skipped" %
                       (self.m_nmemb, self.m_nfail))
    self.m_logger.info("No det-z:       %5d" % self.m_no_detz)
    self.m_logger.info("No attenuation: %5d" % self.m_no_sifoil)
    self.m_logger.info("No wavelength:  %5d" % self.m_no_wavelength)
