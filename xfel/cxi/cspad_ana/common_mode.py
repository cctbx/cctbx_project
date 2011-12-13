# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$

"""Base class for dark subtraction and common-mode correction

XXX Better named cspad_base?
"""

__version__ = "$Revision$"

import cPickle as pickle
import cspad_tbx
import logging
import numpy

from parse_calib         import Section
from parse_calib         import calib2sections
from pypdsdata           import xtc
from scitbx.array_family import flex


class common_mode_correction(object):
  """Dark subtraction and alternate implementation of common mode
  substituting for cspad_tbx.

  Known problems: the algorithm relies on a substantial part of the
  sensor having no signal, which is not the case if a water ring
  crosses the sensor.
  """


  def __init__(self,
               address,
               calib_dir = None,
               common_mode_correction="none",
               photon_threshold=2,
               dark_path = None,
               dark_stddev = None):
    """The common_mode_correction class constructor stores the
    parameters passed from the pyana configuration file in instance
    variables.

    @param address         Address string XXX Que?!
    @param calib_dir       Directory with calibration information
    @param common_mode_correction The type of common mode correction to apply
    @param dark_path       Path to input average dark image
    @param dark_stddev     Path to input standard deviation dark
                           image, required if @p dark_path is given
    """

    self.logger = logging.getLogger(self.__class__.__name__)
    self.logger.setLevel(logging.INFO)

    # This is for messages that are picked up by Nat's monitoring program
    self.stats_logger = logging.getLogger("stats logger")
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    self.stats_logger.addHandler(handler)
    self.stats_logger.removeHandler(self.stats_logger.handlers[0])
    self.stats_logger.setLevel(logging.INFO)

    self.address = cspad_tbx.getOptString(address)
    self.dark_path = cspad_tbx.getOptString(dark_path)
    self.dark_stddev = cspad_tbx.getOptString(dark_stddev)
    self.common_mode_correction = cspad_tbx.getOptString(common_mode_correction)
    self.photon_threshold = cspad_tbx.getOptFloat(photon_threshold)

    self.distance = None
    self.cspad_img = None # The current image - set by self.event()
    self.sifoil = None
    self.sum_common_mode = 0
    self.sumsq_common_mode = 0
    self.wavelength = None # The current wavelength - set by self.event()

    assert self.common_mode_correction in \
        ("gaussian", "mean", "median", "mode", "none")

    # Get and parse metrology.  There is no metrology information for
    # the Sc1 detector, so make it up.
    self.sections = None
    if (self.address == "CxiDs1-0|Cspad-0" and calib_dir is not None):
      self.sections = calib2sections(cspad_tbx.getOptString(calib_dir))
    elif (self.address == "CxiSc1-0|Cspad2x2-0"):
      # The sections are rotated by 90 degrees with respect to the
      # "standing up" convention.
      self.sections = [[Section(90, (185 / 2 + 0,   (2 * 194 + 3) / 2)),
                        Section(90, (185 / 2 + 185, (2 * 194 + 3) / 2))]]
    if (self.sections is None):
      self.logger.error("Failed to load metrology")

    # Load the dark image and ensure it is signed and at least 32 bits
    # wide, since it will be used for differencing.  If a dark image
    # is provided, a standard deviation image is required, and all the
    # ADU scales must match up.
    self.dark_img = None
    if (dark_path is not None):
      assert dark_stddev is not None
      try:
        stream = open(dark_path, "rb")
        dark_dict = pickle.load(stream)
        #assert "ADU_SCALE" not in dark_dict # force use of recalculated dark
        self.dark_img = dark_dict['DATA']
        assert isinstance(self.dark_img, flex.double)
        stream.close()

        stream = open(dark_stddev, "rb")
        dark_stddev_dict = pickle.load(stream)
        self.dark_stddev = dark_stddev_dict['DATA']
        assert isinstance(self.dark_stddev, flex.double)
        stream.close()
      except IOError:
        self.logger.error("Failed to load dark image")
        raise
      self.dark_mask = (self.dark_stddev > 0)


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    # XXX Not needed now that the distance is read in the event?
    env.update(evt)

    self.config = env.getConfig(xtc.TypeId.Type.Id_CspadConfig, self.address)
    if (self.config is None):
      self.logger.error("beginjob(): no config")

    (self.beam_center, self.active_areas) = \
        cspad_tbx.cbcaa(self.config, self.sections)

    self.nfail  = 0
    self.nshots = 0
    self.nmemb = 0


  def common_mode(self, img, stddev, mask):
    """The common_mode() function returns the mode of image stored in
    the array pointed to by @p img.  @p mask must be such that the @p
    stddev at the selected pixels is greater than zero.

    @param img    2D integer array of the image
    @param stddev 2D integer array of the standard deviation of each
                  pixel in @p img
    @param mask   2D Boolean array, @c True if the pixel is to be
                  included, @c False otherwise
    @return       Mode of the image, as a real number
    """

    # Flatten the image and take out inactive pixels XXX because we
    # cannot take means and medians of 2D arrays?
    img_1d = img.as_1d().select(mask.as_1d()).as_double()
    assert img_1d.size() > 0

    if (self.common_mode_correction == "mean"):
      # The common mode is approximated by the mean of the pixels with
      # signal-to-noise ratio less than a given threshold.  XXX Breaks
      # if the selection is empty!
      THRESHOLD_SNR = 2
      img_snr = img_1d / stddev.as_double().as_1d().select(mask.as_1d())
      return (flex.mean(img_1d.select(img_snr < THRESHOLD_SNR)))

    elif (self.common_mode_correction == "median"):
      return (flex.median(img_1d))

    # Identify the common-mode correction as the peak histogram of the
    # histogram of pixel values (the "standard" common-mode correction, as
    # previously implemented in this class).
    hist_min = -40
    hist_max = 40
    n_slots = 100

    hist = flex.histogram(img_1d, hist_min, hist_max, n_slots=n_slots)
    slots = hist.slots()
    i = flex.max_index(slots)
    common_mode = list(hist.slot_infos())[i].center()

    if (self.common_mode_correction == "mode"):
      return (common_mode)

    # Determine the common-mode correction from the peak of a single
    # Gaussian function fitted to the histogram.
    from scitbx.math.curve_fitting import single_gaussian_fit
    x = hist.slot_centers()
    y = slots.as_double()
    fit = single_gaussian_fit(x, y)
    scale, mu, sigma = fit.scale, fit.mu, fit.sigma
    self.logger.debug("fitted gaussian: mu=%.3f, sigma=%.3f" %(mu, sigma))
    mode = common_mode
    common_mode = mu
    if abs(mode-common_mode) > 1000: common_mode = mode # XXX
    self.logger.debug("delta common mode corrections: %.3f" %(mode-common_mode))

    if 0 and abs(mode-common_mode) > 0:
      #if 0 and skew > 0.5:
      # view histogram and fitted gaussian
      from numpy import exp
      from matplotlib import pyplot
      x_all = x
      n, bins, patches = pyplot.hist(section_img.as_1d().as_numpy_array(), bins=n_slots, range=(hist_min, hist_max))
      y_all = scale * flex.exp(-flex.pow2(x_all-mu) / (2 * sigma**2))
      scale = slots[flex.max_index(slots)]
      y_all *= scale/flex.max(y_all)
      pyplot.plot(x_all, y_all)
      pyplot.show()

    return (common_mode)


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    # Increase the event counter, even if this event is to be skipped.
    self.nshots += 1
    if (evt.get("skip_event")):
      return

    distance = cspad_tbx.env_distance(env)
    if (distance is None):
      self.nfail += 1
      self.logger.warn("event(): no distance, shot skipped")
      evt.put(True, "skip_event")
      return
    if (self.distance is not None and self.distance != distance):
      self.logger.warn("event(): distance changed mid-run: % 8.4f -> % 8.4f" %
        (self.distance, distance))
    self.distance = distance

    sifoil = cspad_tbx.env_sifoil(env)
    if (sifoil is None):
      self.nfail += 1
      self.logger.warn("event(): no Si-foil thickness, shot skipped")
      evt.put(True, "skip_event")
      return
    if (self.sifoil is not None and self.sifoil != sifoil):
      self.logger.warn("event(): Si-foil changed mid-run: % 8i -> % 8d" %
        (self.sifoil, sifoil))
    self.sifoil = sifoil

    self.timestamp = cspad_tbx.evt_timestamp(evt)
    if (self.timestamp is None):
      self.nfail += 1
      self.logger.warn("event(): no timestamp, shot skipped")
      evt.put(True, "skip_event")
      return

    self.wavelength = cspad_tbx.evt_wavelength(evt)
    if (self.wavelength is None):
      self.nfail += 1
      self.logger.warn("event(): no wavelength, shot skipped")
      evt.put(True, "skip_event")
      return

    # Early return if the image is already stored in the event.
    # Otherwise, get it from the stream.  XXX It is probably not safe
    # to key the image on self.address, so we should come up with our
    # own namespace.
    self.cspad_img = evt.get(self.address)
    if (self.cspad_img is not None):
      return
    self.cspad_img = cspad_tbx.image(
      self.address, self.config, evt, env, self.sections)
    if (self.cspad_img is None):
      self.nfail += 1
      self.logger.error("event(): no image, shot skipped")
      evt.put(True, "skip_event")
      return

    # Make a signed 32-bit integer flex array of the full detector
    # image, and scale the ADU values appropriately.
    self.cspad_img = flex.int(self.cspad_img.astype(numpy.int32))
    self.cspad_img = self.cspad_img.as_double()

    # If a dark image was provided, subtract it from the image.  There
    # is no point in doing common-mode correction unless the dark
    # image was subtracted.
    if (self.dark_img is not None):
      self.cspad_img -= self.dark_img

      if (self.common_mode_correction != "none"):
        # Mask out inactive pixels prior to common mode correction.
        # Pixels are marked as inactive either due to low ADU values
        # or non-positive standard deviations in dark image.  XXX Make
        # the threshold tunable?
        cspad_mask = self.dark_mask.deep_copy()

        # Extract each active section from the assembled detector
        # image and apply the common mode correction.
        q_mask = self.config.quadMask()
        for q in xrange(len(self.sections)):
          if (not((1 << q) & q_mask)):
            continue

          s_mask = self.config.roiMask(q)
          for s in xrange(len(self.sections[q])):
            # XXX DAQ misconfiguration?  This mask appears not to work
            # reliably for the Sc1 detector.
#            if (not((1 << s) & s_mask)):
#              continue
            corners   = self.sections[q][s].corners()
            i_row     = int(round(min(c[0] for c in corners)))
            i_column  = int(round(min(c[1] for c in corners)))
            n_rows    = int(round(max(c[0] for c in corners))) - i_row
            n_columns = int(round(max(c[1] for c in corners))) - i_column

            section_img    = self.cspad_img.matrix_copy_block(
              i_row  = i_row,  i_column  = i_column,
              n_rows = n_rows, n_columns = n_columns)
            section_mask   = cspad_mask.matrix_copy_block(
              i_row  = i_row,  i_column  = i_column,
              n_rows = n_rows, n_columns = n_columns)
            section_stddev = self.dark_stddev.matrix_copy_block(
              i_row  = i_row,  i_column  = i_column,
              n_rows = n_rows, n_columns = n_columns)

            common_mode = self.common_mode(
              section_img, section_stddev, section_mask)
            self.sum_common_mode += common_mode
            self.sumsq_common_mode += common_mode**2

            # Apply the rounded integer common mode correction to the
            # section, and paste it back into the image.
            self.cspad_img.matrix_paste_block_in_place(
              block    = section_img - common_mode,
              i_row    = i_row,
              i_column = i_column)

    # Store the image in the event.
    evt.put(self.cspad_img, self.address)


  def endjob(self, env):
    """
    @param env Environment object
    """

    if 0 and self.dark_path is not None and self.nmemb > 1:
      print self.sum_common_mode, self.sumsq_common_mode
      self.mean_common_mode = self.sum_common_mode / self.nmemb
      print self.mean_common_mode
      self.stddev_commond_mode = math.sqrt((self.sumsq_common_mode
        - self.sum_common_mode * self.mean_common_mode) / (self.nmemb - 1))

      self.logger.info("mean common mode: %.3f" %self.mean_common_mode)
      self.logger.info("std. dev. common mode: %.3f" %self.stddev_commond_mode)

  def do_sigma_scaling(self):
    # Divide each pixel value by it's dark standard deviation. Since we are led
    # to believe that the standard deviation of a pixel is proportional to the
    # gain of said pixel, this approximates a gain correction.
    assert self.dark_img is not None
    flex_cspad_img = self.cspad_img.as_double()
    flex_cspad_img_sel = flex_cspad_img.as_1d().select(self.dark_mask.as_1d())
    flex_dark_stddev = self.dark_stddev.select(self.dark_mask.as_1d()).as_double()
    assert flex_dark_stddev.count(0) == 0
    flex_cspad_img_sel /= flex_dark_stddev
    flex_cspad_img.as_1d().set_selected(self.dark_mask.as_1d().iselection(), flex_cspad_img_sel)
    self.cspad_img = flex_cspad_img
    if 0: # for debugging
      from matplotlib import pyplot
      hist_min, hist_max = flex.min(flex_cspad_img_sel.as_double()), flex.max(flex_cspad_img_sel.as_double())
      print hist_min, hist_max
      n_slots = 100
      n, bins, patches = pyplot.hist(flex_cspad_img_sel.as_1d().as_numpy_array(), bins=n_slots, range=(hist_min, hist_max))
      pyplot.show()

  def do_photon_counting(self):
    # This only makes sense in combination with some sort of gain correction
    # XXX TODO: count 2, 3, 4, ..., photons
    PHOTON_THRESHOLD = self.photon_threshold
    self.cspad_img.set_selected(self.cspad_img<PHOTON_THRESHOLD, 0)
    self.cspad_img.set_selected(self.cspad_img>=PHOTON_THRESHOLD, 1)
