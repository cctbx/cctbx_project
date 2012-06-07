# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$

"""Base class for dark subtraction and common-mode correction

XXX Better named cspad_base?
"""

__version__ = "$Revision$"

import numpy

from parse_calib         import Section
from parse_calib         import calib2sections
from pypdsdata           import xtc

from libtbx import easy_pickle
from scitbx.array_family import flex
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana.mod_event_info import mod_event_info


class common_mode_correction(mod_event_info):
  """Dark subtraction and alternate implementation of common mode
  substituting for cspad_tbx.

  Known problems: the algorithm relies on a substantial part of the
  sensor having no signal, which is not the case if a water ring
  crosses the sensor.
  """


  def __init__(self,
               address,
               calib_dir=None,
               common_mode_correction="none",
               photon_threshold=None,
               two_photon_threshold=None,
               dark_path=None,
               dark_stddev=None,
               gain_map_path=None,
               cache_image=True,
               roi=None,
               laser_1_status=None,
               laser_4_status=None,
               laser_wait_time=None,
               **kwds):
    """The common_mode_correction class constructor stores the
    parameters passed from the pyana configuration file in instance
    variables.

    @param address         Address string XXX Que?!
    @param calib_dir       Directory with calibration information
    @param common_mode_correction The type of common mode correction to apply
    @param dark_path       Path to input average dark image
    @param dark_stddev     Path to input standard deviation dark
                           image, required if @p dark_path is given
    @param laser_1_status  0 or 1 to indicate that the laser should be off or on respectively
    @param laser_4_status  0 or 1 to indicate that the laser should be off or on respectively
    @param laser_wait_time Length of time in milliseconds to wait after a laser
                           change of status to begin accepting images again.
                           (rejection of images occurs immediately after status
                           change).
    """

    # Cannot use the super().__init__() construct here, because
    # common_mode_correction refers to the argument, and not the
    # class.
    mod_event_info.__init__(self, address=address, **kwds)

    self.dark_path = cspad_tbx.getOptString(dark_path)
    self.dark_stddev = cspad_tbx.getOptString(dark_stddev)
    gain_map_path = cspad_tbx.getOptString(gain_map_path)
    self.common_mode_correction = cspad_tbx.getOptString(common_mode_correction)
    self.photon_threshold = cspad_tbx.getOptFloat(photon_threshold)
    self.two_photon_threshold = cspad_tbx.getOptFloat(two_photon_threshold)
    self.cache_image = cspad_tbx.getOptBool(cache_image)
    self.filter_laser_1_status = cspad_tbx.getOptInteger(laser_1_status)
    self.filter_laser_4_status = cspad_tbx.getOptInteger(laser_4_status)
    if self.filter_laser_1_status is not None:
      self.filter_laser_1_status = bool(self.filter_laser_1_status)
    if self.filter_laser_4_status is not None:
      self.filter_laser_4_status = bool(self.filter_laser_4_status)
    self.filter_laser_wait_time = cspad_tbx.getOptInteger(laser_wait_time)

    self.cspad_img = None # The current image - set by self.event()
    self.sum_common_mode = 0
    self.sumsq_common_mode = 0
    self.roi = cspad_tbx.getOptROI(roi) # used to ignore the signal region in chebyshev fit

    assert self.common_mode_correction in \
        ("gaussian", "mean", "median", "mode", "none", "chebyshev")

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
      raise RuntimeError("Failed to load metrology")

    # Load the dark image and ensure it is signed and at least 32 bits
    # wide, since it will be used for differencing.  If a dark image
    # is provided, a standard deviation image is required, and all the
    # ADU scales must match up.
    self.dark_img = None
    if (dark_path is not None):
      assert dark_stddev is not None
      dark_dict = easy_pickle.load(dark_path)
      #assert "ADU_SCALE" not in dark_dict # force use of recalculated dark
      self.dark_img = dark_dict['DATA']
      assert isinstance(self.dark_img, flex.double)

      dark_stddev_dict = easy_pickle.load(dark_stddev)
      self.dark_stddev = dark_stddev_dict['DATA']
      assert isinstance(self.dark_stddev, flex.double)
      self.dark_mask = (self.dark_stddev > 0)

    self.gain_map = None
    if gain_map_path is not None:
      self.gain_map
      gain_dict = easy_pickle.load(gain_map_path)
      self.gain_map = gain_dict['DATA']
      assert isinstance(self.gain_map, flex.double)


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    # XXX Not needed now that the distance is read in the event?
    env.update(evt)

    # XXX Quick hack to get the config object for the detector at
    # self.address.
    if self.address == 'CxiDs1-0|Cspad-0':
      self.config = env.getConfig(
        xtc.TypeId.Type.Id_CspadConfig, self.address)
    elif self.address == 'CxiSc1-0|Cspad2x2-0':
      self.config = env.getConfig(
        xtc.TypeId.Type.Id_Cspad2x2Config, self.address)

    if self.config is None:
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
    scale, mu, sigma = fit.a, fit.b, fit.c
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

    super(common_mode_correction, self).event(evt, env)
    if (evt.get("skip_event")):
      return

    if self.filter_laser_1_status is not None:
      if (self.laser_1_status.status != self.filter_laser_1_status or
          (self.laser_1_ms_since_change is not None and
           self.laser_1_ms_since_change < self.filter_laser_wait_time)):
        evt.put(True, "skip_event")
        return
    if self.filter_laser_4_status is not None:
      if (self.laser_4_status.status != self.filter_laser_4_status or
          (self.laser_4_ms_since_change is not None and
           self.laser_4_ms_since_change < self.filter_laser_wait_time)):
        evt.put(True, "skip_event")
        return

    # Early return if the full detector image is already stored in the
    # event.  Otherwise, get it from the stream as a double-precision
    # floating-point flex array.  XXX It is probably not safe to key
    # the image on self.address, so we should come up with our own
    # namespace.
    self.cspad_img = evt.get(self.address)
    if (self.cspad_img is not None):
      return
    self.cspad_img = cspad_tbx.image(
      self.address, self.config, evt, env, self.sections)
    if (self.cspad_img is None):
      self.nfail += 1
      self.logger.warn("event(): no image, shot skipped")
      evt.put(True, "skip_event")
      return
    self.cspad_img = flex.double(self.cspad_img.astype(numpy.float64))

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

        if self.roi is not None and self.common_mode_correction == "chebyshev":
          roi_mask = cspad_mask[self.roi[2]:self.roi[3], :]
          roi_mask = flex.bool(roi_mask.accessor(), False)
          cspad_mask.matrix_paste_block_in_place(
            block=roi_mask,
            i_row=self.roi[2],
            i_column=0)

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

            if section_mask.count(True) == 0: continue

            if self.common_mode_correction == "chebyshev":
              assert len(self.sections[q]) == 2
              if s == 0:
                section_imgs = [section_img]
                section_masks = [section_mask]
                i_rows = [i_row]
                i_columns = [i_column]
                continue
              else:
                section_imgs.append(section_img)
                section_masks.append(section_mask)
                i_rows.append(i_row)
                i_columns.append(i_column)

                chebyshev_corrected_imgs = self.chebyshev_common_mode(
                  section_imgs, section_masks)
                for i in range(2):
                  section_imgs[i].as_1d().copy_selected(
                    section_masks[i].as_1d().iselection(),
                    chebyshev_corrected_imgs[i].as_1d())
                  self.cspad_img.matrix_paste_block_in_place(
                    block=section_imgs[i],
                    i_row=i_rows[i],
                    i_column=i_columns[i])

            else:
              common_mode = self.common_mode(
                section_img, section_stddev, section_mask)
              self.sum_common_mode += common_mode
              self.sumsq_common_mode += common_mode**2

              # Apply the common mode correction to the
              # section, and paste it back into the image.
              self.cspad_img.matrix_paste_block_in_place(
                block    = section_img - common_mode,
                i_row    = i_row,
                i_column = i_column)

    if self.gain_map is not None:
      self.cspad_img *= self.gain_map

    if self.cache_image:
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
    assert self.gain_map is None # not appropriate to do sigma scaling and gain correction at the same time!
    flex_cspad_img = self.cspad_img.as_double()
    flex_cspad_img_sel = flex_cspad_img.as_1d().select(self.dark_mask.as_1d())
    flex_dark_stddev = self.dark_stddev.select(self.dark_mask.as_1d()).as_double()
    assert flex_dark_stddev.count(0) == 0
    flex_dark_stddev /= flex.mean(flex_dark_stddev)
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


  def chebyshev_common_mode(self, imgs, masks):
    assert len(imgs) == 2
    assert len(masks) == 2
    corrected_imgs = []
    # first fit the variation along the columns of the detector
    sum_y = flex.double()
    for i, (img, mask) in enumerate(zip(imgs, masks)):
      img -= flex.mean(img)
      masked_img = img.deep_copy()
      masked_img.set_selected(~mask, 0.)
      rows, columns = masked_img.all()
      sum_y.extend(flex.sum(masked_img, axis=1))
    # if the region of interest is across a row, use the noise from the
    # same row on the other section
    midpoint = sum_y.size()//2
    for i in range(0, midpoint):
      if sum_y[i] == 0:
        sum_y[i] == sum_y[i + midpoint]
    for i in range(midpoint, sum_y.size()):
      if sum_y[i] == 0:
        sum_y[i] == sum_y[midpoint - i]
    # we assume that both sections have the same variation
    y_obs = sum_y[:midpoint] + sum_y[midpoint:]
    y_obs /= (2 * columns)
    x_obs = flex.double(range(y_obs.size()))
    w_obs = flex.double(x_obs.size(), 1)
    # don't let the edge pixels influence the fit
    w_obs[0] = 1e16
    w_obs[-1] = 1e16
    y_fitted = self.chebyshev_fit(x_obs, y_obs, w_obs, n_terms=5)

    y_correction = flex.double()
    for i in range(columns):
      y_correction.extend(y_fitted)
    y_correction.reshape(flex.grid(columns, rows))
    y_correction.matrix_transpose_in_place()
    if 0:
      from matplotlib import pyplot
      pyplot.imshow(y_correction.as_numpy_array())
      pyplot.show()

    # now fit the variation along the rows
    n_terms = 10
    for img, mask in zip(imgs, masks):
      img -= y_correction
      masked_img = img.deep_copy()
      masked_img.set_selected(~mask, 0.)
      n_masked_rows = flex.sum(masked_img, axis=1).count(0)
      rows, columns = masked_img.all()
      sum_x = flex.sum(masked_img, axis=0)
      sum_x /= (rows - n_masked_rows)
      assert img.all() == (185, 391) # XXX
      # fit one polynome for both asics
      x_obs = sum_x[:194] + sum_x[197:]
      x_obs /= 2
      y_obs = flex.double(range(x_obs.size()))
      w_obs = flex.double(y_obs.size(), 1)
      # mask out the edges and the gap down the middle from the fit
      w_obs.set_selected(y_obs == 0, 1e16)
      w_obs[0] = 1e16
      w_obs[-1] = 1e16
      x_fitted = self.chebyshev_fit(y_obs, x_obs, w_obs, n_terms=10)
      x_calc = x_fitted.deep_copy()
      x_calc.extend(flex.double([0,0,0])) # The 3 pixel gap between asics
      x_calc.extend(x_fitted)

      correction = flex.double()
      for i in range(rows):
        correction.extend(x_calc)
      correction.reshape(img.accessor())
      zero_pixels_sel = (img == 0)
      img -= correction
      if 0:
        from matplotlib import pyplot
        pyplot.imshow(correction.as_numpy_array())
        pyplot.show()

      img.set_selected(zero_pixels_sel, 0)
      corrected_imgs.append(img)

    return corrected_imgs


  def chebyshev_fit(self, x_obs, y_obs, w_obs, n_terms=None):
    from scitbx.math import chebyshev_polynome
    from scitbx.math import chebyshev_lsq_fit
    if n_terms is None:
      # determining the number of terms takes much, much longer than the fit
      n_terms = chebyshev_lsq_fit.cross_validate_to_determine_number_of_terms(
        x_obs, y_obs, w_obs,
        min_terms=5, max_terms=20,
        n_goes=20, n_free=20)
    self.logger.info("Fitting with %i terms" %n_terms)
    fit = chebyshev_lsq_fit.chebyshev_lsq_fit(n_terms, x_obs, y_obs, w_obs)
    self.logger.info("Least Squares residual: %7.6f" %(fit.f))
    fit_funct = chebyshev_polynome(
      n_terms, fit.low_limit, fit.high_limit, fit.coefs)
    y_fitted = fit_funct.f(x_obs)
    if 0:
      # debugging plots
      from matplotlib import pyplot
      pyplot.clf()
      pyplot.plot(x_obs, y_obs)
      pyplot.plot(x_obs, y_fitted)
      pyplot.draw()
      pyplot.show()
    return y_fitted


  def do_photon_counting(self):
    # This only makes sense in combination with some sort of gain correction
    # XXX TODO: count 2, 3, 4, ..., photons
    PHOTON_THRESHOLD = self.photon_threshold
    TWO_PHOTON_THRESHOLD = self.two_photon_threshold
    if [PHOTON_THRESHOLD, TWO_PHOTON_THRESHOLD].count(None) > 0:
      self.logger.info("Skipping photon counting: photon_threshold is not defined")
      return
    self.cspad_img.set_selected(self.cspad_img<PHOTON_THRESHOLD, 0)
    self.cspad_img.set_selected(self.cspad_img>=TWO_PHOTON_THRESHOLD, 2)
    self.cspad_img.set_selected(self.cspad_img>=PHOTON_THRESHOLD, 1)
    self.logger.debug("zero photon counts: %i" %self.cspad_img.count(0))
    self.logger.debug("one photon counts: %i" %self.cspad_img.count(1))
    self.logger.debug("two photon counts: %i" %self.cspad_img.count(2))
    self.logger.info("No. photons: %i" %flex.sum(self.cspad_img))
    s, ms = self.evt_time
    evt_time = s + ms/1000
    self.stats_logger.info("N_PHOTONS %.3f %s" %(evt_time, flex.sum(self.cspad_img)))
