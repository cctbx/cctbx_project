# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$

"""Hitfinding for CSPad images

XXX
"""

__version__ = "$Revision$"

from scitbx.array_family import flex
from xfel.cxi.cspad_ana.hitfinder_tbx import distl_hitfinder
from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx

# import matplotlib
# matplotlib.use("PDF")


class mod_hitfind(common_mode.common_mode_correction, distl_hitfinder):
  """Class for hitfinding within the pyana framework
  """

  def __init__(self,
               address,
               dispatch               = None,
               integration_dirname    = None,
               integration_basename   = None,
               out_dirname            = None,
               out_basename           = None,
               roi                    = None,
               distl_min_peaks        = None,
               distl_flags            = None,
               threshold              = None,
               xtal_target            = None,
               negate_hits            = False,
               **kwds):
    """The mod_hitfind class constructor stores the parameters passed
    from the pyana configuration file in instance variables.  All
    parameters, except @p address are optional, and hence need not be
    defined in pyana.cfg.

    @param address      Address string XXX Que?!
    @param dispatch     Function to call
    @param out_dirname  Directory portion of output image
    @param out_basename Filename prefix of output image
    @param roi          Region of interest for thresholding, on the
                        form fast_low:fast_high,slow_low:slow_high
    @param threshold    Minimum value in region of interest to pass
    """

    super(mod_hitfind, self).__init__(address=address, **kwds)

    self.m_dispatch             = cspad_tbx.getOptString(dispatch)
    self.m_integration_basename = cspad_tbx.getOptString(integration_basename)
    self.m_integration_dirname  = cspad_tbx.getOptString(integration_dirname)
    self.m_out_basename         = cspad_tbx.getOptString(out_basename)
    self.m_out_dirname          = cspad_tbx.getOptString(out_dirname)
    self.m_distl_min_peaks      = cspad_tbx.getOptInteger(distl_min_peaks)
    self.m_distl_flags          = cspad_tbx.getOptStrings(distl_flags)
    self.m_threshold            = cspad_tbx.getOptInteger(threshold)
    self.m_xtal_target          = cspad_tbx.getOptString(xtal_target)
    self.m_negate_hits          = cspad_tbx.getOptBool(negate_hits)

    # A ROI should not contain any ASIC boundaries, as these are
    # noisy.  Hence circular ROI:s around the beam centre are probably
    # not such a grand idea.
    self.m_roi = cspad_tbx.getOptROI(roi)

    # Verify that dist_min_peaks is either "restrictive" or
    # "permissive", but not both.  ^ is the logical xor operator
    if self.m_distl_min_peaks is not None:
      if (not (('permissive'  in self.m_distl_flags) ^
               ('restrictive' in self.m_distl_flags))):
        raise RuntimeError("""Sorry, with the distl_min_peaks option,
          distl_flags must be set to 'permissive' or 'restrictive'.""")
      if (self.m_roi is not None):
        raise RuntimeError("""Sorry, either specify region of interest
          (roi) or distl_min_peaks, but not both.""")


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_hitfind, self).beginjob(evt, env)
    self.set_up_hitfinder()

  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.
    XXX more?

    Previously, common-mode correction was applied only after initial
    threshold filtering.  Since the common_mode class applies the
    (lengthy) common-mode correction immediately after reading the
    image from the stream, this optimisation is currently not
    (elegantly) doable.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_hitfind, self).event(evt, env)
    if (evt.get("skip_event")):
      return

    # ***** HITFINDING ***** XXX For hitfinding it may be interesting
    # to look at the fraction of subzero pixels in the dark-corrected
    # image.
    if (self.m_threshold is not None):
      # If a threshold value is given it can be applied in one of three ways:
      #    1.  Apply it over the whole image
      if (self.m_roi is None and self.m_distl_min_peaks is None):
        vmax = flex.max(self.cspad_img)
        if (vmax < self.m_threshold):
          if not self.m_negate_hits:
            # Tell downstream modules to skip this event if the threshold was not met.
            evt.put(True, "skip_event")
            return
        elif self.m_negate_hits:
          evt.put(True, "skip_event")
          return

      #    2. Apply threshold over a rectangular region of interest.
      elif (self.m_roi is not None):
        vmax = flex.max(self.cspad_img[self.m_roi[2]:self.m_roi[3],
                                       self.m_roi[0]:self.m_roi[1]])
        if (vmax < self.m_threshold):
          if not self.m_negate_hits:
            evt.put(True, "skip_event")
            return
        elif self.m_negate_hits:
          evt.put(True, "skip_event")
          return

      #    3. Determine the spotfinder spots within the central ASICS, and accept the
      #       image as a hit if there are m_distl_min_peaks exceeding m_threshold.
      #       As a further requirement, the peaks must exceed 2.5 * the 90-percentile
      #       pixel value of the central ASICS.  This filter was added to avoid high-background
      #       false positives.
      elif (self.m_distl_min_peaks is not None):

        peak_heights,outvalue = self.distl_filter(
          self.cspad_img.iround(), # XXX correct?
          self.distance,
          self.timestamp,
          self.wavelength)
        if ('permissive' in self.m_distl_flags):
          number_of_accepted_peaks = (peak_heights > self.m_threshold).count(True)
        else:
          number_of_accepted_peaks = ((peak_heights > self.m_threshold).__and__(outvalue==0)).count(True)

        sec,ms = cspad_tbx.evt_time(evt)
        evt_time = sec + ms/1000
        self.stats_logger.info("BRAGG %.3f %d" %(evt_time, number_of_accepted_peaks))


        if number_of_accepted_peaks < self.m_distl_min_peaks:
          self.logger.info("Subprocess %02d: Spotfinder NO  HIT image #%05d @ %s; %d spots > %d" %(
              env.subprocess(), self.nshots, self.timestamp, number_of_accepted_peaks, self.m_threshold))
          if not self.m_negate_hits:
            evt.put(True, "skip_event")
            return
        else:
          self.logger.info("Subprocess %02d: Spotfinder YES HIT image #%05d @ %s; %d spots > %d" %(
              env.subprocess(), self.nshots, self.timestamp, number_of_accepted_peaks, self.m_threshold))
          if self.m_negate_hits:
            evt.put(True, "skip_event")
            return

    self.logger.info("Subprocess %02d: process image #%05d @ %s" %
                     (env.subprocess(), self.nshots, self.timestamp))

    d = cspad_tbx.dpack(
      active_areas    = self.active_areas,
      beam_center_x   = cspad_tbx.pixel_size * self.beam_center[0],
      beam_center_y   = cspad_tbx.pixel_size * self.beam_center[1],
      data            = self.cspad_img.iround(), # XXX ouch!
      distance        = self.distance,
      timestamp       = self.timestamp,
      sequence_number = self.nshots,
      wavelength      = self.wavelength,
      xtal_target     = self.m_xtal_target)

    if (self.m_dispatch == "index"):
      import sys
      from xfel.cxi.integrate_image_api import integrate_one_image
      info = integrate_one_image(d,
                          integration_dirname  = self.m_integration_dirname,
                          integration_basename = self.m_integration_basename)
      sys.stdout = sys.__stdout__
      sys.stderr = sys.__stderr__
      if (info is None):
        evt.put(True, "skip_event")
        return

    elif (self.m_dispatch == "nop"):
      pass

    elif (self.m_dispatch == "view"): #interactive image viewer

      from cxi_user.xfel_targets import targets
      args = ["indexing.data=dummy",
              "distl.detector_format_version=CXI 5.1",
              ] + targets[self.m_xtal_target]

      from spotfinder.applications.xfel import cxi_phil
      horizons_phil = cxi_phil.cxi_versioned_extract(args)
      horizons_phil.indexing.data = d

      from xfel.cxi import display_spots
      display_spots.parameters.horizons_phil = horizons_phil
      display_spots.wrapper_of_callback().display(horizons_phil.indexing.data)

    elif (self.m_dispatch == "spots"): #interactive spotfinder viewer

      from cxi_user.xfel_targets import targets
      args = ["indexing.data=dummy",
              "distl.detector_format_version=CXI 5.1",
              ] + targets[self.m_xtal_target]

      from spotfinder.applications.xfel import cxi_phil
      horizons_phil = cxi_phil.cxi_versioned_extract(args)
      horizons_phil.indexing.data = d

      from xfel.cxi import display_spots
      display_spots.parameters.horizons_phil = horizons_phil

      from rstbx.new_horizons.index import pre_indexing_validation,pack_names
      pre_indexing_validation(horizons_phil)
      imagefile_arguments = pack_names(horizons_phil)
      horizons_phil.persist.show()
      from spotfinder.applications import signal_strength
      info = signal_strength.run_signal_strength_core(horizons_phil,imagefile_arguments)

      work = display_spots.wrapper_of_callback(info)
      work.display_with_callback(horizons_phil.indexing.data)

    elif (self.m_dispatch == "write_dict"):
      if (self.m_out_dirname  is not None or
          self.m_out_basename is not None):
        cspad_tbx.dwritef(d, self.m_out_dirname, self.m_out_basename)

    # Diagnostic message emitted only when all the processing is done.
    if (env.subprocess() >= 0):
      self.logger.info("Subprocess %02d: accepted #%05d @ %s" %
                       (env.subprocess(), self.nshots, self.timestamp))
    else:
      self.logger.info("Accepted #%05d @ %s" %
                       (self.nshots, self.timestamp))


  def endjob(self, env):
    """The endjob() function logs the number of processed shots.

    @param env Environment object
    """

    super(mod_hitfind, self).endjob(env)
    if (env.subprocess() >= 0):
      self.logger.info("Subprocess %02d: processed %d shots" %
                       (env.subprocess(), self.nshots))
    else:
      self.logger.info("Processed %d shots" % self.nshots)
