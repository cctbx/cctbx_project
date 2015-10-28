# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

"""Hitfinding for CSPad images

XXX
"""
from __future__ import division

__version__ = "$Revision$"

from scitbx.array_family import flex
from xfel.cxi.cspad_ana.hitfinder_tbx import distl_hitfinder
from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import rayonix_tbx
from xfel.cxi.cspad_ana import skip_event_flag
from xfel.detector_formats import detector_format_version as detector_format_function
import getpass

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
               trial_id               = None,
               db_logging             = False,
               sql_buffer_size        = 1,
               db_host                = None,
               db_name                = None,
               db_user                = None,
               db_password            = None,
               **kwds):
    """The mod_hitfind class constructor stores the parameters passed
    from the pyana configuration file in instance variables.  All
    parameters, except @p address are optional, and hence need not be
    defined in pyana.cfg.

    @param address      Full data source address of the DAQ device
    @param dispatch     Function to call
    @param integration_dirname
                        Directory portion of output integration file
                        pathname
    @param integration_basename
                        Filename prefix of output integration file
                        pathname
    @param out_dirname  Directory portion of output image pathname
    @param out_basename Filename prefix of output image pathname
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
    self.m_trial_id             = cspad_tbx.getOptInteger(trial_id)
    self.m_db_logging           = cspad_tbx.getOptBool(db_logging)
    self.m_sql_buffer_size      = cspad_tbx.getOptInteger(sql_buffer_size)
    self.m_db_host              = cspad_tbx.getOptString(db_host)
    self.m_db_name              = cspad_tbx.getOptString(db_name)
    self.m_db_user              = cspad_tbx.getOptString(db_user)
    self.m_db_password          = cspad_tbx.getOptString(db_password)
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

    if self.m_db_logging:
      self.buffered_sql_entries = []
      assert self.m_sql_buffer_size >= 1



  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_hitfind, self).beginjob(evt, env)
    self.set_up_hitfinder()

    if self.m_db_logging:
      from cxi_xdr_xes.cftbx.cspad_ana import db
      self.logger.info("Connecting to db...")
      dbobj = db.dbconnect(self.m_db_host, self.m_db_name, self.m_db_user, self.m_db_password)
      assert dbobj.open
      self.logger.info("Connected.")

      try:
        self.trial = self.m_trial_id # TODO: beat the race condition and use db.get_next_trial_id if
                                      # this is not set or is zero or less
        db.create_tables(dbobj)

      except Exception,e:
        self.logger.info("Couldn't create root tables: %s"%(e))
      dbobj.close()

    """ This doesn't work.  The many threads add the values over and over to the master db :(
    try:
      trial = 1244
      cmd = "SELECT * FROM %s WHERE trial = %%s"%(db.root_table_name)
      count = cursor.execute(cmd, trial)
      self.logger.info("Count is %s"%(count))

      if count < 3:
        cmd = "INSERT INTO %s (trial,experiment,user,datatable) VALUES (%%s,%%s,%%s,'cxi_braggs_front');"%(db.root_table_name)
        #self.logger.info("here!!")
        #self.logger.info(cmd%(123,env.experiment(),getpass.getuser()))
        cursor.execute(cmd, (trial,env.experiment(),getpass.getuser()))

        cmd = "INSERT INTO %s (trial,experiment,user,datatable) VALUES (%%s,%%s,%%s,'cxi_braggs_back');"%(db.root_table_name)
        cursor.execute(cmd, (trial,env.experiment(),getpass.getuser()))

        cmd = "INSERT INTO %s (trial,experiment,user,datatable) VALUES (%%s,%%s,%%s,'cxi_xes');"%(db.root_table_name)
        cursor.execute(cmd, (trial,env.experiment(),getpass.getuser()))

        self.db.commit()
    except Exception,e:
      self.logger.info("Couldn't create root entries: %s"%(e))
    """

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

    # This module only applies to detectors for which a distance is
    # available.
    distance = cspad_tbx.env_distance(self.address, env, self._detz_offset)
    if distance is None:
      self.nfail += 1
      self.logger.warning("event(): no distance, shot skipped")
      evt.put(skip_event_flag(), "skip_event")
      return

    device = cspad_tbx.address_split(self.address)[2]

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
            evt.put(skip_event_flag(), "skip_event")
            return
        elif self.m_negate_hits:
          evt.put(skip_event_flag(), "skip_event")
          return

      #    2. Apply threshold over a rectangular region of interest.
      elif (self.m_roi is not None):
        vmax = flex.max(self.cspad_img[self.m_roi[2]:self.m_roi[3],
                                       self.m_roi[0]:self.m_roi[1]])
        if (vmax < self.m_threshold):
          if not self.m_negate_hits:
            evt.put(skip_event_flag(), "skip_event")
            return
        elif self.m_negate_hits:
          evt.put(skip_event_flag(), "skip_event")
          return

      #    3. Determine the spotfinder spots within the central ASICS, and accept the
      #       image as a hit if there are m_distl_min_peaks exceeding m_threshold.
      #       As a further requirement, the peaks must exceed 2.5 * the 90-percentile
      #       pixel value of the central ASICS.  This filter was added to avoid high-background
      #       false positives.
      elif (self.m_distl_min_peaks is not None):
        if device == 'marccd':
          self.hitfinder_d['BEAM_CENTER_X'] = self.beam_center[0]
          self.hitfinder_d['BEAM_CENTER_Y'] = self.beam_center[1]
        elif device == 'Rayonix':
          self.hitfinder_d['BEAM_CENTER_X'] = self.beam_center[0]
          self.hitfinder_d['BEAM_CENTER_Y'] = self.beam_center[1]

        peak_heights,outvalue = self.distl_filter(
          self.address,
          self.cspad_img.iround(), # XXX correct?
          distance,
          self.timestamp,
          self.wavelength)
        if ('permissive' in self.m_distl_flags):
          number_of_accepted_peaks = (peak_heights > self.m_threshold).count(True)
        else:
          number_of_accepted_peaks = ((peak_heights > self.m_threshold).__and__(outvalue==0)).count(True)

        sec,ms = cspad_tbx.evt_time(evt)
        evt_time = sec + ms/1000
        self.stats_logger.info("BRAGG %.3f %d" %(evt_time, number_of_accepted_peaks))

        skip_event = False
        if number_of_accepted_peaks < self.m_distl_min_peaks:
          self.logger.info("Subprocess %02d: Spotfinder NO  HIT image #%05d @ %s; %d spots > %d" %(
            env.subprocess(), self.nshots, self.timestamp, number_of_accepted_peaks, self.m_threshold))

          if not self.m_negate_hits:
            skip_event = True
        else:
          self.logger.info("Subprocess %02d: Spotfinder YES HIT image #%05d @ %s; %d spots > %d" %(
            env.subprocess(), self.nshots, self.timestamp, number_of_accepted_peaks, self.m_threshold))

          if self.m_negate_hits:
            skip_event = True

        if skip_event:
          if self.m_db_logging:
            # log misses to the database
            self.queue_entry((self.trial, evt.run(), "%.3f"%evt_time, number_of_accepted_peaks, distance,
                              self.sifoil, self.wavelength, False))
          evt.put(skip_event_flag(), "skip_event")
          return
        # the indexer will log this hit when it is ran. Bug: if the spotfinder is ran by itself, this
        # hit will not be logged in the db.
        evt.put(number_of_accepted_peaks, 'sfspots')

    self.logger.info("Subprocess %02d: process image #%05d @ %s" %
                     (env.subprocess(), self.nshots, self.timestamp))

    # See r17537 of mod_average.py.
    if device == 'Cspad':
      pixel_size = cspad_tbx.pixel_size
      saturated_value = cspad_tbx.cspad_saturated_value
    elif device == 'marccd':
      pixel_size = evt.get("marccd_pixel_size")
      saturated_value = evt.get("marccd_saturated_value")
    elif device == 'Rayonix':
      pixel_size = rayonix_tbx.get_rayonix_pixel_size(self.bin_size)
      saturated_value = rayonix_tbx.rayonix_saturated_value

    d = cspad_tbx.dpack(
      active_areas=self.active_areas,
      address=self.address,
      beam_center_x=pixel_size * self.beam_center[0],
      beam_center_y=pixel_size * self.beam_center[1],
      data=self.cspad_img.iround(), # XXX ouch!
      distance=distance,
      pixel_size=pixel_size,
      saturated_value=saturated_value,
      timestamp=self.timestamp,
      wavelength=self.wavelength,
      xtal_target=self.m_xtal_target)

    if (self.m_dispatch == "index"):
      import sys
      from xfel.cxi.integrate_image_api import integrate_one_image
      info = integrate_one_image(d,
                                 integration_dirname  = self.m_integration_dirname,
                                 integration_basename = self.m_integration_basename)
      sys.stdout = sys.__stdout__
      sys.stderr = sys.__stderr__

      indexed = info is not None

      if self.m_db_logging:
        sec,ms = cspad_tbx.evt_time(evt)
        evt_time = sec + ms/1000
        sfspots = evt.get('sfspots')
        if sfspots is None:
          if indexed:
            n_spots = len(info.spotfinder_results.images[info.frames[0]]['spots_total'])
          else:
            n_spots = 0
        else:
          n_spots = sfspots

        self.queue_entry((self.trial, evt.run(), "%.3f"%evt_time, n_spots, distance,
                          self.sifoil, self.wavelength, indexed))

      if (not indexed):
        evt.put(skip_event_flag(), "skip_event")
        return

    elif (self.m_dispatch == "nop"):
      pass

    elif (self.m_dispatch == "view"): #interactive image viewer

      args = ["indexing.data=dummy"]
      detector_format_version = detector_format_function(
        self.address, evt.GetTime())
      if detector_format_version is not None:
        args += ["distl.detector_format_version=%" % detector_format_version]

      from xfel.phil_preferences import load_cxi_phil
      horizons_phil = load_cxi_phil(self.m_xtal_target, args)
      horizons_phil.indexing.data = d

      from xfel.cxi import display_spots
      display_spots.parameters.horizons_phil = horizons_phil
      display_spots.wrapper_of_callback().display(horizons_phil.indexing.data)

    elif (self.m_dispatch == "spots"): #interactive spotfinder viewer

      args = ["indexing.data=dummy"]
      detector_format_version = detector_format_function(
        self.address, evt.GetTime())
      if detector_format_version is not None:
        args += ["distl.detector_format_version=%s" % detector_format_version]

      from xfel.phil_preferences import load_cxi_phil
      horizons_phil = load_cxi_phil(self.m_xtal_target, args)
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
      self.logger.warning(
        "event(): deprecated dispatch 'write_dict', use mod_dump instead")
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

  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

  def endjob(self, obj1, obj2=None):
    """The endjob() function logs the number of processed shots.

    @param evt Event object (psana only)
    @param env Environment object
    """

    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2

    super(mod_hitfind, self).endjob(env)
    if (env.subprocess() >= 0):
      self.logger.info("Subprocess %02d: processed %d shots" %
                       (env.subprocess(), self.nshots))
    else:
      self.logger.info("Processed %d shots" % self.nshots)

    if self.m_db_logging:
      self.commit_entries()

  def queue_entry(self, entry):
    self.buffered_sql_entries.append(entry)
    if len(self.buffered_sql_entries) >= self.m_sql_buffer_size:
      self.commit_entries()

  def commit_entries(self):
    if self.m_sql_buffer_size > 1 and len(self.buffered_sql_entries) > 0:
      from cxi_xdr_xes.cftbx.cspad_ana import db
      dbobj = db.dbconnect()
      cursor = dbobj.cursor()
      cmd = "INSERT INTO %s (trial,run,eventstamp,hitcount,distance,sifoil,wavelength,indexed) VALUES "%(db.table_name)
      comma = ""
      for entry in self.buffered_sql_entries:
        cmd += comma + "(%s,%s,%s,%s,%s,%s,%s,%s)"%entry
        comma = ", "
      cursor.execute(cmd)
      dbobj.commit()
      dbobj.close()
      self.buffered_sql_entries = []
