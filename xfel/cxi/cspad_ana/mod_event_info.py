# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

from __future__ import absolute_import, division, print_function
import logging

from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import skip_event_flag

class mod_event_info(object):
  """Extract basic information from the evt and env objects for each event.
  """


  def __init__(self, address, detz_offset=575, check_beam_status=True, verbose=False, delta_k=0.0, override_energy=None, **kwds):
    """The mod_event_info class constructor stores the
    parameters passed from the pyana configuration file in instance
    variables.

    @param address           Full data source address of the DAQ
                             device
    @param detz_offset       The distance from the interaction region
                             to the back of the detector stage, in mm
    @param check_beam_status Flag used to skip checking the beam
                             parameters
    @param delta_k           Correction to the K value used when calculating
                             wavelength
    @param override_energy   Use this energy instead of what is in the
                             XTC stream
    """
    logging.basicConfig()
    self.logger = logging.getLogger(self.__class__.__name__)
    self.logger.setLevel(logging.INFO)

    # The subclasses accept keyword arguments; warn about any
    # unhandled arguments at the end of the inheritance chain.
    if len(kwds) > 0:
      self.logger.warning("Ignored unknown arguments: " +
                          ", ".join(kwds))

    # This is for messages that are picked up by Nat's monitoring program
    self.stats_logger = logging.getLogger("stats logger")
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    self.stats_logger.addHandler(handler)
    self.stats_logger.removeHandler(self.stats_logger.handlers[0])
    self.stats_logger.setLevel(logging.INFO)

    self._detz_offset = cspad_tbx.getOptFloat(detz_offset)
    self.delta_k = cspad_tbx.getOptFloat(delta_k)
    self.override_energy = cspad_tbx.getOptFloat(override_energy)

    self.address = cspad_tbx.getOptString(address)
    self.verbose = cspad_tbx.getOptBool(verbose)
    self.check_beam_status = cspad_tbx.getOptBool(check_beam_status)
    self.distance = None
    self.sifoil = None
    self.wavelength = None # The current wavelength - set by self.event()
    self.laser_1_status = laser_status(laser_id=1)
    self.laser_4_status = laser_status(laser_id=4)


  def __del__(self):
    logging.shutdown()


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    # XXX Not needed now that the distance is read in the event?
    if hasattr(env, "update"):
      env.update(evt)

    self.nfail  = 0
    self.nshots = 0
    self.nmemb = 0


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    # Increase the event counter, even if this event is to be skipped.
    self.nshots += 1
    if (evt.get("skip_event")):
      return

    sifoil = cspad_tbx.env_sifoil(env)
    if self.check_beam_status and sifoil is None:
      self.nfail += 1
      self.logger.warning("event(): no Si-foil thickness, shot skipped")
      evt.put(skip_event_flag(), "skip_event")
      return
    if (self.sifoil is not None and self.sifoil != sifoil):
      self.logger.warning("event(): Si-foil changed mid-run: % 8i -> % 8d" %
        (self.sifoil, sifoil))
    self.sifoil = sifoil
    if self.verbose: self.logger.info("Si-foil thickness: %i" %sifoil)

    self.evt_time = cspad_tbx.evt_time(evt) # tuple of seconds, milliseconds
    self.timestamp = cspad_tbx.evt_timestamp(self.evt_time) # human readable format
    if (self.timestamp is None):
      self.nfail += 1
      self.logger.warning("event(): no timestamp, shot skipped")
      evt.put(skip_event_flag(), "skip_event")
      return
    if self.verbose: self.logger.info(self.timestamp)

    if self.override_energy is None:
      self.wavelength = cspad_tbx.evt_wavelength(evt, self.delta_k)
      if self.wavelength is None:
        if self.check_beam_status:
          self.nfail += 1
          self.logger.warning("event(): no wavelength, shot skipped")
          evt.put(skip_event_flag(), "skip_event")
          return
        else:
          self.wavelength = 0
    else:
      self.wavelength = 12398.4187/self.override_energy
    if self.verbose: self.logger.info("Wavelength: %.4f" %self.wavelength)

    self.pulse_length = cspad_tbx.evt_pulse_length(evt)
    if self.pulse_length is None:
      if self.check_beam_status:
        self.nfail += 1
        self.logger.warning("event(): no pulse length, shot skipped")
        evt.put(skip_event_flag(), "skip_event")
        return
      else:
        self.pulse_length = 0
    if self.verbose: self.logger.info("Pulse length: %s" %self.pulse_length)

    self.beam_charge = cspad_tbx.evt_beam_charge(evt)
    if self.beam_charge is None:
      if self.check_beam_status:
        self.nfail += 1
        self.logger.warning("event(): no beam charge, shot skipped")
        evt.put(skip_event_flag(), "skip_event")
        return
      else:
        self.beam_charge = 0
    if self.verbose: self.logger.info("Beam charge: %s" %self.beam_charge)

    self.injector_xyz = cspad_tbx.env_injector_xyz(env)
    #if self.injector_xyz is not None:
      #self.logger.info("injector_z: %i" %self.injector_xyz[2].value)

    self.laser_1_status.set_status(cspad_tbx.env_laser_status(env, laser_id=1), self.evt_time)
    self.laser_4_status.set_status(cspad_tbx.env_laser_status(env, laser_id=4), self.evt_time)
    self.laser_1_ms_since_change = self.laser_1_status.ms_since_last_status_change(self.evt_time)
    self.laser_4_ms_since_change = self.laser_4_status.ms_since_last_status_change(self.evt_time)
    if self.verbose:
      if self.laser_1_ms_since_change is not None:
        self.logger.info("ms since laser 1 status change: %i" %self.laser_1_ms_since_change)
      if self.laser_4_ms_since_change is not None:
        self.logger.info("ms since laser 4 status change: %i" %self.laser_4_ms_since_change)
      if self.laser_1_status is not None and self.laser_1_status.status is not None:
        self.logger.info("Laser 1 status: %i" %int(self.laser_1_status.status))
      if self.laser_4_status is not None and self.laser_4_status.status is not None:
        self.logger.info("Laser 4 status: %i" %int(self.laser_4_status.status))


  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

  def endjob(self, obj1, obj2=None):
    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2


class laser_status(object):

  _status = None
  status_change_timestamp = None

  def __init__(self, laser_id, status=None):
    self._status = status

  @property
  def status(self):
    return self._status

  def set_status(self, status, evt_time):
    if self._status is not None and self._status != status:
      self.status_change_timestamp = evt_time
    self._status = status

  def ms_since_last_status_change(self, evt_time):
    if self.status_change_timestamp is not None:
      return (1000 * (evt_time[0] - self.status_change_timestamp[0])
              + (evt_time[1] - self.status_change_timestamp[1]))
