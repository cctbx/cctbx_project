import logging

from pypdsdata import xtc

import libtbx
from xfel.cxi.cspad_ana import cspad_tbx


class mod_event_info(object):
  """Extract basic information from the evt and env objects for each event.
  """


  def __init__(self, address, detz_offset="575", verbose=False):
    """The mod_event_info class constructor stores the
    parameters passed from the pyana configuration file in instance
    variables.

    @param address     Address string XXX Que?!
    @param detz_offset Detector-sample offset in mm, corresponding to
                       longest detector-sample distance
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

    self._detz_offset = cspad_tbx.getOptFloat(detz_offset)

    self.address = cspad_tbx.getOptString(address)
    self.verbose = cspad_tbx.getOptBool(verbose)
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
    env.update(evt)

    self.config = env.getConfig(xtc.TypeId.Type.Id_CspadConfig, self.address)
    if (self.config is None):
      self.logger.error("beginjob(): no config")

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

    distance = cspad_tbx.env_distance(env, self._detz_offset)
    if (distance is None):
      self.nfail += 1
      self.logger.warn("event(): no distance, shot skipped")
      evt.put(True, "skip_event")
      return
    if (self.distance is not None and self.distance != distance):
      self.logger.warn("event(): distance changed mid-run: % 8.4f -> % 8.4f" %
        (self.distance, distance))
    self.distance = distance
    if self.verbose: self.logger.info("Distance: %.4f" %distance)

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
    if self.verbose: self.logger.info("Si-foil thickness: %i" %sifoil)

    self.timestamp = cspad_tbx.evt_timestamp(evt) # human readable format
    self.evt_time = cspad_tbx.evt_time(evt) # tuple of seconds, milliseconds
    if (self.timestamp is None):
      self.nfail += 1
      self.logger.warn("event(): no timestamp, shot skipped")
      evt.put(True, "skip_event")
      return
    if self.verbose: self.logger.info(self.timestamp)

    self.wavelength = cspad_tbx.evt_wavelength(evt)
    if (self.wavelength is None):
      self.nfail += 1
      self.logger.warn("event(): no wavelength, shot skipped")
      evt.put(True, "skip_event")
      return
    if self.verbose: self.logger.info("Wavelength: %.4f" %self.wavelength)

    self.pulse_length = cspad_tbx.evt_pulse_length(evt)
    if (self.pulse_length is None):
      self.nfail += 1
      self.logger.warn("event(): no pulse length, shot skipped")
      evt.put(True, "skip_event")
      return
    if self.verbose: self.logger.info("Pulse length: %s" %self.pulse_length)

    self.beam_charge = cspad_tbx.evt_beam_charge(evt)
    if (self.beam_charge is None):
      self.nfail += 1
      self.logger.warn("event(): no beam charge, shot skipped")
      evt.put(True, "skip_event")
      return
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
      self.logger.info("Laser 1 status: %i" %int(self.laser_1_status.status))
      self.logger.info("Laser 4 status: %i" %int(self.laser_4_status.status))


  def endjob(self, env):
    return


class laser_status(object):

  _status = None
  status_change_timestamp = None

  def __init__(self, laser_id, status=None):
    self._status = status

  class status(libtbx.property):
    def fget(self):
      return self._status

  def set_status(self, status, evt_time):
    if self._status is not None and self._status != status:
      self.status_change_timestamp = evt_time
    self._status = status

  def ms_since_last_status_change(self, evt_time):
    if self.status_change_timestamp is not None:
      return (1000 * (evt_time[0] - self.status_change_timestamp[0])
              + (evt_time[1] - self.status_change_timestamp[1]))
