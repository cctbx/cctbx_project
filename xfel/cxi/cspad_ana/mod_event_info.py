import logging

from pypdsdata           import xtc

from xfel.cxi.cspad_ana import cspad_tbx


class mod_event_info(object):
  """Extract basic information from the evt and env objects for each event.
  """


  def __init__(self, address):
    """The mod_event_info class constructor stores the
    parameters passed from the pyana configuration file in instance
    variables.

    @param address         Address string XXX Que?!
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
    self.distance = None
    self.sifoil = None
    self.wavelength = None # The current wavelength - set by self.event()


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

    self.pulse_length = cspad_tbx.evt_pulse_length(evt)
    if (self.pulse_length is None):
      self.nfail += 1
      self.logger.warn("event(): no pulse length, shot skipped")
      evt.put(True, "skip_event")
      return
    self.logger.info("Pulse length: %s" %self.pulse_length)

    self.beam_charge = cspad_tbx.evt_beam_charge(evt)
    if (self.beam_charge is None):
      self.nfail += 1
      self.logger.warn("event(): no beam charge, shot skipped")
      evt.put(True, "skip_event")
      return
    self.logger.info("Beam charge: %s" %self.beam_charge)

    self.injector_xyz = cspad_tbx.env_injector_xyz(env)
    if self.injector_xyz is not None:
      self.logger.info("injector_z: %i" %self.injector_xyz[2].value)


  def endjob(self, env):
    return
