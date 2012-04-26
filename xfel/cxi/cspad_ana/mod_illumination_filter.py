# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$

"""The mod_illumination_filter module filter events by their
illumination condition.  When illumination conditions matter, laser 4
should always be off, and laser 1 must have had its state changed for
least laser_wait_time ms before the event.  Events with missing
timestamps are always skipped.
"""

__version__ = "$Revision$"

import logging

from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana.mod_event_info import laser_status


class mod_illumination_filter(object):
  def __init__(self, illumination, laser_wait_time="2000"):
    """Initialise laser status, and validate input.  The @p
    illumination parameter is mandatory.

    @param illumination    If @c dark or @c light, only pass dark or
                           light shots, respectively.  If @c other,
                           pass events that are neither light or dark,
                           due to recent state changes to the lasers.
    @param laser_wait_time Number of ms the lasers have to be stable
                           before classifying events as light or dark.
                           Jan F. Kern recommends a value between 1000
                           and 2000 ms.
    """

    self.logger = logging.getLogger(self.__class__.__name__)
    self.logger.setLevel(logging.INFO)

    self._filter = cspad_tbx.getOptString(illumination)
    if (self._filter != "dark" and
        self._filter != "light" and
        self._filter != "other"):
      raise RuntimeError(
        "Parameter illumination must be either "
        "\"light\", \"dark\", or \"other\"")

    self._wait = cspad_tbx.getOptFloat(laser_wait_time)
    if (self._wait is None or self._wait < 0):
      raise RuntimeError(
        "Parameter laser_wait_time must be number >= 0")

    self.laser_1 = laser_status(laser_id=1)
    self.laser_4 = laser_status(laser_id=4)

    self.naccepted = 0
    self.nshots = 0


  def __del__(self):
    logging.shutdown()


  def beginjob(self, evt, env):
    """The beginjob() function initialises the last status change to
    the beginning of the Unix epoch.  XXX Arbitrary--could choose the
    smallest (negative) integers instead, the earliest time possible.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    t = (0, 0)
    self.laser_1.set_status(cspad_tbx.env_laser_status(env, laser_id=1), t)
    self.laser_4.set_status(cspad_tbx.env_laser_status(env, laser_id=4), t)


  def event(self, evt, env):
    """The event() function puts a "skip_event" object with value @c
    True into the event if the shot is to be skipped.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    self.nshots += 1
    if (evt.get("skip_event")):
      return

    # Get time as a (seconds, milliseconds) tuple.
    t = cspad_tbx.evt_time(evt)
    if (t is None):
      self.logger.warn("event(): no timestamp, shot skipped")
      evt.put(True, "skip_event")
      return

    # Update laser status.
    self.laser_1.set_status(cspad_tbx.env_laser_status(env, laser_id=1), t)
    self.laser_4.set_status(cspad_tbx.env_laser_status(env, laser_id=4), t)

    if (self.laser_4.status or
        self.laser_1.ms_since_last_status_change(t) < self._wait or
        self.laser_4.ms_since_last_status_change(t) < self._wait):
      # If laser 4 is on or was switched off less than self._wait ms
      # ago, the shot falls in the "other" category.  If laser 1
      # changed its state less than self._wait ms ago the shot falls
      # in the "other" category.

      if (self._filter != "other"):
        evt.put(True, "skip_event")
        return

    elif (self.laser_1.status):
      # If laser 1 is on the shot falls in the "light" category.
      if (self._filter != "light"):
        evt.put(True, "skip_event")
        return

    elif (not self.laser_1.status):
      # If laser 1 is off the shot falls in the "dark" category.
      if (self._filter != "dark"):
        evt.put(True, "skip_event")
        return

    else:
      # NOTREACHED
      self.logger.error("Could not determine shot category")
      raise RuntimeError("XXX")

    self.naccepted += 1


  def endjob(self, env):
    self.logger.info(
      "Saw %d shots, accepted %d, skipped %d" %
      (self.nshots, self.naccepted, self.nshots - self.naccepted))
