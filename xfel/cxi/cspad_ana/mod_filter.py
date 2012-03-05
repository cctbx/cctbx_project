# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$

"""The mod_filter module extracts timestamps from a file and uses them
to filter events.  A single valid timestamp can be embedded anywhere
in a word, where individual words must be separated by any amount of
white space.  Lines where the first non-white space character is a
hash mark are ignored.

By default, only events whose timestamps match are passed through to
downwind modules.  If @c negate is @c True, mismatching events are
selected instead.  Events with missing timestamps are always skipped.
"""

__version__ = "$Revision$"

import logging
import re

from xfel.cxi.cspad_ana import cspad_tbx


class mod_filter(object):
  def __init__(
    self, timestamps_path=None, timestamps_range=None, negate="False"):
    """Extracts timestamps from the file whose name is the string
    pointed to by @p timestamps_path.

    @param timestamps_path  Path to file containing timestamps
    @param timestamps_range Inclusive timestamp range, XXX Experimental--test and document open ranges.
    @param negate           Select shots not matching any of the
                            timestamps
    """

    self.logger = logging.getLogger(self.__class__.__name__)
    self.logger.setLevel(logging.INFO)

    if (not((timestamps_path is None) ^ (timestamps_range is None))):
      raise RuntimeError(
        "Must specify either timestamps_path or timestamps_range")

    self.negate = cspad_tbx.getOptBool(negate)

    # Regular expression to match timestamp.
    p = re.compile("\d{4}-\d{2}-\d{2}T\d{2}:\d{2}Z\d{2}\.\d{3}")

    if (timestamps_path is not None):
      f = open(timestamps_path, "r")
      self.timestamps_list = []
      for line in f.readlines():
        s = line.strip()
        if (len(s) == 0 or s[0] == "#"):
          continue
        for t in s.split():
          m = p.findall(t)
          if (len(m) == 1):
            self.timestamps_list.append(m[0])
      f.close()
      self.timestamps_range = None
    else:
      try:
        s = timestamps_range.split(",")
        if (len(s) == 1 or len(s) == 2 and len(s[1]) == 0):
          self.timestamps_range = (self._ts2sms(s[0]), None)
        elif (len(s) == 2 and len(s[0]) == 0):
          self.timestamps_range = (None, self._ts2sms(s[1]))
        elif (len(s) == 2):
          self.timestamps_range = (self._ts2sms(s[0]), self._ts2sms(s[1]))
        else:
          raise ValueError()

        # Ensure lower bound is less than upper bound.
        if (self.timestamps_range[0][0] > self.timestamps_range[1][0] or
            self.timestamps_range[0][0] == self.timestamps_range[1][0] and
            self.timestamps_range[0][1] > self.timestamps_range[1][1]):
          raise ValueError()

      except ValueError:
        raise RuntimeError(
          "Failed to parse timestamp range %s" % timestamps_range)
      self.timestamps_list = None

    self.naccepted = 0
    self.nshots = 0


  def __del__(self):
    logging.shutdown()


  def _tir(self, t, r):
    """The _tir() function returns @c True if the timestamp @p t lies
    in the inclusive range given by @p r, and @c False otherwise."""

    if (r[0] is not None):
      if (t[0] < r[0][0] or t[0] == r[0][0] and t[1] < r[0][1]):
        return False

    if (r[1] is not None):
      if (t[0] > r[1][0] or t[0] == r[1][0] and t[1] > r[1][1]):
        return False

    return True


  def _ts2sms(self, timestamp):
    """The _ts2sms() function converts a string representation of a
    timestamp to a two-tuple seconds and milliseconds since the epoch.
    The function raises @c ValueError if @p timestamp cannot be
    interpreted."""
    from calendar import timegm
    from time import strptime

    if (len(timestamp) != 23 or timestamp[19] != "."):
      raise ValueError()

    # int() and strptime() both raise ValueError on error.  See
    # cspad_tbx.evt_timestamp for the format.
    s = timegm(strptime(timestamp[0:19], "%Y-%m-%dT%H:%MZ%S"))
    ms = int(timestamp[20:23])

    return (s, ms)

  def beginjob(self, evt, env):
    pass


  def event(self, evt, env):
    """The event() function puts a "skip_event" object into the event
    if the shot is to be skipped.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    self.nshots += 1
    if (evt.get("skip_event")):
      return

    if (self.timestamps_list is not None):
      t = cspad_tbx.evt_timestamp(evt)
      if (t is None):
        evt.put(True, "skip_event")
        return
      elif (self.negate and t in self.timestamps_list):
        evt.put(True, "skip_event")
        return
      elif (not self.negate and t not in self.timestamps_list):
        evt.put(True, "skip_event")
        return

    else:
      t = cspad_tbx.evt_time(evt)
      if (t is None):
        evt.put(True, "skip_event")
        return
      if (self.negate and self._tir(t, self.timestamps_range)):
        evt.put(True, "skip_event")
        return
      elif (not self.negate and not self._tir(t, self.timestamps_range)):
        evt.put(True, "skip_event")
        return

    #self.logger.info("*** using shot %i ***" %self.nshots)
    self.naccepted += 1


  def endjob(self, env):
    self.logger.info(
      "Saw %d shots, accepted %d, skipped %d" %
      (self.nshots, self.naccepted, self.nshots - self.naccepted))
