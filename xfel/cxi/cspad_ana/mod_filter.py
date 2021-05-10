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
from __future__ import absolute_import, division, print_function
from six.moves import range

__version__ = "$Revision$"

import logging
import re

from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import skip_event_flag

class mod_filter(object):
  def __init__(
    self, timestamps_path=None, timestamps_interval=None, negate="False"):
    """Extracts timestamps from the file whose name is the string
    pointed to by @p timestamps_path.

    @param timestamps_path     Path to file containing timestamps
    @param timestamps_interval Comma-separated inclusive endpoints of
                               a timestamp interval.  The lower or the
                               upper endpoint may be omitted, in which
                               case it will be treated as -infinity or
                               +infinity.
    @param negate              Select shots not matching any of the
                               timestamps
    """

    self.logger = logging.getLogger(self.__class__.__name__)
    self.logger.setLevel(logging.INFO)

    if (not((timestamps_path is None) ^ (timestamps_interval is None))):
      raise RuntimeError(
        "Must specify either timestamps_path or timestamps_interval")

    self.negate = cspad_tbx.getOptBool(negate)

    if (timestamps_path is not None):
      p_old = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}Z\d{2}\.\d{3}")
      p_new = re.compile(r"\d{17}")
      f = open(timestamps_path, "r")
      self.timestamps_list = []
      for line in f.readlines():
        s = line.strip()
        if (len(s) == 0 or s[0] == "#"):
          continue
        for t in s.split():
          m = p_old.findall(t)
          if len(m) == 0:
            m = p_new.findall(t)
          if (len(m) == 1):
            self.timestamps_list.append(self._ts2sms(m[0]))
      f.close()
      self.timestamps_interval = None
    else:
      try:
        s = timestamps_interval.split(",")
        if (len(s) == 1 or len(s) == 2 and len(s[1]) == 0):
          self.timestamps_interval = (self._ts2sms(s[0]), None)
        elif (len(s) == 2 and len(s[0]) == 0):
          self.timestamps_interval = (None, self._ts2sms(s[1]))
        elif (len(s) == 2):
          self.timestamps_interval = (self._ts2sms(s[0]), self._ts2sms(s[1]))
        else:
          raise ValueError()

        # Ensure lower endpoint is earlier than upper endpoint.
        if (self.timestamps_interval[0] is not None and
            self.timestamps_interval[1] is not None and
            self._timestamp_compar(self.timestamps_interval[0],
                                   self.timestamps_interval[1]) > 0):
          raise ValueError()

      except ValueError:
        raise RuntimeError(
          "Failed to parse timestamp interval %s" % timestamps_interval)
      self.timestamps_list = None

    self.naccepted = 0
    self.nshots = 0


  def __del__(self):
    logging.shutdown()


  def _timestamp_compar(self, l, r):
    """The _timestamp_compar() function returns an integer less than,
    equal to, or greater than zero if @p l is considered to be
    respectively less than, equal to or greater than @p r.

    @param l First timestamp, two-tuple of s and ms
    @param r Second timestamp, two-tuple of s and ms
    @return  Negative if @p l < @p r, zero if @p l == @p r, or
             positive if @p l > @p r
    """

    for i in range(2):
      if (l[i] < r[i]):
        return -1
      if (l[i] > r[i]):
        return +1
    return 0


  def _tir(self, t, i):
    """The _tir() function returns @c True if the timestamp @p t lies
    in the inclusive interval given by @p i, and @c False otherwise.
    """

    if (i[0] is not None and self._timestamp_compar(t, i[0]) < 0 or
        i[1] is not None and self._timestamp_compar(t, i[1]) > 0):
      return False
    return True


  def _ts2sms(self, timestamp):
    """The _ts2sms() function converts a string representation of a
    timestamp to a two-tuple of seconds and milliseconds since the
    epoch.  The function raises @c ValueError if @p timestamp cannot
    be interpreted.  The _ts2sms() function is the inverse of
    cspad_tbx.evt_timestamp().

    @param timestamp String representation of human-readable ISO 8601
                     timestamp
    @return          Tuple of the time in seconds and milliseconds
                     since the epoch
    """

    from calendar import timegm
    from time import strptime

    if (len(timestamp) != 17 and (len(timestamp) != 23 or timestamp[19] != ".")):
      raise ValueError()

    if len(timestamp) == 17:
      s = timegm(strptime(timestamp[0:14], "%Y%m%d%H%M%S"))
      ms = int(timestamp[14:])

    else:
      # int() and strptime() both raise ValueError on error.
      s = timegm(strptime(timestamp[0:19], "%Y-%m-%dT%H:%MZ%S"))
      ms = int(timestamp[20:23])

    return (s, ms)


  def beginjob(self, evt, env):
    pass


  def event(self, evt, env):
    """The event() function puts a "skip_event" object with value @c
    True into the event if the shot is to be skipped.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    self.nshots += 1
    if (evt.get("skip_event")):
      return

    if (self.timestamps_list is not None):
      t = cspad_tbx.evt_time(evt)
      if (t is None):
        self.logger.warning("event(): no timestamp, shot skipped. Shot: %s"%self.nshots)
        evt.put(skip_event_flag(), "skip_event")
        return
      elif (self.negate and t in self.timestamps_list):
        evt.put(skip_event_flag(), "skip_event")
        return
      elif (not self.negate and t not in self.timestamps_list):
        evt.put(skip_event_flag(), "skip_event")
        return

    else:
      t = cspad_tbx.evt_time(evt)
      if (t is None):
        self.logger.warning("event(): no timestamp, shot skipped. Shot: %s"%self.nshots)
        evt.put(skip_event_flag(), "skip_event")
        return
      if (self.negate and self._tir(t, self.timestamps_interval)):
        evt.put(skip_event_flag(), "skip_event")
        return
      elif (not self.negate and not self._tir(t, self.timestamps_interval)):
        evt.put(skip_event_flag(), "skip_event")
        return

    self.logger.info("event(): event accepted. Shot: %s, TS: %s"%(self.nshots,cspad_tbx.evt_timestamp(t)))
    self.naccepted += 1

  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

  def endjob(self, obj1, obj2=None):
    """
    @param evt Event object (psana only)
    @param env Environment object
    """

    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2

    self.logger.info(
      "Saw %d shots, accepted %d, skipped %d" %
      (self.nshots, self.naccepted, self.nshots - self.naccepted))
