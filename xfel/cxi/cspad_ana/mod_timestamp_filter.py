# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$

"""The mod_timestamp_filter module extracts timestamps from a file and
uses them to filter events.  Individual timestamps must be separated
by any amount of white space.  Lines where the first non-white space
character is a hash mark are ignored.
"""

__version__ = "$Revision$"

import logging

from xfel.cxi.cspad_ana import cspad_tbx


class mod_timestamp_filter(object):
  """XXX"""

  def __init__(self, timestamps):
    """Extracts timestamps from the file whose name is the string
    pointed to by @p timestamps.  XXX Could do more clever matching
    here?

    @param timestamps Path to file containing timestamps
    """

    self.logger = logging.getLogger(self.__class__.__name__)
    self.logger.setLevel(logging.INFO)

    stream = open(cspad_tbx.getOptString(timestamps), "r")
    self.timestamps = []
    for line in stream.readlines():
      s = line.strip()
      if (len(s) == 0 or s[0] == "#"):
        continue
      self.timestamps += s.split()
    stream.close()


  def __del__(self):
    logging.shutdown()


  def beginjob(self, evt, env):
    pass


  def event(self, evt, env):
    """The event() function puts a "skip_event" object into the event
    unless the timestamp of the shot is in @c self.timestamps.  XXX
    Implement negation?

    @param evt Event data object, a configure object
    @param env Environment object
    """

    timestamp = cspad_tbx.evt_timestamp(evt)
    if (timestamp is None or timestamp not in self.timestamps):
      evt.put(True, "skip_event")
      return
    self.logger.info("Accepted %s" % timestamp)
    return


  def endjob(self, env):
    pass
