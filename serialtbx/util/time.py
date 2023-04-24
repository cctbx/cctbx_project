from __future__ import division

import math, time

def now_s_ms():
  """The now_s_ms() function returns the time since
  midnight, 1 January 1970 UTC (Unix time) to millisecond precision.

  @return    Unix time as a tuple of seconds and milliseconds
  """

  t = time.time()
  s = int(math.floor(t))
  return (s, int(round((t - s) * 1000)))

def timestamp(t=None):
  """The timestamp() function returns a string representation of
  an extended human-readable ISO 8601 timestamp.  If @p t is @c None
  the current time is used.

  @param t Tuple of the time in seconds and milliseconds
  @return  Human-readable ISO 8601 timestamp in string representation
  """

  if t is None:
    t = now_s_ms()
  return time.strftime("%Y-%m-%dT%H:%MZ%S", time.gmtime(t[0])) + \
      (".%03d" % t[1])
