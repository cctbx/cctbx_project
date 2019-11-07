from __future__ import absolute_import, division, print_function

import time, math

def now():
  return "%02d:%02d:%02d" % (time.localtime().tm_hour, time.localtime().tm_min, time.localtime().tm_sec)

def duration(t1, t2):
  diff = t2 - t1
  seconds = int(math.floor(diff))
  frac = diff - seconds
  hh = seconds // 3600
  mm = seconds // 60
  if hh > 0:
    mm = mm % 60
  ss = seconds % 60
  return "%02dh %02dm %fs" % (hh, mm, ss + frac)

