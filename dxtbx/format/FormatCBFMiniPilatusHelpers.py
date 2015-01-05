from __future__ import division
#!/usr/bin/env python
# FormatCBFMiniPilatus.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Helpers for FormatCBFMiniPilatus...

import time

def get_pilatus_timestamp(timestamp_string):
  import calendar
  timestamp, milliseconds = timestamp_string.split('.')

  for format in ['%Y-%b-%dT%H:%M:%S',
                 '%Y-%m-%dT%H:%M:%S',
                 '%Y/%b/%d %H:%M:%S']:

    try:
      struct_time = time.strptime(timestamp, format)
      return calendar.timegm(struct_time) + float('0.' + milliseconds)

    except: # intentional
      pass

  raise RuntimeError, 'timestamp %s not recognised' % timestamp
