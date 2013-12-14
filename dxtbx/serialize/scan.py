#!/usr/bin/env python
#
# dxtbx.serialize.scan.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

def to_dict(scan):
  ''' Convert the scan model to a dictionary

  Params:
      scan The scan model

  Returns:
      A dictionary of the parameters

  '''
  if scan == None:
    return None
  return scan.to_dict()

def from_dict(d, t=None):
  ''' Convert the dictionary to a scan model

  Params:
      d The dictionary of parameters
      t The template dictionary to use

  Returns:
      The scan model

  '''
  from dxtbx.model import Scan
  from scitbx.array_family import flex # import dependency

  # If None, return None
  if d == None:
    if t == None: return None
    else: return from_dict(t, None)
  elif t != None:
    d = dict(t.items() + d.items())
  if not isinstance(d['exposure_time'], list):
    d['exposure_time'] = [d['exposure_time']]

  # Create the model from the dictionary
  return Scan.from_dict(d)
