from __future__ import division
#!/usr/bin/env python
#
# dxtbx.serialize.detector.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def to_dict(detector):
  ''' Convert the detector model to a dictionary

  Params:
      detector The detector model

  Returns:
      A dictionary of the parameters

  '''
  if detector == None:
    return None
  return detector.to_dict()

def from_dict(d, t=None):
  ''' Convert the dictionary to a detector model

  Params:
      d The dictionary of parameters
      t The template dictionary to use

  Returns:
      The detector model

  '''
  from dxtbx.model import Detector, HierarchicalDetector
  from dxtbx.array_family import flex # import dependency

  # If None, return None
  if d == None:
    if t == None: return None
    else: return from_dict(t, None)
  elif t != None:
    if isinstance(d, list):
      d = { 'panels' : d }
    d2 = dict(t.items() + d.items())
  else:
    if isinstance(d, list):
      d = { 'panels' : d }

  # Create the model from the dictionary
  if "hierarchy" in d:
    return HierarchicalDetector.from_dict(d)
  else:
    return Detector.from_dict(d)
