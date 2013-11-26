from __future__ import division
#!/usr/bin/env python
#
# dxtbx.serialize.goniometer.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

def to_dict(gonio):
  ''' Convert the goniometer model to a dictionary

  Params:
      gonio The goniometer model

  Returns:
      A dictionary of the parameters

  '''
  if gonio == None:
    return None
  return gonio.to_dict()

def from_dict(d, t=None):
  ''' Convert the dictionary to a goniometer model

  Params:
      d The dictionary of parameters
      t The template dictionary to use

  Returns:
      The goniometer model

  '''
  from dxtbx.model import Goniometer

  # If None, return None
  if d == None:
    if t == None: return None
    else: return from_dict(t, None)
  elif t != None:
    d = dict(t.items() + d.items())

  # Create the model from the dictionary
  return Goniometer.from_dict(d)
