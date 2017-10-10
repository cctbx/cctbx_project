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
from __future__ import absolute_import, division
from builtins import map
from functools import reduce

def tuple_almost_equal(a, b, eps=1e-7):
  ''' Check if the tuples are equal. '''
  return reduce(lambda a,b: a and b, list(map(lambda a,b: abs(a - b) < eps, a, b)))
