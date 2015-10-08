from __future__ import division
#!/usr/bin/env python
# test_goniometer.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Tests for the goniometer class.

import math

from dxtbx.model.goniometer import Goniometer
from dxtbx.model.goniometer import goniometer_factory

def compare_tuples(a, b, tol = 1.0e-6):

  assert(len(a) == len(b))

  for j in range(len(a)):
    if math.fabs(b[j] - a[j]) > tol:
      return False

  return True

def test_goniometer():
  '''A test class for the goniometer class.'''

  from libtbx import easy_pickle

  axis = (1, 0, 0)
  fixed = (1, 0, 0, 0, 1, 0, 0, 0, 1)

  xg = Goniometer(axis, fixed)

  assert(len(xg.get_rotation_axis()) == 3)
  assert(len(xg.get_fixed_rotation()) == 9)

  assert(compare_tuples(xg.get_rotation_axis(), axis))
  assert(compare_tuples(xg.get_fixed_rotation(), fixed))

  single = goniometer_factory.single_axis()

  assert(len(single.get_rotation_axis()) == 3)
  assert(len(single.get_fixed_rotation()) == 9)

  assert(compare_tuples(single.get_rotation_axis(), axis))
  assert(compare_tuples(single.get_fixed_rotation(), fixed))

  kappa = goniometer_factory.kappa(50.0, 0.0, 0.0, 0.0, '-y', 'omega')

  assert(len(kappa.get_rotation_axis()) == 3)
  assert(len(kappa.get_fixed_rotation()) == 9)
  assert(compare_tuples(kappa.get_rotation_axis(), axis))
  assert(compare_tuples(kappa.get_fixed_rotation(), fixed))

  kappa = goniometer_factory.kappa(50.0, 0.0, 0.0, 0.0, '-y', 'omega')

  assert(len(kappa.get_rotation_axis()) == 3)
  assert(len(kappa.get_fixed_rotation()) == 9)

  assert(compare_tuples(kappa.get_rotation_axis(), axis))
  assert(compare_tuples(kappa.get_fixed_rotation(), fixed))

  kappa = goniometer_factory.kappa(50.0, 0.0, 0.0, 0.0, '-y', 'phi')

  assert(len(kappa.get_rotation_axis()) == 3)
  assert(len(kappa.get_fixed_rotation()) == 9)

  assert(compare_tuples(kappa.get_rotation_axis(), axis))
  assert(compare_tuples(kappa.get_fixed_rotation(), fixed))

  kappa = goniometer_factory.kappa(50.0, 0.0, 30.0, 0.0, '-y', 'omega')

  assert(len(kappa.get_rotation_axis()) == 3)
  assert(len(kappa.get_fixed_rotation()) == 9)

  assert(compare_tuples(kappa.get_rotation_axis(), axis))
  assert(not compare_tuples(kappa.get_fixed_rotation(), fixed))

  import libtbx.load_env
  import os

  dxtbx_dir = libtbx.env.dist_path('dxtbx')

  image = os.path.join(dxtbx_dir, 'tests', 'phi_scan_001.cbf')
  cbf = goniometer_factory.imgCIF(image)

  kappa = goniometer_factory.kappa(50.0, -10.0, 30.0, 0.0, '-y', 'phi')

  s = easy_pickle.dumps(kappa)
  kappa2 = easy_pickle.loads(s)
  assert kappa == kappa2

  image = os.path.join(dxtbx_dir, 'tests', 'omega_scan.cbf')
  cbf = goniometer_factory.imgCIF(image)

  kappa = goniometer_factory.kappa(50.0, -10.0, 30.0, 20.0, '-y', 'omega')

  s = easy_pickle.dumps(kappa)
  kappa2 = easy_pickle.loads(s)
  assert kappa == kappa2

  print 'OK'

if __name__ == '__main__':

  test_goniometer()
