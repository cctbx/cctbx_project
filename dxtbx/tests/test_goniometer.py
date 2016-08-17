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

from dxtbx.model.goniometer import Goniometer, MultiAxisGoniometer
from dxtbx.model.goniometer import goniometer_factory
from libtbx import easy_pickle
from libtbx.test_utils import Exception_expected

def compare_tuples(a, b, tol = 1.0e-6):

  assert(len(a) == len(b))

  for j in range(len(a)):
    if math.fabs(b[j] - a[j]) > tol:
      return False

  return True

def test_goniometer():
  '''A test class for the goniometer class.'''


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

def test_multi_axis_goniometer():
  from libtbx.test_utils import approx_equal
  from scitbx.array_family import flex

  alpha = 50
  omega = -10
  kappa = 30
  phi = 20
  direction = '-y'

  from dxtbx.model.goniometer import KappaGoniometer
  kappa_omega_scan = KappaGoniometer(
    alpha, omega, kappa, phi, direction, 'omega')
  axes = (kappa_omega_scan.get_phi_axis(), kappa_omega_scan.get_kappa_axis(),
          kappa_omega_scan.get_omega_axis())
  angles = (kappa_omega_scan.get_phi_angle(), kappa_omega_scan.get_kappa_angle(),
            kappa_omega_scan.get_omega_angle())

  # First test a kappa goniometer with omega as the scan axis
  axes = flex.vec3_double(axes)
  angles = flex.double(angles)
  names = flex.std_string(('phi', 'kappa', 'omega'))
  scan_axis = 2

  multi_axis_omega_scan = goniometer_factory.multi_axis(
    axes, angles, names, scan_axis)
  assert approx_equal(multi_axis_omega_scan.get_fixed_rotation(), kappa_omega_scan.get_fixed_rotation())
  assert approx_equal(multi_axis_omega_scan.get_rotation_axis(), kappa_omega_scan.get_rotation_axis())

  recycle_omega = MultiAxisGoniometer.from_dict(multi_axis_omega_scan.to_dict())
  assert approx_equal(recycle_omega.get_axes(), multi_axis_omega_scan.get_axes())
  assert approx_equal(recycle_omega.get_angles(), multi_axis_omega_scan.get_angles())
  assert recycle_omega.get_scan_axis() == multi_axis_omega_scan.get_scan_axis()

  # Now test a kappa goniometer with phi as the scan axis
  kappa_phi_scan = KappaGoniometer(
    alpha, omega, kappa, phi, direction, 'phi')

  scan_axis = 0
  multi_axis_phi_scan = goniometer_factory.multi_axis(
    axes, angles, names, scan_axis)
  assert approx_equal(multi_axis_phi_scan.get_fixed_rotation(), kappa_phi_scan.get_fixed_rotation())
  from scitbx import matrix
  assert approx_equal(matrix.sqr(multi_axis_phi_scan.get_setting_rotation()) * multi_axis_phi_scan.get_rotation_axis(), kappa_phi_scan.get_rotation_axis())

  recycle_phi = MultiAxisGoniometer.from_dict(multi_axis_phi_scan.to_dict())
  assert approx_equal(recycle_phi.get_axes(), multi_axis_phi_scan.get_axes())
  assert approx_equal(recycle_phi.get_angles(), multi_axis_phi_scan.get_angles())
  assert recycle_phi.get_scan_axis() == multi_axis_phi_scan.get_scan_axis()

  s = easy_pickle.dumps(multi_axis_phi_scan)
  recycle = easy_pickle.loads(s)
  assert recycle == multi_axis_phi_scan

  assert approx_equal(recycle.get_axes(), multi_axis_phi_scan.get_axes())
  assert approx_equal(recycle.get_angles(), multi_axis_phi_scan.get_angles())
  assert recycle.get_scan_axis() == multi_axis_phi_scan.get_scan_axis()

  # Check exception is raised if scan axis is out range
  try: goniometer_factory.multi_axis(axes, angles, names, 3)
  except RuntimeError, e: pass
  else: raise Exception_expected

  # Single axis is just a special case of a multi axis goniometer
  single_axis = goniometer_factory.multi_axis(
    flex.vec3_double(((1,0,0),)), flex.double((0,)), flex.std_string(('PHI',)), 0)
  assert single_axis.get_fixed_rotation() == (1,0,0,0,1,0,0,0,1)
  assert single_axis.get_setting_rotation() == (1,0,0,0,1,0,0,0,1)
  assert single_axis.get_rotation_axis() == (1,0,0)

  print 'OK'

if __name__ == '__main__':

  test_goniometer()
  test_multi_axis_goniometer()
