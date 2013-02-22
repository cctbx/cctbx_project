#!/usr/bin/env python
# test_goniometer.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Tests for the goniometer class.

import math
import os
import sys

from dxtbx.model.goniometer import goniometer
from dxtbx.model.goniometer import goniometer_factory

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

    xg = goniometer(axis, fixed)

    assert(len(xg.get_axis()) == 3)
    assert(len(xg.get_fixed()) == 9)

    assert(compare_tuples(xg.get_axis(), axis))
    assert(compare_tuples(xg.get_fixed(), fixed))

    single = goniometer_factory.single_axis()

    assert(len(single.get_axis()) == 3)
    assert(len(single.get_fixed()) == 9)

    assert(compare_tuples(single.get_axis(), axis))
    assert(compare_tuples(single.get_fixed(), fixed))

    kappa = goniometer_factory.kappa(50.0, 0.0, 0.0, 0.0, '-y', 'omega')

    assert(len(single.get_axis()) == 3)
    assert(len(single.get_fixed()) == 9)

    assert(compare_tuples(kappa.get_axis(), axis))
    assert(compare_tuples(kappa.get_fixed(), fixed))

    kappa = goniometer_factory.kappa(50.0, 0.0, 0.0, 0.0, '-y', 'omega')

    assert(len(single.get_axis()) == 3)
    assert(len(single.get_fixed()) == 9)

    assert(compare_tuples(kappa.get_axis(), axis))
    assert(compare_tuples(kappa.get_fixed(), fixed))

    kappa = goniometer_factory.kappa(50.0, 0.0, 0.0, 0.0, '-y', 'phi')

    assert(len(single.get_axis()) == 3)
    assert(len(single.get_fixed()) == 9)

    assert(compare_tuples(kappa.get_axis(), axis))
    assert(compare_tuples(kappa.get_fixed(), fixed))

    kappa = goniometer_factory.kappa(50.0, 0.0, 30.0, 0.0, '-y', 'omega')

    assert(len(single.get_axis()) == 3)
    assert(len(single.get_fixed()) == 9)

    assert(compare_tuples(kappa.get_axis(), axis))
    assert(not compare_tuples(kappa.get_fixed(), fixed))

    cbf = goniometer_factory.imgCIF('phi_scan_001.cbf')

    print cbf

    kappa = goniometer_factory.kappa(50.0, -10.0, 30.0, 0.0, '-y', 'phi')

    print kappa

    cbf = goniometer_factory.imgCIF('omega_scan.cbf')

    print cbf

    kappa = goniometer_factory.kappa(50.0, -10.0, 30.0, 20.0, '-y', 'omega')

    print kappa


if __name__ == '__main__':

    test_goniometer()
