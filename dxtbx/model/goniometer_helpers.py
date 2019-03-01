from __future__ import absolute_import, division

#!/usr/bin/env python
# goniometer_helpers.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Helper functions for goniometer

from scitbx import matrix
from scitbx.math import r3_rotation_axis_and_angle_from_matrix

# this function no longer appears to be needed - the cbflib call now gives the
# correct rotation axis


def cbf_gonio_to_effective_axis_fixed(cbf_gonio):
    axis = matrix.col(cbf_gonio.get_rotation_axis())
    start, increment = cbf_gonio.get_rotation_range()

    # xia2-56 handle gracefully reverse turning goniometers - this assumes
    # that the angles will be correctly inverted in the scan factory
    if increment < 0:
        start = -start
        increment = -increment
        axis = -axis

    x = cbf_gonio.rotate_vector(0.0, 1, 0, 0)
    y = cbf_gonio.rotate_vector(0.0, 0, 1, 0)
    z = cbf_gonio.rotate_vector(0.0, 0, 0, 1)
    R = matrix.rec(x + y + z, (3, 3)).transpose()
    S = axis.axis_and_angle_as_r3_rotation_matrix(start, deg=True)

    fixed = S.inverse() * R
    return axis, fixed


def cbf_gonio_to_effective_axis_fixed_old(cbf_gonio):
    """Given a cbf goniometer handle, first determine the real rotation
  axis, then determine the fixed component of rotation which is rotated
  about this axis."""

    # First construct the real rotation axis, as the difference in rotating
    # the identity matrix at the end of the scan and the beginning.

    x = cbf_gonio.rotate_vector(0.0, 1, 0, 0)
    y = cbf_gonio.rotate_vector(0.0, 0, 1, 0)
    z = cbf_gonio.rotate_vector(0.0, 0, 0, 1)

    R = matrix.rec(x + y + z, (3, 3)).transpose()

    x1 = cbf_gonio.rotate_vector(1.0, 1, 0, 0)
    y1 = cbf_gonio.rotate_vector(1.0, 0, 1, 0)
    z1 = cbf_gonio.rotate_vector(1.0, 0, 0, 1)

    R1 = matrix.rec(x1 + y1 + z1, (3, 3)).transpose()

    RA = R1 * R.inverse()

    rot = r3_rotation_axis_and_angle_from_matrix(RA)

    # Then, given this, determine the component of the scan which is fixed -
    # which will need to be with respect to the unrotated axis. N.B. this
    # will not be unique, but should be correct modulo a free rotation about
    # the shifted axis.

    start = cbf_gonio.get_rotation_range()[0]

    # want positive rotations => if negative invert axis
    axis = matrix.col(rot.axis)
    angle = rot.angle()
    if angle < 0:
        axis = -1 * axis
        # common sense would suggest in here that if the angle is -ve should
        # be made +ve - works OK for omega scans but not phi scans, probably
        # incomplete goniometer definition problem...
        # start = -start

    S = axis.axis_and_angle_as_r3_rotation_matrix(start, deg=True)

    return axis, S.inverse() * R
