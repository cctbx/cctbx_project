from __future__ import division
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

def cbf_gonio_to_effective_axis_fixed(cbf_gonio):
  '''Given a cbf goniometer handle, first determine the real rotation
  axis, then determine the fixed component of rotation which is rotated
  about this axis.'''

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

  axis = matrix.col(rot.axis)

  S = axis.axis_and_angle_as_r3_rotation_matrix(start, deg = True)

  return axis, S.inverse() * R
