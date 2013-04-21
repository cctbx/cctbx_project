# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id: metrology.py 299 2012-04-04 18:46:56Z hattne $

"""The metrology module does conversion of this and that.  XXX Need to
describe the keys of the dictionaries here.  XXX Just like in OpenGL,
all matrices are to be pre-multiplied (with column vectors)?
"""
from __future__ import division

from scitbx import matrix


def _transform(o, t):
  """The _transform() function returns the transformation matrices in
  homogeneous coordinates between the parent and child frames.  The
  forward transform maps coordinates in the parent frame to the child
  frame, and the backward transform provides the inverse.  The last
  row of the product of any two homogeneous transformation matrices is
  always (0, 0, 0, 1).

  @param o Orientation of child w.r.t. parent, as a unit quaternion
  @param t Translation of child w.r.t. parent
  @return  Two-tuple of the forward and backward transformation
           matrices
  """

  Rb = o.unit_quaternion_as_r3_rotation_matrix()
  tb = t
  Tb = matrix.sqr((Rb(0, 0), Rb(0, 1), Rb(0, 2), tb(0, 0),
                   Rb(1, 0), Rb(1, 1), Rb(1, 2), tb(1, 0),
                   Rb(2, 0), Rb(2, 1), Rb(2, 2), tb(2, 0),
                   0,        0,        0,        1))
  Rf = Rb.transpose()
  tf = -Rf * t
  Tf = matrix.sqr((Rf(0, 0), Rf(0, 1), Rf(0, 2), tf(0, 0),
                   Rf(1, 0), Rf(1, 1), Rf(1, 2), tf(1, 0),
                   Rf(2, 0), Rf(2, 1), Rf(2, 2), tf(2, 0),
                   0,        0,        0,        1))
  return (Tf, Tb)


def get_projection_matrix(dim_pixel, dim_readout):
  """The get_projection_matrix() function returns the projection
  matrices in homogeneous coordinates between readout frame and
  fractional row and column indices.  The forward transform maps
  coordinates in the readout frame, i.e. [x, y, z, 1], to readout
  indices, i.e. [i, j, 1], where i and j are the row and column
  indices, respectively.  The backward transform provides the inverse.
  The get_projection_matrix() function assumes that coordinates in the
  readout frame are pixel centers.  XXX Bad name!  Better
  readout2metric(), readout_projection_matrices()?

  @note If coordinates in the readout frame are the top, left corner
        of the pixel, pass dim_readout=(readout_width_in_pixels + 1,
        readout_height_in_pixels + 1).

  @param dim_pixel   Two-tuple of pixel width and height, in meters
  @param dim_readout Two-tuple of readout width and height, in pixels
  @return            Two-tuple of the forward and backward projection
                     matrices
  """

  Pb = matrix.rec(
    elems=(0, +dim_pixel[0], dim_pixel[0] * (1 - dim_readout[0]) / 2,
           -dim_pixel[1], 0, dim_pixel[1] * (dim_readout[1] - 1) / 2,
           0, 0, 0,
           0, 0, 1),
    n=(4, 3))
  Pf = matrix.rec(
    elems=(0, -1 / dim_pixel[1], 0, (dim_readout[1] - 1) / 2,
           +1 / dim_pixel[0], 0, 0, (dim_readout[0] - 1) / 2,
           0, 0, 0, 1),
    n=(3, 4))
  return (Pf, Pb)


def metrology_as_transformation_matrices(params):
  """XXX What if the child is invisible from the parent (pointing away
  from the parent)?  XXX What if two elements are overlapping?
  Implement bounding boxes (only needed for transformation/projection
  from eye to ASIC).  XXX Projections implemented elsewhere.

  XXX Should probably be something like "phil_to_matrix_dict" and put
  the metrology2phil() function here, too!

  XXX It was not such a great idea to include the final projection
  here--it made it impossible to recover z!

  XXX Normalisation of orientation may be necessary to guarantee that
  rotation matrices are orthogonal if the serialisation truncated the
  numbers.

  @param params Pure Python object, extracted from a metrology Phil
                object
  @return       Dictionary of homogeneous transformation matrices
  """

  d = params.detector
  o_d = matrix.col(d.orientation)
  t_d = matrix.col(d.translation)
  T_d = _transform(o_d, t_d)

  matrices = {(0,): T_d}

  for p in d.panel:
    o_p = matrix.col(p.orientation).normalize()
    t_p = matrix.col(p.translation)
    T_p = _transform(o_p, t_p)
    T_p = (T_p[0] * T_d[0], T_d[1] * T_p[1])

    matrices[(d.serial, p.serial)] = T_p

    for s in p.sensor:
      o_s = matrix.col(s.orientation).normalize()
      t_s = matrix.col(s.translation)
      T_s = _transform(o_s, t_s)
      T_s = (T_s[0] * T_p[0], T_p[1] * T_s[1])

      matrices[(d.serial, p.serial, s.serial)] = T_s

      for a in s.asic:
        o_a = matrix.col(a.orientation).normalize()
        t_a = matrix.col(a.translation)
        T_a = _transform(o_a, t_a)
        T_a = (T_a[0] * T_s[0], T_s[1] * T_a[1])

        matrices[(d.serial, p.serial, s.serial, a.serial)] = (T_a[0], T_a[1])

  return matrices
