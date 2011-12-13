# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# XXX Module summary here
#
# XXX Zap excessive documentation
#
# $Id$

# XXX Move imports into functions?
import glob
import math
import numpy
import os
import re


class Section(object):
  """Class for a section, or a pair of ASIC:s, or a two-by-one.  A
  quadrant (and to some extent a section) is really an imaginary
  object--the detector only reads out the ASIC:s.  The metrology
  should be accurate to +/- one pixel, usually better but occasionally
  worse.  XXX Is it a section, a two-by-one, or a sensor?

  This class adopts a matrix-oriented coordinate system.  The origin
  is in the top left corner, the first coordinate increases downwards,
  the second coordinate increases to the right, and the third
  coordinate increases towards the viewer.  In this right-handed
  coordinate system, a rotation by a positive angle is
  counter-clockwise in the plane of the two first coordinates.

  @note There are a few hard coded numbers throughout this class.
  """


  # The size of a quadrant, in pixels.  Mikhail S. Dubrovin
  # empirically found that 850-by-850 pixels are enough to accommodate
  # any possibly section misalignment.
  q_size = (850, 850)


  def __init__(self, angle = 0, center = (0, 0)):
    """By default, an untransformed section is standing up, with its
    long side along the first (vertical) coordinate.  The two
    194-by-185 pixel ASIC:s of the section are separated by a
    three-pixel gap.

    XXX This is the metrology convention.  The data in the XTC stream
    is rotated by 90 degrees (transposed).
    """

    self.angle  = angle
    self.center = center
    self.size   = (2 * 194 + 3, 185)


  def corners(self, right = True):
    """The corners() function returns an array of the four corners of
    the section, in counter-clockwise order.  Each vertex is a
    two-dimensional array of the plane components.  XXX Maybe better
    named corners_section()?

    @param right @c True to restrict rotations to right angles
    @return      Coordinates of the four corners, in counter-clockwise
                 order
    """

    # The coordinates of the corners of the untransformed section, in
    # counter-clockwise order starting at the upper, left corner.
    coords = [[-0.5 * self.size[0], -0.5 * self.size[1]],
              [-0.5 * self.size[0], +0.5 * self.size[1]],
              [+0.5 * self.size[0], +0.5 * self.size[1]],
              [+0.5 * self.size[0], -0.5 * self.size[1]]]

    # Determine the cosine and sine of the rotation angle, rounded to
    # a multiple of 90 degrees if appropriate.
    if (right):
      a = math.radians(90.0 * round(self.angle / 90.0))
    else:
      a = math.radians(self.angle)
    c = math.cos(a)
    s = math.sin(a)

    # Apply plane rotation, and translation.
    for i in xrange(len(coords)):
      p         = coords[i]
      coords[i] = [c * p[0] - s * p[1] + self.center[0],
                   s * p[0] + c * p[1] + self.center[1]]
    return (coords)


  def corners_asic(self):
    """The corners_asic() function returns a list of pixel indices
    which position the two ASIC:s on the detector, in order from top
    to bottom by the "standing up" convention.  Each list of indices
    gives the coordinates of the top, left and bottom, right corners
    of the section's ASIC:s in "spotfinder format".
    """
    a = int(round(self.angle / 90.0)) % 4
    c = self.corners(True)

    if (a == 0):
      # The section is "standing up", and the top left corner is given
      # by the first corner.
      ul0 = [int(round(c[a][0])), int(round(c[a][1]))]
      ul1 = [ul0[0] + 194 + 3, ul0[1]]
      dlr = [194, 185]
    elif (a == 2):
      # The section is "standing up", and the top left corner is given
      # by the third corner.
      ul1 = [int(round(c[a][0])), int(round(c[a][1]))]
      ul0 = [ul1[0] + 194 + 3, ul1[1]]
      dlr = [194, 185]
    elif (a == 1):
      # The section is "laying down", and the top left corner is given
      # by the second corner.
      ul0 = [int(round(c[a][0])), int(round(c[a][1]))]
      ul1 = [ul0[0], ul0[1] + 194 + 3]
      dlr = [185, 194]
    elif (a == 3):
      # The section is "laying down", and the top left corner is given
      # by the forth corner.
      ul1 = [int(round(c[a][0])), int(round(c[a][1]))]
      ul0 = [ul1[0], ul1[1] + 194 + 3]
      dlr = [185, 194]

    coords = [
      [ul0[0], ul0[1], ul0[0] + dlr[0], ul0[1] + dlr[1]],
      [ul1[0], ul1[1], ul1[0] + dlr[0], ul1[1] + dlr[1]]]
    return (coords)


  def qrotate(self, angle):
    """ The qrotate() function rotates the section counter-clockwise
    by @p angle degrees around the centre of its quadrant.  The
    rotation angle is rounded to an integer multiple of 90 degrees
    prior to transformation.  Rotation around the quadrant centre
    changes the location and the orientation of the section.

    @param angle Rotation angle, in degrees
    """

    q = int(round(angle / 90.0)) % 4
    a = 90.0 * q

    if (q == 0):
      pass
    elif (q == 1):
      self.srotate(a)
      self.center = (self.q_size[1] - self.center[1],
                     0              + self.center[0])
    elif (q == 2):
      self.srotate(a)
      self.center = (self.q_size[0] - self.center[0],
                     self.q_size[1] - self.center[1])
    elif (q == 3):
      self.srotate(a)
      self.center = (0              + self.center[1],
                     self.q_size[0] - self.center[0])


  def srotate(self, angle):
    """The srotate() function rotates the section counter-clockwise by
    @p angle degrees around its centre.  Rotation within the quadrant
    only changes the orientation of the section.

    @param angle Rotation angle, in degrees
    """

    self.angle = self.angle + angle


  def translate(self, displacement):
    """The translate() function displaces the section.

    @param displacement Two-dimensional array of the additive
                        displacement
    """

    self.center = (self.center[0] + displacement[0],
                   self.center[1] + displacement[1])


def fread_vector(stream):
  """The fread_vector() function reads a vector from the stream
  pointed to by @p stream and returns it as a numpy array.

  @param stream Stream object
  @return       Tensor as numpy array
  """

  return (numpy.array(
      [float(d) for d in re.split("\s+", stream.readline()) if len(d) > 0]))


def fread_matrix(stream):
  """The fread_matrix() function reads a vector or matrix from the
  stream pointed to by @p stream and returns it as a numpy array.

  @param stream Stream object
  @return       Tensor as numpy array
  """

  A = fread_vector(stream)
  while (True):
    v = fread_vector(stream)
    if (v.shape[0] == 0):
      return (A)
    A = numpy.vstack((A, v))


def fread_tensor3(stream):
  """The fread_tensor3() function reads a tensor of rank no greater
  than 3 from the stream pointed to by @p stream and returns it as a
  numpy array.

  @param stream Stream object
  @return       Tensor as numpy array
  """

  T = fread_matrix(stream)
  while (True):
    A = fread_matrix(stream)
    if (len(A.shape) < 2 or A.shape[0] == 0 or A.shape[1] == 0):
      return (T)
    T = numpy.dstack((T, A))


def calib2tensor3(dirname, component):
  """The calib2tensor3() function reads the latest calibration tensor
  for @p component from @p dirname.  Any obsoleted calibration data is
  ignored.  The function returns the tensor as a numpy array on
  successful completion.

  @param dirname   Directory with calibration information
  @param component Kind of calibration data sought
  @return          Tensor as numpy array
  """

  basename = "*-end.data"
  path     = glob.glob(os.path.join(dirname, component, basename))[-1]
  stream   = open(path, "r")
  T        = fread_tensor3(stream)
  stream.close()
  return (T)


def calib2sections(dirname):
  """The calib2sections() function reads calibration information
  stored in the directory tree beneath @p dirname and returns a
  two-dimensional array of Section objects.  The first index in the
  returned array identifies the quadrant, and the second index
  identifies the section within the quadrant.

  @param dirname Directory with calibration information
  @return        Section objects
  """

  # Get centres of the sections, and apply corrections.
  s_cen = calib2tensor3(dirname, "center") \
      +   calib2tensor3(dirname, "center_corr")

  # Get the rotation of sections, and apply corrections.  Note that
  # sections 0, 1 and 4, 5 are antiparallel!
  s_rot = calib2tensor3(dirname, "rotation") \
      +   calib2tensor3(dirname, "tilt")

  # Get the margin, gap, and shift adjustments of the sections within
  # each quadrant.
  s_mgs = calib2tensor3(dirname, "marg_gap_shift")

  # Get offsets of the quadrants, and apply corrections.
  q_off = calib2tensor3(dirname, "offset") \
      +   calib2tensor3(dirname, "offset_corr")

  # Get rotation of the quadrants, and apply corrections.
  q_rot = calib2tensor3(dirname, "quad_rotation") \
      +   calib2tensor3(dirname, "quad_tilt")

  # The third coordinate is ignored for now, even though optical
  # measurement gives a variation in Z up to 0.6 mm.
  sections = []
  for q in xrange(s_cen.shape[0]):
    sections.append([])
    for s in xrange(s_cen.shape[1]):
      sec = Section()
      sec.translate((s_mgs[0, 0] + s_cen[q, s, 0],
                     s_mgs[1, 0] + s_cen[q, s, 1]))
      sec.srotate(s_rot[q, s])
      sec.qrotate(q_rot[q])
      sec.translate((s_mgs[0, 1] + q_off[0, q],
                     s_mgs[1, 1] + q_off[1, q]))

      # XXX I still don't understand this bit!
      if (q == 0):
        sec.translate((-s_mgs[0, 2] + s_mgs[0, 3],
                       -s_mgs[1, 2] - s_mgs[1, 3]))
      elif (q == 1):
        sec.translate((-s_mgs[0, 2] - s_mgs[0, 3],
                       +s_mgs[1, 2] - s_mgs[1, 3]))
      elif (q == 2):
        sec.translate((+s_mgs[0, 2] - s_mgs[0, 3],
                       +s_mgs[1, 2] + s_mgs[1, 3]))
      elif (q == 3):
        sec.translate((+s_mgs[0, 2] + s_mgs[0, 3],
                       -s_mgs[1, 2] + s_mgs[1, 3]))

      sections[q].append(sec)
  return (sections)
