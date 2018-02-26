from __future__ import absolute_import, division
#!/usr/bin/env python
# detector_helpers.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Helpers for the detector class...

import math
from scitbx import matrix

def read_xds_xparm(xds_xparm_file):
  '''Parse the XDS XPARM file, which contains a description of the detector
  and experimental geometry, to a dictionary.'''

  data = map(float, open(xds_xparm_file, 'r').read().split())

  assert(len(data) == 42)

  starting_frame = int(data[0])
  phi_start, phi_width = data[1:3]
  axis = data[3:6]

  wavelength = data[6]
  beam = data[7:10]

  nx, ny = map(int, data[10:12])
  px, py = data[12:14]

  distance = data[14]
  ox, oy = data[15:17]

  x, y = data[17:20], data[20:23]
  normal = data[23:26]

  spacegroup = int(data[26])
  cell = data[27:33]

  a, b, c = data[33:36], data[36:39], data[39:42]

  results = {
      'starting_frame':starting_frame,
      'phi_start':phi_start, 'phi_width':phi_width,
      'axis':axis, 'wavelength':wavelength, 'beam':beam,
      'nx':nx, 'ny':ny, 'px':px, 'py':py, 'distance':distance,
      'ox':ox, 'oy':oy, 'x':x, 'y':y, 'normal':normal,
      'spacegroup':spacegroup, 'cell':cell, 'a':a, 'b':b, 'c':c
      }

  return results

def compute_frame_rotation(original, final):
  '''Compute reference frame rotation to rotate from the original frame
  given by original = (x, y, z) to the to reference frame given by
  final = (_x, _y, _z). Returns M where M.x = _x etc.'''

  x, y, z = original
  _x, _y, _z = final

  O = matrix.sqr(x.elems + y.elems + z.elems).transpose()
  assert((O.determinant() - 1.0) < 1.0e-7)

  F = matrix.sqr(_x.elems + _y.elems + _z.elems).transpose()
  assert((F.determinant() - 1.0) < 1.0e-7)

  # #1 rotate about x ^ (1, 0, 0) - if they are not coincident,
  # rotate about _x ^ _y if they are colinear but in opposite
  # directions

  if _x.angle(x) % math.pi:
    _ra_x = _x.cross(x)
    _a_x = _x.angle(x)
  elif math.fabs(_x.angle(x) - math.pi) < 1.0e-7:
    _ra_x = _x.cross(_y)
    _a_x = math.pi
  else:
    _ra_x = _x
    _a_x = 0.0

  _m_x = _ra_x.axis_and_angle_as_r3_rotation_matrix(- _a_x)

  # then rotate z to _z by rotating about _x (which is now coincident
  # with x)

  _ra_z = _x
  _a_z = _z.angle(_m_x * z)
  _m_z = _ra_z.axis_and_angle_as_r3_rotation_matrix(- _a_z)

  _m = _m_z * _m_x

  assert(math.fabs(_m.determinant() - 1.0) < 1.0e-7)

  return _m

def find_undefined_value(cbf_handle):
  '''Given a cbf handle, get the value for the undefined pixel.'''

  cbf_handle.find_category('array_intensities')
  cbf_handle.find_column('undefined_value')
  return cbf_handle.get_doublevalue()

class detector_helper_sensors:
  '''A helper class which allows enumeration of detector sensor technologies
  which should help in identifying specific detectors when needed. These are
  currently limited to IMAGE_PLATE CCD PAD.'''

  SENSOR_CCD = 'SENSOR_CCD'
  SENSOR_PAD = 'SENSOR_PAD'
  SENSOR_IMAGE_PLATE = 'SENSOR_IMAGE_PLATE'
  SENSOR_UNKNOWN = 'SENSOR_UNKNOWN'

  @staticmethod
  def check_sensor(sensor_type):
    if sensor_type in [detector_helper_sensors.SENSOR_CCD,
                       detector_helper_sensors.SENSOR_PAD,
                       detector_helper_sensors.SENSOR_IMAGE_PLATE,
                       detector_helper_sensors.SENSOR_UNKNOWN]:
      return True
    return False

  @staticmethod
  def all():
    return [detector_helper_sensors.SENSOR_CCD,
            detector_helper_sensors.SENSOR_PAD,
            detector_helper_sensors.SENSOR_IMAGE_PLATE]

def set_slow_fast_beam_centre_mm(detector, beam, beam_centre, panel_id=None):
  """ detector and beam are dxtbx objects,
      beam_centre is a tuple of (slow, fast) mm coordinates.
      supports 2-theta offset detectors, assumes correct centre provided
      for 2-theta=0
  """
  beam_s, beam_f = beam_centre

  # Ensure panel_id is set
  us0 = matrix.col(beam.get_unit_s0())
  if panel_id is None:
    panel_id = detector.get_panel_intersection(us0)
    if panel_id < 0:
      panel_id = detector.get_panel_intersection(-us0)
    if panel_id < 0:
      panel_id = 0

  # Get data from the chosen panel
  panel = detector[panel_id]
  f = matrix.col(panel.get_fast_axis())
  s = matrix.col(panel.get_slow_axis())
  n = matrix.col(panel.get_normal())
  o = matrix.col(panel.get_origin())

  # Attempt to find the axis an angle of an applied 2theta shift
  cos_angle = n.cos_angle(us0)
  if (cos_angle < 0):
    axi = us0.cross(-n)
    ang = us0.angle(-n)
  else:
    axi = us0.cross(n)
    ang = us0.angle(n)

  # Assume a 2theta offset if obliquity >= 5 deg
  two_theta = abs(ang) >= 5.0 * math.pi / 180.0

  # Undo 2theta shift
  if two_theta:
    R = axi.axis_and_angle_as_r3_rotation_matrix(ang)
    Rinv = R.inverse()
    try:
      h = detector.hierarchy()
      h.set_frame(fast_axis=Rinv*matrix.col(h.get_fast_axis()),
                  slow_axis=Rinv*matrix.col(h.get_slow_axis()),
                  origin=Rinv*matrix.col(h.get_origin()))
    except AttributeError:
      for p in detector:
        p.set_frame(fast_axis=Rinv*matrix.col(p.get_fast_axis()),
                    slow_axis=Rinv*matrix.col(p.get_slow_axis()),
                    origin=Rinv*matrix.col(p.get_origin()))

  # Lab coord of desired beam centre
  if us0.accute_angle(n, deg=True) > 89.9:
    raise RuntimeError("Beam is in the plane of the detector panel")
  beam_dist = panel.get_directed_distance() / us0.dot(n)
  beam_centre_lab = beam_dist * us0

  # Lab coord of the current position where we want the beam centre
  intersection_lab = matrix.col(panel.get_lab_coord((beam_f, beam_s)))

  # If the detector has a hierarchy, just update the root note
  try:
    h = detector.hierarchy()
    translation = beam_centre_lab - intersection_lab
    new_origin = matrix.col(h.get_origin()) + translation
    h.set_frame(fast_axis=h.get_fast_axis(),
                slow_axis=h.get_slow_axis(),
                origin=new_origin)
  except AttributeError:
    # No hierarchy, update each panel instead by finding the offset of
    # its origin from the current position of the desired beam centre. Use
    # this to reposition the panel origin wrt the final beam centre
    for p in detector:
      origin = matrix.col(p.get_origin())
      offset = origin - intersection_lab
      new_origin = beam_centre_lab + offset
      p.set_frame(fast_axis=p.get_fast_axis(),
                  slow_axis=p.get_slow_axis(),
                  origin=new_origin)

  # sanity check to make sure we have got the new beam centre correct
  new_beam_centre = detector[panel_id].get_bidirectional_ray_intersection(us0)
  assert (matrix.col(new_beam_centre) -
          matrix.col((beam_f, beam_s))).length() < 1e-4

  # Re-apply 2theta shift if required
  if two_theta:
    try:
      h = detector.hierarchy()
      h.set_frame(fast_axis=R*matrix.col(h.get_fast_axis()),
                  slow_axis=R*matrix.col(h.get_slow_axis()),
                  origin=R*matrix.col(h.get_origin()))
    except AttributeError:
      for p in detector:
        p.set_frame(fast_axis=R*matrix.col(p.get_fast_axis()),
                    slow_axis=R*matrix.col(p.get_slow_axis()),
                    origin=R*matrix.col(p.get_origin()))

  return

def set_mosflm_beam_centre(detector, beam, mosflm_beam_centre):
  """ detector and beam are dxtbx objects,
      mosflm_beam_centre is a tuple of mm coordinates.
      supports 2-theta offset detectors, assumes correct centre provided
      for 2-theta=0
  """
  return set_slow_fast_beam_centre_mm(detector, beam, mosflm_beam_centre)

def set_detector_distance(detector, distance):
  '''
  Set detector origin from distance along normal

  '''
  from scitbx import matrix
  assert len(detector) == 1
  normal = matrix.col(detector[0].get_normal())
  origin = matrix.col(detector[0].get_origin())
  d = origin.dot(normal)
  l = origin - d * normal
  origin = distance * normal + l
  fast_axis = detector[0].get_fast_axis()
  slow_axis = detector[0].get_slow_axis()
  detector[0].set_frame(fast_axis, slow_axis, origin)

