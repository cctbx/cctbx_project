#!/usr/bin/env python
# FormatSMVRigakuSaturnSN09040159.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for Rigaku Saturn images, for
# the instrument at CSHL, SN 09040159.

from __future__ import division

from scitbx import matrix

from dxtbx.format.FormatSMVRigakuSaturn import FormatSMVRigakuSaturn

class FormatSMVRigakuSaturn09040159(FormatSMVRigakuSaturn):
  '''A class for reading SMV format Rigaku Saturn images, and correctly
  constructing a model for the experiment from this. This is the custom
  version for instrument with serial number 09040159 at CSHL.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like a Rigaku Saturn SMV format image,
    i.e. we can make sense of it. Essentially that will be if it contains
    all of the keys we are looking for. Also checks the serial number.'''

    size, header = FormatSMVRigakuSaturn.get_smv_header(image_file)

    detector_prefix = header['DETECTOR_NAMES'].split()[0].strip()
    try:
      serial_number = header['%sSERIAL_NUMBER' % detector_prefix]
    except KeyError:
      return False

    if serial_number != '09040159':
      return False

    return True

  def _detector(self):
    '''Return a model for the detector, allowing for two-theta offsets
    and the detector position. This will be rather more complex... and
    overloads the definition for the general Rigaku Saturn detector by
    rotating the detector by -90 degrees about the beam.'''

    detector_name = self._header_dictionary[
        'DETECTOR_NAMES'].split()[0].strip()

    detector_axes = map(float, self._header_dictionary[
        '%sDETECTOR_VECTORS' % detector_name].split())

    R = matrix.col((0, 0, 1)).axis_and_angle_as_r3_rotation_matrix(
        -90, deg = True)

    detector_fast = R * matrix.col(tuple(detector_axes[:3]))
    detector_slow = R * matrix.col(tuple(detector_axes[3:]))

    beam_pixels = map(float, self._header_dictionary[
        '%sSPATIAL_DISTORTION_INFO' % detector_name].split()[:2])
    pixel_size = map(float, self._header_dictionary[
        '%sSPATIAL_DISTORTION_INFO' % detector_name].split()[2:])
    image_size = map(int, self._header_dictionary[
        '%sDETECTOR_DIMENSIONS' % detector_name].split())

    detector_origin = - (beam_pixels[0] * pixel_size[0] * detector_fast + \
                         beam_pixels[1] * pixel_size[1] * detector_slow)

    gonio_axes = map(float, self._header_dictionary[
        '%sGONIO_VECTORS' % detector_name].split())
    gonio_values = map(float, self._header_dictionary[
        '%sGONIO_VALUES' % detector_name].split())
    gonio_units = self._header_dictionary[
        '%sGONIO_UNITS' % detector_name].split()
    gonio_num_axes = int(self._header_dictionary[
        '%sGONIO_NUM_VALUES' % detector_name])

    rotations = []
    translations = []

    for j, unit in enumerate(gonio_units):
      axis = matrix.col(gonio_axes[3 * j:3 * (j + 1)])
      if unit == 'deg':
        rotations.append(axis.axis_and_angle_as_r3_rotation_matrix(
            gonio_values[j], deg = True))
        translations.append(matrix.col((0.0, 0.0, 0.0)))
      elif unit == 'mm':
        rotations.append(matrix.sqr((1.0, 0.0, 0.0,
                                     0.0, 1.0, 0.0,
                                     0.0, 0.0, 1.0)))
        translations.append(gonio_values[j] * axis)
      else:
        raise RuntimeError, 'unknown axis unit %s' % unit

    rotations.reverse()
    translations.reverse()

    for j in range(gonio_num_axes):
      detector_fast = rotations[j] * detector_fast
      detector_slow = rotations[j] * detector_slow
      detector_origin = rotations[j] * detector_origin
      detector_origin = translations[j] + detector_origin

    overload = int(self._header_dictionary['SATURATED_VALUE'])
    underload = 0

    return self._detector_factory.complex(
        'CCD', detector_origin.elems, detector_fast.elems,
        detector_slow.elems, pixel_size, image_size, (underload, overload))

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatSMVRigakuSaturnSN07400090.understand(arg)
