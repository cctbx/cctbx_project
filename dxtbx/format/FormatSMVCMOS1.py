#!/usr/bin/env python
# FormatSMVCMOS1.py
#   Copyright (C) 2015 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for CMOS1 images, from ALS 4.2.2

from __future__ import division

from dxtbx.format.FormatSMV import FormatSMV

class FormatSMVCMOS1(FormatSMV):
  '''A class for reading SMV format CMOS1 images. 'Abstract' class.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like a CMOS1 d*TREK SMV format image,
    i.e. we can make sense of it. Essentially that will be if it contains
    all of the keys we are looking for.'''

    size, header = FormatSMV.get_smv_header(image_file)

    wanted_header_items = [
        'DETECTOR_NUMBER', 'DETECTOR_NAMES',
        'BYTE_ORDER', 'DIM', 'SIZE1', 'SIZE2', 'Data_type'
    ]

    for header_item in wanted_header_items:
      if not header_item in header:
        return False

    detector_prefixes = header['DETECTOR_NAMES'].split()

    if len(detector_prefixes) != 1:
      return False

    detector_prefix = detector_prefixes[0]

    more_wanted_header_items = [
        'DETECTOR_DIMENSIONS', 'DETECTOR_SIZE', 'DETECTOR_VECTORS',
        'GONIO_NAMES', 'GONIO_UNITS', 'GONIO_VALUES', 'GONIO_VECTORS'
        ]

    for header_item in more_wanted_header_items:
      if not '%s%s' % (detector_prefix, header_item) in header:
        return False

    det_desc = '%sDETECTOR_DESCRIPTION' % detector_prefix
    if 'CMOS-1' in header.get(det_desc, ''):
        return True

    return False

  def __init__(self, image_file):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    assert(self.understand(image_file))

    FormatSMV.__init__(self, image_file)

    detector_prefixes = self._header_dictionary['DETECTOR_NAMES'].split()
    self._prefix = detector_prefixes[0]
    return

  def _start(self):
    FormatSMV._start(self)
    self._header_size = int(self._header_dictionary['HEADER_BYTES'])

  def _goniometer(self):
    axis = tuple(map(float, self._header_dictionary[
        'ROTATION_VECTOR'].split()))

    return self._goniometer_factory.known_axis(axis)

  def _detector(self):
    from scitbx import matrix
    detector_name = self._header_dictionary[
        'DETECTOR_NAMES'].split()[0].strip()

    detector_axes = map(float, self._header_dictionary[
        '%sDETECTOR_VECTORS' % detector_name].split())

    detector_fast = matrix.col(tuple(detector_axes[:3]))
    detector_slow = matrix.col(tuple(detector_axes[3:]))

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

    overload = int(float(self._header_dictionary['SATURATED_VALUE']))
    underload = 0

    return self._detector_factory.complex(
        'CCD', detector_origin.elems, detector_fast.elems,
        detector_slow.elems, pixel_size, image_size, (underload, overload))

  def _beam(self):
    beam_direction = map(float, self._header_dictionary[
        'SOURCE_VECTORS'].split()[:3])

    polarization = map(float, self._header_dictionary[
        'SOURCE_POLARZ'].split())

    p_fraction = polarization[0]
    p_plane = polarization[1:]

    wavelength = float(self._header_dictionary['WAVELENGTH'])

    return self._beam_factory.complex(
        beam_direction, p_fraction, p_plane, wavelength)

  def _scan(self):
    import calendar
    import time

    rotation = map(float, self._header_dictionary['ROTATION'].split())

    format = self._scan_factory.format('SMV')

    time_record = self._header_dictionary['DATE']
    date_record, ms = time_record.split('.')
    epoch = calendar.timegm(time.strptime(date_record, '%a %b %d %Y %H:%M:%S'))
    epoch += 0.001 * int(ms)

    exposure_time = rotation[3]
    osc_start = rotation[0]
    osc_range = rotation[2]

    return self._scan_factory.single(
        self._image_file, format, exposure_time,
        osc_start, osc_range, epoch)

  def get_raw_data(self):
    '''Get the pixel intensities (i.e. read the image and return as a
       flex array of integers.)'''

    # currently have no non-little-endian machines...

    from boost.python import streambuf
    from dxtbx import read_uint16
    from scitbx.array_family import flex
    assert(len(self.get_detector()) == 1)
    size = self.get_detector()[0].get_image_size()
    f = FormatSMV.open_file(self._image_file)
    f.read(self._header_size)
    raw_data = read_uint16(streambuf(f), int(size[0] * size[1]))
    image_size = self.get_detector()[0].get_image_size()
    raw_data.reshape(flex.grid(image_size[1], image_size[0]))

    return raw_data


if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatSMVCMOS1.understand(arg)
