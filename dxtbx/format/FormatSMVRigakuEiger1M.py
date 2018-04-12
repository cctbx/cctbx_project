#!/usr/bin/env python
# FormatSMVRigakuEiger1M.py
#   Copyright (C) 2013 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the SMV image reader for Rigaku Pilatus 200L images.
# Inherits from FormatSMVRigaku. Be aware: this is completely unrelated
# to the HDF5 Eiger format.

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatSMVRigaku import FormatSMVRigaku

class FormatSMVRigakuEiger1M(FormatSMVRigaku):
  '''A class for reading SMV format Rigaku Pilatus 200L images.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like a Rigaku d*TREK SMVRigaku format image,
    i.e. we can make sense of it. Essentially that will be if it contains
    all of the keys we are looking for.'''

    size, header = FormatSMVRigaku.get_smv_header(image_file)

    if not 'DETECTOR_NAMES' in header:
      return False

    name = header['DETECTOR_NAMES']

    if not '%sDETECTOR_DESCRIPTION' % name in header:
      return False

    return 'Eiger1M' in header['%sDETECTOR_DESCRIPTION' % name]

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatSMVRigaku.__init__(self, image_file, **kwargs)

  def _start(self):

    FormatSMVRigaku._start(self)

  def _goniometer(self):
    '''Initialize the structure for the goniometer - this will need to
    correctly compose the axes given in the image header. In this case
    this is made rather straightforward as the image header has the
    calculated rotation axis stored in it. We could work from the
    rest of the header and construct a goniometer model.'''

    axis = tuple(map(float, self._header_dictionary[
        'ROTATION_VECTOR'].split()))

    return self._goniometer_factory.known_axis(axis)

  def _detector(self):
    '''Return a model for the detector, allowing for two-theta offsets
    and the detector position. This will be rather more complex...'''

    from scitbx import matrix

    detector_name = self._header_dictionary[
        'DETECTOR_NAMES'].split()[0].strip()

    detector_axes = map(float, self._header_dictionary[
        '%sDETECTOR_VECTORS' % detector_name].split())

    fast = matrix.col(tuple(detector_axes[:3]))
    slow = matrix.col(tuple(detector_axes[3:]))

    distortion = map(int, self._header_dictionary[
      '%sSPATIAL_DISTORTION_VECTORS' % detector_name].split())

    # multiply through by the distortion to get the true detector fast, slow

    detector_fast, detector_slow = distortion[0] * fast + distortion[1] * slow, \
      distortion[2] * fast + distortion[3] * slow

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
        raise RuntimeError('unknown axis unit %s' % unit)

    rotations.reverse()
    translations.reverse()

    for j in range(gonio_num_axes):
      detector_fast = rotations[j] * detector_fast
      detector_slow = rotations[j] * detector_slow
      detector_origin = rotations[j] * detector_origin
      detector_origin = translations[j] + detector_origin

    overload = int(float(self._header_dictionary['SATURATED_VALUE']))
    underload = -1

    return self._detector_factory.complex(
        'PAD', detector_origin.elems, detector_fast.elems,
        detector_slow.elems, pixel_size, image_size, (underload, overload))

  def _beam(self):
    '''Return a simple model for the beam.'''

    beam_direction = map(float, self._header_dictionary[
        'SOURCE_VECTORS'].split()[:3])

    polarization = map(float, self._header_dictionary[
        'SOURCE_POLARZ'].split())

    p_fraction = polarization[0]
    p_plane = polarization[1:]

    wavelength = float(self._header_dictionary['SCAN_WAVELENGTH'])

    return self._beam_factory.complex(
        beam_direction, p_fraction, p_plane, wavelength)

  def _scan(self):
    '''Return the scan information for this image.'''

    import time

    rotation = map(float, self._header_dictionary['ROTATION'].split())

    format = self._scan_factory.format('SMV')

    epoch = time.mktime(time.strptime(self._header_dictionary[
        'DTREK_DATE_TIME'], '%d-%b-%Y %H:%M:%S'))

    exposure_time = rotation[3]
    osc_start = rotation[0]
    osc_range = rotation[2]

    return self._scan_factory.single(
        self._image_file, format, exposure_time,
        osc_start, osc_range, epoch)

  def get_raw_data(self):
    '''Read the data - assuming it is streaming 4-byte unsigned ints following the
    header...'''

    from boost.python import streambuf
    from dxtbx import read_int32
    from scitbx.array_family import flex
    assert(len(self.get_detector()) == 1)
    size = self.get_detector()[0].get_image_size()
    f = self.open_file(self._image_file)
    f.read(int(self._header_dictionary['HEADER_BYTES']))
    raw_data = read_int32(streambuf(f), int(size[0] * size[1]))
    raw_data.reshape(flex.grid(size[1], size[0]))

    return raw_data


if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatSMVRigakuEiger1M.understand(arg))
