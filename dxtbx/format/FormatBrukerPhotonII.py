#!/usr/bin/env python
# FormatBrukerPhotonII.py
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatBruker import FormatBruker
from scitbx import matrix

class FormatBrukerPhotonII(FormatBruker):

  @staticmethod
  def understand(image_file):

    try:
      header_lines = FormatBruker.read_header_lines(image_file)
    except IOError:
      return False

    header_dic = {}

    for l in header_lines:
      k_v = l.split(':', 1)
      if len(k_v) == 1: continue
      k, v = [v.strip() for v in k_v]
      if k in header_dic:
        header_dic[k] = header_dic[k] + "\n" + v
      else:
        header_dic[k] = v

    dettype = header_dic.get('DETTYPE')
    if dettype is None: return False
    if not dettype.startswith('CMOS-PHOTONII'): return False

    return True

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment. Easy from Rigaku Saturn images as
    they contain everything pretty much we need...'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    self._image_file = image_file
    FormatBruker.__init__(self, image_file, **kwargs)

  def _start(self):

    try:
      header_lines = FormatBruker.read_header_lines(self._image_file)
    except IOError:
      return False

    self.header_dict = {}
    for l in header_lines:
      k, v = [v.strip() for v in l.split(':', 1)]
      if k in self.header_dict:
        self.header_dict[k] = self.header_dict[k] + "\n" + v
      else:
        self.header_dict[k] = v

    from iotbx.detectors.bruker import BrukerImage
    self.detectorbase = BrukerImage(self._image_file)

    return

  def _goniometer(self):
    # goniometer angles in ANGLES are 2-theta, omega, phi, chi (FIXED)
    # AXIS indexes into this list to define the scan axis (in FORTRAN counting)
    # START and RANGE define the start and step size for each image
    # assume omega is 1,0,0 axis; chi about 0,0,1 at datum

    from scitbx import matrix
    angles = map(float, self.header_dict['ANGLES'].split())

    beam = matrix.col((0, 0, 1))
    phi = matrix.col((1, 0, 0)).rotate(-beam, angles[3], deg=True)

    if self.header_dict['AXIS'][0] == '2':
      # OMEGA scan
      axis = (-1, 0, 0)
      incr = float(self.header_dict['INCREME'])
      if incr < 0:
        axis = (1, 0, 0)
      fixed = phi.axis_and_angle_as_r3_rotation_matrix(angles[2], deg=True)
      return self._goniometer_factory.make_goniometer(axis, fixed.elems)
    else:
      # PHI scan
      assert(self.header_dict['AXIS'][0] == '3')
      omega = matrix.col((1, 0, 0))
      axis = phi.rotate(omega, angles[1], deg=True)
      return self._goniometer_factory.known_axis(axis.elems)

  def _detector(self):
    # goniometer angles in ANGLES are 2-theta, omega, phi, chi (FIXED)
    two_theta = float(self.header_dict['ANGLES'].split()[0])

    overload = float(self.header_dict['MAXIMUM'].split()[0])
    underload = float(self.header_dict['MINIMUM'].split()[0])

    # This is mysterious. The images are 16 bit, but the header claims
    # maximum count per pixel of >160'000. The product sheet talks of
    # 'adaptive oversampling' without explaining this. For now, just reset
    # the overload value
    overload = 65535 - 1

    fast = matrix.col((1, 0, 0))
    slow = matrix.col((0, 1, 0))
    beam = matrix.col((0, 0, 1))
    pixel_mm = 5.0 / float(self.header_dict['DETTYPE'].split()[1])
    beam_pixel = map(float, self.header_dict['CENTER'].split()[:2])
    distance_mm = 10.0 * float(self.header_dict['DISTANC'].split()[1])
    origin = - distance_mm * beam - fast * pixel_mm * beam_pixel[1] - \
      slow * pixel_mm * beam_pixel[0]
    origin = origin.rotate(-fast, two_theta, deg = True)
    slow = slow.rotate(-fast, two_theta, deg = True)
    pixel_size = pixel_mm, pixel_mm
    # ncols is nfast, nrows is nslow
    image_size = int(self.header_dict['NCOLS'].split()[0]), \
      int(self.header_dict['NROWS'].split()[0])

    # Not a CCD, but is an integrating detector
    return self._detector_factory.complex(
        'CCD', origin.elems, fast.elems, slow.elems, pixel_size, image_size,
      (underload, overload))

  def _beam(self):
    wavelength = float(self.header_dict['WAVELEN'].split()[0])

    return self._beam_factory.simple(wavelength)

  def _scan(self):

    start = float(self.header_dict['START'].split()[0])
    incr = float(self.header_dict['INCREME'].split()[0])
    if incr < 0:
      start *= -1
      incr *= -1

    return self._scan_factory.single(
      filename = self._image_file,
      format = "BrukerCCD",
      exposure_times = 1,
      osc_start = start,
      osc_width = incr,
      epoch = None)

  def get_raw_data(self):
    '''Get the pixel intensities (i.e. read the image and return as a
    flex array of integers.)'''

    from boost.python import streambuf
    from dxtbx import read_uint16, read_uint16_bs, is_big_endian
    from scitbx.array_family import flex
    f = self.open_file(self._image_file, 'rb')
    header_size = int(self.header_dict['HDRBLKS']) * 512
    f.read(header_size)

    # 16 bits per pixel
    assert int(self.header_dict['NPIXELB'].split()[0]) == 2

    nrows = int(self.header_dict['NROWS'].split()[0])
    ncols = int(self.header_dict['NCOLS'].split()[0])

    if is_big_endian():
      raw_data = read_uint16_bs(streambuf(f), nrows*ncols)
    else:
      raw_data = read_uint16(streambuf(f), nrows*ncols)

    image_size = (nrows, ncols)
    raw_data.reshape(flex.grid(image_size[0], image_size[1]))

    return raw_data

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatBrukerPhotonII.understand(arg))
