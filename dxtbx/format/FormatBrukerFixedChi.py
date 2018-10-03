#!/usr/bin/env python
# FormatBrukerFixedChi.py
#   Copyright (C) 2016 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatBruker import FormatBruker
from scitbx import matrix

class FormatBrukerFixedChi(FormatBruker):

  @staticmethod
  def understand(image_file):

    hdr = FormatBruker.read_header_lines(image_file)
    hdr_dic = FormatBruker.parse_header(hdr)

    if 'FixedChiStage' not in hdr_dic['MODEL']: return False
    if 'PHOTONII' in hdr_dic['DETTYPE']: return False

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
    self.header_dict = { }
    header_text = open(self._image_file).read().split('......')[0]
    for j in range(0, len(header_text), 80):
      record = header_text[j:j + 80]
      if record.startswith('CFR:'):
        continue
      key = record[:8]
      assert(key[-1] == ':')
      values = record[8:]
      self.header_dict[key.replace(':', '').strip()] = [
        v.strip() for v in values.split()]
    from iotbx.detectors.bruker import BrukerImage
    self.detectorbase = BrukerImage(self._image_file)

  def _goniometer(self):
    # goniometer angles in ANGLES are 2-theta, omega, phi, chi (FIXED)
    # AXIS indexes into this list to define the scan axis (in FORTRAN counting)
    # START and RANGE define the start and step size for each image
    # assume omega is 1,0,0 axis; chi about 0,0,1 at datum

    from scitbx import matrix
    angles = map(float, self.header_dict['ANGLES'])

    beam = matrix.col((0, 0, 1))
    phi = matrix.col((1, 0, 0)).rotate(-beam, angles[3], deg=True)

    if self.header_dict['AXIS'][0] == '2':
      # OMEGA scan
      axis = (-1, 0, 0)
      incr = float(self.header_dict['INCREME'][0])
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
    two_theta = float(self.header_dict['ANGLES'][0])
    overload = 60000
    underload = 0

    fast = matrix.col((1, 0, 0))
    slow = matrix.col((0, 1, 0))
    beam = matrix.col((0, 0, 1))
    pixel_mm = 5.0 / float(self.header_dict['DETTYPE'][1])
    beam_pixel = map(float, self.header_dict['CENTER'][:2])
    distance_mm = 10.0 * float(self.header_dict['DISTANC'][1])
    origin = - distance_mm * beam - fast * pixel_mm * beam_pixel[1] - \
      slow * pixel_mm * beam_pixel[0]
    origin = origin.rotate(-fast, two_theta, deg = True)
    slow = slow.rotate(-fast, two_theta, deg = True)
    pixel_size = pixel_mm, pixel_mm
    image_size = int(self.header_dict['NROWS'][0]), \
      int(self.header_dict['NCOLS'][0])

    return self._detector_factory.complex(
        'CCD', origin.elems, fast.elems, slow.elems, pixel_size, image_size,
      (underload, overload))

  def _beam(self):
    wavelength = float(self.header_dict['WAVELEN'][0])

    return self._beam_factory.simple(wavelength)

  def _scan(self):

    start = float(self.header_dict['START'][0])
    incr = float(self.header_dict['INCREME'][0])
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

  # FIXME implement get_raw_data including LINEAR factor (multiplier) =>
  # explains the observed GAIN around 10

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatBrukerFixedChi.understand(arg))
