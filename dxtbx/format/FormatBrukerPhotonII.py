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

    # The Photon II format can't currently use BrukerImage, see
    # https://github.com/cctbx/cctbx_project/issues/65
    #from iotbx.detectors.bruker import BrukerImage
    #self.detectorbase = BrukerImage(self._image_file)

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

    overload = float(self.header_dict['CCDPARM'].split()[-1])
    underload = -1

    fast = matrix.col((1, 0, 0))
    slow = matrix.col((0, 1, 0))
    beam = matrix.col((0, 0, 1))
    pixel_mm = 5.0 / float(self.header_dict['DETTYPE'].split()[1])
    beam_pixel = map(float, self.header_dict['CENTER'].split()[:-3:-1])
    distance_mm = 10.0 * float(self.header_dict['DISTANC'].split()[1])
    origin = - distance_mm * beam - fast * pixel_mm * beam_pixel[1] - \
      slow * pixel_mm * beam_pixel[0]
    # 2theta rotation appears to be around the slow axis
    origin = origin.rotate(slow, two_theta, deg = True)
    fast = fast.rotate(slow, two_theta, deg = True)
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
    from dxtbx import read_uint8, read_uint16, read_uint16_bs, read_uint32, read_uint32_bs
    from dxtbx import is_big_endian
    from scitbx.array_family import flex
    f = self.open_file(self._image_file, 'rb')
    header_size = int(self.header_dict['HDRBLKS']) * 512
    f.read(header_size)

    if is_big_endian():
      read_2b = read_uint16_bs
      read_4b = read_uint32_bs
    else:
      read_2b = read_uint16
      read_4b = read_uint32

    # NPIXELB stores the number of bytes/pixel for the data and the underflow
    # table. We expect 1 byte for underflows and either 2 or 1 byte per pixel
    # for the data
    npixelb = [int(e) for e in self.header_dict['NPIXELB'].split()]
    assert npixelb[1] == 1

    if npixelb[0] == 1:
      read_data = read_uint8
    elif npixelb[0] == 2:
      read_data = read_2b
    else:
      from dxtbx import IncorrectFormatError
      raise IncorrectFormatError("{0} bytes per pixel is not supported".format(
        npixelb[0]))

    nrows = int(self.header_dict['NROWS'].split()[0])
    ncols = int(self.header_dict['NCOLS'].split()[0])

    raw_data = read_data(streambuf(f), nrows*ncols)

    image_size = (nrows, ncols)
    raw_data.reshape(flex.grid(*image_size))

    (num_underflows,
     num_2b_overflows,
     num_4b_overflows) = [int(e) for e in self.header_dict['NOVERFL'].split()]

    # read underflows
    if num_underflows > 0:
      # stored values are padded to a multiple of 16 bytes
      nbytes = num_underflows + 15 & ~(15)
      underflow_vals = read_uint8(streambuf(f), nbytes)[:num_underflows]
    else:
      underflow_vals = None

    # handle 2 byte overflows
    if num_2b_overflows > 0:
      # stored values are padded to a multiple of 16 bytes
      nbytes = num_2b_overflows * 2 + 15 & ~(15)
      overflow_vals = read_2b(streambuf(f), nbytes // 2)[:num_2b_overflows]
      overflow = flex.int(nrows * ncols, 0)
      sel = (raw_data == 255).as_1d()
      overflow.set_selected(sel, overflow_vals - 255)
      overflow.reshape(flex.grid(*image_size))
      raw_data += overflow

    # handle 4 byte overflows
    if num_4b_overflows > 0:
      # stored values are padded to a multiple of 16 bytes
      nbytes = num_4b_overflows * 4 + 15 & ~(15)
      overflow_vals = read_4b(streambuf(f), nbytes // 4)[:num_4b_overflows]
      overflow = flex.int(nrows * ncols, 0)
      sel = (raw_data == 65535).as_1d()
      overflow.set_selected(sel, overflow_vals - 65535)
      overflow.reshape(flex.grid(*image_size))
      raw_data += overflow

    # handle underflows
    if underflow_vals is not None:
      sel = (raw_data == 0).as_1d()
      underflow = flex.int(nrows * ncols, 0)
      underflow.set_selected(sel, underflow_vals)
      underflow.reshape(flex.grid(*image_size))
      raw_data += underflow

    # handle baseline. num_underflows == -1 means no baseline subtraction. See
    # https://github.com/cctbx/cctbx_project/files/1262952/BISFrameFileFormats.zip
    if num_underflows != -1:
      num_exposures = [int(e) for e in self.header_dict['NEXP'].split()]
      baseline = num_exposures[2]
      raw_data += baseline

    return raw_data

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatBrukerPhotonII.understand(arg))
