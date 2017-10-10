#!/usr/bin/env python
# FormatSMVTimePix_SU.py
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""Implementation of a format class to specifically recognise images
from a Timepix-based detector, installed on a JEOL-2100 electron microscope at
Stockholm University, which have been preprocessed by Wei Wan's RED software,
to produce SMV files. See https://doi.org/10.1107/S0021889813027714 for RED.
The detector itself has a 2x2 array of Timepix modules"""

from __future__ import absolute_import, division
from __future__ import print_function
import time
from dxtbx.format.FormatSMV import FormatSMV
from dxtbx.model.detector import Detector

class FormatSMVTimePix_SU(FormatSMV):
  '''A class for reading SMV format images from the RED software, where the
  main identifying feature is BEAMLINE=TimePix_SU. This is intended to
  construct a model for an electron diffraction experiment using a microscope
  at Stockholm University.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an ADSC SMV format image, i.e. we
    can make sense of it. Essentially that will be if it contains all of
    the keys we are looking for and not some we are not (i.e. that belong
    to a Rigaku Saturn.)'''

    size, header = FormatSMV.get_smv_header(image_file)

    # only recognise TimePix_SU
    if header.get('BEAMLINE') != 'TimePix_SU': return False

    # check the header contains the things we're going to use
    wanted_header_items = ['BEAM_CENTER_X', 'BEAM_CENTER_Y',
                           'DISTANCE', 'WAVELENGTH', 'PIXEL_SIZE',
                           'OSC_START', 'OSC_RANGE', 'PHI', 'SIZE1', 'SIZE2',
                           'BYTE_ORDER', 'DETECTOR_SN']
    for header_item in wanted_header_items:
      if not header_item in header:
        return False

    # check the pixel size is 55 microns
    if not float(header['PIXEL_SIZE']) == 0.055: return False

    # check there are 512*512 pixels
    if not (header['SIZE1']) == '512': return False
    if not (header['SIZE2']) == '512': return False

    return True

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatSMV.__init__(self, image_file, **kwargs)

    return

  def _start(self):

    FormatSMV._start(self)

  def detectorbase_start(self):
    if not hasattr(self, "detectorbase") or self.detectorbase is None:
      from iotbx.detectors import SMVImage
      self.detectorbase = SMVImage(self._image_file)
      self.detectorbase.open_file = self.open_file
      self.detectorbase.readHeader()

  def _goniometer(self):
    '''Return a model for a simple single-axis goniometer. For this beamline
    this should be close to the provided values and neither aligned with the
    detector slow or fast axes.'''

    return self._goniometer_factory.known_axis((-0.755, -0.656, 0.0))

  def _detector(self):
    '''4 panel detector, 55 micron pixels except for pixels at the outer
    edge of each chip, which are 165 microns wide.'''
    from scitbx import matrix

    # expect 55 mu pixels here, but make it general anyway
    pixel_size = tuple([float(self._header_dictionary['PIXEL_SIZE'])] * 2)
    image_size = (int(self._header_dictionary['SIZE1']),
                  int(self._header_dictionary['SIZE2']))
    panel_size = tuple([int(e/2) for e in image_size])

    # outer pixels have three times the width
    panel_size_mm = (pixel_size[0] * 3 + (panel_size[0] - 2) * pixel_size[0],
                     pixel_size[1] * 3 + (panel_size[1] - 2) * pixel_size[1])
    image_size_mm = (panel_size_mm[0] * 2, panel_size_mm[1] * 2)
    trusted_range = (-1, 65535)
    material = 'Si'
    thickness = 0.3 # assume 300 mu thick

    # Initialise detector frame
    fast = matrix.col((1.0, 0.0, 0.0))
    slow = matrix.col((0.0, -1.0, 0.0))
    beam_centre = (float(self._header_dictionary['BEAM_CENTER_X']),
                   float(self._header_dictionary['BEAM_CENTER_Y']))

    bx_px, by_px = beam_centre
    # the beam centre is in pixels. We want to convert to mm, taking the
    # different size of outer pixels into account. Use this local function
    # to do that
    def px_to_mm(px, px_size_1d, panel_size_1d):
      mm = 0
      if px > 1: # add first outer pixel
        mm += px_size_1d * 3
      else: # or fraction of first outer pixel
        mm += px * px_size_1d * 3
        return mm

      if px > panel_size_1d - 1: # add full panel of inner pixels
        mm += (panel_size_1d - 2) * px_size_1d
      else: # or fraction of inner pixels
        mm += (px - 1) * px_size_1d
        return mm

      if px > panel_size_1d: # add second outer pixel
        mm += px_size_1d * 3
      else: # or fraction of second outer pixel
        mm += (px - (panel_size_1d - 1)) * px_size_1d * 3
        return mm

      if px > panel_size_1d + 1: # add first outer pixel of second panel
        mm += px_size_1d * 3
      else: # or fraction of first outer pixel of second panel
        mm += (px - panel_size_1d) * px_size_1d * 3
        return mm

      if px > (2 * panel_size_1d - 1): # add second full panel of inner pixels
        mm += (panel_size_1d - 2) * px_size_1d
        # plus remaining fraction of the second outer pixel
        mm += (px - (2 * panel_size_1d - 1)) * px_size_1d * 3
      else: # or fraction of inner pixels of the second panel
        mm += (px - panel_size_1d - 1) * px_size_1d
      return mm

    bx_mm = px_to_mm(bx_px, pixel_size[0], panel_size[0])
    by_mm = px_to_mm(by_px, pixel_size[1], panel_size[1])
    beam_centre_mm = (bx_mm, by_mm)

    # the beam centre is defined from the origin along fast, slow. To determine
    # the lab frame origin we place the beam centre down the -z axis
    cntr = matrix.col((0.0, 0.0, -1 * float(self._header_dictionary['DISTANCE'])))
    orig = cntr - bx_mm * fast - by_mm * slow

    d = Detector()

    root = d.hierarchy()
    root.set_local_frame(
      fast.elems,
      slow.elems,
      orig.elems)

    self.coords = {}
    panel_idx = 0

    # set panel extent in pixel numbers and x, y mm shifts. Note that the
    # outer pixels are 0.165 mm in size. These are excluded from the panel
    # extents.
    pnl_data = []
    pnl_data.append({'xmin':1, 'ymin':1,
                     'xmax':255, 'ymax':255,
                     'xmin_mm': 1 * 0.165,
                     'ymin_mm': 1 * 0.165})
    pnl_data.append({'xmin':257, 'ymin':1,
                     'xmax':511, 'ymax':255,
                     'xmin_mm': 3 * 0.165 + (511 - 257) * pixel_size[0],
                     'ymin_mm': 1 * 0.165})
    pnl_data.append({'xmin':1, 'ymin':257,
                     'xmax':255, 'ymax':511,
                     'xmin_mm': 1 * 0.165,
                     'ymin_mm': 3 * 0.165 + (511 - 257) * pixel_size[1]})
    pnl_data.append({'xmin':257, 'ymin':257,
                     'xmax':511, 'ymax':511,
                     'xmin_mm': 3 * 0.165 + (511 - 257) * pixel_size[0],
                     'ymin_mm': 3 * 0.165 + (511 - 257) * pixel_size[1]})

    # redefine fast, slow for the local frame
    fast = matrix.col((1.0, 0.0, 0.0))
    slow = matrix.col((0.0, 1.0, 0.0))

    for ipanel, pd in enumerate(pnl_data):
      xmin = pd['xmin']
      xmax = pd['xmax']
      ymin = pd['ymin']
      ymax = pd['ymax']
      xmin_mm = pd['xmin_mm']
      ymin_mm = pd['ymin_mm']

      origin_panel = fast * xmin_mm + slow * ymin_mm

      panel_name = "Panel%d" % panel_idx
      panel_idx += 1

      p = d.add_panel()
      p.set_type('SENSOR_PAD')
      p.set_name(panel_name)
      p.set_raw_image_offset((xmin, ymin))
      p.set_image_size((xmax-xmin, ymax-ymin))
      p.set_trusted_range(trusted_range)
      p.set_pixel_size((pixel_size[0], pixel_size[1]))
      p.set_thickness(thickness)
      p.set_material('Si')
      #p.set_mu(mu)
      #p.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, t0))
      p.set_local_frame(
        fast.elems,
        slow.elems,
        origin_panel.elems)
      p.set_raw_image_offset((xmin, ymin))
      self.coords[panel_name] = (xmin, ymin, xmax, ymax)

    return d

  def _beam(self):
    '''Return an unpolarized beam model.'''

    wavelength = float(self._header_dictionary['WAVELENGTH'])

    return self._beam_factory.make_polarized_beam(
        sample_to_source=(0.0, 0.0, 1.0),
        wavelength=wavelength,
        polarization=(0, 1, 0),
        polarization_fraction=0.5)

  def _scan(self):
    '''Return the scan information for this image.'''
    import calendar

    format = self._scan_factory.format('SMV')
    exposure_time = float(self._header_dictionary['TIME'])
    epoch = None

    # PST, PDT timezones not recognised by default...

    epoch = 0
    try:
      date_str = self._header_dictionary['DATE']
      date_str = date_str.replace('PST', '').replace('PDT', '')
    except KeyError as e:
      date_str = ''
    for format_string in ['%a %b %d %H:%M:%S %Y', '%a %b %d %H:%M:%S %Z %Y']:
      try:
        epoch = calendar.timegm(time.strptime(date_str, format_string))
        break
      except ValueError as e:
        pass

    # assert(epoch)
    osc_start = float(self._header_dictionary['OSC_START'])
    osc_range = float(self._header_dictionary['OSC_RANGE'])

    return self._scan_factory.single(
        self._image_file, format, exposure_time,
        osc_start, osc_range, epoch)

  def get_raw_data(self):
    '''Get the pixel intensities (i.e. read the image and return as a
    flex array of integers.)'''

    from boost.python import streambuf
    from dxtbx import read_uint16, read_uint16_bs, is_big_endian
    from scitbx.array_family import flex
    f = self.open_file(self._image_file, 'rb')
    f.read(self._header_size)

    if self._header_dictionary['BYTE_ORDER'] == 'big_endian':
      big_endian = True
    else:
      big_endian = False

    if big_endian == is_big_endian():
      raw_data = read_uint16(streambuf(f), int(512 * 512))
    else:
      raw_data = read_uint16_bs(streambuf(f), int(512 * 512))

    image_size = (512,512)
    raw_data.reshape(flex.grid(image_size[1], image_size[0]))

    # split into separate panels
    self._raw_data = []
    d = self.get_detector()
    for panel in d:
      xmin, ymin, xmax, ymax = self.coords[panel.get_name()]
      self._raw_data.append(raw_data[ymin:ymax,xmin:xmax])

    return tuple(self._raw_data)

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatSMVTimePix_SU.understand(arg))
