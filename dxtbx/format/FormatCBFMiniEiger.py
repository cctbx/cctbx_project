#!/usr/bin/env python
# FormatCBFMiniEiger.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Eiger images. Inherits from
# FormatCBFMini.

from __future__ import absolute_import, division
from __future__ import print_function

from builtins import map
from builtins import zip
from builtins import range
from dxtbx.format.FormatCBFMini import FormatCBFMini
from dxtbx.format.FormatCBFMiniPilatusHelpers import \
     get_pilatus_timestamp
from dxtbx.model import ParallaxCorrectedPxMmStrategy
import os

if 'DXTBX_OVERLOAD_SCALE' in os.environ:
  dxtbx_overload_scale = float(os.environ['DXTBX_OVERLOAD_SCALE'])
else:
  dxtbx_overload_scale = 1

def determine_eiger_mask(xdetector):
  '''Return an appropriate pixel mask for a Eiger detector.'''

  size = xdetector[0].get_image_size()

  # Hardcoded module size and gap size
  module_size_fast, module_size_slow = (1030, 514)
  gap_size_fast, gap_size_slow = (10, 37)

  # Edge dead areas not included, only gaps between modules matter
  n_fast, remainder = divmod(size[0], module_size_fast)
  assert (n_fast-1) * gap_size_fast == remainder

  n_slow, remainder = divmod(size[1], module_size_slow)
  assert (n_slow-1) * gap_size_slow == remainder

  # Specify the dead areas between the modules, i.e. the rows and columns
  # where there are no active pixels
  mask = []
  for i_fast in range(n_fast-1):
    mask.append([
      (i_fast+1) * module_size_fast + i_fast * gap_size_fast + 1,
      (i_fast+1) * module_size_fast + i_fast * gap_size_fast + gap_size_fast,
      1, size[1]
    ])
  for i_slow in range(n_slow-1):
    mask.append([
      1, size[0],
      (i_slow+1) * module_size_slow + i_slow * gap_size_slow + 1,
      (i_slow+1) * module_size_slow + i_slow * gap_size_slow + gap_size_slow
    ])

  return mask

class FormatCBFMiniEiger(FormatCBFMini):
  '''A class for reading mini CBF format Eiger images, and correctly
  constructing a model for the experiment from this.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Eiger mini CBF format image,
    i.e. we can make sense of it.'''

    header = FormatCBFMini.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# detector' in record.lower() and \
             'eiger' in record.lower():
        return True

    return False

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatCBFMini.__init__(self, image_file, **kwargs)

    self._raw_data = None
    return

  def _start(self):
    FormatCBFMini._start(self)

  def _goniometer(self):
    if 'Phi' in self._cif_header_dictionary:
      phi_value = float(self._cif_header_dictionary['Phi'].split()[0])

    return self._goniometer_factory.single_axis()

  def _detector(self):
    distance = float(
        self._cif_header_dictionary['Detector_distance'].split()[0])

    beam_xy = self._cif_header_dictionary['Beam_xy'].replace(
        '(', '').replace(')', '').replace(',', '').split()[:2]

    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    beam_x, beam_y = list(map(float, beam_xy))

    pixel_xy = self._cif_header_dictionary['Pixel_size'].replace(
        'm', '').replace('x', '').split()

    pixel_x, pixel_y = list(map(float, pixel_xy))


    # extract sensor thickness, assume < 10mm; default to 450 microns unfound
    try:
        thickness = float(
          self._cif_header_dictionary['Silicon'].split()[-2]) * 1000.0
        material = 'Si'
        if thickness > 10:
            thickness = thickness / 1000.0
    except KeyError:
        thickness = 0.450
        material = 'Si'

    nx = int(
        self._cif_header_dictionary['X-Binary-Size-Fastest-Dimension'])
    ny = int(
        self._cif_header_dictionary['X-Binary-Size-Second-Dimension'])

    if 'Count_cutoff' in self._cif_header_dictionary:
      overload = int(self._cif_header_dictionary['Count_cutoff'].split()[0])
    else:
      # missing from data transformed with GPhL converter - dials#376
      overload = 100000000
    underload = -1

    try:
        identifier =  self._cif_header_dictionary['Detector']
    except KeyError:
        identifier = 'Unknown Eiger'

    from cctbx.eltbx import attenuation_coefficient
    table = attenuation_coefficient.get_table("Si")
    mu = table.mu_at_angstrom(wavelength) / 10.0
    t0 = thickness

    detector = self._detector_factory.simple(
        'PAD', distance * 1000.0, (beam_x * pixel_x * 1000.0,
                                   beam_y * pixel_y * 1000.0), '+x', '-y',
        (1000 * pixel_x, 1000 * pixel_y),
        (nx, ny), (underload, overload), [],
        ParallaxCorrectedPxMmStrategy(mu, t0))

    for f0, s0, f1, s1 in determine_eiger_mask(detector):
      detector[0].add_mask(f0, s0, f1, s1)

    for panel in detector:
        panel.set_thickness(thickness)
        panel.set_material(material)
        panel.set_identifier(identifier)
        panel.set_mu(mu)

    return detector

  def _beam(self):
    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    beam = self._beam_factory.simple(wavelength)

    try:
      flux = float(self._cif_header_dictionary['Flux'].split()[0])
      beam.set_flux(flux)
    except KeyError:
      pass

    try:
      transmission = float(self._cif_header_dictionary['Transmission'].split()[0])
      beam.set_transmission(transmission)
    except KeyError:
      pass

    return beam

  def _scan(self):
    format = self._scan_factory.format('CBF')

    exposure_time = float(
        self._cif_header_dictionary['Exposure_period'].split()[0])
    osc_start = float(
        self._cif_header_dictionary['Start_angle'].split()[0])
    osc_range = float(
        self._cif_header_dictionary['Angle_increment'].split()[0])

    if "timestamp" in self._cif_header_dictionary:
      timestamp = get_pilatus_timestamp(
      self._cif_header_dictionary['timestamp'])
    else:
      timestamp = 0.0

    return self._scan_factory.single(
        self._image_file, format, exposure_time,
        osc_start, osc_range, timestamp)

  def read_cbf_image(self, cbf_image):
    from cbflib_adaptbx import uncompress
    import binascii
    from scitbx.array_family import flex

    start_tag = binascii.unhexlify('0c1a04d5')

    data = self.open_file(cbf_image, 'rb').read()
    data_offset = data.find(start_tag) + 4
    cbf_header = data[:data_offset - 4]

    fast = 0
    slow = 0
    length = 0

    for record in cbf_header.split('\n'):
      if 'X-Binary-Size-Fastest-Dimension' in record:
        fast = int(record.split()[-1])
      elif 'X-Binary-Size-Second-Dimension' in record:
        slow = int(record.split()[-1])
      elif 'X-Binary-Number-of-Elements' in record:
        length = int(record.split()[-1])
      elif 'X-Binary-Size:' in record:
        size = int(record.split()[-1])

    assert(length == fast * slow)

    pixel_values = uncompress(packed = data[data_offset:data_offset + size],
                              fast = fast, slow = slow)

    return pixel_values

  def get_raw_data(self):
    if self._raw_data is None:
      data = self.read_cbf_image(self._image_file)
      self._raw_data = data

    return self._raw_data

  def get_mask(self, goniometer=None):
    from scitbx.array_family import flex
    detector = self.get_detector()
    mask = [flex.bool(flex.grid(reversed(p.get_image_size())), True)
            for p in detector]
    for i, p in enumerate(detector):
      untrusted_regions = p.get_mask()
      for j, (f0, s0, f1, s1) in enumerate(untrusted_regions):
        sub_array = flex.bool(flex.grid(s1-s0+1, f1-f0+1), False)
        mask[i].matrix_paste_block_in_place(sub_array, s0-1, f0-1)

    if len(detector) == 1:
      raw_data = [self.get_raw_data()]
    else:
      raw_data = self.get_raw_data()
      assert(len(raw_data) ==  len(detector))
    trusted_mask = [
      p.get_trusted_range_mask(im) for im, p in zip(raw_data, detector)]

    # returns merged untrusted pixels and active areas using bitwise AND
    # (pixels are accepted if they are inside of the active areas AND
    # inside of the trusted range)
    return tuple([m & tm for m, tm in zip(mask, trusted_mask)])

  def detectorbase_start(self):
    from iotbx.detectors.eiger_minicbf import EigerCBFImage
    self.detectorbase = EigerCBFImage(self._image_file)
    self.detectorbase.readHeader()

  def get_vendortype(self):
    from dxtbx.format.FormatPilatusHelpers import get_vendortype_eiger as gv
    return gv(self.get_detector())

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatCBFMiniEiger.understand(arg))
