#!/usr/bin/env python
# FormatCBFMiniEiger.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Eiger images. Inherits from
# FormatCBFMini.

from __future__ import division

from dxtbx.format.FormatCBFMini import FormatCBFMini
from dxtbx.format.FormatCBFMiniPilatusHelpers import \
     get_pilatus_timestamp
from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask
from dxtbx.model import ParallaxCorrectedPxMmStrategy

class FormatCBFMiniEiger(FormatCBFMini):
  '''A class for reading mini CBF format Eiger images, and correctly
  constructing a model for the experiment from this.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Eiger mini CBF format image,
    i.e. we can make sense of it.'''

    header = FormatCBFMini.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# Detector' in record and \
             'EIGER' in record:
        return True

    return False

  def __init__(self, image_file):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    assert(self.understand(image_file))

    FormatCBFMini.__init__(self, image_file)

    return

  def _start(self):
    FormatCBFMini._start(self)
    try:
      from iotbx.detectors.pilatus_minicbf import PilatusImage
      self.detectorbase = PilatusImage(self._image_file)
      self.detectorbase.readHeader()
    except KeyError, e:
      pass

  def _goniometer(self):
    return self._goniometer_factory.single_axis_reverse()

  def _detector(self):
    distance = float(
        self._cif_header_dictionary['Detector_distance'].split()[0])

    beam_xy = self._cif_header_dictionary['Beam_xy'].replace(
        '(', '').replace(')', '').replace(',', '').split()[:2]

    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    beam_x, beam_y = map(float, beam_xy)

    pixel_xy = self._cif_header_dictionary['Pixel_size'].replace(
        'm', '').replace('x', '').split()

    pixel_x, pixel_y = map(float, pixel_xy)

    thickness = float(
      self._cif_header_dictionary['Silicon'].split()[-2]) * 1000.0

    nx = int(
        self._cif_header_dictionary['X-Binary-Size-Fastest-Dimension'])
    ny = int(
        self._cif_header_dictionary['X-Binary-Size-Second-Dimension'])

    overload = int(
        self._cif_header_dictionary['Count_cutoff'].split()[0])
    underload = -1

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

    for f0, s0, f1, s1 in [(0, 513, 1030, 550)]:
      detector[0].add_mask(f0, s0, f1, s1)

    return detector

  def _beam(self):
    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    return self._beam_factory.simple(wavelength)

  def _scan(self):
    format = self._scan_factory.format('CBF')

    exposure_time = float(
        self._cif_header_dictionary['Exposure_period'].split()[0])
    osc_start = float(
        self._cif_header_dictionary['Start_angle'].split()[0])
    osc_range = float(
        self._cif_header_dictionary['Angle_increment'].split()[0])

    timestamp = get_pilatus_timestamp(
        self._cif_header_dictionary['timestamp'])

    return self._scan_factory.single(
        self._image_file, format, exposure_time,
        osc_start, osc_range, timestamp)

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFMiniEiger.understand(arg)
