#!/usr/bin/env python
# FormatCBFFullPilatus.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Pilatus implementation of fullCBF format, for use with Dectris detectors.

from __future__ import division

#import pycbf

from dxtbx.format.FormatCBFFullPilatus import FormatCBFFullPilatus
#from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask

class FormatCBFFullPilatusCdTe300K(FormatCBFFullPilatus):
  '''An image reading class for full CBF format images from Pilatus
  detectors.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an CBF format image, i.e. we can
    make sense of it.'''

    header = FormatCBFFullPilatus.get_cbf_header(image_file)

    for record in header.split('\n'):
      if 'CdTe sensor' in record:
        return True
      if '_array_data.data' in record:
        break

    return False

  def __init__(self, image_file):
    '''Initialise the image structure from the given file.'''

    assert(self.understand(image_file))

    FormatCBFFullPilatus.__init__(self, image_file)

    return

  def _detector(self):
    '''Return a working detector instance, with added mask regions.'''

    cif_header = FormatCBFFullPilatus.get_cbf_header(self._image_file)

    self._cif_header_dictionary = { }

    for record in cif_header.split('\n'):
      if not '#' in record[:1]:
        continue

      if len(record[1:].split()) <= 2 and record.count(':') == 2:
        self._cif_header_dictionary['timestamp'] = record[1:].strip()
        continue

      tokens = record.replace('=', '').replace(':', '').split()[1:]

      self._cif_header_dictionary[tokens[0]] = ' '.join(tokens[1:])

    for record in self._mime_header.split('\n'):
      if not record.strip():
        continue
      token, value = record.split(':')
      self._cif_header_dictionary[token.strip()] = value.strip()


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

    nx = int(
        self._cif_header_dictionary['X-Binary-Size-Fastest-Dimension'])
    ny = int(
        self._cif_header_dictionary['X-Binary-Size-Second-Dimension'])

    overload = int(
        self._cif_header_dictionary['Count_cutoff'].split()[0])
    underload = -1

    two_theta = float(self._cif_header_dictionary['Detector_2theta'].split()[0])

    # hard code "correct" values for these
    distance = 0.11338
    beam_x = 222.0
    beam_y = 383.0

    detector = self._detector_factory.two_theta(
        'PAD', distance * 1000.0, (beam_x * pixel_x * 1000.0,
                                   beam_y * pixel_y * 1000.0), '+x', '-y',
                                   '+x', two_theta,
        (1000 * pixel_x, 1000 * pixel_y),
        (nx, ny), (underload, overload), [], None)

    import re
    m = re.search('^#\s*(\S+)\ssensor, thickness\s*([0-9.]+)\s*m\s*$', \
                  self._cif_header, re.MULTILINE)

    if m:
      # header gives thickness in metres, we store mm
      thickness = float(m.group(2)) * 1000
      material = m.group(1)

      if material == 'Silicon':
        material = 'Si'

      for panel in detector:
        panel.set_thickness(thickness)
        panel.set_material(material)

      try:
        # a header only CBF file will not have a beam object
        beam = self._beam()

      except:
        pass

      if beam:
        # attenuation coefficient depends on the beam wavelength
        wavelength = beam.get_wavelength()

        from cctbx.eltbx import attenuation_coefficient
        from dxtbx.model import ParallaxCorrectedPxMmStrategy
        # this will fail for undefined composite materials (ie all except CdTe)
        table = attenuation_coefficient.get_table(material)
        # mu_at_angstrom returns cm^-1, but need mu in mm^-1
        mu = table.mu_at_angstrom(wavelength) / 10.0

        for panel in detector:
          panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, thickness))

    return detector

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFFullPilatus.understand(arg)
