#!/usr/bin/env python
# FormatCBFFullPilatus.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Pilatus implementation of fullCBF format, for use with Dectris detectors.

from __future__ import division

import pycbf

from dxtbx.format.FormatCBFFull import FormatCBFFull
from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask

class FormatCBFFullPilatus(FormatCBFFull):
  '''An image reading class for full CBF format images from Pilatus
  detectors.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an CBF format image, i.e. we can
    make sense of it.'''

    header = FormatCBFFull.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '_array_data.header_convention' in record and \
             'PILATUS' in record:
        return True

    return False

  def __init__(self, image_file):
    '''Initialise the image structure from the given file.'''

    assert(self.understand(image_file))

    FormatCBFFull.__init__(self, image_file)

    return

  def _start(self):
    '''Open the image file as a cbf file handle, and keep this somewhere
    safe.'''

    FormatCBFFull._start(self)

    self._cbf_handle = pycbf.cbf_handle_struct()
    self._cbf_handle.read_widefile(self._image_file, pycbf.MSG_DIGEST)

    return

  def _detector(self):
    '''Return a working detector instance, with added mask regions.'''

    detector = self._detector_factory.imgCIF_H(self._cbf_handle,
                                               'PAD')

    # FIXME Set proper overload
    detector[0].set_trusted_range((-1, 1000000))
    for f0, s0, f1, s1 in determine_pilatus_mask(detector):
      detector[0].add_mask(f0, s0, f1, s1)

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
