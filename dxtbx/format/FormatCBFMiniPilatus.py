#!/usr/bin/env python
# FormatCBFMiniPilatus.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Pilatus images. Inherits from
# FormatCBFMini.

from __future__ import division, print_function

import os

from dxtbx.format.FormatCBFMini import FormatCBFMini
from dxtbx.format.FormatCBFMiniPilatusHelpers import get_pilatus_timestamp
from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask

class FormatCBFMiniPilatus(FormatCBFMini):
  '''A class for reading mini CBF format Pilatus images, and correctly
  constructing a model for the experiment from this.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Pilatus mini CBF format image,
    i.e. we can make sense of it.'''

    if 'ENABLE_PHOTON_FACTORY_TWO_EIGER' in os.environ:
      return False

    header = FormatCBFMini.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# Detector' in record:
        if 'EIGER' in record.upper():
          return False
        if 'TIMEPIX' in record.upper():
          return False

    for record in header.split('\n'):
      if '_array_data.header_convention' in record and \
             'PILATUS' in record:
        return True
      if '_array_data.header_convention' in record and \
             'SLS' in record:
        return True
      if '_array_data.header_convention' in record and \
             '?' in record:
        return True
      if '# Detector' in record and \
             'PILATUS' in record:  #CBFlib v0.8.0 allowed
        return True

    return False

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatCBFMini.__init__(self, image_file, **kwargs)

    return

  def _start(self):
    FormatCBFMini._start(self)

  def _detector(self):
    '''Return a model for a simple detector, presuming no one has
    one of these on a two-theta stage. Assert that the beam centre is
    provided in the Mosflm coordinate frame.'''
    detector = FormatCBFMini._detector(self)

    for f0, f1, s0, s1 in determine_pilatus_mask(detector):
      detector[0].add_mask(f0-1, s0-1, f1, s1)

    return detector

  def _scan(self):
    '''Return the scan information for this image.'''

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

  def get_mask(self, goniometer=None):
    from scitbx.array_family import flex
    detector = self.get_detector()
    mask = [flex.bool(flex.grid(reversed(p.get_image_size())), True)
             for p in detector]
    for i, p in enumerate(detector):
      untrusted_regions = p.get_mask()
      for j, (f0, s0, f1, s1) in enumerate(untrusted_regions):
        sub_array = flex.bool(flex.grid(s1-s0, f1-f0), False)
        mask[i].matrix_paste_block_in_place(sub_array, s0, f0)

    if len(detector) == 1:
      raw_data = [self.get_raw_data()]
    else:
      raw_data = self.get_raw_data()
      assert(len(raw_data) ==  len(detector))
    trusted_mask = [
      p.get_trusted_range_mask(im) for im, p in zip(raw_data, detector)]

    # returns merged untrusted pixels and active areas using bitwise AND (pixels are accepted
    # if they are inside of the active areas AND inside of the trusted range)
    return tuple([m & tm for m, tm in zip(mask, trusted_mask)])

  def get_vendortype(self):
    from dxtbx.format.FormatPilatusHelpers import get_vendortype as gv
    return gv(self.get_detector())

  @staticmethod
  def as_file(detector,beam,gonio,scan,data,path,
              header_convention="PILATUS_1.2",det_type="PILATUS3 6M"):
    FormatCBFMini.as_file(detector,beam,gonio,scan,data,path,
              header_convention,det_type)
if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatCBFMiniPilatus.understand(arg))
