#!/usr/bin/env python
# FormatSMVHamamatsuSPring8BL32XU.py
#
#  Copyright (C) (2015) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division

from dxtbx.format.FormatSMVHamamatsu import FormatSMVHamamatsu

class FormatSMVHamamatsuSPring8BL32XU(FormatSMVHamamatsu):

  @staticmethod
  def understand(image_file):

    size, header = FormatSMVHamamatsu.get_smv_header(image_file)

    wanted_header_items = ['DETECTOR_NAME']

    for header_item in wanted_header_items:
      if not header_item in header:
        return 0

    return header["DETECTOR_NAME"] == "Hamamatsu C10158DK"

  def _start(self):

    FormatSMVHamamatsu._start(self)
    from iotbx.detectors.hamamatsu import HamamatsuImage
    self.detectorbase = HamamatsuImage(self._image_file)
    self.detectorbase.readHeader()

  def _goniometer(self):
    '''Return a model for a simple single-axis reversed direction goniometer.'''

    return self._goniometer_factory.single_axis_reverse()
