#!/usr/bin/env python
#
# FormatNexus.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dxtbx.format.FormatHDF5 import FormatHDF5
from dxtbx.model import Beam # import dependency
from dxtbx.model import Detector # import dependency
from dxtbx.model import Goniometer # import dependency
from dxtbx.model import Scan # import dependency
from dxtbx.format.nexus import is_nexus_file
from dxtbx.format.nexus import NXmxReader
from dxtbx.format.nexus import BeamFactory
from dxtbx.format.nexus import DetectorFactory
from dxtbx.format.nexus import GoniometerFactory
from dxtbx.format.nexus import ScanFactory
from dxtbx.format.nexus import DataFactory

class FormatNexus(FormatHDF5):

  def __init__(self, image_file):
    assert(self.understand(image_file))
    FormatHDF5.__init__(self, image_file)

  @staticmethod
  def understand(image_file):
    try:
      is_nexus = is_nexus_file(image_file)
    except IOError, e:
      return False
    return is_nexus

  def _start(self):

    # Read the file structure
    self._reader = NXmxReader(self._image_file)

    # Only support 1 set of models at the moment
    assert len(self._reader.entries) == 1, \
      "Currently only supports 1 NXmx entry"
    assert len(self._reader.entries[0].data) == 1, \
      "Currently only supports 1 NXdata"
    assert len(self._reader.entries[0].instruments) == 1, \
      "Currently only supports 1 NXinstrument"
    assert len(self._reader.entries[0].samples) == 1, \
      "Currently only supports 1 NXsample"
    assert len(self._reader.entries[0].instruments[0].detectors) == 1, \
      "Currently only supports 1 NXdetector"
    assert len(self._reader.entries[0].instruments[0].detectors[0].modules) == 1, \
      "Currently only supports 1 NXdetector_module"
    assert len(self._reader.entries[0].samples[0].beams) == 1, \
      "Currently only supports 1 NXbeam"

    # Get the NXmx model objects
    entry = self._reader.entries[0]
    instrument = entry.instruments[0]
    detector = instrument.detectors[0]
    sample = entry.samples[0]
    beam = sample.beams[0]
    data = entry.data[0]

    # Construct the models
    self._beam_model = BeamFactory(beam).model
    self._detector_model = DetectorFactory(detector, self._beam_model).model
    self._goniometer_model = GoniometerFactory(sample).model
    self._scan_model = ScanFactory(sample).model
    self._raw_data = DataFactory(data).model

  def _end(self):
    return

  def _goniometer(self):
    return self._goniometer_model

  def _detector(self):
    return self._detector_model

  def _beam(self):
    return self._beam_model

  def _scan(self):
    return self._scan_model

  def get_goniometer(self, index=None):
    return self._goniometer()

  def get_detector(self, index=None):
    return self._detector()

  def get_beam(self, index=None):
    return self._beam()

  def get_scan(self, index=None):
    if index is None:
      return self._scan()
    return self._scan()[index]

  def get_raw_data(self, index):
    return self._raw_data[index]

  def get_num_images(self):
    return self._scan().get_num_images()

  def get_image_file(self, index=None):
    return self._image_file

  def get_detectorbase(self, index=None):
    return None
