#!/usr/bin/env python
# FormatCBFMiniPilatusDLS12M.py
#
#   Copyright (C) 2014 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Pilatus images, for the P12M-DLS

from __future__ import division

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.model import ParallaxCorrectedPxMmStrategy
from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask

def read_cbf_image(cbf_image):
  from scitbx.array_family import flex
  from cbflib_adaptbx import uncompress, compress
  import binascii

  start_tag = binascii.unhexlify('0c1a04d5')

  data = open(cbf_image, 'rb').read()
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

class FormatCBFMiniPilatusDLS12M(FormatCBFMiniPilatus):

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Pilatus mini CBF format image,
    i.e. we can make sense of it.'''

    header = FormatCBFMiniPilatus.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# Detector' in record and \
             'PILATUS' in record and 'S/N 120-0100' in header:
        return True

    return False

  def __init__(self, image_file):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    assert(self.understand(image_file))

    FormatCBFMiniPilatus.__init__(self, image_file)

    self._raw_data = None

    return

  def _detector(self):

    # module positions from detector blueprints - modelling at the moment as
    # 24 modules, each consisting of 5 sensors (the latter is ignored)

    from dxtbx.model.detector import HierarchicalDetector
    from scitbx import matrix
    import math

    detector = HierarchicalDetector()
    root = detector.hierarchy()
    root.set_frame(
      (1, 0, 0),
      (0, 1, 0),
      (0, 0, -250))

    x = matrix.col((1, 0, 0))
    y = matrix.col((0, 1, 0))
    z = matrix.col((0, 0, 1))

    obs_beam_y = 2587
    ideal_beam_y = 2594
    beam_shift_y = 0.172 * (2594 - 2587)

    for j in range(24):

      angle = math.pi * (-12.2 + 0.5 * 7.903 + j * (7.903 + 0.441)) / 180.0
      fast = matrix.col((-1, 0, 0))
      slow = matrix.col((0, math.sin(angle), - math.cos(angle)))
      normal = fast.cross(slow)
      # from observation of beam image on panel 12-down 3-across @ 1117,2587
      # FIXME this should be determined from the drawings & then applied as
      # a bulk shift after the detector is constructed
      origin = 250.0 * normal - 192.3 * fast - 16.8 * slow + 250 * z + \
          beam_shift_y * y
      p = detector.add_panel()

      # OBS! you need to set the panel to a root before set local frame...
      root.add_panel(p)
      p.set_name('row-%02d' % j)
      p.set_image_size((2463, 195))
      p.set_trusted_range((-1, 1000000))
      p.set_pixel_size((0.172, 0.172))
      p.set_local_frame(
        fast.elems,
        slow.elems,
        origin.elems)

    return detector

  def get_raw_data(self, index=None):
    if self._raw_data is None:
      data = read_cbf_image(self._image_file)
      self._raw_data = []

      for j, panel in enumerate(self.get_detector()):
        shift_y = 195 + 17
        xmin, ymin, xmax, ymax = 0, j * shift_y, 2463, j * shift_y + 195
        self._raw_data.append(data[ymin:ymax,xmin:xmax])

    if index is not None:
      return self._raw_data[index]

    return self._raw_data[0]

  def _goniometer(self):
    '''Return a model for a simple single-axis goniometer. This should
    probably be checked against the image header.'''

    return self._goniometer_factory.single_axis_reverse()

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFMiniPilatusDLS12M.understand(arg)
