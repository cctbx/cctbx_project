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

__mask = None

# if group_rows == True, then interpret data as 24 panels, where each row
# of 5 panels is grouped as one "panel"
# elif group_rows == False, then interpret data as 120 panels,
# 24 rows * 5 columns - bodge through ENV variable in 1st cut...

import os
if 'P12M_120_PANEL' in os.environ:
  group_rows = False
else:
  group_rows = True

def read_mask():
  global __mask
  if not __mask:
    import os
    import cPickle as pickle
    from scitbx.array_family import flex
    source_dir = os.path.split(__file__)[0]
    mask_file = os.path.join(source_dir, 'FormatCBFMiniPilatusDLS12M.pickle')
    __mask = pickle.load(open(mask_file, 'rb'))
  return __mask

def read_cbf_image(cbf_image):
  from cbflib_adaptbx import uncompress
  import binascii
  from scitbx.array_family import flex

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


  isel = read_mask()
  pixel_values.as_1d().set_selected(isel, -2)

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

    x = matrix.col((1, 0, 0))
    y = matrix.col((0, 1, 0))
    z = matrix.col((0, 0, 1))

    obs_beam_y = 2587
    ideal_beam_y = 2594
    beam_shift_y = 0.172 * (2594 - 2587)

    distance = float(
        self._cif_header_dictionary['Detector_distance'].split()[0]) * 1000.0

    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    thickness = float(
      self._cif_header_dictionary['Silicon'].split()[2]) * 1000.0

    detector = HierarchicalDetector()
    root = detector.hierarchy()
    root.set_frame(
      (1, 0, 0),
      (0, 1, 0),
      (0, 0, - (250 + distance)))

    from cctbx.eltbx import attenuation_coefficient
    table = attenuation_coefficient.get_table("Si")
    mu = table.mu_at_angstrom(wavelength) / 10.0
    t0 = thickness
    px_mm = ParallaxCorrectedPxMmStrategy(mu, t0)

    for j in range(24):

      angle = math.pi * (-12.2 + 0.5 * 7.903 + j * (7.903 + 0.441)) / 180.0
      fast = matrix.col((-1, 0, 0))
      slow = matrix.col((0, math.sin(angle), - math.cos(angle)))
      normal = fast.cross(slow)
      # for longer wavelength data sets move 192.3 below to 184.9
      if wavelength < 1.128:
        off_x = 191.9
      else:
        off_x = 184.9

      if group_rows:
        origin = 250.0 * normal - off_x * fast - 16.8 * slow + 250 * z + \
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
        p.set_px_mm_strategy(px_mm)

      else:
        row_origin = 250.0 * normal - off_x * fast - 16.8 * slow + 250 * z + \
            beam_shift_y * y

        for i in range(5):
          p = detector.add_panel()
          origin = row_origin + i * (487+7) * 0.172 * fast

          # OBS! you need to set the panel to a root before set local frame...
          root.add_panel(p)
          p.set_name('row-%02d' % j)
          p.set_image_size((487, 195))
          p.set_trusted_range((-1, 1000000))
          p.set_pixel_size((0.172, 0.172))
          p.set_local_frame(
            fast.elems,
            slow.elems,
            origin.elems)
          p.set_px_mm_strategy(px_mm)

    return detector

  def get_raw_data(self, index=None):
    if self._raw_data is None:
      data = read_cbf_image(self._image_file)
      self._raw_data = []

      if group_rows:
        for j, panel in enumerate(self.get_detector()):
          shift_y = 195 + 17
          xmin, ymin, xmax, ymax = 0, j * shift_y, 2463, j * shift_y + 195
          self._raw_data.append(data[ymin:ymax,xmin:xmax])
      else:
        for j in range(24):
          for i in range(5):
            shift_y = 195 + 17
            shift_x = 487 + 7
            xmin, ymin, xmax, ymax \
              = i * shift_x, j * shift_y, i * shift_x + 487, j * shift_y + 195
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
