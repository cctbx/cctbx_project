#!/usr/bin/env python
# FormatCBFFullPilatusDLSI19.py
#
#   Copyright (C) 2016 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Pilatus implementation of fullCBF format, for use with Dectris detectors.

from __future__ import division

from dxtbx.format.FormatCBFFullPilatus import FormatCBFFullPilatus
from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask

class FormatCBFFullPilatusDLSI19(FormatCBFFullPilatus):
  '''An image reading class for full CBF format images from Pilatus
  detectors. For DLS I19'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an CBF format image, i.e. we can
    make sense of it.'''

    return False

    header = FormatCBFFullPilatus.get_cbf_header(image_file)

    for record in header.split('\n'):
      if 'DLS_I19 synchrotron' in record:
        return True

    return False

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file.'''

    assert(self.understand(image_file))

    FormatCBFFullPilatus.__init__(self, image_file, **kwargs)

    self._raw_data = None
    self._cif_header_dictionary = { }

    return

  def _start(self):
    '''Open the image file as a cbf file handle, and keep this somewhere
    safe.'''

    FormatCBFFullPilatus._start(self)

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

    return

  def _beam(self):
    '''Return a working beam instance. Override polarization to be 0.999.'''

    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    return self._beam_factory.make_polarized_beam(
            sample_to_source=(0.000499745,-0.000706688,1),
            wavelength=wavelength,
            polarization=(0,1,0),
            polarization_fraction=0.999)

  def _detector(self):
    '''Return a working detector instance, with added mask regions.'''

    from scitbx import matrix
    from cctbx.eltbx import attenuation_coefficient
    from dxtbx.model import ParallaxCorrectedPxMmStrategy

    cbf_handle = self._get_cbf_handle()
    cbf_detector = cbf_handle.construct_detector(0)

    fast = matrix.col((0.00116567, 0.999989, 0.00446934))
    slow = matrix.col((0.999999, -0.00116139, -0.000959859))
    origin = matrix.col(tuple(cbf_detector.get_pixel_coordinates(0, 0)))
    centre = matrix.col((0.0005, 1.8756, 1.7050))
    two_theta = float(self._cif_header_dictionary['Detector_2theta'].split()[0])

    cbf_detector.__swig_destroy__(cbf_detector)
    del(cbf_detector)

    distance = float(
        self._cif_header_dictionary['Detector_distance'].split()[0])

    # FIXME distance offset

    beam_xy = self._cif_header_dictionary['Beam_xy'].replace(
        '(', '').replace(')', '').replace(',', '').split()[:2]

    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    beam_x, beam_y = map(float, beam_xy)

    pixel_xy = self._cif_header_dictionary['Pixel_size'].replace(
        'm', '').replace('x', '').split()

    pixel_x, pixel_y = map(float, pixel_xy)

    # FIXME construct correct origin from beam x, y pixel position


    # FIXME apply correct two-theta shift

    thickness = float(
      self._cif_header_dictionary['Silicon'].split()[2]) * 1000.0

    nx = int(
        self._cif_header_dictionary['X-Binary-Size-Fastest-Dimension'])
    ny = int(
        self._cif_header_dictionary['X-Binary-Size-Second-Dimension'])

    overload = int(self._cif_header_dictionary['Count_cutoff'].split()[0])
    underload = -1

    # take into consideration here the thickness of the sensor also the
    # wavelength of the radiation (which we have in the same file...)
    from cctbx.eltbx import attenuation_coefficient
    table = attenuation_coefficient.get_table("Si")
    mu = table.mu_at_angstrom(wavelength) / 10.0
    t0 = thickness

    detector = self._detector_factory.make_detector(
      'PAD', fast.elems, slow.elems, origin.elems,
      (1000 * pixel_x, 1000 * pixel_y),
      (nx, ny), trusted_range=(underload, overload),
      px_mm=ParallaxCorrectedPxMmStrategy(mu, t0))

    for f0, s0, f1, s1 in determine_pilatus_mask(detector):
      detector[0].add_mask(f0, s0, f1, s1)

    detector[0].set_thickness(thickness)
    detector[0].set_material('Si')
    detector[0].set_mu(mu)

    return detector

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
    # (pixels are accepted if they are inside of the active areas AND inside
    # of the trusted range)
    return tuple([m & tm for m, tm in zip(mask, trusted_mask)])

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFFullPilatus.understand(arg)
