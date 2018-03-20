#!/usr/bin/env python
# FormatCBFFullBrukerPhotonII.py
#
# Bruker Photon II implementation of fullCBF format.

from __future__ import absolute_import, division, print_function

from dxtbx.format.FormatCBFFull import FormatCBFFull

class FormatCBFFullBrukerPhotonII(FormatCBFFull):
  '''An image reading class for full CBF format images from Bruker Photon II
  detectors.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an CBF format image, i.e. we can
    make sense of it.'''

    header = FormatCBFFull.get_cbf_header(image_file)

    recognise = 0
    records = header.split('\n')

    for i, record in enumerate(records):
      if '_diffrn_measurement.diffrn_id' in record and 'BRUKER' in record:
        recognise += 1

      if '_array_element_size.size' in record:
        try:
          size1 = float(records[i+1].split()[-1])
          size2 = float(records[i+2].split()[-1])
        except ValueError:
          return False
        if size1 == 0.000135 and size2 == 0.000135:
          recognise += 1

    if recognise == 2: return True
    return False

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatCBFFull.__init__(self, image_file, **kwargs)

    self._raw_data = None

    return

  def _start(self):
    '''Open the image file as a cbf file handle, and keep this somewhere
    safe.'''

    FormatCBFFull._start(self)

    return

  def _goniometer(self):
    '''Return a default goniometer instance.'''

    return self._goniometer_factory.single_axis()

  def _beam(self):
    '''Return a default beam instance.'''

    wavelength = self._cbf_handle.get_wavelength()
    beam = self._beam_factory.simple(wavelength)
    return beam

  def _detector(self):
    '''Return a default detector instance.'''

    detector = self._detector_factory.imgCIF_H(self._get_cbf_handle(),
                                               'PAD')

    # no idea if this is right
    gain = self._cbf_handle.get_gain(0)[0]
    for panel in detector:
      panel.set_gain(gain)
      panel.set_trusted_range((-1,65535))

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

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print(FormatCBFFullBrukerPhotonII.understand(arg))
