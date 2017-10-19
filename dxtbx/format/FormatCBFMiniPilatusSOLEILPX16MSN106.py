#!/usr/bin/env python
# FormatCBFMiniPilatusSOLEILPX16MSN106.py
#   Copyright (C) 2014 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Set up for Soleil PX1, with full kappa goniometer

from __future__ import absolute_import, division

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
#from dxtbx.model import ParallaxCorrectedPxMmStrategy
#from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask

def read_cbf_image_as_double(cbf_image):
  from cbflib_adaptbx import uncompress
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

  return pixel_values.as_double()

class FormatCBFMiniPilatusSOLEILPX16MSN106(FormatCBFMiniPilatus):
  '''A class for reading mini CBF format Pilatus images for 6M SN 106 @
  Synchrotrol Soleil PX1.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Pilatus mini CBF format image,
    i.e. we can make sense of it.'''

    header = FormatCBFMiniPilatus.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# Detector' in record and \
             'PILATUS' in record and 'S/N 60-0106, Soleil' in record:
        return True

    return False

  def __init__(self, image_file, **kwargs):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)

    FormatCBFMiniPilatus.__init__(self, image_file, **kwargs)

    self._raw_data = None

    return

  def _goniometer(self):
    '''Construct a goniometer from the records in the mini CBF header.'''

    if ('Alpha' in self._cif_header_dictionary and
        'Kappa' in self._cif_header_dictionary):
      # Kappa
      alpha = float(self._cif_header_dictionary['Alpha'].split()[0])
      omega = float(self._cif_header_dictionary['Chi'].split()[0])
      kappa = float(self._cif_header_dictionary['Kappa'].split()[0])
      phi = float(self._cif_header_dictionary['Phi'].split()[0])

      axis = self._cif_header_dictionary['Oscillation_axis']

      scanaxis  = {'OMEGA':'Omega', 'PHI':'Phi'}

      assert axis in scanaxis

      # this is the direction the arm points in at datum
      direction = '+z'

      return self._goniometer_factory.make_kappa_goniometer(
        alpha, omega, kappa, phi, direction, scanaxis[axis])

    else:
      # Smargon
      from scitbx.array_family import flex

      phi = float(self._cif_header_dictionary['Phi'].split()[0])
      chi = float(self._cif_header_dictionary['Chi'].split()[0])
      omega = float(self._cif_header_dictionary['Omega'].split()[0])

      names = flex.std_string(("PHI", "CHI", "OMEGA"))
      axes = flex.vec3_double(((1,0,0), (0,0,-1), (1,0,0)))
      angles = flex.double((phi, chi, omega))

      axis = self._cif_header_dictionary['Oscillation_axis'].upper()
      assert axis in names, axis
      scan_axis = flex.first_index(names, axis)

      return self._goniometer_factory.make_multi_axis_goniometer(
        axes, angles, names, scan_axis)

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFMiniPilatusSOLEILPX16MSN106.understand(arg)
