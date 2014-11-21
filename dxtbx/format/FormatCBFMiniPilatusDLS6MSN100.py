#!/usr/bin/env python
# FormatCBFMiniPilatusDLS6MSN100.py
#   Copyright (C) 2014 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# An implementation of the CBF image reader for Pilatus images, from the Pilatus
# 6M SN 100 currently on Diamond I04.

from __future__ import division

from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus
from dxtbx.model import ParallaxCorrectedPxMmStrategy
from dxtbx.format.FormatPilatusHelpers import determine_pilatus_mask

class FormatCBFMiniPilatusDLS6MSN100(FormatCBFMiniPilatus):
  '''A class for reading mini CBF format Pilatus images for 6M SN 100 @ DLS.'''

  @staticmethod
  def understand(image_file):
    '''Check to see if this looks like an Pilatus mini CBF format image,
    i.e. we can make sense of it.'''

    header = FormatCBFMiniPilatus.get_cbf_header(image_file)

    for record in header.split('\n'):
      if '# Detector' in record and \
             'PILATUS' in record and 'S/N 60-0100 Diamond' in header:
        return True

    return False

  def __init__(self, image_file):
    '''Initialise the image structure from the given file, including a
    proper model of the experiment.'''

    # assert(self.understand(image_file))

    FormatCBFMiniPilatus.__init__(self, image_file)

    return

  # FIXME this beamline has a kappa goniometer so should really be supporting this
  # in here...

  def _detector(self):
    '''Detector model, allowing for small offsets in the positions of 60
    detector modules'''

    # Module positional offsets in x, y, in pixels - for the moment ignoring the
    # rotational offsets as these are not well defined. To be honest these
    # positional offsets are also not well defined as I do not know how they
    # should be applied...

    x = {
      (0, 0):  -0.477546,   (0, 1):   0.130578,   (0, 2):   0.045041,
      (0, 3):  -0.439872,   (0, 4):  -0.382077,   (1, 0):   0.087405,
      (1, 1):   0.249597,   (1, 2):   0.184265,   (1, 3):   0.158342,
      (1, 4):   0.025225,   (2, 0):  -0.179892,   (2, 1):  -0.010974,
      (2, 2):  -0.139207,   (2, 3):   0.282851,   (2, 4):  -0.442219,
      (3, 0):  -0.185027,   (3, 1):   0.218601,   (3, 2):   0.092585,
      (3, 3):    0.35862,   (3, 4):   -0.29161,   (4, 0):   0.145368,
      (4, 1):   0.609289,   (4, 2):   0.396265,   (4, 3):    0.41625,
      (4, 4):    0.07152,   (5, 0):   0.247142,   (5, 1):   0.046563,
      (5, 2):   0.248714,   (5, 3):  -0.044628,   (5, 4):  -0.391509,
      (6, 0):   0.516643,   (6, 1):   0.358453,   (6, 2):   0.069219,
      (6, 3):   0.095861,   (6, 4):  -0.167403,   (7, 0):  -0.381352,
      (7, 1):   -0.35338,   (7, 2):   0.348656,   (7, 3):   0.024543,
      (7, 4):   0.328706,   (8, 0):   0.150886,   (8, 1):   0.244987,
      (8, 2):  -0.102911,   (8, 3):    0.16633,   (8, 4):   0.386622,
      (9, 0):   0.037924,   (9, 1):   0.314392,   (9, 2):   0.238818,
      (9, 3):   0.815028,   (9, 4):  -0.048818,  (10, 0):  -0.670524,
      (10, 1):  -0.304119,  (10, 2):   0.252284,  (10, 3):   -0.05485,
      (10, 4):  -0.355264,  (11, 0):  -0.404947,  (11, 1):  -0.020622,
      (11, 2):   0.648473,  (11, 3):  -0.277175,  (11, 4):  -0.711951
    }

    y = {
      (0, 0):  -0.494797,   (0, 1):  -0.212976,   (0, 2):   0.085351,
      (0, 3):    0.35494,   (0, 4):   0.571189,   (1, 0):  -0.421708,
      (1, 1):   0.061914,   (1, 2):   0.238996,   (1, 3):   0.146692,
      (1, 4):   0.407145,   (2, 0):  -0.313212,   (2, 1):  -0.225025,
      (2, 2):   0.031613,   (2, 3):  -0.047839,   (2, 4):    0.42716,
      (3, 0):  -0.361193,   (3, 1):   0.057663,   (3, 2):   0.022357,
      (3, 3):   0.062717,   (3, 4):   0.150611,   (4, 0):   0.035511,
      (4, 1):  -0.271567,   (4, 2):   0.007761,   (4, 3):  -0.124021,
      (4, 4):   0.093017,   (5, 0):  -0.238897,   (5, 1):  -0.179724,
      (5, 2):  -0.113608,   (5, 3):   0.017841,   (5, 4):  -0.012933,
      (6, 0):  -0.166337,   (6, 1):  -0.272922,   (6, 2):  -0.194665,
      (6, 3):  -0.058535,   (6, 4):  -0.405404,   (7, 0):  -0.318824,
      (7, 1):  -0.311276,   (7, 2):  -0.205223,   (7, 3):  -0.292664,
      (7, 4):  -0.474762,   (8, 0):  -0.039504,   (8, 1):  -0.239887,
      (8, 2):  -0.343485,   (8, 3):  -0.459429,   (8, 4):  -0.426901,
      (9, 0):  -0.187805,   (9, 1):   0.282727,   (9, 2):  -0.601164,
      (9, 3):  -0.467605,   (9, 4):  -0.589271,  (10, 0):   0.028311,
      (10, 1):  -0.391571,  (10, 2):  -0.463112,  (10, 3):  -0.358092,
      (10, 4):  -0.285396,  (11, 0):    0.01863,  (11, 1):  -0.380099,
      (11, 2):  -0.234953,  (11, 3):  -0.593992,  (11, 4):  -0.801247
    }

    distance = float(
        self._cif_header_dictionary['Detector_distance'].split()[0])

    beam_xy = self._cif_header_dictionary['Beam_xy'].replace(
        '(', '').replace(')', '').replace(',', '').split()[:2]

    beam_x, beam_y = map(float, beam_xy)

    wavelength = float(
        self._cif_header_dictionary['Wavelength'].split()[0])

    pixel_xy = self._cif_header_dictionary['Pixel_size'].replace(
        'm', '').replace('x', '').split()

    pixel_x, pixel_y = map(float, pixel_xy)

    thickness = float(
      self._cif_header_dictionary['Silicon'].split()[2]) * 1000.0

    nx = int(
        self._cif_header_dictionary['X-Binary-Size-Fastest-Dimension'])
    ny = int(
        self._cif_header_dictionary['X-Binary-Size-Second-Dimension'])

    overload = int(
        self._cif_header_dictionary['Count_cutoff'].split()[0])
    underload = -1

    # take into consideration here the thickness of the sensor also the
    # wavelength of the radiation (which we have in the same file...)
    from cctbx.eltbx import attenuation_coefficient
    table = attenuation_coefficient.get_table("Si")
    mu = table.mu_at_angstrom(wavelength) / 10.0
    t0 = thickness

    # FIXME would also be very nice to be able to take into account the
    # misalignment of the individual modules given the calibration information...

    detector = self._detector_factory.simple(
        'PAD', distance * 1000.0, (beam_x * pixel_x * 1000.0,
                                   beam_y * pixel_y * 1000.0), '+x', '-y',
        (1000 * pixel_x, 1000 * pixel_y),
        (nx, ny), (underload, overload), [],
        ParallaxCorrectedPxMmStrategy(mu, t0))

    for f0, s0, f1, s1 in determine_pilatus_mask(detector):
      detector[0].add_mask(f0, s0, f1, s1)

    return detector

if __name__ == '__main__':

  import sys

  for arg in sys.argv[1:]:
    print FormatCBFMiniPilatusDLS6MSN100.understand(arg)
