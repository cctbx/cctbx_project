#!/usr/bin/env python
# FormatRAXIS.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Implementation of an ImageFormat class to read RAXIS format image which,
# since it is only really used for one brand of detector, will be the only
# class in the pile.
#
# Following taken from:
#
# http://www.rigaku.com/downloads/software/readimage.html
#
# Header Structure:
#
# First 256 bytes
#
# 52 characters: name (10) version (10) xname (20) system (12)
# 6 floats: a, b, c, alpha, beta, gamma
# 12 characters: spacegroup
# float: mosaic
# 80 characters: notes
# 84 characters: padding
#
# Next 256 bytes:
#
# 36 characters: date (12) username (20) target (4)
# float wavelength
# 20 character: typo of mono
# float: mono 2theta (degrees)
# 24 characters: collimator settings (20) filter type (4)
# 3 floats: distance, voltage, current
# 92 characters: focus stuff (12) xray chunder (80)
# long int: ipshape - 0 flat (good) 1 cylinder (bad)
# float weiss (1?)
# 56 characters: padding
#
# Next 256 bytes:
#
# 8 bytes (xtal axes)
# 3 float - phi0, phistart, phiend
# long int frame no
# float exposure time (minutes!)
# 6 floats - beam x, beam z omega chi two theta mu
# 204 characters - padding
#
# Next 256 bytes:
#
# 2 long ints - nx, nz
# 2 float - pixel size
# 4 long, 3 float (not really useful)
# 20 char host (10) ip type (10)
# 3 long - scan order horizontal, vertical, back / front:
#
#  long  dr_x;    /* horizontal scanning code: 0=left->right, 1=>right->left */
#  long  dr_z;    /* vertical scanning code: 0=down->up, 1=up->down */
#  long  drxz;    /* front/back scanning code: 0=front, 1=back */
#
# 2 useless floats
# 2 long - a magic number to show validity of this section, and number of axes
# 15 floats - up to 5 goniometer axis vectors
# 5 floats - up to 5 start angles
# 5 floats - up to 5 end angles
# 5 floats - up to 5 offset values
# long - which axis is scan axis
# 40 characters - gonio axis names, space or comma separated
#
# Then some more chunder follows - however I don't think it contains anything
# useful. So need to read first 1K of the image header.

from __future__ import division

import time
import datetime
import struct
import math

from dxtbx.format.Format import Format

class FormatRAXIS(Format):
  '''A class to support the RAXIS detector format from Rigaku.'''

  @staticmethod
  def understand(image_file):
    '''See if this looks like an RAXIS format image - clue is first
    5 letters of file should be RAXIS.'''

    if Format.open_file(image_file).read(5) == 'RAXIS':
      return True

    return False

  def __init__(self, image_file):
    assert(self.understand(image_file))

    Format.__init__(self, image_file)

    return

  def _start(self):
    self._header_bytes = Format.open_file(self._image_file).read(1024)

    if self._header_bytes[812:822].strip() in ['SGI', 'IRIS']:
      self._f = '>f'
      self._i = '>i'
    else:
      self._f = '<f'
      self._i = '<i'

    from iotbx.detectors.raxis import RAXISImage
    self.detectorbase = RAXISImage(self._image_file)
    self.detectorbase.readHeader()

  def _goniometer(self):
    '''Return a model for the goniometer from the values stored in the
    header. Assumes same reference frame as is used for the Saturn
    instrument.'''

    # OK for the moment I am not completely clear on how omega, chi and
    # phi are defined so for the moment assert (i) that only one rotation
    # axis has a non-zero offset (ii) that this is the scan axis and
    # (iii) that this is aligned with (1, 0, 0) in the file... N.B.
    # with these scanners it is (I believe) always the case that the
    # (1, 0, 0) direction points downwards from the sample to the
    # goniometer.

    i = self._i
    f = self._f
    header = self._header_bytes

    # check magic number to see if this section is valid
    if struct.unpack(i, header[852:856]) == 1:
      n_axes = struct.unpack(i, header[856:860])[0]
    else:
      n_axes = 0
    scan_axis = struct.unpack(i, header[980:984])[0]

    for j in range(n_axes):

      axis_x = struct.unpack(f, header[860 + j * 12:864 + j * 12])[0]
      axis_y = struct.unpack(f, header[864 + j * 12:868 + j * 12])[0]
      axis_z = struct.unpack(f, header[868 + j * 12:872 + j * 12])[0]

      axis_start = struct.unpack(f, header[920 + j * 4:924 + j * 4])[0]
      axis_end = struct.unpack(f, header[940 + j * 4:944 + j * 4])[0]
      axis_offset = struct.unpack(f, header[960 + j * 4:964 + j * 4])[0]

      if j == scan_axis:
        assert(math.fabs(axis_x - 1) < 0.001)
        assert(math.fabs(axis_y) < 0.001)
        assert(math.fabs(axis_z) < 0.001)
      else:
        assert(math.fabs(axis_start % 180.0) < 0.001)
        assert(math.fabs(axis_end % 180.0) < 0.001)

    return self._goniometer_factory.single_axis()

  def _detector(self):
    '''Return a model for the detector as defined in the image header,
    with the additional knowledge about how things are arranged i.e. that
    the principle rotation axis vector points from the sample downwards.'''

    # As above will have to make a bunch of assumptions about how the
    # detector is set up, namely that the rotation axis points downwards,
    # the fast axis points along and the slow axis points up. These will
    # (as well as I can understand) be asserted. Also assuming that the
    # two-theta axis is along (1, 0, 0) which is all that makes sense.

    i = self._i
    f = self._f
    header = self._header_bytes

    det_h = struct.unpack(i, header[832:836])[0]
    det_v = struct.unpack(i, header[836:840])[0]
    det_f = struct.unpack(i, header[840:844])[0]

    assert(det_h == 0)
    assert(det_v == 0)
    assert(det_f == 0)

    nx = struct.unpack(i, header[768:772])[0]
    ny = struct.unpack(i, header[772:776])[0]

    dx = struct.unpack(f, header[776:780])[0]
    dy = struct.unpack(f, header[780:784])[0]

    distance = struct.unpack(f, header[344:348])[0]
    two_theta = struct.unpack(f, header[556:560])[0]

    beam_x = struct.unpack(f, header[540:544])[0]
    beam_y = struct.unpack(f, header[544:548])[0]

    beam = (beam_x * dx, beam_y * dy)

    return self._detector_factory.two_theta(
        'IMAGE_PLATE', distance, beam, '+y', '-x', '+x', two_theta,
        (dx, dy), (nx, ny), (0, 1000000), [])

  def _beam(self):
    '''Return a simple model for the beam. This assumes a laboratory source
    which has an unpolarized beam.'''

    wavelength = struct.unpack(self._f, self._header_bytes[292:296])[0]

    return self._beam_factory.complex((0.0, 0.0, 1.0), 0.5,
                                      (0.0, 1.0, 0.0), wavelength)

  def _scan(self):
    '''Return the scan information for this image.'''
    import calendar
    i = self._i
    f = self._f
    header = self._header_bytes

    # We expect an invalid goniometer section, indicated by wrong magic number
    magic_no = struct.unpack(">i", header[852:856])

    format = self._scan_factory.format('RAXIS')
    exposure_time = struct.unpack(f, header[536:540])[0]

    y, m, d = map(int, header[256:268].strip().split('-'))

    epoch = calendar.timegm(datetime.datetime(y, m, d, 0, 0, 0).timetuple())

    if magic_no == 1:
      s = struct.unpack(i, header[980:984])[0]
      osc_start = struct.unpack(f, header[920 + s * 4:924 + s * 4])[0]
      osc_end = struct.unpack(f, header[940 + s * 4:944 + s * 4])[0]

    else:
      osc_dat = struct.unpack(f, header[520:524])[0]
      osc_start = osc_dat + struct.unpack(f, header[524:528])[0]
      osc_end = osc_dat + struct.unpack(f, header[528:532])[0]

    osc_range = osc_end - osc_start

    return self._scan_factory.single(
        self._image_file, format, exposure_time,
        osc_start, osc_range, epoch)
