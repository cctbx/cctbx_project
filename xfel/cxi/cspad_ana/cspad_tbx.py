# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

"""Toolbox for images from the Cornell SLAC Pixel Array Detector
(CSpad).

XXX Better named cspad_common?

XXX Read out detector temperature (see Hart et al., 2012)?
"""
from __future__ import absolute_import, division, print_function
from six.moves import range

import math
import numpy
import os
import time

from libtbx import easy_pickle
from scitbx.array_family import flex
from xfel.cxi.cspad_ana.parse_calib import Section
import six
from six.moves import zip
from six.moves import map

__version__ = "$Revision$"


# The CAMP and CSpad counters are both 14 bits wide (Strüder et al
# 2010; Philipp et al., 2007), which means the physical limit is 2**14 - 1.
# However, in practice, when the pixels are in the low gain mode, after
# correcting by a gain value of around 6.87, the pixels tend to saturate
# around 90000. See xpp experiment xppe0314, run 184 as evidence.
cspad_saturated_value = 90000

# The dark average for the CSPAD detector is around 1100-1500. A pixel
# histogram of a minimum projection of an uncorrected (raw) light run shows
# a mostly flat tail up to ~800 ADU with a few bumps in the tail which
# represent true underloads. Assume a dark average of 1200 ADU. After dark
# subtraction, 800 - 1200 gives a minimum trusted value of -400. Reject
# pixels less than this.
cspad_min_trusted_value = -400

# As long as the mask value is outside of the trusted range, the pixel should
# be ignored by any downstream software.
cspad_mask_value = -100000

# The side length of a square quadrant from the old XtcExplorer code.
# XXX This should be obsoleted!
npix_quad = 850

# The pixel size in mm.  The pixel size is fixed and square, with side
# length of 110 µm (Philipp et al., 2007).  XXX Should really clarify
# this with Sol and Chris.
#
# XXX Andor: 13.5 µm square, CAMP: 75 µm, square (Strüder et al.,
# 2010)
pixel_size = 110e-3

# origin of section in quad coordinate system.  x-position
# correspond to column number.  XXX Note/reference the source!
# XXX This should be obsoleted!
xpos_sec2x1 = [[ 414,  626,    0,    0,  213,    1,  418,  419],  # 2:5 were not measured
               [ 421,  634,    0,    0,  213,    1,  424,  425],
               [ 417,  630,    0,    1,  212,    0,  425,  426],
               [ 416,  630,    0,    0,  213,    1,  420,  421]] # 2:5 were not measured
# y-position correspond to maxrows - row number
ypos_sec2x1 = [[   0,    0,  214,    1,  425,  425,  615,  402],  # 2:5 were not measured
               [   0,    0,  214,    1,  425,  425,  615,  402],
               [   0,    0,  215,    3,  431,  431,  616,  403],
               [   0,    0,  214,    1,  425,  425,  615,  403]] # 2:5 were not measured


def address_split(address, env=None):
  """The address_split() function splits an address into its four
  components.  Address strings are on the form
  detector-detectorID|device-deviceID, where the detectors must be in
  dir(xtc.DetInfo.Detector) and device must be in
  (xtc.DetInfo.Device).
  @param address Full data source address of the DAQ device
  @param env     Optional env to dereference an alias into an address
  @return        Four-tuple of detector name, detector ID, device, and
                 device ID
  """

  import re

  # pyana
  m = re.match(
    '^(?P<det>\S+)\-(?P<det_id>\d+)\|(?P<dev>\S+)\-(?P<dev_id>\d+)$', address)
  if m is not None:
    return (m.group('det'), m.group('det_id'), m.group('dev'), m.group('dev_id'))

  # psana
  m = re.match(
    '^(?P<det>\S+)\.(?P<det_id>\d+)\:(?P<dev>\S+)\.(?P<dev_id>\d+)$', address)
  if m is not None:
    return (m.group('det'), m.group('det_id'), m.group('dev'), m.group('dev_id'))

  # psana DetInfo string
  m = re.match(
    '^DetInfo\((?P<det>\S+)\.(?P<det_id>\d+)\:(?P<dev>\S+)\.(?P<dev_id>\d+)\)$', address)
  if m is not None:
    return (m.group('det'), m.group('det_id'), m.group('dev'), m.group('dev_id'))

  if env is not None:
    # Try to see if this is a detector alias, and if so, dereference it. Code from psana's Detector/PyDetector.py
    amap = env.aliasMap()
    alias_src = amap.src(address) # string --> DAQ-style psana.Src

    # if it is an alias, look up the full name
    if amap.alias(alias_src) != '':         # alias found
      address = str(alias_src)
      return address_split(address)

  return (None, None, None, None)


def cbcaa(config, sections):
  """The cbcaa() function uses on-disk calibration data to estimate
  the beam centre and the active detector areas.  The beam centre is
  approximated as the average of the four ASIC corners closest to the
  detector centre.  That is the first corner of the section 1 in every
  quadrant.  Note that first corner index is vertical coordinate,
  second index is the horizontal coordinate.  XXX Construct the active
  areas in "spotfinder format", i.e. opposing corners.  XXX This is a
  really bad function name!  XXX The beam centre may be extracted from
  the ebeam object?

  @param config   XXX
  @param sections XXX Directory with calibration information
  @return         Tuple of 2D beam centre, and active areas in
                  "spotfinder format"
  """

  aa = flex.int()
  if (sections is None):
    # The active areas of the detector, (UL_slow, UL_fast, LR_slow,
    # LR_fast) A two-by-one is 185-by-392 pixels with a 4-pixel gap.
    # An ASIC is 185-by-194 pixels.  XXX Still need to sort out the
    # inclusive/exclusive detail.  Upper-left corner is inclusive,
    # lower-right corner is exclusive.  XXX Should subtract one from
    # x, y on lower-right corner and verify with zoom.  XXX All this
    # should probably go by now
    for q in range(4): # loop over quadrants
      for i in range(8): # loop over two-by-one:s XXX variable name!

        # Skip this two-by-one if it is missing.
        if not (config.roiMask(q) & 0x1 << i):
          continue

        # XXX Note the different coordinate systems in use here!
        xpos =       xpos_sec2x1[q][i] # x-value of lower, left corner
        ypos = 850 - ypos_sec2x1[q][i] # y-value of lower, left corner

        if (i == 0 or i == 1 or i == 4 or i == 5):
          UL1_x = xpos
          UL2_x = xpos
          UL1_y = ypos - 194 - 4 - 194
          UL2_y = ypos - 194

          LR1_x = UL1_x + 185
          LR2_x = UL2_x + 185
          LR1_y = UL1_y + 194
          LR2_y = UL2_y + 194

        elif (i == 2 or i == 3 or i == 6 or i == 7):
          UL1_x = xpos
          UL2_x = xpos + 194 + 4
          UL1_y = ypos - 185
          UL2_y = ypos - 185

          LR1_x = UL1_x + 194
          LR2_x = UL2_x + 194
          LR1_y = UL1_y + 185
          LR2_y = UL2_y + 185

        # Quadrant rotations, counter-clockwise.  Zeroth quadrant
        # needs no special action.
        if (q == 0):
          pass

        elif (q == 1):
          UL1_x, UL1_y               = 850 + 850 - UL1_y, UL1_x
          LR1_x, LR1_y               = 850 + 850 - LR1_y, LR1_x

          UL2_x, UL2_y               = 850 + 850 - UL2_y, UL2_x
          LR2_x, LR2_y               = 850 + 850 - LR2_y, LR2_x

          UL1_x, LR1_x               = LR1_x, UL1_x
          UL2_x, LR2_x               = LR2_x, UL2_x

        elif (q == 2):
          UL1_x, UL1_y               = 850 + 850 - UL1_x, 850 + 850 - UL1_y
          LR1_x, LR1_y               = 850 + 850 - LR1_x, 850 + 850 - LR1_y

          UL2_x, UL2_y               = 850 + 850 - UL2_x, 850 + 850 - UL2_y
          LR2_x, LR2_y               = 850 + 850 - LR2_x, 850 + 850 - LR2_y

          UL1_x, UL1_y, LR1_x, LR1_y = LR1_x, LR1_y, UL1_x, UL1_y
          UL2_x, UL2_y, LR2_x, LR2_y = LR2_x, LR2_y, UL2_x, UL2_y

        elif (q == 3):
          UL1_x, UL1_y               = UL1_y, 850 + 850 - UL1_x
          LR1_x, LR1_y               = LR1_y, 850 + 850 - LR1_x

          UL2_x, UL2_y               = UL2_y, 850 + 850 - UL2_x
          LR2_x, LR2_y               = LR2_y, 850 + 850 - LR2_x

          UL1_y, LR1_y               = LR1_y, UL1_y
          UL2_y, LR2_y               = LR2_y, UL2_y

        # This is row-major matrix layout; FAST <=> x, SLOW <=> y.
        aa.extend(flex.int([UL1_y, UL1_x, LR1_y, LR1_x]))
        aa.extend(flex.int([UL2_y, UL2_x, LR2_y, LR2_x]))

    # The beam centre is estimated as the centre of the image.
    return ([npix_quad, npix_quad], aa)

  # Old way of computing beam center, phased out 05/19/15
  #bc     = [0, 0]

  # XXX Make up a quadrant mask for the emission detector.  Needs to
  # be checked!
  if len(sections) <= 1:
    q_mask = 1
  else:
    q_mask = config.quadMask()

  for q in range(len(sections)):
    if (not((1 << q) & q_mask)):
      continue

    # Old way of computing beam center, phased out 05/19/15
    #corner = sections[q][1].corners(True)[0]
    #bc     = [bc[0] + corner[1] / len(sections),
    #          bc[1] + corner[0] / len(sections)]

    # XXX Make up section mask for the emission detector.  Needs to be
    # checked!
    try:
      import _pdsdata
      types = _pdsdata.cspad2x2.ConfigV1, _pdsdata.cspad2x2.ConfigV2
    except ImportError:
      import psana
      types = psana.CsPad2x2.ConfigV1, psana.CsPad2x2.ConfigV2
    if len(sections) == 1 and type(config) in types:
      s_mask = config.roiMask()
    else:
      s_mask = config.roiMask(q)
    for s in range(len(sections[q])):
      if (not((1 << s) & s_mask)):
        continue
      c = sections[q][s].corners_asic()
      aa.extend(flex.int(c[0]))
      aa.extend(flex.int(c[1]))

  # The beam center was defined above as the center of the innermost 4 sensors. Recently,
  # that center has drifted too much from the true image center (Spring 2015). So, here we
  # use the true image center instead.
  return [882.5,882.5], aa


def CsPad2x2Image(data, config, sections):
  """The CsPad2x2Image() function assembles a two-dimensional image
  from the Sc1 detector readout in @p data.

  @param data     Detector readout from XTC stream
  @param config   XXX
  @param sections XXX Directory with calibration information
  @return         Assembled detector image
  """

  assert (data.shape[2] == 2)

  det  = numpy.zeros((2 * 185, 2 * 194 + 3))

  # XXX config.sections is now a function returning a list?  Since the
  # masking was disabled in December commenting out this bit does not
  # cause any further breakage XXX Does this still work for runs 4 and
  # 5?
#  s = config.sections
#  mask = map(s, range(2))

  # For this detector, the quadrant index is always zero.
  q_idx = 0
  for s in range(2):
    # XXX DAQ misconfiguration?  This mask appears not to work
    # reliably for the Sc1 detector.
#    if (s not in mask[q_idx]):
#      continue
    asics  = numpy.vsplit(numpy.rot90(data[:, :, s], -1), 2)
    gap    = numpy.zeros((3, 185), dtype = data.dtype)
    s_data = numpy.vstack((asics[0], gap, asics[1]))

    angle  = sections[q_idx][s].angle
    center = sections[q_idx][s].center
    rplace(det, s_data, angle, center)
  return (det)

def evt_get_quads(address, evt, env):
  try:
    # pyana
    quads = evt.getCsPadQuads(address, env)
  except AttributeError:
    # psana
    from psana import Source, CsPad
    src = Source(address)
    cspad = evt.get(CsPad.DataV2, src)
    if cspad is None:
      return None
    quads = [cspad.quads(i) for i in range(cspad.quads_shape()[0])]
  return quads

def CsPadDetector(address, evt, env, sections, right=True, quads=None):
  """The CsPadDetector() function assembles a two-dimensional image
  from the Ds1 detector readout in @p data3d and the calibration
  information in @p sections.  XXX General question: do
  variable/function names make sense?

  @param address  Full data source address of the DAQ device
  @param evt      Event data object, a configure object
  @param env      Environment object
  @param sections XXX Directory with calibration information
  @param right    @c True to restrict rotations to right angles
  @return         Assembled detector image
  """

  device = address_split(address)[2]
  if device is None or device != 'Cspad':
    return None

  # Get a current configure object for the detector
  config = getConfig(address, env)
  if config is None:
    return None

  # For consistency, one could/should verify that len(quads) is equal
  # to len(sections).
  if quads is None:
    quads = evt_get_quads(address, evt, env)

  if quads is None or len(quads) != len(sections):
    return None

  # This is from Mikhail S. Dubrovin's
  # HDF5Explorer/src/ConfigCSpad.py, which uses a detector size of
  # 1765-by-1765 pixels.
  extra_space = (1765 - 2 * Section.q_size[0],
                 1765 - 2 * Section.q_size[1])

  # Start out with a blank image of the detector.  This assumes that
  # the type of the first section in the first quadrant is identical
  # to the type of all the other sections.
  det = numpy.zeros((2 * Section.q_size[0] + extra_space[0],
                     2 * Section.q_size[1] + extra_space[1]),
                    dtype=quads[0].data()[0].dtype)

  ### need to swap the quadrants for data collected mid October, 2013
  evttime = time.gmtime(evt_time(evt)[0])
  swap = evttime.tm_year == 2013 and evttime.tm_mon == 10 and evttime.tm_mday >= 20 and evttime.tm_mday <= 25

  for quad in quads:
    q_data = quad.data()
    q_idx = quad.quad()
    if swap:
      q_idx = [0,1,3,2][q_idx]
    try:
      # pyana
      # example: if the third sensor (2x1) is disabled, q_mask = [0,1,3,4,5,6,7]
      q_mask = config.sections(q_idx)
    except AttributeError:
      # psana
      # as above, using config.roiMask, a bitstring where the ith bit is true if the ith sensor is active. x << y means bitwise shift
      # x, y times, and & is the bitwise AND operator
      q_mask = [i for i in range(len(sections[q_idx])) if 1 << i & config.roiMask(q_idx)]

    # For consistency, assert that there is data for each unmasked
    # section.
    assert len(q_data) == len(q_mask)
    for (s_data, s_idx) in zip(q_data, q_mask):
      # Rotate the section from the XTC-stream by -90 degrees to
      # conform to the "standing up" convention used by the
      # calibration data, and insert a 3-pixel gap between the ASIC:s.
      # This requires the horizontal dimension of the unrotated
      # section to be even.
      assert s_data.shape[1] % 2 == 0
      asics = numpy.vsplit(numpy.rot90(s_data, -1), 2)
      gap = numpy.zeros((3, s_data.shape[0]), dtype=s_data.dtype)
      s_data = numpy.vstack((asics[0], gap, asics[1]))

      # Place the section in the detector image, either by forcing
      # rotation to right angles or by interpolating.
      angle = sections[q_idx][s_idx].angle
      center = sections[q_idx][s_idx].center
      if right:
        rplace(det, s_data, angle, center)
      else:
        iplace(det, s_data, angle, center)
  return det


def CsPadElement(data3d, qn, config):
  """Construct one image for each quadrant, each with 8 sections from
  a data3d = 3 x 2*194 x 185 data array.  This function was originally
  written by Ingrid Ofte for pyana's XtcExplorer module.  XXX
  Documentation!
  """

  # If any sections are missing, insert zeros.
  mask = list(map(config.sections, range(4)))
  if (len(data3d) < 8):
    zsec = numpy.zeros((185, 388), dtype = data3d.dtype)
    for i in range(8) :
      if (i not in mask[qn]):
        data3d = numpy.insert(data3d, i, zsec, axis = 0)

  pairs = []
  for i in range(8) :
    # Insert gap between ASIC:s in the 2x1.
    asics = numpy.hsplit(data3d[i], 2)
    gap   = numpy.zeros((185, 4), dtype = data3d.dtype)
    pair  = numpy.hstack((asics[0], gap, asics[1]))

    # Sections 2,3 and 6,7 are as is.  The others need some rotation,
    # implemented as a matrix transposition here.
    if (i == 0 or i == 1):
      pair = pair[:, ::-1].T
    if (i == 4 or i == 5):
      pair = pair[::-1, :].T
    pairs.append(pair)

  # Make the array for this quadrant, and insert the 2x1 sections.
  quadrant = numpy.zeros((npix_quad, npix_quad), dtype = data3d.dtype)
  for sec in range(8):
    nrows, ncols = pairs[sec].shape

    # x,y  in quadrant coordinate system
    xpos = xpos_sec2x1[qn][sec]
    ypos = ypos_sec2x1[qn][sec]
    colp = xpos
    rowp = npix_quad - ypos
    quadrant[(rowp - nrows):rowp, colp:(colp + ncols)] = \
        pairs[sec][0:nrows, 0:ncols]

  # Finally, rotate the quadrant as needed.
  if (qn > 0):
    quadrant = numpy.rot90(quadrant, 4 - qn)
  return quadrant

def dpack(active_areas=None,
          address=None,
          beam_center_x=None,
          beam_center_y=None,
          ccd_image_saturation=None,
          data=None,
          distance=None,
          pixel_size=pixel_size,
          saturated_value=None,
          timestamp=None,
          wavelength=None,
          xtal_target=None,
          min_trusted_value=None):
  """XXX Check completeness.  Should fill in sensible defaults."""

  # Must have data.
  if data is None:
    return None

  # Create a time stamp of the current time if none was supplied.
  if timestamp is None:
    timestamp = evt_timestamp()

  # For unknown historical reasons, the dictionary must contain both
  # CCD_IMAGE_SATURATION and SATURATED_VALUE items.
  if ccd_image_saturation is None:
    if saturated_value is None:
      ccd_image_saturation = cspad_saturated_value
    else:
      ccd_image_saturation = saturated_value
  if saturated_value is None:
    saturated_value = ccd_image_saturation

  # Use a minimum value if provided for the pixel range
  if min_trusted_value is None:
    min_trusted_value = cspad_min_trusted_value

  # By default, the beam center is the center of the image.  The slow
  # (vertical) and fast (horizontal) axes correspond to x and y,
  # respectively.
  if beam_center_x is None:
    beam_center_x = pixel_size * data.focus()[1] / 2
  if beam_center_y is None:
    beam_center_y = pixel_size * data.focus()[0] / 2

  # By default, the entire detector image is an active area.  There is
  # no sensible default for distance nor wavelength.  XXX But setting
  # wavelength to zero may be disastrous?
  if active_areas is None:
    # XXX Verify order with non-square detector
    active_areas = flex.int((0, 0, data.focus()[0], data.focus()[1]))
  if distance is None:
    distance = 0
  if wavelength is None:
    wavelength = 0

  # The size must match the image dimensions.  The length along the
  # slow (vertical) axis is SIZE1, the length along the fast
  # (horizontal) axis is SIZE2.
  return {'ACTIVE_AREAS': active_areas,
          'BEAM_CENTER_X': beam_center_x,
          'BEAM_CENTER_Y': beam_center_y,
          'CCD_IMAGE_SATURATION': ccd_image_saturation,
          'DATA': data,
          'DETECTOR_ADDRESS': address,
          'DISTANCE': distance,
          'PIXEL_SIZE': pixel_size,
          'SATURATED_VALUE': saturated_value,
          'MIN_TRUSTED_VALUE': min_trusted_value,
          'SIZE1': data.focus()[0],
          'SIZE2': data.focus()[1],
          'TIMESTAMP': timestamp,
          'SEQUENCE_NUMBER': 0, # XXX Deprecated
          'WAVELENGTH': wavelength,
          'xtal_target': xtal_target}


def hdf5pack(hdf5_file,
             active_areas=None,
             address=None,
             attenuation=None,
             beam_center_x=None,
             beam_center_y=None,
             ccd_image_saturation=None,
             data=None,
             distance=None,
             pixel_size=None,
             pulse_length=None,
             saturated_value=None,
             timestamp=None,
             wavelength=None,
             xtal_target=None):
  """Similar but far from identical to the HDF5 output from CASS.  XXX
  Poor diagnostics--we don't know if it failed or not.

  @note Does not include the deprecated SEQUENCE_NUMBER attribute.
        While some redundant items are written in order to keep the
        HDF5 synchronised to the pickle format, neither SIZE1 nor
        SIZE2 are included.
  """

  # Need this because we cannot write None values to the HDF5 file.
  if address is None:
    address = repr(None)
  if attenuation is None:
    attenuation = 0
  if xtal_target is None:
    xtal_target = repr(None)
  if pixel_size is None:
    pixel_size = globals()['pixel_size'] # XXX CSpad-specific!
  if pulse_length is None:
    pulse_length = 0

  d = dpack(address=address,
            active_areas=active_areas,
            beam_center_x=beam_center_x,
            beam_center_y=beam_center_y,
            ccd_image_saturation=ccd_image_saturation,
            data=data,
            distance=distance,
            pixel_size=pixel_size,
            saturated_value=saturated_value,
            timestamp=timestamp,
            wavelength=wavelength,
            xtal_target=xtal_target)
  if d is None:
    return

  grp_event = hdf5_file.create_group(d['TIMESTAMP'])
  grp_detector = grp_event.create_group(address)
  for (key, value) in six.iteritems(d):
    if key == 'ACTIVE_AREAS':
      grp_detector.create_dataset(key, data=value.as_numpy_array())
    elif key == 'DATA':
      # Compress the image data with gzip at the default level (4).
      # CASS seems to use maximum compression level (9), which gives a
      # moderate decrease in file size at the price of much longer
      # running time.
      grp_detector.create_dataset(
        key, compression='gzip', data=value.as_numpy_array())
    else:
      grp_event.create_dataset(key, data=[value])
  grp_event.create_dataset('ATTENUATION', data=[attenuation])
  grp_event.create_dataset('PULSE_LENGTH', data=[pulse_length])

def write_tiff(d, dirname=None, basename=None):
  """The write an image tiff.  Basic implementation no frills, no metadata
  """

  if basename is None:
    basename = ""
  if dirname is None:
    dirname = "."
  if not os.path.isdir(dirname):
    os.makedirs(dirname)

  # The output path should not contain any funny characters which may
  # not work in all environments.  This constructs a sequence number à
  # la evt_seqno() from the dictionary's timestamp.
  t = d['TIMESTAMP']
  s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]

  path = os.path.join(dirname, basename + s + '.tiff')

  #assure that the 2-byte data are within the unsigned limits
  selecthi = d["DATA"]>65535
  d["DATA"].set_selected(selecthi,0)
  selectlo = d["DATA"]<0
  d["DATA"].set_selected(selectlo,0)

  idata = d["DATA"].as_numpy_array()
  idata =  idata.astype("uint16")
  import cv2 # psdm install should have this extension
  cv2.imwrite(path,idata)
  return path

def dwritef(d, dirname=None, basename=None):
  """The dwritef() function pickles the dictionary pointed to by @p d
  to the file whose directory and filename portions are pointed to by
  @p dirname and @p basename, respectively.  The directory at @p
  dirname, as well as any intermediate directories, are recursively
  created if they do not already exist.  The name of the written file
  is the concatenation of the @p basename parameter and a sequence
  number derived from the timestamp in the dictionary, @p d.

  @param d        Dictionary, as created by e.g. dpack()
  @param dirname  Directory portion of output file
  @param basename Filename prefix of output file
  @return         Path of output file
  """

  if basename is None:
    basename = ""
  if dirname is None:
    dirname = "."
  if not os.path.isdir(dirname):
    os.makedirs(dirname)

  # The output path should not contain any funny characters which may
  # not work in all environments.  This constructs a sequence number à
  # la evt_seqno() from the dictionary's timestamp.
  t = d['TIMESTAMP']
  s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]

  # XXX Several non-pyana tools rely on the .pickle extension.  Fix
  # those before migrating to .pkl.
  path = os.path.join(dirname, basename + s + '.pickle')
  easy_pickle.dump(path, d)
  return path


def dwritef2(obj, path):
  """The dwritef2() function writes the object @p obj to the Python
  pickle file whose path is pointed to by @p path.  Non-existent
  directories of @p path are created as necessary.

  @param obj  Object to write, as created by e.g. dpack()
  @param path Path of output file
  @return     Path of output file
  """

  dirname = os.path.dirname(path)
  if dirname is not "" and not os.path.isdir(dirname):
    os.makedirs(dirname)

  easy_pickle.dump(path, obj)
  return path


def pathsubst(format_string, evt, env, **kwargs):
  """The pathsubst() function provides variable substitution and value
  formatting as described in PEP 3101.  The function returns a copy of
  the input string, @p format_string, with field names replaced by
  their appropriate values as determined by either @p evt, @p env, or
  the user-supplied keyworded arguments, @p kwargs.

  chunk:      Chunk number or -1 if unknown.

  epoch:      Time of the event, in number of seconds since midnight,
              1 January 1970 UTC (Unix time), to millisecond
              precision.

  experiment: Experiment name, or empty string if unknown.

  expNum:     Experiment number or -1 if unknown.

  instrument: Instrument name, or empty string if unknown.

  iso8601:    The time of the event as an extended human-readable ISO
              8601 timestamp, to millisecond precision, or the empty
              string if unknown.  Not suitable for file names, because
              it contains characters that do not play well with
              certain file systems (e.g. NTFS).

  jobName:    Job name.

  jobNameSub: Combination of job name and subprocess index as a string
              which is unique for all subprocesses in a job.

  run:        Run number or -1 if unknown.

  seqno:      Sequence number or -1 if unknown.

  stream:     Stream number or -1 if unknown.

  subprocess: Subprocess number.  This is a non-negative integer in
              the range [0, nproc) when multiprocessing, or -1 for a
              single-process job.

  user:       The "login name" of the user.

  In addition to the standard conversion flags, the pathsubst()
  function implements the <code>!u</code> and <code>!l</code> flags
  for conversion to upper- and lower-case strings, respectively.

  Literal braces can be escaped by doubling, i.e. <code>{</code> is
  written <code>{{</code>, and <code>}</code> as <code>}}</code>.

  @note Chunk number, expNum, run number, and stream number are
        determined from the input XTC file name.  If a file does not
        adhere to the standard format, it may not be possible to
        determine these quantities.

  @note String substitution requires PSDM pyana version 0.10.3 or
        greater.

  @param format_string String containing replacement fields
  @param evt           Event data object, a configure object
  @param env           Environment object
  @param kwargs        User-supplied replacements, on the form
                       <code>field_name=value</code>
  @return               Copy of @p format_string, with replacement
                       fields substituted by their appropriate values
  """

  from getpass import getuser
  from string import Formatter

  class CaseFormatter(Formatter):
    def convert_field(self, value, conversion):
      # Extends the stock Formatter class with lower() and upper()
      # conversion types.

      if conversion == 'l':
        return str(value).lower()
      elif conversion == 'u':
        return str(value).upper()
      return super(CaseFormatter, self).convert_field(value, conversion)

    def get_value(self, key, args, kwargs_local):
      # The get_value() function sequentially applies user-supplied
      # and standard substitutions, and implements suitable defaults
      # in case a field name evaluates to None.  XXX Emit a warning
      # when this happens?
      if key in kwargs:
        return kwargs[key]

      value = super(CaseFormatter, self).get_value(key, args, kwargs_local)
      if value is None:
        if key == 'chunk':
          return -1
        elif key == 'expNum':
          return -1
        elif key == 'iso8601':
          return ''
        elif key == 'run':
          return -1
        elif key == 'seqno':
          return -1
        elif key == 'stream':
          return -1
      return value

  t = evt_time(evt)
  if t is not None:
    epoch = t[0] + t[1] / 1000
  else:
    epoch = None
  fmt = CaseFormatter()

  try:
    # psana
    expNum = env.expNum()
  except AttributeError:
    # pyana
    expNum = evt.expNum()

  try:
    # pyana
    chunk = evt.chunk()
  except AttributeError:
    # not supported in psana
    chunk = None

  try:
    # pyana
    stream = evt.stream()
  except AttributeError:
    # not supported in psana
    stream = None

  # If chunk or stream numbers cannot be determined, which may happen
  # if the XTC file has a non-standard name, evt.chunk() and
  # evt.stream() will return None.
  return fmt.format(format_string,
                    chunk=chunk,
                    epoch=epoch,
                    experiment=env.experiment(),
                    expNum=expNum,
                    instrument=env.instrument(),
                    iso8601=evt_timestamp(t),
                    jobName=env.jobName(),
                    jobNameSub=env.jobNameSub(),
                    run=evt.run(),
                    seqno=int(evt_seqno(evt)),
                    stream=stream,
                    subprocess=env.subprocess(),
                    user=getuser())

def get_ebeam(evt):
  try:
    # pyana
    ebeam = evt.getEBeam()
  except AttributeError as e:
    from psana import Source, Bld
    src = Source('BldInfo(EBeam)')
    ebeam = evt.get(Bld.BldDataEBeamV6, src)
    if ebeam is None:
      ebeam = evt.get(Bld.BldDataEBeamV5, src)
    if ebeam is None:
      ebeam = evt.get(Bld.BldDataEBeamV4, src)
    if ebeam is None:
      ebeam = evt.get(Bld.BldDataEBeamV3, src)
    if ebeam is None:
      ebeam = evt.get(Bld.BldDataEBeamV2, src)
    if ebeam is None:
      ebeam = evt.get(Bld.BldDataEBeamV1, src)
    if ebeam is None:
      ebeam = evt.get(Bld.BldDataEBeamV0, src)
    if ebeam is None:
      ebeam = evt.get(Bld.BldDataEBeam, src) # recent version of psana will return a V7 event or higher if this type is asked for

  return ebeam

def env_laser_status(env, laser_id):
  """The return value is a bool that indicates whether the laser in
  question was on for that particular shot.  Bear in mind that sample
  hit by the laser will only encounter the X-rays some time after,
  depending on the flow rate.
  """

  if env is not None:
    pv_in = env.epicsStore().value('CXI:LAS:SHT:%02i:IN' % laser_id)
    pv_out = env.epicsStore().value('CXI:LAS:SHT:%02i:OUT' % laser_id)

    if pv_in is None or pv_out is None:
      return

    if hasattr(pv_in, "values"):
      if len(pv_in.values) != 1:
        return
      laser_off = pv_in.values[0]
    else:
      laser_off = pv_in

    if hasattr(pv_out, "values"):
      if len(pv_out.values) != 1:
        return
      laser_on = pv_out.values[0]
    else:
      laser_on = pv_out

    if laser_on and laser_off:
      # According to LCLS staff, this means the laser is not plugged in
      return False

    return bool(laser_on)


def env_injector_xyz(env):
  """Returns the coordinates of the sample injector.  XXX units unknown?"""
  if env is not None:
    return tuple([
      env.epicsStore().value("CXI:USR:MZM:0%i:ENCPOSITIONGET" %(i+1))
                             for i in range(3)])


def env_detz(address, env):
  """The env_detz() function returns the position of the detector with
  the given address string on the z-axis in mm.  The zero-point is as
  far away as possible from the sample, and values decrease as the
  detector is moved towards the sample.
  @param address Full data source address of the DAQ device
  @param env     Environment object
  @return        Detector z-position, in mm
  """

  if env is not None:
    detector = address_split(address, env)[0]
    if detector is None:
      return None
    elif detector == 'CxiDs1':
      pv = env.epicsStore().value('CXI:DS1:MMS:06.RBV')
      if pv is None:
        # Even though potentially unsafe, fall back on the commanded
        # value if the corresponding read-back value cannot be read.
        # According to Sébastien Boutet, this particular motor has not
        # caused any problem in the past.
        pv = env.epicsStore().value('CXI:DS1:MMS:06')
      if pv is None:
        # Try the other detector. These are sometimes inconsistent
        pv = env.epicsStore().value('CXI:DS2:MMS:06.RBV')
    elif detector == 'CxiDsd' or detector == 'CxiDs2':
      # XXX Note inconsistency in naming: Dsd vs Ds2!
      pv = env.epicsStore().value('CXI:DS2:MMS:06.RBV')
      if pv is None:
        # Try the other detector. These are sometimes inconsistent
        pv = env.epicsStore().value('CXI:DS1:MMS:06.RBV')
    elif detector == 'XppGon':
      # There is no distance recorded for the XPP's CSPAD on the robot
      # arm.  Always return zero to allow the distance to be set using
      # the offset.
      return 0
    elif detector == 'XppEndstation' or \
         detector == 'MfxEndstation':
      # There is no distance recorded for the XPP's or MFX's Rayonix
      # on the robot arm.  Always return zero to allow the distance to
      # be set using the offset.
      return 0
    else:
      return None

    if pv is None:
      return None

    if hasattr(pv, "values"):
      if len(pv.values) == 1:
        return pv.values[0]
      else:
        return None
    return pv

  return None


def env_distance(address, env, offset):
  """The env_distance() function returns the distance between the
  sample and the detector with the given address string in mm.  The
  distance between the sample and the the detector's zero-point can
  vary by an inch or more between different LCLS runs.  According to
  Sébastien Boutet the offset should be stable to within ±0.5 mm
  during a normal experiment.

  @param address Full data source address of the DAQ device
  @param env     Environment object
  @param offset  Detector-sample offset in mm, corresponding to
                 longest detector-sample distance
  @return        Detector-sample distance, in mm
  """

  detz = env_detz(address, env)
  if detz is not None:
    return detz + offset
  return None


def env_sifoil(env):
  """The env_sifoil() function returns the total thickness of Si-foil,
  in um, that attenuates the beam.  According to an e-mail from Garth
  Williams, the centres of the attenuators are in the beam at around 0
  mm, and leave the beam at something like -7 mm.  The "out" position
  is at approximately -15 mm.

  @param env Environment object
  @return    Total thickness of attenuating Si-foil
  """

  if (env is None):
    return (None)

  # the pv name (? XXX) and the length of Si-foil it corresponds to
  # XXX static?
  dia = { "XRT:DIA:MMS:02.RBV":    20,
          "XRT:DIA:MMS:03.RBV":    40,
          "XRT:DIA:MMS:04.RBV":    80,
          "XRT:DIA:MMS:05.RBV":   160,
          "XRT:DIA:MMS:06.RBV":   320,
          "XRT:DIA:MMS:07.RBV":   640,
          "XRT:DIA:MMS:08.RBV":  1280,
          "XRT:DIA:MMS:09.RBV":  2560,
          "XRT:DIA:MMS:10.RBV":  5120,
          "XRT:DIA:MMS:11.RBV": 10240 }

  si_tot = 0
  for pvname, si_len in six.iteritems(dia):
    pv = env.epicsStore().value(pvname)

    # XXX Why is this an EpicsPvTime object?  The absorption
    # coefficient of Si is E-18 * n_{0} * lambda^2, (for lambda >= 5
    # um, Schroder, D. K., R. N. Thomos, and J. C. Swartz, IEEE
    # Trans. Electron. Dev. ED-25, 2(1978) 254-261).  See also
    # http://henke.lbl.gov/optical_constants/filter2.html

    #print "For ", pvname, " got ", pv, " and ", pv.values[0]
    if pv is not None: # and pv.units          == "mm"
      if hasattr(pv, "values"):
        # pyana
        if len(pv.values) == 1 and abs(pv.values[0]) <  7:
          si_tot += si_len
      else:
        # psana
        if abs(pv) < 7:
          si_tot += si_len

  return (si_tot)


def env_wavelength_sxr(evt, env):
  """The env_wavelength_sxr() function returns the wavelength in
  Ångström of the environment pointed to by @p env at the time of the
  event @p evt.  The function returns a positive value or @c None if
  no wavelength is available for the event.  See Heimann et al. (2011)
  Rev. Sci. Instrum. 82, 093104.

  @note The wavelength in eV is 12398.4187 divided by the value
        returned from env_wavelength_sxr().

  @param evt Event data object, a configure object
  @param env Environment object
  @return    Wavelength, in Ångström
  """

  from calendar import timegm
  from time import strptime

  if evt is None or env is None:
    return None

  t = evt.getTime()
  if t is None:
    return None

  es = env.epicsStore()
  if es is None:
    return None

  # Note that the monochromator coefficients could change from day to
  # day.  Unless specific values for the requested time are available,
  # attempt to retrieve them from EPICS.
  #
  # The compiler could recognize that strptime() and timegm() are pure
  # and reduce the test expression to an integer comparison.
  f = '%Y-%m-%d, %H:%M %Z'
  s = t.seconds()
  if s is None:
    return None
  elif s < timegm(strptime('2012-11-12, 17:00 UTC', f)):
    return None
  elif s < timegm(strptime('2012-11-17, 17:00 UTC', f)):
    abc = [+3.65920, -0.76851, +0.02105]
  elif s < timegm(strptime('2012-11-20, 17:00 UTC', f)):
    abc = [+4.18190, -0.77650, +0.01020]

  if 'abc' not in locals():
    pv = []
    for name in ['SXR:IOC:POLY:POLY:Lambda:O1:G3:A',
                 'SXR:IOC:POLY:POLY:Lambda:O1:G3:B',
                 'SXR:IOC:POLY:POLY:Lambda:O1:G3:C']:
      pv.append(es.value(name))
      if pv[-1] is None or len(pv[-1].values) != 1:
        return None
      pv[-1] = pv[-1].values[0]
      if pv[-1] is None:
        return None
    abc = [pv[i] for i in range(3)]

  # Get the grating motor position from EPICS.
  pv = es.value('SXR:MON:MMS:06.RBV')
  if pv is not None and len(pv.values) == 1:
    x = pv.values[0]
    e = 10 * (abc[0] + abc[1] * x + abc[2] * x**2)
    if e > 0:
      return e
  return None


def evt_pulse_energy(evt):
  """The evt_pulse_energy() function returns the energy, or the
  intensity, of the pulse in arbitrary units.  The returned value
  should be proportional to the number of photons in the pulse, and
  may be negative due to noise.

  @note An absolute, but less accurate, estimate of the number of
        photons in the pulse may be obtained from the gas monitor
        detector's fMilliJoulesPerPulse value.

  @param evt Event data object, a configure object
  @return    Pulse intensity, in arbitrary units
  """

  from pypdsdata.xtc import TypeId

  if evt is None:
    return None

  gmd = evt.get(key=TypeId.Type.Id_GMD)
  if hasattr(gmd, 'fRelativeEnergyPerPulse') and evt.expNum() == 208:
    # Note that for L632 (experiment number 208)
    # fRelativeEnergyPerPulse actually gives the negated value
    # sought.  Details are given in Moeller, S. (2012) "GMD Look
    # up Sheet for variable names in the DAQ (BLD) versus the C++
    # code".
    return -gmd.fRelativeEnergyPerPulse

  elif hasattr(gmd, 'fCorrectedSumPerPulse'):
    # This relatively pressure-independent quantity in arbitrary
    # units is preferable.  It is also known as
    # SXR:GMD:BLD:CumSumAllPeaks.
    return gmd.fCorrectedSumPerPulse
  return None


def evt_pulse_length(evt):
  """The evt_pulse_length() function returns the pulse length in fs.
  It is calculated as the ratio of the charge (in nC) and the peak
  current (in A).

  @param evt Event data object, a configure object
  @return    Pulse length, in fs
  """

  if (evt is not None):
    ebeam = get_ebeam(evt)

    if ebeam is None:
      return

    try:
      if ebeam.fEbeamPkCurrBC2 > 0:
        return 1e6 * ebeam.fEbeamCharge / ebeam.fEbeamPkCurrBC2
    except AttributeError:
      if ebeam.ebeamPkCurrBC2() > 0:
        return 1e6 * ebeam.ebeamCharge() / ebeam.ebeamPkCurrBC2()
  return None


def evt_repetition_rate(evt, address='*'):
  """The evt_repetition_rate() function returns the repetition rate of
  the instrument in Hz.  See
  https://confluence.slac.stanford.edu/display/PCDS/EVR+Event+Codes

  @param evt     Event data object, a configure object
  @param address Data source address of the DAQ device
  @return        Integer repetition rate, in Hz
  """

  evr = evt.getEvrData(address)
  if evr is not None:
    event_code_map = [120, 60, 30, 10, 5, 1]
    for i in range(evr.numFifoEvents() - 1, -1, -1):
      # Search for the last repetition rate event code.
      j = evr.fifoEvent(i).EventCode
      if j >= 40 and j <= 45:
        # These are the NO BEAM event codes.
        return event_code_map[j - 40]
      if j >= 140 and j <= 145:
        # These are the undocumented BEAM event codes.
        return event_code_map[j - 140]
  return None


def evt_beam_charge(evt):
  """The evt_beam_charge() function returns the charge of the pulse in
  nC.

  @param evt Event data object, a configure object
  @return    Pulse charge, in nC
  """

  if evt is not None:
    ebeam = get_ebeam(evt)

    if ebeam is None:
      return
    try:
      ebeam = evt.getEBeam()
      return ebeam.fEbeamCharge
    except AttributeError:
      return ebeam.ebeamCharge()
  return None


def evt_seqno(evt=None):
  """The evt_seqno() function returns string representation of a
  sequence number.  If @p evt is not @c None the return value reflects
  the time at which @p evt occurred, otherwise the current time is
  used.  If @p evt does not contain a time, evt_seqno() returns @c
  None.  XXX Should probably return an integer type instead?

  @param evt Event data object, a configure object
  @return    String representation of sequence number
  """

  t = evt_time(evt=evt)
  if t is None:
    return None
  return time.strftime("%Y%m%d%H%M%S", time.gmtime(t[0])) + ("%03d" % t[1])


def evt_time(evt=None):
  """The evt_time() function returns the time of the event @p evt since
  midnight, 1 January 1970 UTC (Unix time) to millisecond precision.
  If @p evt does not contain a time, evt_time() returns @c None.  If
  @p evt is @c None the return value reflects current time is used.

  @note Millisecond precision is sufficient, because at 120 Hz, shots
        are taken at 8.3 ms intervals.

  @param evt Event data object, a configure object
  @return    Unix time as a tuple of seconds and milliseconds
  """

  if evt is None:
    t = time.time()
    s = int(math.floor(t))
    return (s, int(round((t - s) * 1000)))

  if hasattr(evt, "getTime"):
    t = evt.getTime()
    if t is None:
      return None
    return (t.seconds(), t.nanoseconds() // 1000000)
  else:
    from psana import EventId
    id = evt.get(EventId)
    return (id.time()[0], id.time()[1] // 1000000)


def evt_timestamp(t=None):
  """The evt_timestamp() function returns a string representation of
  an extended human-readable ISO 8601 timestamp.  If @p t is @c None
  the current time is used.  The function returns @c None on failure.

  @param t Tuple of the time in seconds and milliseconds
  @return  Human-readable ISO 8601 timestamp in string representation
  """

  if t is None:
    t = evt_time(evt=None)
    if t is None:
      return None
  return time.strftime("%Y-%m-%dT%H:%MZ%S", time.gmtime(t[0])) + \
      (".%03d" % t[1])


def evt_wavelength(evt, delta_k=0):
  """The evt_wavelength() function returns the wavelength in Ångström
  of the event pointed to by @p evt.  From Margaritondo & Rebernik
  Ribic (2011): the dimensionless relativistic γ-factor is derived
  from beam energy in MeV and the electron rest mass, K is a
  dimensionless "undulator parameter", and L is the macroscopic
  undulator period in Ångström.  See also
  https://people.eecs.berkeley.edu/~attwood/srms/2007/Lec10.pdf

  @param evt     Event data object, a configure object
  @param delta_k Optional K-value correction
  @return        Wavelength, in Ångström
  """

  if evt is not None:
    ebeam = get_ebeam(evt)

    if hasattr(ebeam, 'fEbeamPhotonEnergy') and ebeam.fEbeamPhotonEnergy > 0:
      # pyana
      return 12398.4187 / ebeam.fEbeamPhotonEnergy
    if hasattr(ebeam, 'ebeamPhotonEnergy') and ebeam.ebeamPhotonEnergy() > 0:
      # psana
      return 12398.4187 / ebeam.ebeamPhotonEnergy()

    if hasattr(ebeam, 'fEbeamL3Energy') and ebeam.fEbeamL3Energy > 0:
      # pyana
      gamma = ebeam.fEbeamL3Energy / 0.510998910
    elif hasattr(ebeam, 'ebeamL3Energy') and ebeam.ebeamL3Energy() > 0:
      # psana
      gamma = ebeam.ebeamL3Energy() / 0.510998910
    else:
      return None
    K = 3.5 + delta_k
    L = 3.0e8
    return L / (2 * gamma**2) * (1 + K**2 / 2)
  return None

def old_address_to_new_address(address):
  """ Change between old and new style detector addresses.
  I.E. CxiDs1-0|Cspad-0 becomes CxiDs1.0:Cspad.0
  @param address detector address to convert
  """
  return address.replace('-','.').replace('|',':')

def getConfig(address, env):
  """ Given a detector address, find the config object in an env object
  that goes with it.
  @param address detector address
  @param env environment object to search"""

  if hasattr(env, 'configStore'):
    good_key = None
    address = old_address_to_new_address(address)
    for key in env.configStore().keys():
      if address in str(key.src()) and key.type() is not None:
        good_key = key
        break
    if good_key is None:
      return None
    return env.configStore().get(good_key.type(),good_key.src())
  else:
    # Try the pyana method for older data
    from pypdsdata.xtc import TypeId
    return env.getConfig(TypeId.Type.Id_CspadConfig, address)

def getOptBool(s):
  if s is None or s == "None": return False
  elif isinstance(s, bool):
    return s
  s = s.strip().lower()
  return s == "true"

def getOptEvalOrString(s) :
  """Allow python code macros in the pyana configuration file, e.g.
  dark_path   = "/location_of_darks/r%%04d/Ds1-avg.pickle"%%(max([{True:dark,False:0}[3 > dark] for dark in [1,2,6,9,12,14,17,19]]))
  """
  possible_string = getOptString(s)
  try:
    eval_string = eval(possible_string,{},{})
    return eval_string
  except (SyntaxError, TypeError):
    return possible_string

def getOptString(s) :
  """XXX Return the string, strip of any white space (make sure there
  are no newline characters here).  This function was originally
  written by Ingrid Ofte for pyana's XtcExplorer module.
  """

  if (s is None):
    return (None)

  s = s.strip()
  if (s == "" or s == "No" or s == "None"):
    return (None)
  return (s)


def getOptStrings(s, default=None) :
  """XXX Return a list of strings.  This function was originally
  written by Ingrid Ofte for pyana's XtcExplorer module.
  """
  if (s is None):
    return default

  # strip off any leading or trailing whitespace
  s = s.strip()

  # make sure there are no newline characters here
  s = s.split("\n")
  s = " ".join(s)

  # make a list
  l = s.split()

  if (len(l) == 0 or (len(l) == 1 and (s == "" or s == "No" or s == "None"))):
    return ([])

  # all other cases:
  return (l)


def getOptInteger(s):
  """XXX Return a single integer.  This function was originally
  written by Ingrid Ofte for pyana's XtcExplorer module.  XXX What if
  conversion fails?
  """

  if (s is None or s == "" or s == "None"):
    return None
  return (int(s))

def getOptFloat(s):
  """Return a single float.
  """

  if (s is None or s == "" or s == "None"):
    return None
  return (float(s))

def getOptROI(s):
  """Return a tuple of the region of interest.
     Format: roi = fast_low:fast_high,slow_low:slow_high
  """
  roi_str    = getOptString(s)
  if (roi_str is not None and roi_str != ""):
    ivl        = roi_str.strip("()").split(",")
    ivl_x      = ivl[0].split(":")
    ivl_y      = ivl[1].split(":")
    roi = [ivl_x[0], ivl_x[1], ivl_y[0], ivl_y[1]]
    for i in range(4):
      if roi[i] == "": roi[i] = None
      else: roi[i] = int(roi[i])
    return tuple(roi)


def image(address, config, evt, env, sections=None):
  """Assemble the uint16 detector image, and sum it up as int32.  Sum
  the image of squared intensities as uint64.  XXX Documentation! XXX
  Would be nice to get rid of the constant string names.  XXX Better
  named evt_image()?

  @param address  Full data source address of the DAQ device
  @param config   XXX This should go--get up-to-date object on the fly!
  @param evt      Event data object, a configure object
  @param env      Environment object
  @param sections XXX
  @return         XXX
  """

  device = address_split(address)[2]
  if device is None:
    return None

  elif device == 'Andor':
    # XXX There is no proper getter for Andor frames yet, and
    # evt.getFrameValue(address) does not appear to work.
    from pypdsdata.xtc import TypeId
    value = evt.get(TypeId.Type.Id_AndorFrame, address)
    if value is not None:
      img = value.data()
      return img

  elif device == 'Cspad':
    if sections is not None:
      return CsPadDetector(address, evt, env, sections)
    else:
      # XXX This is obsolete code, provided for backwards
      # compatibility with the days before detector metrology was
      # used.
      assert False # sections always required now as of Sep 1 2014
      quads = evt.getCsPadQuads(address, env)
      qimages = numpy.empty((4, npix_quad, npix_quad), dtype='uint16')
      for q in quads:
        qimages[q.quad()] = CsPadElement(q.data(), q.quad(), config)
      return numpy.vstack((numpy.hstack((qimages[0], qimages[1])),
                           numpy.hstack((qimages[3], qimages[2]))))

  elif device == 'Cspad2x2':
    from pypdsdata.xtc import TypeId
    quads = evt.get(TypeId.Type.Id_Cspad2x2Element, address)
    if quads is not None:
      return CsPad2x2Image(quads.data(), config, sections)

  elif device == 'pnCCD':
    value = evt.getPnCcdValue(address, env)
    if value is not None:
      # Returns the image data as a numpy 1024-by-1024 uint16 array
      # XXX Should be split up into tiles (halves) to allow metrology
      # to be adjusted?  Will require a sections parameter!
      img = value.data()

      # Deal with overflows.  XXX This might be dependent on the
      # particular version of pyana.  CASS ignores the two most
      # significant bits, which is different from what is done below,
      # but Lutz Foucar says they do contain data which could be used.
      img[img > 2**14 - 1] = 2**14 - 1
      return img
  return None

def image_xpp(address, evt, env, aa, quads = None):
  """Assemble the uint16 detector image, see also
  cspad_tbx.CsPadDetector().  XXX Documentation! XXX Would be nice to
  get rid of the constant string names.  XXX Better named evt_image()?

  @param address Full data source address of the DAQ device
  @param evt     Event data object, a configure object
  @param env     Environment object
  @param aa      Active areas, in lieu of full metrology object
  @param quads   Data, if None get it from the event
  @return        XXX
  """

  if address != 'XppGon-0|Cspad-0':
    return None

  # Get a current configure object for the detector
  config = getConfig(address, env)
  if config is None:
    return None

  if quads is None:
    # For consistency, one could/should verify that len(quads) is equal
    # to len(sections).
    quads = evt_get_quads(address, evt, env)
    if quads is None or len(quads) != len(aa) // (8 * 2 * 4):
      return None

  # Start out with a blank image of the detector.  Mikhail
  # S. Dubrovin's HDF5Explorer/src/ConfigCSpad.py uses a detector
  # size of 1765-by-1765 pixels.  This assumes that the type of the
  # first section in the first quadrant is identical to the type of
  # all the other sections.
  det = numpy.zeros((1765, 1765), dtype=quads[0].data()[0].dtype)

  for quad in quads:
    q_data = quad.data()
    q_idx = quad.quad()
    try:
      # pyana
      # example: if the third sensor (2x1) is disabled, q_mask = [0,1,3,4,5,6,7]
      q_mask = config.sections(q_idx)
    except AttributeError:
      # psana
      # as above, using config.roiMask, a bitstring where the ith bit is true if the ith sensor is active. x << y means bitwise shift
      # x, y times, and & is the bitwise AND operator
      q_mask = [i for i in range(config.numSect()//config.numQuads()) if 1 << i & config.roiMask(q_idx)]

    # For consistency, one could/should verify that len(q_data) is
    # equal to len(sections[q_idx]).
    assert len(q_data) == len(q_mask)
    for (s_data, s_idx) in zip(q_data, q_mask):
      # Rotate the "lying down" sensor readout from the XTC stream by
      # an integer multiple of 90 degrees to match the orientation on
      # the detector.  This assumes that the horizontal dimension of
      # the unrotated sensor is even.  Note that the XPP CSPAD is
      # rotated by 180 degrees with respect to the optical metrology
      # measurements.
      assert s_data.shape[1] % 2 == 0
      if   q_idx == 0 and s_idx in [2, 3, 6, 7] or \
           q_idx == 1 and s_idx in [0, 1]       or \
           q_idx == 3 and s_idx in [4, 5]:
        asics = numpy.hsplit(numpy.rot90(s_data, 0 + 2), 2)
        asics.reverse()
      elif q_idx == 0 and s_idx in [0, 1]       or \
           q_idx == 2 and s_idx in [4, 5]       or \
           q_idx == 3 and s_idx in [2, 3, 6, 7]:
        asics = numpy.vsplit(numpy.rot90(s_data, 1 + 2), 2)
      elif q_idx == 1 and s_idx in [4, 5]       or \
           q_idx == 2 and s_idx in [2, 3, 6, 7] or \
           q_idx == 3 and s_idx in [0, 1]:
        asics = numpy.hsplit(numpy.rot90(s_data, 2 + 2), 2)
      elif q_idx == 0 and s_idx in [4, 5]       or \
           q_idx == 1 and s_idx in [2, 3, 6, 7] or \
           q_idx == 2 and s_idx in [0, 1]:
        asics = numpy.vsplit(numpy.rot90(s_data, 3 + 2), 2)
        asics.reverse()
      else:
        # NOTREACHED
        return None

      # Use the active areas to place the two ASICS on the
      # destination detector image.
      for a_idx in range(len(asics)):
        aa_idx = q_idx * (8 * 2 * 4) + s_idx * (2 * 4) + a_idx * 4
        det[aa[aa_idx + 0]:aa[aa_idx + 2],
            aa[aa_idx + 1]:aa[aa_idx + 3]] = asics[a_idx]

  return det


def iplace(dst, src, angle, center):
    """The iplace() function places @p src in @p dst centred on @p
    center after rotating it by @p angle degrees counter-clockwise.
    The source image is mapped onto the destination image by bilinear
    interpolation.  While this may introduce interpolation artifacts
    it is significantly simpler than many other interpolation
    methods--and bog slow.

    @p dst    Destination image
    @p src    Source image
    @p angle  Rotation angle, in degrees
    @p center Centre of @p src in @p dst, after rotation
    """

    a = math.radians(angle)
    c = math.cos(a)
    s = math.sin(a)

    # Find the origin-centred bounding box of the rotated source
    # image.  Due to the symmetry of a rectangle, the extrema can be
    # determined by the transformed coordinates of two adjacent
    # corners.
    hsize = [0.5 * max(abs(c * src.shape[0] - s * src.shape[1]),
                       abs(c * src.shape[0] + s * src.shape[1])),
             0.5 * max(abs(s * src.shape[0] + c * src.shape[1]),
                       abs(s * src.shape[0] - c * src.shape[1]))]
    xlim  = [int(math.floor(-hsize[0])),
             int(math.ceil( +hsize[0])) + 1]
    ylim  = [int(math.floor(-hsize[1])),
             int(math.ceil( +hsize[1])) + 1]

    # For each pixel in the bounding box, determine the real-valued
    # components in coordinate system of the untransformed source
    # image, (xp, yp).  Then do bilinear interpolation based on the
    # four pixels with integer coordinates around (xp, yp).
    for x in range(xlim[0], xlim[1]):
        for y in range(ylim[0], ylim[1]):
            xp =  c * x + s * y + 0.5 * src.shape[0]
            yp = -s * x + c * y + 0.5 * src.shape[1]
            if (xp >= 0 and math.ceil(xp) < src.shape[0] and
                yp >= 0 and math.ceil(yp) < src.shape[1]):

                xi =[int(math.floor(xp)), int(math.ceil(xp))]
                yi =[int(math.floor(yp)), int(math.ceil(yp))]

                xf = xp - xi[0]
                yf = yp - yi[0]

                dst[int(round(x + center[0])),
                    int(round(y + center[1]))] =              \
                    src[xi[0], yi[0]] * (1 - xf) * (1 - yf) + \
                    src[xi[1], yi[0]] * xf       * (1 - yf) + \
                    src[xi[0], yi[1]] * (1 - xf) * yf       + \
                    src[xi[1], yi[1]] * xf       * yf


def rplace(dst, src, angle, center):
    """The rplace() function places @p src in @p dst centred on @p
    centre after rotating it by @p angle degrees counter-clockwise.
    The rotation angle is rounded to the nearest integer multiple of
    90 degrees before transformation.

    @p dst    Destination image
    @p src    Source image
    @p angle  Rotation angle, in degrees
    @p center Centre of @p src in @p dst, after rotation
    """

    # Rotate the source image, and determine the upper, left corner of
    # its location in the destination image.
    rot = numpy.rot90(src, int(round(angle / 90.0)) % 4)
    ulc = [int(round(center[0] - 0.5 * rot.shape[0])),
           int(round(center[1] - 0.5 * rot.shape[1]))]

    dst[ulc[0]:(ulc[0] + rot.shape[0]),
        ulc[1]:(ulc[1] + rot.shape[1])] = rot

# For the moment, the XPP CSPAD detector's metrology is stored here
# as a series of active areas
_xpp_active_areas = {
  'XPP 7.1': { # metrology recorded 1/24/13 and processed by flatfile.py
    'active_areas': flex.int([
           865, 1121, 1059, 1306, 1062, 1121, 1256, 1306,
           864,  909, 1058, 1094, 1061,  909, 1255, 1094,
          1083, 1534, 1268, 1728, 1083, 1337, 1268, 1531,
           871, 1538, 1056, 1732,  871, 1341, 1056, 1535,
          1495, 1326, 1689, 1511, 1298, 1326, 1492, 1511,
          1496, 1539, 1690, 1724, 1299, 1539, 1493, 1724,
          1482, 1105, 1667, 1299, 1482,  908, 1667, 1102,
          1270, 1107, 1455, 1301, 1270,  910, 1455, 1104,
          1123,  706, 1308,  900, 1123,  509, 1308,  703,
           910,  706, 1095,  900,  910,  509, 1095,  703,
          1535,  498, 1729,  683, 1338,  498, 1532,  683,
          1534,  711, 1728,  896, 1337,  711, 1531,  896,
          1324,   77, 1509,  271, 1324,  274, 1509,  468,
          1537,   75, 1722,  269, 1537,  272, 1722,  466,
          1104,   97, 1298,  282,  907,   97, 1101,  282,
          1105,  310, 1299,  495,  908,  310, 1102,  495,
           706,  457,  900,  642,  509,  457,  703,  642,
           705,  669,  899,  854,  508,  669,  702,  854,
           496,   36,  681,  230,  496,  233,  681,  427,
           709,   38,  894,  232,  709,  235,  894,  429,
            77,  256,  271,  441,  274,  256,  468,  441,
            77,   44,  271,  229,  274,   44,  468,  229,
            98,  467,  283,  661,   98,  664,  283,  858,
           311,  467,  496,  661,  311,  664,  496,  858,
           457,  867,  642, 1061,  457, 1064,  642, 1258,
           670,  865,  855, 1059,  670, 1062,  855, 1256,
            37, 1084,  231, 1269,  234, 1084,  428, 1269,
            37,  871,  231, 1056,  234,  871,  428, 1056,
           256, 1495,  441, 1689,  256, 1298,  441, 1492,
            43, 1497,  228, 1691,   43, 1300,  228, 1494,
           469, 1481,  663, 1666,  666, 1481,  860, 1666,
           467, 1269,  661, 1454,  664, 1269,  858, 1454]),
    'rotations' : flex.int([
                   3,3,3,3,2,2,2,2,1,1,1,1,2,2,2,2,
                   2,2,2,2,1,1,1,1,0,0,0,0,1,1,1,1,
                   1,1,1,1,0,0,0,0,3,3,3,3,0,0,0,0,
                   0,0,0,0,3,3,3,3,2,2,2,2,3,3,3,3
                  ])
  },
  'XPP 7.2': { # metrology recorded 1/29/13 and processed by flatfile.py
    'active_areas': flex.int([
           868, 1122, 1062, 1307, 1065, 1122, 1259, 1307,
           868,  910, 1062, 1095, 1065,  910, 1259, 1095,
          1087, 1534, 1272, 1728, 1087, 1337, 1272, 1531,
           874, 1536, 1059, 1730,  874, 1339, 1059, 1533,
          1497, 1328, 1691, 1513, 1300, 1328, 1494, 1513,
          1499, 1541, 1693, 1726, 1302, 1541, 1496, 1726,
          1483, 1105, 1668, 1299, 1483,  908, 1668, 1102,
          1271, 1106, 1456, 1300, 1271,  909, 1456, 1103,
          1122,  705, 1307,  899, 1122,  508, 1307,  702,
           909,  705, 1094,  899,  909,  508, 1094,  702,
          1534,  497, 1728,  682, 1337,  497, 1531,  682,
          1533,  710, 1727,  895, 1336,  710, 1530,  895,
          1323,   76, 1508,  270, 1323,  273, 1508,  467,
          1536,   75, 1721,  269, 1536,  272, 1721,  466,
          1103,   97, 1297,  282,  906,   97, 1100,  282,
          1103,  310, 1297,  495,  906,  310, 1100,  495,
           705,  456,  899,  641,  508,  456,  702,  641,
           704,  669,  898,  854,  507,  669,  701,  854,
           495,   35,  680,  229,  495,  232,  680,  426,
           707,   38,  892,  232,  707,  235,  892,  429,
            75,  256,  269,  441,  272,  256,  466,  441,
            75,   43,  269,  228,  272,   43,  466,  228,
            97,  467,  282,  661,   97,  664,  282,  858,
           310,  466,  495,  660,  310,  663,  495,  857,
           456,  866,  641, 1060,  456, 1063,  641, 1257,
           669,  865,  854, 1059,  669, 1062,  854, 1256,
            36, 1084,  230, 1269,  233, 1084,  427, 1269,
            35,  870,  229, 1055,  232,  870,  426, 1055,
           254, 1494,  439, 1688,  254, 1297,  439, 1491,
            42, 1496,  227, 1690,   42, 1299,  227, 1493,
           468, 1481,  662, 1666,  665, 1481,  859, 1666,
           465, 1268,  659, 1453,  662, 1268,  856, 1453]),
    'rotations' : flex.int([
                   3,3,3,3,2,2,2,2,1,1,1,1,2,2,2,2,
                   2,2,2,2,1,1,1,1,0,0,0,0,1,1,1,1,
                   1,1,1,1,0,0,0,0,3,3,3,3,0,0,0,0,
                   0,0,0,0,3,3,3,3,2,2,2,2,3,3,3,3
                  ])
  },
  'XPP 8.1': { # metrology recorded 10/09/13 and processed by flatfile.py
    'active_areas': flex.int([
           863, 1118, 1057, 1303, 1060, 1118, 1254, 1303,
           865,  913, 1059, 1098, 1062,  913, 1256, 1098,
          1070, 1532, 1255, 1726, 1070, 1335, 1255, 1529,
           863, 1532, 1048, 1726,  863, 1335, 1048, 1529,
          1484, 1335, 1678, 1520, 1287, 1335, 1481, 1520,
          1484, 1543, 1678, 1728, 1287, 1543, 1481, 1728,
          1475, 1110, 1660, 1304, 1475,  913, 1660, 1107,
          1268, 1109, 1453, 1303, 1268,  912, 1453, 1106,
          1119,  707, 1304,  901, 1119,  510, 1304,  704,
           912,  707, 1097,  901,  912,  510, 1097,  704,
          1533,  506, 1727,  691, 1336,  506, 1530,  691,
          1533,  715, 1727,  900, 1336,  715, 1530,  900,
          1334,   84, 1519,  278, 1334,  281, 1519,  475,
          1541,   85, 1726,  279, 1541,  282, 1726,  476,
          1108,  103, 1302,  288,  911,  103, 1105,  288,
          1108,  311, 1302,  496,  911,  311, 1105,  496,
           706,  460,  900,  645,  509,  460,  703,  645,
           706,  666,  900,  851,  509,  666,  703,  851,
           507,   38,  692,  232,  507,  235,  692,  429,
           713,   38,  898,  232,  713,  235,  898,  429,
            82,  241,  276,  426,  279,  241,  473,  426,
            82,   37,  276,  222,  279,   37,  473,  222,
           103,  459,  288,  653,  103,  656,  288,  850,
           310,  460,  495,  654,  310,  657,  495,  851,
           460,  862,  645, 1056,  460, 1059,  645, 1253,
           666,  863,  851, 1057,  666, 1060,  851, 1254,
            38, 1070,  232, 1255,  235, 1070,  429, 1255,
            38,  864,  232, 1049,  235,  864,  429, 1049,
           242, 1484,  427, 1678,  242, 1287,  427, 1481,
            37, 1484,  222, 1678,   37, 1287,  222, 1481,
           458, 1475,  652, 1660,  655, 1475,  849, 1660,
           459, 1267,  653, 1452,  656, 1267,  850, 1452]),
    'rotations' : flex.int([
                   3,3,3,3,2,2,2,2,1,1,1,1,2,2,2,2,
                   2,2,2,2,1,1,1,1,0,0,0,0,1,1,1,1,
                   1,1,1,1,0,0,0,0,3,3,3,3,0,0,0,0,
                   0,0,0,0,3,3,3,3,2,2,2,2,3,3,3,3
                  ])
  },

  #SOME BIG ISSUES REMAIN WITH Sacla.MPCCD.8tile format
  # Evidently the data from Takanori 22 Sep 2015 already has slight rotation
  # applied to the MPCCD modules, as the data rectangles displayed in cctbx.image_viewer are tilted
  # This is inconsistent with the expectation that npy.py should get the raw data, not preprocessed.

  'Sacla.MPCCD.8tile': { # as given by Takanori 22 Sep 2015
    'active_areas': flex.int([
           112,  189,  622, 1212,  647,  188, 1156, 1212,
          1180,  140, 1691, 1163, 1714,  140, 2226, 1163,
           159, 1231,  671, 2254,  694, 1230, 1206, 2253,
          1229, 1180, 1740, 2203, 1762, 1180, 2274, 2202,
     ]),
    'rotations' : flex.int([
                   0,0,0,0,0,0,0,0,
                  ])
  },

}
_xpp_active_areas['XPP 11.1'] = _xpp_active_areas['XPP 9.1'] = _xpp_active_areas['XPP 8.1']
xpp_active_areas = _xpp_active_areas
