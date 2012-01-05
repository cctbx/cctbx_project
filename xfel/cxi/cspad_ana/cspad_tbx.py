# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$

"""Toolbox for images from the Cornell SLAC Pixel Array Detector
(CSPad).

XXX Better named cspad_common?
"""

__version__ = "$Revision$"

import math
import numpy
import os
import time

from pypdsdata import xtc

from libtbx import easy_pickle
from scitbx.array_family import flex
from xfel.cxi.cspad_ana.parse_calib import Section


# The CSPad counters are 14 bits wide (Philipp et al., 2007).  XXX
# Capitalise these constants.
dynamic_range = 2**14 - 1

# The side length of a square quadrant from the old XtcExplorer code.
# XXX This should be obsoleted!
npix_quad = 850

# The pixel size in mm.  The pixel size is fixed and square at 110 um
# by 110 um (Philipp et al., 2007).  XXX UTF-8 comments?
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


def cbcaa(config, sections):
  """The cbcaa() function uses on-disk calibration data to estimate
  the beam centre and the active detector areas.  The beam centre is
  approximated as the average of the four ASIC corners closest to the
  detector centre.  That is the first corner of the section 1 in every
  quadrant.  Note that first corner index is vertical coordinate,
  second index is the horizontal coordinate.  XXX Construct the active
  areas in "spotfinder format", i.e. opposing corners.  XXX This is a
  really bad function name!

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
    for q in xrange(4): # loop over quadrants
      for i in xrange(8): # loop over two-by-one:s XXX variable name!

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

  bc     = [0, 0]
  q_mask = config.quadMask()
  for q in xrange(len(sections)):
    if (not((1 << q) & q_mask)):
      continue

    corner = sections[q][1].corners(True)[0]
    bc     = [bc[0] + corner[1] / len(sections),
              bc[1] + corner[0] / len(sections)]

    s_mask = config.roiMask(q)
    for s in xrange(len(sections[q])):
      if (not((1 << s) & s_mask)):
        continue
      c = sections[q][s].corners_asic()
      aa.extend(flex.int(c[0]))
      aa.extend(flex.int(c[1]))

  return (bc, aa)


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
  mask = map(config.sections, xrange(2))

  # For this detector, the quadrant index is always zero.
  q_idx = 0
  for s in xrange(2):
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


def CsPadDetector(data3d, config, sections, right = True):
  """The CsPadDetector() function assembles a two-dimensional image
  from the Ds1 detector readout in @p data3d and the calibration
  information in @p sections.  XXX General question: do
  variable/function names make sense?

  @param data3d   Detector readout from XTC stream
  @param config   XXX
  @param sections XXX Directory with calibration information
  @param right    @c True to restrict rotations to right angles
  @return         Assembled detector image
  """

  # For consistency, one could/should verify that len(data3d) is equal
  # to len(sections).
  assert (len(data3d) == len(sections))

  # This is from Mikhail S. Dubrovin's
  # HDF5Explorer/src/ConfigCSpad.py, which uses a detector size of
  # 1765-by-1765 pixels.
  extra_space = (1765 - 2 * Section.q_size[0],
                 1765 - 2 * Section.q_size[1])

  # Start out with a blank image of the detector.  This assumes that
  # the type of the first section in the first quadrant is identical
  # to the type of all the other sections.
  det  = numpy.zeros((2 * Section.q_size[0] + extra_space[0],
                      2 * Section.q_size[1] + extra_space[1]),
                     dtype = data3d[0].data()[0].dtype)
  mask = map(config.sections, xrange(4))

  for q in xrange(len(data3d)):
    q_data = data3d[q].data()
    q_idx  = data3d[q].quad()

    # For consistency, one could/should verify that len(q_data) is
    # equal to len(sections[q_idx]).
    assert (len(q_data) == len(sections[q_idx]))
    for s in xrange(len(q_data)):
      if (s not in mask[q_idx]):
        continue

      # Rotate the section from the XTC-stream by -90 degrees to match
      # the "standing up" convention used by the calibration data, and
      # insert a 3-pixel gap between the ASIC:s.  This assumes that
      # the horizontal dimension of the unrotated section is even.
      asics  = numpy.vsplit(numpy.rot90(q_data[s], -1), 2)
      gap    = numpy.zeros((3, 185), dtype = q_data[s].dtype)
      s_data = numpy.vstack((asics[0], gap, asics[1]))

      # Place the section in the detector image, either by forcing
      # rotation by right angles or by interpolating.
      angle  = sections[q_idx][s].angle
      center = sections[q_idx][s].center
      if (right):
        rplace(det, s_data, angle, center)
      else:
        iplace(det, s_data, angle, center)
  return (det)


def CsPadElement(data3d, qn, config):
  """Construct one image for each quadrant, each with 8 sections from
  a data3d = 3 x 2*194 x 185 data array.  This function was originally
  written by Ingrid Ofte for pyana's XtcExplorer module.  XXX
  Documentation!
  """

  # If any sections are missing, insert zeros.
  mask = map(config.sections, xrange(4))
  if (len(data3d) < 8):
    zsec = numpy.zeros((185, 388), dtype = data3d.dtype)
    for i in xrange(8) :
      if (i not in mask[qn]):
        data3d = numpy.insert(data3d, i, zsec, axis = 0)

  pairs = []
  for i in xrange(8) :
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
  for sec in xrange(8):
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


def dpack(active_areas         = None,
          beam_center_x        = None,
          beam_center_y        = None,
          ccd_image_saturation = None,
          data                 = None,
          distance             = None,
          pixel_size           = pixel_size,
          saturated_value      = None,
          size1                = None,
          size2                = None,
          timestamp            = None,
          sequence_number      = 0,
          wavelength           = None,
          xtal_target          = None):
  """XXX Check completeness.  Should fill in sensible defaults.  Must
  have data."""

  if (data is None):
    return (None)

  # Create a time stamp of the current time if none was supplied.
  if (timestamp is None):
    timestamp = evt_timestamp()

  # For unknown historical reasons, the dictionary must contain both
  # CCD_IMAGE_SATURATION and SATURATED_VALUE items.
  if (ccd_image_saturation is None):
    ccd_image_saturation = dynamic_range
  if (saturated_value is None):
    saturated_value = dynamic_range

  # The size must match the image dimensions.  The length along the
  # slow axis is size1, the length along the fast axis is size2.
  if (size1 is None):
    size1 = data.focus()[0]
  if (size2 is None):
    size2 = data.focus()[1]

  # The slow (vertical) axis is x, while size2 is the fast
  # (horizontal) y-axis.
  if (beam_center_x is None):
    beam_center_x = pixel_size * size2
  if (beam_center_y is None):
    beam_center_y = pixel_size * size1

  # There is no sensible default for the active areas, distance, nor
  # wavelength.  XXX But setting wavelength to zero may be disastrous?
  if (active_areas is None):
    active_areas = flex.int()
  if (distance is None):
    distance = 0
  if (wavelength is None):
    wavelength = 0

  return (dict(ACTIVE_AREAS         = active_areas,
               BEAM_CENTER_X        = beam_center_x,
               BEAM_CENTER_Y        = beam_center_y,
               CCD_IMAGE_SATURATION = ccd_image_saturation,
               DATA                 = data,
               DISTANCE             = distance,
               PIXEL_SIZE           = pixel_size,
               SATURATED_VALUE      = saturated_value,
               SIZE1                = size1,
               SIZE2                = size2,
               TIMESTAMP            = timestamp,
               SEQUENCE_NUMBER      = sequence_number,
               WAVELENGTH           = wavelength,
               xtal_target          = xtal_target))


def dwritef(d, dirname = None, basename = None):
  """The dwritef() function pickles the dictionary pointed to by @p d
  to the file whose directory and filename portions are pointed to by
  @p dirname and @p basename, respectively.  The directory at @p
  dirname, as well as any intermediate directories, are recursively
  created if they do not already exist.  The name of the written file
  is the concatenation of the @p basename parameter, and the timestamp
  and sequence number from the dictionary, @p d.  In order to be able
  to read the file later on, the file name must be of the form
  path_NNN.ext, where path and ext are arbitrary leaders and trailers,
  and NNN is a zero-padded integer number of arbitrary length.

  @param d        Dictionary, as created by e.g. dpack()
  @param dirname  Directory portion of output file
  @param basename Filename prefix of output file
  """

  if (basename is None):
    basename = ""
  if (dirname is None):
    dirname = "."
  if (not os.path.isdir(dirname)):
    os.makedirs(dirname)

  path   = os.path.join(dirname, basename       \
                          +      d['TIMESTAMP'] \
                          +      ("_%05d.pickle" % d['SEQUENCE_NUMBER']))
  easy_pickle.dump(path, d)
  return (path)


def env_detz(env):
  """The env_detz() function returns the position of the Ds1 detector
  on the z-axis in mm.  The zero-point is as far away as possible from
  the sample, and values decrease as the detector is moved towards the
  sample.  According to Marvin Seibert, the current configuration
  allows for approximately 50 cm of in-vacuum motion.

  @param env Environment object
  @return    Detector z-position, in mm
  """

  if (env is not None):
    pv = env.epicsStore().value("CXI:DS1:MMS:06.RBV")
    if (pv is None):
      # XXX Should probably return None if the RBV cannot be read.
      pv = env.epicsStore().value("CXI:DS1:MMS:06")
    if (pv is not None and len(pv.values) == 1):
      return (pv.values[0])
  return (None)


def env_distance(env):
  """The env_distance() function returns distance from the sample to
  the detector in mm.  According to an e-mail from Garth Williams, a
  value of -500 mm corresponds to a sample-detector distance of about
  73 mm.

  @param env Environment object
  @return    Detector-sample distance, in mm
  """

  detz = env_detz(env)
  if (detz is not None):
    #return (detz + (500 + 73))
    return (detz + 578) # Experimentally determined for run 22
  return (None)


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
  for pvname, si_len in dia.iteritems():
    pv = env.epicsStore().value(pvname)

    # XXX Why is this an EpicsPvTime object?  The absorption
    # coefficient of Si is E-18 * n_{0} * lambda^2, (for lambda >= 5
    # um, Schroder, D. K., R. N. Thomos, and J. C. Swartz, IEEE
    # Trans. Electron. Dev. ED-25, 2(1978) 254-261)

    #print "For ", pvname, " got ", pv, " and ", pv.values[0]
    if (pv is not None # and pv.units          == "mm"
        and                len(pv.values)    == 1
        and                abs(pv.values[0]) <  7):
      si_tot += si_len
  return (si_tot)

def evt_time(evt = None):
  """The evt_time() function returns a tuple of the time in seconds
  and milliseconds.  The return value reflects the time at which @p
  evt occurred.  If @p evt is @c None or does not contain a time, the
  current time is used.

  @param evt Event data object, a configure object
  @return    Tuple of the time in seconds and milliseconds
  """

  if (evt is not None):
    t = evt.getTime()
    if (t is not None):
      s  = t.seconds()
      ms = t.nanoseconds() // 1000000
      return (s, ms)

  t  = time.time()
  s  = int(math.floor(t))
  ms = int(round((t - s) * 1000))
  return (s, ms)


def evt_timestamp(evt = None):
  """The evt_timestamp() function returns a string representation of
  an extended human-readable ISO 8601 timestamp.  The return value
  reflects the time at which @p evt occurred.  If @p evt is @c None or
  does not contain a time, the current time is used.  Millisecond
  accuracy is sufficient, because at 120 Hz, shots are taken at 8.3 ms
  intervals.  XXX Note that the string contains colons which may look
  funny in the Finder.

  @param evt Event data object, a configure object
  @return    Human-readable ISO 8601 timestamp in string
             representation
  """

  s, ms = evt_time(evt=evt)
  timestamp = time.strftime("%Y-%m-%dT%H:%MZ%S", time.gmtime(s)) \
      +       (".%03d" % ms)
  return (timestamp)


def evt_wavelength(evt):
  """The evt_wavelength() function returns the wavelength in
  Aangstroem of the event pointed to by @p evt.  From Margaritondo &
  Rebernik Ribic (2011): the dimensionless relativistic gamma-factor
  is derived from beam energy in MeV and the electron rest mass, K is
  a dimensionless "undulator parameter", and L is the macroscopic
  undulator period in Aangstroem (XXX).  See also
  http://ast.coe.berkeley.edu/srms/2007/Lec10.pdf.

  @param evt Event data object, a configure object
  @return    Wavelength, in Aangstroem
  """

  if (evt is not None):
    # XXX EBeam object (of type bld.BldDataEBeam or bld.BldDataEBeamV0)
    ebeam = evt.getEBeam()
    if (ebeam is not None):
      gamma = ebeam.fEbeamL3Energy / 0.510998910
      K     = 3.5
      L     = 3.0e8
      return (L / (2 * gamma**2) * (1 + K**2 / 2))
  return (None)

def getOptBool(s):
  if s is None: return False
  elif isinstance(s, bool):
    return s
  s = s.strip().lower()
  return s == "true"

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

  if (s is None or s == ""):
    return None
  return (int(s))

def getOptFloat(s):
  """Return a single float.
  """

  if (s is None or s == ""):
    return None
  return (float(s))

def getOptROI(s):
  """Return a tuple ofthe region of interest.
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


def image(address, config, evt, env, sections):
  """Assemble the uint16 detector image, and sum it up as int32.  Sum
  the image of squared intensities as uint64.  XXX Documentation! XXX
  Would be nice to get rid of the constant string names.  XXX Better
  named evt_image()?

  @param address  Address string XXX Que?!
  @param config   XXX
  @param evt      Event data object, a configure object
  @param env      Environment object
  @param sections XXX
  @return         XXX
  """

  if (address == "CxiDs1-0|Cspad-0"):
    quads = evt.getCsPadQuads(address, env)
    if (quads is not None):
      if (sections is not None):
        return (CsPadDetector(quads, config, sections))
      else:
        # XXX This is obsolete code, provided for backwards
        # compatibility with the days before detector metrology was
        # used.
        qimages = numpy.empty((4, npix_quad, npix_quad), dtype = "uint16")
        for q in quads:
          qimages[q.quad()] = CsPadElement(q.data(), q.quad(), config)
        return (numpy.vstack((numpy.hstack((qimages[0], qimages[1])),
                              numpy.hstack((qimages[3], qimages[2])))))

  elif (address == "CxiSc1-0|Cspad2x2-0"):
    quads = evt.get(xtc.TypeId.Type.Id_Cspad2x2Element, address)
    if (quads is not None):
      return (CsPad2x2Image(quads.data(), config, sections))
  return (None)


def image_central(address, config, evt, env):
  """The image_central() function returns a list of the raw data
  arrays from the sections closest to the detector centre for each
  quadrant.

  @param address Address string XXX Que?!
  @param config  XXX
  @param evt     Event data object, a configure object
  @param env     Environment object
  @return        List of numpy arrays of the central sections
  """
  if (address != "CxiDs1-0|Cspad-0"):
    return (None)

  quads = evt.getCsPadQuads(address, env)
  if (quads is None):
    return (None)
  mask = map(config.sections, xrange(4))

  # The section closest to the detector centre has index 1.
  l = []
  s = 1
  for q in quads:
    if (s not in mask[q.quad()]):
      continue
    l.append(q.data()[s])
  return (l)


def image_central2(address, config, evt, env):
  """Low-levelish implementation of image_central().
  """

  # Get the first XTC object matching the address and the type.
  xtcObj = evt.findFirstXtc(address = address,
                            typeId  = xtc.TypeId.Type.Id_CspadElement)
  if (not xtcObj):
    return None

  l = []
  p = xtcObj.payload()
  s = 1
  for i in xrange(config.numQuads()):
    l.append(p.data(config)[s])
    p = p.next(config)
  return (l)


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
    for x in xrange(xlim[0], xlim[1]):
        for y in xrange(ylim[0], ylim[1]):
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
