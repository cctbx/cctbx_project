# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

from __future__ import absolute_import, division, print_function
from six.moves import range

import sys

from libtbx import easy_pickle
from libtbx import phil
from optparse import OptionParser
from scitbx import matrix
from xfel.cxi.cspad_ana.parse_calib import calib2sections


def metrology2phil(calib_dir, verbose):
  """XXX Should really review the SLAC progress since last time
  around!  XXX Note that this is all SLAC-specific (as is the whole
  thing, I guess).
  """
  # XXX Can this fail?  How?
  sections = calib2sections(calib_dir)
  if (sections is None):
    return (None)

  return sections2phil(sections, verbose)

def sections2phil(sections, verbose):
  from xfel.cftbx.detector.metrology import master_phil

  # Properties of CSPad pixels (Philipp et al., 2007).  The counters
  # are 14 bits wide, and the pixels are square with a side length of
  # 110 um.  Because cspad_tbx depends on pyana, it may fail to import
  # in which case a hardcoded fallback is provided.
  try:
    from xfel.cxi.cspad_ana.cspad_tbx import cspad_saturated_value as sv
    from xfel.cxi.cspad_ana.cspad_tbx import pixel_size as ps
    saturated_value = sv
    pixel_size = ps * 1e-3
  except ImportError:
    saturated_value = 90000
    pixel_size = 110e-6

  # Build the Phil object.  XXX Should have det-z?  Probably not,
  # because then it aint't just metrology anymore.  In fact, under
  # orthographic projection the whole translation/orientation thing
  # can be scrapped for the detector.  XXX look up include,
  # include_scope for phil XXX Hardcoded address for now.
  address = "CxiDs1-0|Cspad-0"
  metrology_str = "detector { serial = %d\n" % 0
  metrology_str += "label = %s\n" % address

  # The center of the detector is defined as the average of all the
  # sections it contains.  The first coordinate of the matrix-oriented
  # coordinate system of the Section class maps to -y, the second to
  # x, and the third to z.  While the detector still sits on the
  # origin (zero translation), the mapping from the origin in Section
  # coordinates is still needed.
  t_d = [0, 0, 0]
  nmemb = 0
  for p in range(len(sections)):
    for s in range(len(sections[p])):
      t_d[0] += +sections[p][s].center[1]
      t_d[1] += -sections[p][s].center[0]
      t_d[2] += 0
      nmemb += 1
  if (nmemb > 0):
    for i in range(3):
      t_d[i] /= nmemb
  metrology_str += "translation = 0, 0, 0\n"

  o_d = matrix.col([0, 0, 1]).axis_and_angle_as_unit_quaternion(
    angle=0, deg=True)
  metrology_str += "orientation = %s, %s, %s, %s\n" % \
      tuple(repr(c) for c in o_d)

  for p in range(len(sections)):
    # Loop over quadrants (panels).  XXX Translation of panels is
    # wrongly set to to centre of sections in the quadrant
    # w.r.t. center of the detector.
    metrology_str += "panel { serial = %d\n" % p
    metrology_str += "translation = "
    t_p = [0, 0, 0]
    for s in range(len(sections[p])):
      t_p[0] += +sections[p][s].center[1] - t_d[0]
      t_p[1] += -sections[p][s].center[0] - t_d[1]
      t_p[2] +=  0                        - t_d[2]
    for i in range(3):
      t_p[i] /= len(sections[p])
      metrology_str += "%s" % repr(t_p[i] * pixel_size)
      if (i < 2):
        metrology_str += ", "
      else:
        metrology_str += "\n"

    # XXX Orientation wrongly set to zero.
    o_p = matrix.col([0, 0, 1]).axis_and_angle_as_unit_quaternion(
      angle=0, deg=True)
    metrology_str += "orientation = %s, %s, %s, %s\n" % \
        tuple(repr(c) for c in o_p)

    for s in range(len(sections[p])):
      # Loop over sensors (sections or two-by-one:s).  Note that
      # sensors are rotated by -90 degrees in the SLAC metrology
      # convention w.r.t. their appearance in the XTC stream.
      s_t = [0, 0, 0]
      s_t[0] = +sections[p][s].center[1] - t_p[0] - t_d[0]
      s_t[1] = -sections[p][s].center[0] - t_p[1] - t_d[1]
      s_t[2] =  0                        - t_p[2] - t_d[2]

      metrology_str += "sensor { serial = %d\n" % s
      metrology_str += "translation = "
      for i in range(3):
        metrology_str += "%s" % repr(s_t[i] * pixel_size)
        if (i < 2):
          metrology_str += ", "
        else:
          metrology_str += "\n"

      s_o = matrix.col([0, 0, 1]).axis_and_angle_as_unit_quaternion(
        angle=sections[p][s].angle - 90, deg=True)
      metrology_str += "orientation = %s, %s, %s, %s\n" % \
          tuple(repr(c) for c in s_o)

      for a in range(2):
        # The ASIC:s of the CSPad are 185 rows by 194 columns.  The
        # ASIC:s are horizontally aligned within a section, with a
        # three-column gap between them.
        metrology_str += "asic { serial = %d\n" % a

        if (a == 0):
          metrology_str += "translation = %s, %s, %s\n" % (
            repr(-(194 + 3) / 2 * pixel_size), "0", "0") # XXX hardcoded!
        else:
          metrology_str += "translation = %s, %s, %s\n" % (
            repr(+(194 + 3) / 2 * pixel_size), "0", "0") # XXX hardcoded!
        metrology_str += "orientation = 1, 0, 0, 0\n"
        metrology_str += "pixel_size = %s, %s\n" % (
          repr(pixel_size), repr(pixel_size))
        metrology_str += "dimension = %d, %d\n" % (194, 185) # XXX hardcoded!
        metrology_str += "saturation = %s\n" % repr(float(saturated_value))
        metrology_str += "}\n"
      metrology_str += "}\n"
    metrology_str += "}\n"
  metrology_str += "}\n"

  metrology_phil = master_phil.fetch(
    sources=[phil.parse(metrology_str)])
  return (metrology_phil)


# Run with "metrology2phil.py -o metrology.pkl
# ../metrology/CSPad/run4/CxiDs1.0:Cspad.0".  XXX
# http://docs.python.org/library/optparse.html XXX Could support
# different output (and eventually input) options: plain-text phil and
# XML.
if (__name__ == "__main__"):
  parser = OptionParser()
  parser.add_option("-o", "--output",
                    dest="path_out",
                    help= "Write output to FILE",
                    metavar="FILE")
  parser.add_option("-v", "--verbose",
                    action="store_true",
                    default=False,
                    dest="verbose",
                    help="Print more information about progress")

  # XXX Requires exactly one argument, so output to stdout?
  (options, args) = parser.parse_args()
  if (len(args) < 1):
    sys.stderr.write("Usage: XXX [-v] [-o FILE] directory\n")

  phil = metrology2phil(args[0], options.verbose)
  if (phil is None):
    sys.exit(1)
  if (options.verbose):
    phil.show()

  # Output to file, if given, otherwise use stdout.
  if (options.path_out is not None):
    easy_pickle.dump(options.path_out, phil)
  else:
    sys.stdout.write(phil)
  sys.exit(0)
