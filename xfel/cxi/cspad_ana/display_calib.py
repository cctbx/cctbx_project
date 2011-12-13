#! /usr/bin/python
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# XXX Jiffy summary here
#
# $Id$

import numpy
import sys

import matplotlib
from   matplotlib.collections import PatchCollection
import matplotlib.patches     as     mpatches
import matplotlib.pyplot      as     plt

from optparse import OptionParser

from parse_calib import calib2sections


# XXX http://www.artima.com/weblogs/viewpost.jsp?thread=4829
def display_calib(dirname, right, verbose):
  """XXX Docstring, in fact revise all the documentation

  @param dirname Directory with calibration information
  @param right   @c True to restrict rotations to right angles
  @param verbose @c True to print ASIC coordinates
  """

  fig     = plt.figure(figsize = (10, 10))
  ax      = plt.axes([0, 0, 1, 1])
  plt.axis([0, 1765, 1765, 0])

  colours  = []
  patches  = []
  sections = calib2sections(dirname)
  for q in xrange(len(sections)):
    for s in xrange(len(sections[q])):

      # Get the vertices of the section s in quadrant q, and round
      # rotation angles to integer multiples of 90 degrees by default.
      # Change from matrix-coordinate system to screen coordinate
      # system, where origin is in the top left corner, the first
      # coordinate increases to the right, and the second coordinate
      # increases downwards.  Ensure that the eight sections within
      # the quadrants are coloured consistently.
      vertices = sections[q][s].corners(right)
      for i in xrange(len(vertices)):
        vertices[i] = [vertices[i][1], vertices[i][0]]

      art = mpatches.Circle(vertices[0], 10)
      patches.append(art)
      colours.append(s)

      polygon = mpatches.Polygon(vertices)
      patches.append(polygon)
      colours.append(s)

      plt.text(sections[q][s].center[1], sections[q][s].center[0],
               "(%1d, %1d)" % (q, s),
               family = "sans-serif",
               size   = 14,
               ha     = "center",
               va     = "center")

      # Assuming that rotations are integer multiples of 90 degrees,
      # print the ASIC coordinates in "spotfinder" format, ordered by
      # quadrant, section, and ASIC.  XXX This only makes sense for
      # right = True.
      if (verbose):
        vertices = sections[q][s].corners_asic()
        print "(%4d, %4d, %4d, %4d)" % \
            (vertices[0][0], vertices[0][1], vertices[0][2], vertices[0][3])
        print "(%4d, %4d, %4d, %4d)" % \
            (vertices[1][0], vertices[1][1], vertices[1][2], vertices[1][3])

  collection = PatchCollection(patches, cmap = matplotlib.cm.jet, alpha = 0.4)
  collection.set_array(numpy.array(colours))
  ax.add_collection(collection)
  ax.set_xticks([])
  ax.set_yticks([])
  plt.show()

  return (0)


# Run with "display_calib.py
# /reg/d/ana11/cxi/data/CSPAD-metrology/run4/CxiDs1.0:Cspad.0".  XXX
# http://docs.python.org/library/optparse.html
if (__name__ == "__main__"):
  parser = OptionParser()
  parser.add_option("-r", "--rotate",
                    action  = "store_false",
                    default = True,
                    dest    = "right",
                    help    = "Allow sections to rotate by arbitrary angles")
  parser.add_option("-v", "--verbose",
                    action  = "store_true",
                    default = False,
                    dest    = "verbose",
                    help    = "Print ordered list of diagonal ASIC corners")

  (options, args) = parser.parse_args()
  for arg in args:
    ret = display_calib(arg, options.right, options.verbose)
    if (ret != 0):
      sys.exit(ret)
