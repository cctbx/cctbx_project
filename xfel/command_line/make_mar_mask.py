from __future__ import division
from __future__ import print_function
from six.moves import range
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.make_mar_mask
#
# $Id:
#
# This code dumps a flex array containing masked out pixels given certain specifications
# The result is an image with all pixels which are valid for use set to 0, and
# those that are invalid set to -2 by default, or the value of the option passed in
# to mask_pix_val.
#

from scitbx.array_family import flex
import sys
from libtbx import easy_pickle
from xfel.cxi.cspad_ana.cspad_tbx import dpack, dwritef2

def point_in_polygon(point, poly):
  """ Determine if a point is inside a given polygon or not.  Polygon is a list of (x,y) pairs.
  Code adapted from a dials polygon clipping test algorithm"""
  if len(poly) < 3: return False

  inside = False
  for i in range(len(poly)):
    j = (i+1) % len(poly)
    if (((poly[i][1] > point[1]) != (poly[j][1] > point[1])) and
      (point[0] < (poly[j][0] - poly[i][0]) * (point[1] - poly[i][1]) /
                  (poly[j][1] - poly[i][1]) + poly[i][0])):
      inside = not inside
  return inside

def point_inside_circle(x,y,center_x,center_y,radius):
  """Determine if a given point (x,y) is inside a circle whose center is at
  (center_x,center_y) with radius x."""
  return (x-center_x)**2 + (y - center_y)**2 < radius**2

def run(argv=None):
  import libtbx.option_parser

  if (argv is None):
    argv = sys.argv

  command_line = (libtbx.option_parser.option_parser(
    usage="%s [-v] [-p poly_mask] [-c circle_mask]  [-x mask_pix_val] [-o output] -W mask_width -H mask_height" % libtbx.env.dispatcher_name)
                  .option(None, "--verbose", "-v",
                          action="store_true",
                          default=False,
                          dest="verbose",
                          help="Print more information about progress")
                  .option(None, "--poly_mask", "-p",
                          type="string",
                          default=None,
                          dest="poly_mask",
                          help="Polygon to mask out.  Comma-seperated string of xy pairs.")
                  .option(None, "--circle_mask", "-c",
                          type="string",
                          default=None,
                          dest="circle_mask",
                          help="Circle to mask out.  Comma-seperated string of x, y, and radius.")
                  .option(None, "--mask_pix_val", "-x",
                          type="int",
                          default=-2,
                          dest="mask_pix_val",
                          help="Value for masked out pixels")
                  .option(None, "--mask_width", "-W",
                          type="int",
                          default=None,
                          dest="mask_width",
                          help="Width of output mask")
                   .option(None, "--mask_height", "-H",
                          type="int",
                          default=None,
                          dest="mask_height",
                          help="Height of output mask")
                   .option(None, "--output", "-o",
                          type="string",
                          default="mask.pickle",
                          dest="destpath",
                          help="Output file path, should be *.pickle")
                   .option(None, "--pixel_size", "-s",
                          type="float",
                          default=None,
                          dest="pixel_size",
                          help="Pixel size for detector")
                  ).process(args=argv[1:])

  # Must have width and height set
  if command_line.options.mask_height is None or command_line.options.mask_width is None:
    command_line.parser.print_usage(file=sys.stderr)
    return

  poly_mask = None
  if not command_line.options.poly_mask == None:
    poly_mask = []
    poly_mask_tmp = command_line.options.poly_mask.split(",")
    if len(poly_mask_tmp) % 2 != 0:
      command_line.parser.print_usage(file=sys.stderr)
      return
    odd = True
    for item in poly_mask_tmp:
      try:
        if odd:
          poly_mask.append(int(item))
        else:
          poly_mask[-1] = (poly_mask[-1],int(item))
      except ValueError:
        command_line.parser.print_usage(file=sys.stderr)
        return
      odd = not odd

  circle_mask = None
  if command_line.options.circle_mask is not None:
    circle_mask_tmp = command_line.options.circle_mask.split(",")
    if len(circle_mask_tmp) != 3:
      command_line.parser.print_usage(file=sys.stderr)
      return
    try:
      circle_mask = (int(circle_mask_tmp[0]),int(circle_mask_tmp[1]),int(circle_mask_tmp[2]))
    except ValueError:
      command_line.parser.print_usage(file=sys.stderr)
      return

  mask = flex.int(flex.grid(command_line.options.mask_width,
                             command_line.options.mask_height))

  if poly_mask is not None or circle_mask is not None:
    minx = miny = 0
    maxx = mask.focus()[0]
    maxy = mask.focus()[1]
    if poly_mask is not None:
      minx = min([x[0] for x in poly_mask])
      miny = min([y[1] for y in poly_mask])
      maxx = max([x[0] for x in poly_mask])
      maxy = max([y[1] for y in poly_mask])
    if circle_mask is not None:
      circle_x, circle_y, radius = circle_mask

      if circle_x - radius < minx: minx = circle_x - radius
      if circle_y - radius < miny: miny = circle_y - radius
      if circle_x + radius > maxx: maxx = circle_x + radius
      if circle_y + radius > maxy: maxy = circle_y + radius

    sel = mask == command_line.options.mask_pix_val
    for j in range(miny, maxy):
      for i in range(minx, maxx):
        idx = j * mask.focus()[0] + i
        if not sel[idx]:
          if poly_mask is not None and point_in_polygon((i,j),poly_mask):
            sel[idx] = True
          elif circle_mask is not None and point_inside_circle(i,j,circle_x,circle_y,radius):
            sel[idx] = True
    mask.set_selected(sel,command_line.options.mask_pix_val)

  masked_out = len(mask.as_1d().select((mask == command_line.options.mask_pix_val).as_1d()))

  print("Masked out %d pixels out of %d (%.2f%%)"% \
    (masked_out,len(mask),(masked_out)*100/(len(mask))))

  easy_pickle.dump(command_line.options.destpath, mask)

  d = dpack(
    active_areas=[0,0,command_line.options.mask_width,command_line.options.mask_height],
    address=None,
    beam_center_x=None,
    beam_center_y=None,
    data=mask,
    distance=None,
    timestamp=None,
    wavelength=1,
    xtal_target=None,
    pixel_size=command_line.options.pixel_size,
    saturated_value=None)

  dwritef2(d, command_line.options.destpath)

if (__name__ == "__main__"):
  sys.exit(run())
