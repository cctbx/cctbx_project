from __future__ import absolute_import, division, print_function
from six.moves import range
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.make_mask
#
# $Id: make_mask.py 411 2013-10-16 22:17:45Z aaron $
#
# This code reads three cctbx.xfel pickle format images and builds a mask from
# them.  The first image should be an average from a dark run, the second the
# standard deviation from that run.  The third image should be a maximum projection
# from a run with the beam on.
#
# The result is an image with all pixels which are valid for use set to 0, and
# those that are invalid set to -2 by default, or the value of the option passed in
# to mask_pix_val.
#

from dxtbx.format.Registry import Registry
from xfel.cxi.cspad_ana.cspad_tbx import dpack, dwritef2
from scitbx.array_family import flex
import sys

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
    usage="%s [-v] [-p poly_mask] [-c circle_mask] [-a avg_max] [-s stddev_max] [-m maxproj_min] [-x mask_pix_val] [-o output] avg_path stddev_path max_path" % libtbx.env.dispatcher_name)
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
                  .option(None, "--avg_max", "-a",
                          type="float",
                          default=2000.0,
                          dest="avg_max",
                          help="Maximum ADU that pixels in the average image are allowed to have before masked out")
                  .option(None, "--stddev_max", "-s",
                          type="float",
                          default=10.0,
                          dest="stddev_max",
                          help="Maximum ADU that pixels in the standard deviation image are allowed to have before masked out")
                  .option(None, "--maxproj_min", "-m",
                          type="float",
                          default=300.0,
                          dest="maxproj_min",
                          help="Minimum ADU that pixels in the maximum projection image are allowed to have before masked out")
                  .option(None, "--mask_pix_val", "-x",
                          type="int",
                          default=-2,
                          dest="mask_pix_val",
                          help="Value for masked out pixels")
                  .option(None, "--detector_format_version", "-d",
                          type="string",
                          default=None,
                          dest="detector_format_version",
                          help="detector format version string")
                  .option(None, "--output", "-o",
                          type="string",
                          default="mask_.pickle",
                          dest="destpath",
                          help="output file path, should be *.pickle")
                  ).process(args=argv[1:])

  # Must have exactly three remaining arguments.
  paths = command_line.args
  if (len(paths) != 3):
    command_line.parser.print_usage(file=sys.stderr)
    return

  if command_line.options.detector_format_version is None:
    address = timestamp = None
  else:
    from xfel.cxi.cspad_ana.cspad_tbx import evt_timestamp
    from iotbx.detectors.cspad_detector_formats import address_and_timestamp_from_detector_format_version
    address, timestamp = address_and_timestamp_from_detector_format_version(command_line.options.detector_format_version)
    timestamp = evt_timestamp((timestamp,0))

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

  avg_path    = paths[0]
  stddev_path = paths[1]
  max_path    = paths[2]

  # load the three images
  format_class = Registry.find(avg_path)
  avg_f = format_class(avg_path)
  avg_i = avg_f.get_detectorbase()
  avg_d = avg_i.get_raw_data()

  stddev_f = format_class(stddev_path)
  stddev_i = stddev_f.get_detectorbase()
  stddev_d = stddev_i.get_raw_data()

  max_f = format_class(max_path)
  max_i = max_f.get_detectorbase()
  max_d = max_i.get_raw_data()

  # first find all the pixels in the average that are less than zero or greater
  # than a cutoff and set them to the masking value
  avg_d.set_selected((avg_d <= 0) | (avg_d > command_line.options.avg_max), command_line.options.mask_pix_val)

  # set all the rest of the pixels to zero.  They will be accepted
  avg_d.set_selected(avg_d != command_line.options.mask_pix_val, 0)

  # mask out the overly noisy or flat pixels
  avg_d.set_selected(stddev_d <= 0, command_line.options.mask_pix_val)
  avg_d.set_selected(stddev_d >= command_line.options.stddev_max, command_line.options.mask_pix_val)

  # these are the non-bonded pixels
  avg_d.set_selected(max_d < command_line.options.maxproj_min, command_line.options.mask_pix_val)

  # calculate the beam center
  panel = avg_f.get_detector()[0]
  bcx, bcy = panel.get_beam_centre(avg_f.get_beam().get_s0())

  if poly_mask is not None or circle_mask is not None:
    minx = miny = 0
    maxx = avg_d.focus()[0]
    maxy = avg_d.focus()[1]
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

    sel = avg_d == command_line.options.mask_pix_val
    for j in range(miny, maxy):
      for i in range(minx, maxx):
        idx = j * avg_d.focus()[0] + i
        if not sel[idx]:
          if poly_mask is not None and point_in_polygon((i,j),poly_mask):
            sel[idx] = True
          elif circle_mask is not None and point_inside_circle(i,j,circle_x,circle_y,radius):
            sel[idx] = True
    avg_d.set_selected(sel,command_line.options.mask_pix_val)

  # have to re-layout the data to match how it was stored originally
  shifted_int_data_old = avg_d
  shifted_int_data_new = shifted_int_data_old.__class__(
    flex.grid(shifted_int_data_old.focus()))
  shifted_int_data_new += command_line.options.mask_pix_val

  phil = avg_i.horizons_phil_cache
  manager = avg_i.get_tile_manager(phil)

  for i,shift in enumerate(manager.effective_translations()):
    shift_slow = shift[0]
    shift_fast = shift[1]

    ur_slow = phil.distl.detector_tiling[4 * i + 0] + shift_slow
    ur_fast = phil.distl.detector_tiling[4 * i + 1] + shift_fast
    ll_slow = phil.distl.detector_tiling[4 * i + 2] + shift_slow
    ll_fast = phil.distl.detector_tiling[4 * i + 3] + shift_fast

    #print "Shifting tile at (%d, %d) by (%d, %d)" % (ur_slow-shift_slow, ur_fast-shift_fast, -shift_slow, -shift_fast)

    shifted_int_data_new.matrix_paste_block_in_place(
      block = shifted_int_data_old.matrix_copy_block(
        i_row=ur_slow,i_column=ur_fast,
        n_rows=ll_slow-ur_slow, n_columns=ll_fast-ur_fast),
      i_row = ur_slow - shift_slow,
      i_column = ur_fast - shift_fast
    )

  d = dpack(
    active_areas=avg_i.parameters['ACTIVE_AREAS'],
    address=address,
    beam_center_x=bcx,
    beam_center_y=bcy,
    data=shifted_int_data_new,
    distance=avg_i.distance,
    timestamp=timestamp,
    wavelength=avg_i.wavelength,
    xtal_target=None,
    pixel_size=avg_i.pixel_size,
    saturated_value=avg_i.saturation)

  dwritef2(d, command_line.options.destpath)

  #the minimum number of pixels to mask out cooresponding to the interstitial regions for the CS-PAD
  min_count  = 818265 # (1765 * 1765) - (194 * 185 * 64)
  masked_out = len(avg_d.as_1d().select((avg_d == command_line.options.mask_pix_val).as_1d()))
  assert masked_out >= min_count

  print("Masked out %d pixels out of %d (%.2f%%)"% \
    (masked_out-min_count,len(avg_d)-min_count,(masked_out-min_count)*100/(len(avg_d)-min_count)))

if (__name__ == "__main__"):
  sys.exit(run())
