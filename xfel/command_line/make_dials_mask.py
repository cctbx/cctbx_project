from __future__ import division
from __future__ import print_function
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.make_dials_mask
#
# This code reads three images and builds a mask from
# them.  The first image should be an average from a dark run, the second the
# standard deviation from that run.  The third image should be a maximum projection
# from a run with the beam on.
#
# The result is a dials-style pickle file containing a list of flex.bool objects.
#

from dxtbx.format.Registry import Registry
from scitbx.array_family import flex
import sys
from libtbx import easy_pickle

def run(argv=None):
  import libtbx.option_parser

  if (argv is None):
    argv = sys.argv

  command_line = (libtbx.option_parser.option_parser(
    usage="%s [-v] [-a avg_max] [-s stddev_max] [-m maxproj_min] [-o output] [-b border] avg_path stddev_path max_path" % libtbx.env.dispatcher_name)
                  .option(None, "--verbose", "-v",
                          action="store_true",
                          default=False,
                          dest="verbose",
                          help="Print more information about progress")
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
                  .option(None, "--output", "-o",
                          type="string",
                          default="mask.pickle",
                          dest="destpath",
                          help="output file path, should be *.pickle")
                  .option(None, "--border", "-b",
                          type="int",
                          default=0,
                          dest="border",
                          help="border width in pixels to mask out of each tile")
                  ).process(args=argv[1:])

  # Must have exactly three remaining arguments.
  paths = command_line.args
  if (len(paths) != 3):
    command_line.parser.print_usage(file=sys.stderr)
    return

  avg_path    = paths[0]
  stddev_path = paths[1]
  max_path    = paths[2]

  # load the three images
  format_class = Registry.find(avg_path)
  avg_f = format_class(avg_path)
  avg_d = avg_f.get_raw_data()
  if not isinstance(avg_d, tuple):
    avg_d = (avg_d,)

  stddev_f = format_class(stddev_path)
  stddev_d = stddev_f.get_raw_data()
  if not isinstance(stddev_d, tuple):
    stddev_d = (stddev_d,)

  max_f = format_class(max_path)
  max_d = max_f.get_raw_data()
  if not isinstance(max_d, tuple):
    max_d = (max_d,)

  mask = [flex.bool(flex.grid(p.focus()), True) for p in avg_d]

  for mask_p, avg_p, stddev_p, max_p in zip(mask, avg_d, stddev_d, max_d):
    # first find all the pixels in the average that are less than zero or greater
    # than a cutoff and set them to the masking value
    mask_p &= avg_p > 0
    mask_p &= avg_p <= command_line.options.avg_max

    # mask out the overly noisy or flat pixels
    mask_p &= stddev_p > 0
    mask_p &= stddev_p <= command_line.options.stddev_max # cxi.make_mask uses <

    # these are the non-bonded pixels
    mask_p &= max_p >= command_line.options.maxproj_min

    # Add a border around the image
    if command_line.options.border > 0:
      border = command_line.options.border
      height, width = mask_p.all()
      borderx = flex.bool(flex.grid(border, width), False)
      bordery = flex.bool(flex.grid(height, border), False)
      mask_p[0:border,:] = borderx
      mask_p[-border:,:] = borderx
      mask_p[:,0:border] = bordery
      mask_p[:,-border:] = bordery

  easy_pickle.dump(command_line.options.destpath, tuple(mask))

  masked_out = sum([len(mask_p.as_1d().select((~mask_p).as_1d())) for mask_p in mask])
  total = sum([len(mask_p) for mask_p in mask])

  print("Masked out %d pixels out of %d (%.2f%%)"% \
    (masked_out,total,masked_out*100/total))

if (__name__ == "__main__"):
  sys.exit(run())
