# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.cspad_average
#
# $Id

from __future__ import absolute_import, division, print_function
from six.moves import range

import sys, copy
from six.moves import zip

def run(argv=None):
  """Compute mean, standard deviation, and maximum projection images
  from a set of CSPAD cbf images given on the command line.

  @param argv Command line argument list
  @return     @c 0 on successful termination, @c 1 on error, and @c 2
              for command line syntax errors
  """

  import libtbx.load_env

  from libtbx import option_parser
  from scitbx.array_family import flex
  from dxtbx.format.Registry import Registry
  from xfel.cftbx.detector.cspad_cbf_tbx import cbf_file_to_basis_dict, write_cspad_cbf
#  from xfel.cxi.cspad_ana import cspad_tbx
#  from iotbx.detectors.cspad_detector_formats import reverse_timestamp

  if argv is None:
    argv = sys.argv
  command_line = (option_parser.option_parser(
    usage="%s [-v] [-a PATH] [-m PATH] [-s PATH] " \
    "image1 image2 [image3 ...]" % libtbx.env.dispatcher_name)
                  .option(None, "--average-path", "-a",
                          type="string",
                          default=None,
                          dest="avg_path",
                          metavar="PATH",
                          help="Write average image to PATH")
                  .option(None, "--maximum-path", "-m",
                          type="string",
                          default=None,
                          dest="max_path",
                          metavar="PATH",
                          help="Write maximum projection image to PATH")
                  .option(None, "--stddev-path", "-s",
                          type="string",
                          default=None,
                          dest="stddev_path",
                          metavar="PATH",
                          help="Write standard deviation image to PATH")
                  .option(None, "--verbose", "-v",
                          action="store_true",
                          default=False,
                          dest="verbose",
                          help="Print more information about progress")
                  ).process(args=argv[1:])

  # Note that it is not an error to omit the output paths, because
  # certain statistics could still be printed, e.g. with the verbose
  # option.
  paths = command_line.args
  if len(paths) == 0:
    command_line.parser.print_usage(file=sys.stderr)
    return 2

  # Loop over all images and accumulate statistics.
  nfail = 0
  nmemb = 0
  for path in paths:
    if command_line.options.verbose:
      sys.stdout.write("Processing %s...\n" % path)

    try:
      # Promote the image to double-precision floating point type.
      # All real-valued flex arrays have the as_double() function.
      # Warn if the header items across the set of images do not match
      # up.  Note that discrepancies regarding the image size are
      # fatal.
      if not 'reader' in locals():
        reader = Registry.find(path)
      img = reader(path)
      if 'detector' in locals():
        test_detector = img.get_detector()
        if len(test_detector) != len(detector):
          sys.stderr.write("Detectors do not have the same number of panels\n")
          return 1
        for t, d in zip(test_detector, detector):
          if t.get_image_size() != d.get_image_size():
            sys.stderr.write("Panel sizes do not match\n")
            return 1
          if t.get_pixel_size() != d.get_pixel_size():
            sys.stderr.write("Pixel sizes do not match\n")
            return 1
          if t.get_d_matrix() != d.get_d_matrix():
            sys.stderr.write("Detector panels are not all in the same location. The average will use the positions of the first image.\n")
        detector = test_detector
      else:
        detector = img.get_detector()

      data = [img.get_raw_data()[i].as_1d().as_double() for i in range(len(detector))]
      wavelength = img.get_beam().get_wavelength()
      distance = flex.mean(flex.double([d.get_directed_distance() for d in detector]))

    except Exception:
      nfail += 1
      continue

    # The sum-of-squares image is accumulated using long integers, as
    # this delays the point where overflow occurs.  But really, this
    # is just a band-aid...
    if nmemb == 0:
      max_img = copy.deepcopy(data)
      sum_distance = distance
      sum_img = copy.deepcopy(data)
      ssq_img = [flex.pow2(d) for d in data]
      sum_wavelength = wavelength
      metro = cbf_file_to_basis_dict(path)

    else:
      sel = [(d > max_d).as_1d() for d, max_d in zip(data, max_img)]
      for d, max_d, s in zip(data, max_img, sel): max_d.set_selected(s, d.select(s))

      sum_distance += distance
      for d, sum_d in zip(data, sum_img): sum_d += d
      for d, ssq_d in zip(data, ssq_img): ssq_d += flex.pow2(d)
      sum_wavelength += wavelength

    nmemb += 1

  # Early exit if no statistics were accumulated.
  if command_line.options.verbose:
    sys.stderr.write("Processed %d images (%d failed)\n" % (nmemb, nfail))
  if nmemb == 0:
    return 0

  # Calculate averages for measures where other statistics do not make
  # sense.  Note that avg_img is required for stddev_img.
  avg_img = [sum_d.as_double() / nmemb for sum_d in sum_img]
  avg_distance = sum_distance / nmemb
  avg_wavelength = sum_wavelength / nmemb

  def make_tiles(data, detector):
    """
    Assemble a tiles dictionary as required by write_cspad_cbf, consisting of 4 arrays of shape 8x185x388.
    Assumes the order in the data array matches the order of the enumerated detector panels.
    """
    assert len(data) == 64
    tiles = {}
    s, f = 185, 194

    for q_id in range(4):
      tiles[0,q_id] = flex.double((flex.grid(s*8, f*2)))
      for s_id in range(8):
        for a_id in range(2):
          asic_idx = (q_id*16) + (s_id*2) + a_id
          asic = data[asic_idx]
          asic.reshape(flex.grid((s, f)))

          tiles[0, q_id].matrix_paste_block_in_place(asic, s_id*s, a_id*f)
      tiles[0, q_id].reshape(flex.grid((8, s, f*2)))

    return tiles


  # Output the average image, maximum projection image, and standard
  # deviation image, if requested.
  if command_line.options.avg_path is not None:
    tiles = make_tiles(avg_img, detector)
    write_cspad_cbf(tiles, metro, 'cbf', None, command_line.options.avg_path, avg_wavelength, avg_distance)

  if command_line.options.max_path is not None:
    tiles = make_tiles(max_img, detector)
    write_cspad_cbf(tiles, metro, 'cbf', None, command_line.options.max_path, avg_wavelength, avg_distance)

  if command_line.options.stddev_path is not None:
    stddev_img = [ssq_d.as_double() - sum_d.as_double() * avg_d for ssq_d, sum_d, avg_d in zip(ssq_img, sum_img, avg_img)]

    # Accumulating floating-point numbers introduces errors, which may
    # cause negative variances.  Since a two-pass approach is
    # unacceptable, the standard deviation is clamped at zero.
    for stddev_d in stddev_img:
      stddev_d.set_selected(stddev_d < 0, 0)

    if nmemb == 1:
      stddev_img = [flex.sqrt(stddev_d) for stddev_d in stddev_img]
    else:
      stddev_img = [flex.sqrt(stddev_d / (nmemb - 1)) for stddev_d in stddev_img]

    tiles = make_tiles(stddev_img, detector)
    write_cspad_cbf(tiles, metro, 'cbf', None, command_line.options.stddev_path, avg_wavelength, avg_distance)

  return 0


if __name__ == '__main__':
  sys.exit(run())
