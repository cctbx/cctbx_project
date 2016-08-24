# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME dxtbx.image_average
# LIBTBX_SET_DISPATCHER_NAME cxi.image_average
#

from __future__ import division

import sys


def run(argv=None):
  """Compute mean, standard deviation, and maximum projection images
  from a set of images given on the command line.

  @param argv Command line argument list
  @return     @c 0 on successful termination, @c 1 on error, and @c 2
              for command line syntax errors
  """

  import libtbx.load_env

  from libtbx import easy_pickle, option_parser
  from scitbx.array_family import flex
  from xfel.cxi.cspad_ana import cspad_tbx
  from iotbx.detectors.cspad_detector_formats import reverse_timestamp

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
      d = easy_pickle.load(path)
      distance = d['DISTANCE']
      img = d['DATA'].as_1d().as_double()
      wavelength = d['WAVELENGTH']
      time_tuple = reverse_timestamp(d['TIMESTAMP'])

      # Warn if the header items across the set of images do not match
      # up.  Note that discrepancies regarding the image size are
      # fatal.
      if 'active_areas' in locals():
        if (active_areas != d['ACTIVE_AREAS']).count(True) != 0:
          sys.stderr.write("Active areas do not match\n")
      else:
        active_areas = d['ACTIVE_AREAS']

      if 'beam_center' in locals():
        if beam_center != (d['BEAM_CENTER_X'], d['BEAM_CENTER_Y']):
          sys.stderr.write("Beam centers do not match\n")
      else:
        beam_center = (d['BEAM_CENTER_X'], d['BEAM_CENTER_Y'])

      if 'detector_address' in locals():
        if detector_address != d['DETECTOR_ADDRESS']:
          sys.stderr.write("Detector addresses do not match\n")
      else:
        detector_address = d['DETECTOR_ADDRESS']

      if 'saturated_value' in locals():
        if saturated_value != d['SATURATED_VALUE']:
          sys.stderr.write("Saturated values do not match\n")
      else:
        saturated_value = d['SATURATED_VALUE']

      if 'size' in locals():
        if size != (d['SIZE1'], d['SIZE2']):
          sys.stderr.write("Image sizes do not match\n")
          return 1
      else:
        size = (d['SIZE1'], d['SIZE2'])
      if size != d['DATA'].focus():
        sys.stderr.write("Image size does not match pixel array\n")
        return 1

      if 'pixel_size' in locals():
        if pixel_size != d['PIXEL_SIZE']:
          sys.stderr.write("Pixel sizes do not match\n")
          return 1
      else:
        if 'PIXEL_SIZE' in d:
          pixel_size = d['PIXEL_SIZE']
        else:
          pixel_size = None


    except Exception:
      try:
        # Fall back on reading the image with dxtbx, and shoehorn the
        # extracted information into what would have been found in a
        # pickle file.  XXX This code assumes a monolithic detector!

        from dxtbx.format.Registry import Registry

        format_class = Registry.find(path)
        i = format_class(path)

        beam = i.get_beam()
        assert len(i.get_detector()) == 1
        detector = i.get_detector()[0]

        beam_center = detector.get_beam_centre(beam.get_s0())
        detector_address = format_class.__name__
        distance = detector.get_distance()
        img = i.get_raw_data().as_1d().as_double()
        pixel_size = 0.5 * sum(detector.get_pixel_size())
        saturated_value = int(round(detector.get_trusted_range()[1]))
        size = detector.get_image_size()
        time_tuple = (i.get_scan().get_epochs()[0], 0)
        wavelength = beam.get_wavelength()

        active_areas = flex.int((0, 0, size[0], size[1]))


      except Exception:
        nfail += 1
        continue


    # See also event() in xfel.cxi.cspad_ana.average_tbx.  Record the
    # base time as the timestamp of the first image.
    #
    # The sum-of-squares image is accumulated using long integers, as
    # this delays the point where overflow occurs.  But really, this
    # is just a band-aid...
    if nmemb == 0:
      max_img = img.deep_copy()
      sum_distance = distance
      sum_img = img.deep_copy()
      ssq_img = flex.pow2(img)
      sum_wavelength = wavelength
      sum_time = (0, 0)
      time_base = time_tuple

    else:
      sel = (img > max_img).as_1d()
      max_img.set_selected(sel, img.select(sel))

      sum_distance += distance
      sum_img += img
      ssq_img += flex.pow2(img)
      sum_wavelength += wavelength
      sum_time = (sum_time[0] + (time_tuple[0] - time_base[0]),
                  sum_time[1] + (time_tuple[1] - time_base[1]))

    nmemb += 1

  # Early exit if no statistics were accumulated.
  if command_line.options.verbose:
    sys.stderr.write("Processed %d images (%d failed)\n" % (nmemb, nfail))
  if nmemb == 0:
    return 0

  # Calculate averages for measures where other statistics do not make
  # sense.  Note that avg_img is required for stddev_img.
  avg_img = sum_img.as_double() / nmemb
  avg_distance = sum_distance / nmemb
  avg_timestamp = cspad_tbx.evt_timestamp(
    (time_base[0] + int(round(sum_time[0] / nmemb)),
     time_base[1] + int(round(sum_time[1] / nmemb))))
  avg_wavelength = sum_wavelength / nmemb

  # Output the average image, maximum projection image, and standard
  # deviation image, if requested.
  if command_line.options.avg_path is not None:
    avg_img.resize(flex.grid(size[0], size[1]))
    d = cspad_tbx.dpack(
      active_areas=active_areas,
      address=detector_address,
      beam_center_x=beam_center[0],
      beam_center_y=beam_center[1],
      data=avg_img,
      distance=avg_distance,
      pixel_size=pixel_size,
      saturated_value=saturated_value,
      timestamp=avg_timestamp,
      wavelength=avg_wavelength)
    easy_pickle.dump(command_line.options.avg_path, d)

  if command_line.options.max_path is not None:
    max_img.resize(flex.grid(size[0], size[1]))
    d = cspad_tbx.dpack(
      active_areas=active_areas,
      address=detector_address,
      beam_center_x=beam_center[0],
      beam_center_y=beam_center[1],
      data=max_img,
      distance=avg_distance,
      pixel_size=pixel_size,
      saturated_value=saturated_value,
      timestamp=avg_timestamp,
      wavelength=avg_wavelength)
    easy_pickle.dump(command_line.options.max_path, d)

  if command_line.options.stddev_path is not None:
    stddev_img = ssq_img.as_double() - sum_img.as_double() * avg_img

    # Accumulating floating-point numbers introduces errors, which may
    # cause negative variances.  Since a two-pass approach is
    # unacceptable, the standard deviation is clamped at zero.
    stddev_img.set_selected(stddev_img < 0, 0)
    if nmemb == 1:
      stddev_img = flex.sqrt(stddev_img)
    else:
      stddev_img = flex.sqrt(stddev_img / (nmemb - 1))

    stddev_img.resize(flex.grid(size[0], size[1]))
    d = cspad_tbx.dpack(
      active_areas=active_areas,
      address=detector_address,
      beam_center_x=beam_center[0],
      beam_center_y=beam_center[1],
      data=stddev_img,
      distance=avg_distance,
      pixel_size=pixel_size,
      saturated_value=saturated_value,
      timestamp=avg_timestamp,
      wavelength=avg_wavelength)
    easy_pickle.dump(command_line.options.stddev_path, d)

  return 0


if __name__ == '__main__':
  sys.exit(run())
