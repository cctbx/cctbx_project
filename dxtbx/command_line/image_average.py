# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME dxtbx.image_average
# LIBTBX_SET_DISPATCHER_NAME cxi.image_average
#

from __future__ import absolute_import, division

import sys, dxtbx
from scitbx.array_family import flex

"""
Average single-panel images of any format. Handles many individual images or single container files.
"""

def splitit(l, n):
  """ Utitlity funciton to evenly split a list. Handles edge cases.
  There is probably a 1-liner list comprehension to do this, but it would be super gnarly.
  @param l list to split (not a generator)
  @param n number of chunks to split the list into
  @return list of n lists
  """
  # if list is shorter then n, split into a list of lists, each with one entry
  if len(l) < n:
    n = len(l)
  s = len(l) // (n) # each chunk will either be of size s or size s+1
  m = len(l) % (s*n) # remainder after n chunks of size s = how many chunks of size s+1 are needed
  r = [] # result
  p = 0 # pointer
  for i in xrange(n):
    if i < m:
      r.append(l[p:p+s+1])
      p += s+1
    else:
      r.append(l[p:p+s])
      p += s
  return r

class image_worker(object):
  """ Class to compute running sums while reading image data """
  # Deriving class should implement __init__, load and read

  def __call__(self, subset):
    """ Worker function for multiprocessing """
    nfail = 0
    nmemb = 0

    self.load()

    for item in subset:
      try:
        #XXX This code assumes a monolithic detector!
        beam_center, detector_address, distance, img, pixel_size, saturated_value, size, wavelength, active_areas = \
          self.read(item)
      except Exception, e:
        print str(e)
        nfail += 1
        continue

      # The sum-of-squares image is accumulated using long integers, as
      # this delays the point where overflow occurs.  But really, this
      # is just a band-aid...
      if nmemb == 0:
        max_img = img.deep_copy()
        sum_distance = distance
        sum_img = img.deep_copy()
        ssq_img = flex.pow2(img)
        sum_wavelength = wavelength

      else:
        sel = (img > max_img).as_1d()
        max_img.set_selected(sel, img.select(sel))

        sum_distance += distance
        sum_img += img
        ssq_img += flex.pow2(img)
        sum_wavelength += wavelength

      nmemb += 1
    return nfail, nmemb, max_img, sum_distance, sum_img, ssq_img, sum_wavelength, size, active_areas, detector_address, beam_center, pixel_size, saturated_value

class multi_image_worker(image_worker):
  """ Class for reading container files """
  def __init__(self, command_line, path, imageset):
    self.path = path
    self.command_line = command_line

    self.imageset = imageset

  def load(self):
    """ Called by seperate process during multiprocessing """

    if self.command_line.options.nproc > 1:
      # Need to re-open the file if HDF5 as HDF5 file handles can't be pickled during multiprocessing
      self.imageset.reader().nullify_format_instance()

  def read(self, n):
    """ Read image at postion n"""
    if self.command_line.options.verbose:
      print "Processing %s: %d" % (self.path, n)

    beam = self.imageset.get_beam(n)
    assert len(self.imageset.get_detector(n)) == 1
    detector = self.imageset.get_detector(n)[0]

    beam_center = detector.get_beam_centre(beam.get_s0())
    detector_address = type(self.imageset.reader()).__name__
    distance = detector.get_distance()
    img = self.imageset[n][0].as_1d().as_double()
    pixel_size = 0.5 * sum(detector.get_pixel_size())
    saturated_value = int(round(detector.get_trusted_range()[1]))
    size = detector.get_image_size()
    wavelength = beam.get_wavelength()

    active_areas = flex.int((0, 0, size[0], size[1]))

    return beam_center, detector_address, distance, img, pixel_size, saturated_value, size, wavelength, active_areas

class single_image_worker(image_worker):
  """ Class for averaging single images from individual files """
  def __init__(self, command_line):
    self.command_line = command_line

  def load(self):
    pass # no additional test needed when loading individual files

  def read(self, path):
    if self.command_line.options.verbose:
      print "Processing %s" % path

    from dxtbx.format.Registry import Registry
    from dxtbx.format.FormatMultiImage import FormatMultiImage
    format_class = Registry.find(path)
    assert not issubclass(format_class, FormatMultiImage), "Average container files seperately"
    img_instance = format_class(path)

    beam = img_instance.get_beam()
    assert len(img_instance.get_detector()) == 1
    detector = img_instance.get_detector()[0]

    beam_center = detector.get_beam_centre(beam.get_s0())
    detector_address = format_class.__name__
    distance = detector.get_distance()
    img = img_instance.get_raw_data().as_1d().as_double()
    pixel_size = 0.5 * sum(detector.get_pixel_size())
    saturated_value = int(round(detector.get_trusted_range()[1]))
    size = detector.get_image_size()
    wavelength = beam.get_wavelength()

    active_areas = flex.int((0, 0, size[0], size[1]))
    return beam_center, detector_address, distance, img, pixel_size, saturated_value, size, wavelength, active_areas

def run(argv=None):
  """Compute mean, standard deviation, and maximum projection images
  from a set of images given on the command line.

  @param argv Command line argument list
  @return     @c 0 on successful termination, @c 1 on error, and @c 2
              for command line syntax errors
  """
  import libtbx.load_env

  from libtbx import easy_pickle, option_parser
  from xfel.cxi.cspad_ana import cspad_tbx

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
                  .option(None, "--nproc", "-n",
                          type="int",
                          default=1,
                          dest="nproc",
                          help="Number of processors")
                  .option(None, "--num-images-max", "-N",
                          type="int",
                          default=None,
                          dest="num_images_max",
                          help="Maximum number of frames to average")
                  ).process(args=argv[1:])

  # Note that it is not an error to omit the output paths, because
  # certain statistics could still be printed, e.g. with the verbose
  # option.
  paths = command_line.args
  if len(paths) == 0:
    command_line.parser.print_usage(file=sys.stderr)
    return 2

  if len(paths) == 1:
    # test if the iamge is a multi-image
    from dxtbx.datablock import DataBlockFactory
    datablocks = DataBlockFactory.from_filenames([paths[0]])
    assert len(datablocks) == 1
    datablock = datablocks[0]
    imagesets = datablock.extract_imagesets()
    assert len(imagesets) == 1
    imageset = imagesets[0]
    if not imageset.reader().is_single_file_reader():
      from libtbx.utils import Usage
      raise Usage("Supply more than one image")

    worker = multi_image_worker(command_line, paths[0], imageset)
    if command_line.options.num_images_max is not None and command_line.options.num_images_max < len(imageset):
      iterable = range(command_line.options.num_images_max)
    else:
      iterable = range(len(imageset))
  else:
    # Multiple images provided
    worker = single_image_worker(command_line)
    if command_line.options.num_images_max is not None and command_line.options.num_images_max < len(paths):
      iterable = paths[:command_line.options.num_images_max]
    else:
      iterable = paths
  if command_line.options.nproc > 1:
    iterable = splitit(iterable, command_line.options.nproc)

  from libtbx import easy_mp
  if command_line.options.nproc == 1:
    results = [worker(iterable)]
  else:
    results = easy_mp.parallel_map(func=worker,
                                   iterable=iterable,
                                   processes=command_line.options.nproc)

  nfail = 0
  nmemb = 0
  for i, (r_nfail, r_nmemb, r_max_img, r_sum_distance, r_sum_img, r_ssq_img, r_sum_wavelength, size, active_areas, detector_address, beam_center, pixel_size, saturated_value) in enumerate(results):
    nfail += r_nfail
    nmemb += r_nmemb
    if i == 0:
      max_img = r_max_img
      sum_distance = r_sum_distance
      sum_img = r_sum_img
      ssq_img = r_ssq_img
      sum_wavelength = r_sum_wavelength
    else:
      sel = (r_max_img > max_img).as_1d()
      max_img.set_selected(sel, r_max_img.select(sel))

      sum_distance += r_sum_distance
      sum_img += r_sum_img
      ssq_img += r_ssq_img
      sum_wavelength += r_sum_wavelength

  # Early exit if no statistics were accumulated.
  if command_line.options.verbose:
    sys.stderr.write("Processed %d images (%d failed)\n" % (nmemb, nfail))
  if nmemb == 0:
    return 0

  # Calculate averages for measures where other statistics do not make
  # sense.  Note that avg_img is required for stddev_img.
  avg_img = sum_img.as_double() / nmemb
  avg_distance = sum_distance / nmemb
  avg_wavelength = sum_wavelength / nmemb

  # Output the average image, maximum projection image, and standard
  # deviation image, if requested.
  if command_line.options.avg_path is not None:
    avg_img.resize(flex.grid(size[1], size[0]))
    d = cspad_tbx.dpack(
      active_areas=active_areas,
      address=detector_address,
      beam_center_x=beam_center[0],
      beam_center_y=beam_center[1],
      data=avg_img,
      distance=avg_distance,
      pixel_size=pixel_size,
      saturated_value=saturated_value,
      wavelength=avg_wavelength)
    easy_pickle.dump(command_line.options.avg_path, d)

  if command_line.options.max_path is not None:
    max_img.resize(flex.grid(size[1], size[0]))
    d = cspad_tbx.dpack(
      active_areas=active_areas,
      address=detector_address,
      beam_center_x=beam_center[0],
      beam_center_y=beam_center[1],
      data=max_img,
      distance=avg_distance,
      pixel_size=pixel_size,
      saturated_value=saturated_value,
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

    stddev_img.resize(flex.grid(size[1], size[0]))
    d = cspad_tbx.dpack(
      active_areas=active_areas,
      address=detector_address,
      beam_center_x=beam_center[0],
      beam_center_y=beam_center[1],
      data=stddev_img,
      distance=avg_distance,
      pixel_size=pixel_size,
      saturated_value=saturated_value,
      wavelength=avg_wavelength)
    easy_pickle.dump(command_line.options.stddev_path, d)

  return 0


if __name__ == '__main__':
  sys.exit(run())
