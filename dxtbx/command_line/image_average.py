# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME dxtbx.image_average
# LIBTBX_SET_DISPATCHER_NAME cxi.image_average
#

from __future__ import absolute_import, division, print_function

import sys, copy

import dxtbx
from scitbx.array_family import flex
from dxtbx.format.cbf_writer import FullCBFWriter

"""
Average images of any dxtbx-supported format. Handles many individual images or single container files.
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
        img, distance, wavelength = self.read(item)
      except Exception as e:
        print(str(e))
        nfail += 1
        continue

      assert isinstance(img, tuple)

      # The sum-of-squares image is accumulated using long integers, as
      # this delays the point where overflow occurs.  But really, this
      # is just a band-aid...
      if nmemb == 0:
        max_img = list(copy.deepcopy(img))
        sum_distance = distance
        sum_img = list(copy.deepcopy(img))
        ssq_img = [flex.pow2(p) for p in img]
        sum_wavelength = wavelength

      else:
        for p in xrange(len(img)):
          sel = (img[p] > max_img[p]).as_1d()
          max_img[p].set_selected(sel, img[p].select(sel))

          sum_img[p] += img[p]
          ssq_img[p] += flex.pow2(img[p])
        sum_distance += distance
        sum_wavelength += wavelength

      nmemb += 1
    return nfail, nmemb, max_img, sum_distance, sum_img, ssq_img, sum_wavelength

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
      print("Processing %s: %d" % (self.path, n))

    beam = self.imageset.get_beam(n)
    detector = self.imageset.get_detector(n)

    image_data = self.imageset[n]
    if not isinstance(image_data, tuple):
      image_data = (image_data,)
    img = tuple([image_data[i].as_1d().as_double() for i in xrange(len(detector))])
    wavelength = beam.get_wavelength()

    return img, detector.hierarchy().get_distance(), wavelength

class single_image_worker(image_worker):
  """ Class for averaging single images from individual files """
  def __init__(self, command_line):
    self.command_line = command_line

  def load(self):
    pass # no additional test needed when loading individual files

  def read(self, path):
    if self.command_line.options.verbose:
      print("Processing %s" % path)

    from dxtbx.format.Registry import Registry
    from dxtbx.format.FormatMultiImage import FormatMultiImage
    format_class = Registry.find(path)
    assert not issubclass(format_class, FormatMultiImage), "Average container files seperately"
    img_instance = format_class(path)

    beam = img_instance.get_beam()
    detector = img_instance.get_detector()
    image_data = img_instance.get_raw_data()
    if not isinstance(image_data, tuple):
      image_data = (image_data,)
    img = tuple([image_data[i].as_1d().as_double() for i in xrange(len(detector))])
    wavelength = beam.get_wavelength()

    return img, detector.hierarchy().get_distance(), wavelength

def run(argv=None):
  """Compute mean, standard deviation, and maximum projection images
  from a set of images given on the command line.

  @param argv Command line argument list
  @return     @c 0 on successful termination, @c 1 on error, and @c 2
              for command line syntax errors
  """
  import libtbx.load_env

  from libtbx import option_parser

  if argv is None:
    argv = sys.argv
  command_line = (option_parser.option_parser(
    usage="%s [-v] [-a PATH] [-m PATH] [-s PATH] " \
    "image1 image2 [image3 ...]" % libtbx.env.dispatcher_name)
                  .option(None, "--average-path", "-a",
                          type="string",
                          default="avg.cbf",
                          dest="avg_path",
                          metavar="PATH",
                          help="Write average image to PATH")
                  .option(None, "--maximum-path", "-m",
                          type="string",
                          default="max.cbf",
                          dest="max_path",
                          metavar="PATH",
                          help="Write maximum projection image to PATH")
                  .option(None, "--stddev-path", "-s",
                          type="string",
                          default="stddev.cbf",
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
                  .option(None, "--skip-images", "-S",
                          type="int",
                          default=None,
                          dest="skip_images",
                          help="Number of images to skip at the start of the dataset")
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
    iterable = range(len(imageset))
  else:
    # Multiple images provided
    worker = single_image_worker(command_line)
    iterable = paths

  if command_line.options.skip_images is not None:
    if command_line.options.skip_images >= len(iterable):
      from libtbx.utils import Usage
      raise Usage("Skipping all the images")
    iterable = iterable[command_line.options.skip_images:]
  if command_line.options.num_images_max is not None and command_line.options.num_images_max < iterable:
    iterable = iterable[:command_line.options.num_images_max]
  assert len(iterable) >= 2, "Need more than one image to average"

  if len(paths) > 1:
    from dxtbx.datablock import DataBlockFactory
    datablocks = DataBlockFactory.from_filenames([iterable[0]])
    assert len(datablocks) == 1
    datablock = datablocks[0]
    imagesets = datablock.extract_imagesets()
    assert len(imagesets) == 1
    imageset = imagesets[0]

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
  for i, (r_nfail, r_nmemb, r_max_img, r_sum_distance, r_sum_img, r_ssq_img, r_sum_wavelength) in enumerate(results):
    nfail += r_nfail
    nmemb += r_nmemb
    if i == 0:
      max_img = r_max_img
      sum_distance = r_sum_distance
      sum_img = r_sum_img
      ssq_img = r_ssq_img
      sum_wavelength = r_sum_wavelength
    else:
      for p in xrange(len(sum_img)):
        sel = (r_max_img[p] > max_img[p]).as_1d()
        max_img[p].set_selected(sel, r_max_img[p].select(sel))

        sum_img[p] += r_sum_img[p]
        ssq_img[p] += r_ssq_img[p]

      sum_distance += r_sum_distance
      sum_wavelength += r_sum_wavelength

  # Early exit if no statistics were accumulated.
  if command_line.options.verbose:
    sys.stderr.write("Processed %d images (%d failed)\n" % (nmemb, nfail))
  if nmemb == 0:
    return 0

  # Calculate averages for measures where other statistics do not make
  # sense.  Note that avg_img is required for stddev_img.
  avg_img = tuple([sum_img[p].as_double() / nmemb for p in xrange(len(sum_img))])
  avg_distance = sum_distance / nmemb
  avg_wavelength = sum_wavelength / nmemb

  detector = imageset.get_detector()
  h = detector.hierarchy()
  origin = h.get_local_origin()
  h.set_local_frame(h.get_local_fast_axis(), h.get_local_slow_axis(), (origin[0], origin[1], -avg_distance))
  imageset.get_beam().set_wavelength(avg_wavelength)

  # Output the average image, maximum projection image, and standard
  # deviation image, if requested.
  if command_line.options.avg_path is not None:
    for p in xrange(len(detector)):
      fast, slow = detector[p].get_image_size()
      avg_img[p].resize(flex.grid(slow, fast))

    writer = FullCBFWriter(imageset = imageset)
    cbf = writer.get_cbf_handle(header_only = True)
    writer.add_data_to_cbf(cbf, data = avg_img)
    writer.write_cbf(command_line.options.avg_path, cbf = cbf)

  if command_line.options.max_path is not None:
    for p in xrange(len(detector)):
      fast, slow = detector[p].get_image_size()
      max_img[p].resize(flex.grid(slow, fast))
    max_img = tuple(max_img)

    writer = FullCBFWriter(imageset = imageset)
    cbf = writer.get_cbf_handle(header_only = True)
    writer.add_data_to_cbf(cbf, data = max_img)
    writer.write_cbf(command_line.options.max_path, cbf = cbf)

  if command_line.options.stddev_path is not None:
    stddev_img = []
    for p in xrange(len(detector)):
      stddev_img.append(ssq_img[p].as_double() - sum_img[p].as_double() * avg_img[p])

      # Accumulating floating-point numbers introduces errors, which may
      # cause negative variances.  Since a two-pass approach is
      # unacceptable, the standard deviation is clamped at zero.
      stddev_img[p].set_selected(stddev_img[p] < 0, 0)
      if nmemb == 1:
        stddev_img[p] = flex.sqrt(stddev_img[p])
      else:
        stddev_img[p] = flex.sqrt(stddev_img[p] / (nmemb - 1))

      fast, slow = detector[p].get_image_size()
      stddev_img[p].resize(flex.grid(slow, fast))
    stddev_img = tuple(stddev_img)

    writer = FullCBFWriter(imageset = imageset)
    cbf = writer.get_cbf_handle(header_only = True)
    writer.add_data_to_cbf(cbf, data = stddev_img)
    writer.write_cbf(command_line.options.stddev_path, cbf = cbf)

  return 0

if __name__ == '__main__':
  sys.exit(run())
