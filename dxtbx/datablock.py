#!/usr/bin/env python
#
# datablock.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class ImageRecord(object):
  ''' Record for storing image metadata. '''
  def __init__(self, mtime=None, beam=None, detector=None,
               goniometer=None, scan=None, template=None, index=None):
    self.mtime = mtime
    self.beam = beam
    self.detector = detector
    self.goniometer = goniometer
    self.scan = scan
    self.template = template
    self.index = index

  def clone(self, rhs):
    self.mtime = rhs.mtime
    self.beam = rhs.beam
    self.detector = rhs.detector
    self.goniometer = rhs.goniometer
    self.scan = rhs.scan
    self.template = rhs.template
    self.index = rhs.index

  def __eq__(self, rhs):
    return (self.mtime == rhs.mtime and
            self.beam == rhs.beam and
            self.detector == rhs.detector and
            self.goniometer == rhs.goniometer and
            self.scan == rhs.scan and
            self.template == rhs.template and
            self.index == rhs.index)

  def __ne__(self, rhs):
    return not self.__eq__(rhs)


class DataBlock(object):
  ''' High level container for blocks of sweeps and imagesets. '''

  def __init__(self, filenames=None, format_class=None):
    ''' Instantiate from filenames or a format class.

    If no format is given then the format is found from the first
    filename. All files must use the same format, otherwise an
    exception will be raised.

    Params:
        filenames The filenames to use
        format_class The format class to use

    '''
    from dxtbx.format.Registry import Registry
    from collections import OrderedDict

    # Try to get a format class
    if format_class == None:
      if filenames is not None and len(filenames) > 0:
        format_class = Registry.find(filenames[0])
    self._format_class = format_class

    # Load any images
    self._images = OrderedDict()
    if filenames is not None:
      for f in filenames:
        self.append(f)

  def format_class(self):
    ''' Return the format class. '''
    return self._format_class

  def filenames(self):
    ''' Return the list of filenames. '''
    return list(self._images.iterkeys())

  def metadata(self):
    ''' Return the list of filenames and meta data. '''
    return list(self._images.iteritems())

  def extract_all(self):
    ''' Extract all the images as an image set. '''
    from dxtbx.imageset2 import ImageSetFactory
    return ImageSetFactory.make_imageset(
        self._images.keys(), self._format_class)

  def extract_stills(self):
    ''' Extract all the stills as an image set. '''
    from dxtbx.imageset2 import ImageSetFactory
    stills = [k for k, r in self._images.iteritems() if r.template == None]
    return ImageSetFactory.make_imageset(stills, self._format_class)

  def extract_sweeps(self):
    ''' Extract all the sweeps from the block. '''
    from dxtbx.imageset2 import ImageSetFactory
    from itertools import groupby
    from dxtbx.format.FormatStill import FormatStill

    # Check we're not using stills
    assert(not issubclass(self._format_class, FormatStill))
    sweeps = []

    # Get consecutive groups of scans (which should correspond to sweeps)
    groups = groupby(self._images.itervalues(), key=lambda r: r.scan)
    for scan, records in groups:
      records = list(records)
      if scan is not None and len(records) > 0:
        templates = [r.template for r in records]
        indices = [r.index for r in records]
        if templates[0] is not None:
          assert(len(indices) > 0)
          assert(templates.count(templates[0]) == len(templates))
          assert(all(j == i+1 for i, j in zip(indices[:-1], indices[1:])))
          sweep = ImageSetFactory.make_sweep(
              templates[0],
              indices,
              self._format_class,
              records[0].beam,
              records[0].detector,
              records[0].goniometer,
              records[0].scan)
          sweeps.append(sweep)

    # Return the list of sweeps
    return sweeps

  def append(self, filename):
    ''' Add another image to the block. The image must use the same
    format class, otherwise an exception will be raised. '''
    from os.path import abspath

    # Check the image is not already here and can be understood
    filename = abspath(filename)
    if filename in self._images:
      raise RuntimeError('%s already in data block' % filename)
    if not self.understand(filename):
      raise RuntimeError('cannot understand file: %s' % filename)

    # Get the meta data
    if len(self._images) > 0:
      last = self._images[next(reversed(self._images))]
    else:
      last = None
    self._images[filename] = self._create_record(filename, last)


  def understand(self, filename):
    ''' Check if the data block format understands the given file. This
    function checks the method resolution order for the format class and
    calls the understand methods of each parent down to the bottom level.
    Just calling the format class understand method directly can result
    in problems. This is really a workaround for a bug in the way that
    the format understand method works. '''
    from dxtbx.format.Format import Format
    mro = self._format_class.mro()[::-1]
    if len(mro) <= 2 or mro[0] != object or mro[1] != Format:
      return False
    for m in mro[2:]:
      if m.understand(filename) == False:
        return False
    return True

  def _create_record(self, filename, last=None):
    ''' Get any information about the image we can and create a record. '''
    from dxtbx.sweep_filenames import template_regex
    from os.path import getmtime

    # Read the image
    fmt = self._format_class(filename)

    # Get the meta data from the format
    try: b = fmt.get_beam()
    except Exception: b = None
    try: d = fmt.get_detector()
    except Exception: d = None
    try: g = fmt.get_goniometer()
    except Exception: g = None
    try: s = fmt.get_scan()
    except Exception: s = None

    # Get the template and index if possible
    if s is not None and abs(s.get_oscillation()[1]) > 0.0:
      template, index = template_regex(filename)
    else:
      template, index = None, None

    # Test against last item in list
    if last is not None:
      if last.beam == b:
        b = last.beam
      if last.detector == d:
        d = last.detector
      if last.goniometer == g:
        g = last.goniometer
      try:
        if last.template == template and last.index + 1 == index:
          last.scan += s
          s = last.scan
      except Exception:
        pass

    # Create the record and return
    return ImageRecord(
        mtime=getmtime(filename),
        beam=b, detector=d, goniometer=g, scan=s,
        template=template, index=index)

  def __len__(self):
    ''' The number of images. '''
    return len(self._images)

  def __eq__(self, rhs):
    ''' Check if two blocks are the same. '''
    return (self._format_class == rhs._format_class
            and self._images == rhs._images)

  def __ne__(self, rhs):
    ''' Check if two blocks are not equal. '''
    return not self.__eq__(rhs)

  def __getstate__(self):
    ''' Return the dict for pickling. '''
    return self.__dict__

  def __setstate__(self, state):
    ''' Update the dict from pickling. On reload, update any records for
    images that have changed since the class was pickled. '''
    from os.path import getmtime
    self.__dict__.update(state)
    last = None
    for filename, image in self._images.iteritems():
      if getmtime(filename) > image.mtime:
        image.clone(self._create_record(filename, last))
      last = image


class DataBlockFactory(object):
  ''' Class for creating DataBlock instances'''

  @staticmethod
  def from_filenames(filenames, verbose=False):
    ''' Create a list of data blocks from a list of filenames. '''
    return DataBlockFactory.create_list(filenames, verbose)

  @staticmethod
  def create_list(filenames, verbose=False):
    ''' Create a list of data blocks from a list of filenames. '''
    datablock_list = []
    for f in filenames:
      if verbose: print 'Loading file: %s' % f
      try:
        datablock_list[-1].append(f)
      except Exception:
        datablock_list.append(DataBlockFactory.create_single([f]))
    return datablock_list

  @staticmethod
  def create_single(filenames, verbose=False):
    ''' Create a single data blocks from a list of filenames. '''

    # Ensure we have a list of images
    if len(filenames) < 1:
      raise RuntimeError('Need at least 1 image to create a data block')

    # Create the datablock
    return DataBlock(filenames)


if __name__ == '__main__':

  import sys

  # Get the data blocks from the input files
  # We've set verbose to print out files as they're tested.
  datablock_list = DataBlockFactory.from_filenames(sys.argv[1:], verbose=True)

  # Loop through the data blocks
  for i, datablock in enumerate(datablock_list):

    # Extract any sweeps
    sweeps = datablock.extract_sweeps()
    stills = datablock.extract_stills()

    print "-" * 80
    print "DataBlock %d" % i
    print "  format: %s" % str(datablock.format_class())
    print "  num images: %d" % len(datablock)
    print "  num sweeps: %d" % len(sweeps)
    print "  num stills: %d" % len(stills)

    # Loop through all the sweeps
    for j, sweep in enumerate(sweeps):
      print ""
      print "Sweep %d" % j
      print "  length %d" % len(sweep)
      print sweep.get_beam()
      print sweep.get_goniometer()
      print sweep.get_detector()
      print sweep.get_scan()
