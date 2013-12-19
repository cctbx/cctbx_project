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
  def __init__(self, beam=None, detector=None, goniometer=None, scan=None,
               template=None, index=None):
    self.beam = beam
    self.detector = detector
    self.goniometer = goniometer
    self.scan = scan
    self.template = template
    self.index = index

  def clone(self, rhs):
    self.beam = rhs.beam
    self.detector = rhs.detector
    self.goniometer = rhs.goniometer
    self.scan = rhs.scan
    self.template = rhs.template
    self.index = rhs.index

  def __eq__(self, rhs):
    return (self.beam == rhs.beam and
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
    filenames, records = zip(*[(k, r) for k, r in self._images.iteritems() if
                               r.template == None])
    if len(filenames) == 0:
      return None
    imageset = ImageSetFactory.make_imageset(filenames, self._format_class)
    for i, r in enumerate(records):
      imageset.set_beam(r.beam, i)
      imageset.set_detector(r.detector, i)
      imageset.set_goniometer(r.goniometer, i)
      imageset.set_scan(r.scan, i)
    return imageset

  def extract_sweeps(self):
    ''' Extract all the sweeps from the block. '''
    from dxtbx.imageset2 import ImageSetFactory
    from dxtbx.format.FormatStill import FormatStill

    # Check we're not using stills
    sweeps = []
    if issubclass(self._format_class, FormatStill):
      return sweeps

    # Get consecutive groups of scans (which should correspond to sweeps)
    for group in self.iter_sweep_groups():
      filenames, records = zip(*group)
      indices = [r.index for r in records]
      sweep = ImageSetFactory.make_sweep(
        records[0].template, indices, self._format_class,
        records[0].beam, records[0].detector,
        records[0].goniometer, records[0].scan)
      sweeps.append(sweep)

    # Return the list of sweeps
    return sweeps

  def append(self, filename):
    ''' Add another image to the block. The image must use the same
    format class, otherwise an exception will be raised. '''
    from os.path import abspath

    # Recursively check whether any of the children understand
    # image_file, in which case they are preferred over the parent
    # format.
    def better_understood_by_child(format_class, image_file):
      for child in format_class._children:
        if child.understand(image_file):
          return True
      return False

    # Check the image is not already here and can be understood and is not
    # better understood by a child format class
    filename = abspath(filename)
    if filename in self._images:
      raise RuntimeError('%s already in data block' % filename)
    if not self.understand(filename):
      raise RuntimeError('cannot understand file: %s' % filename)
    if better_understood_by_child(self.format_class(), filename):
      raise RuntimeError('%s is better understood by a child format' % filename)

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
    the format understand method works. Furthermore, the FormatStill
    class does not have an understand method so we have to check that
    the "understand" method is present in the class dictionary before
    we actually do the call. '''
    from dxtbx.format.Format import Format
    mro = self._format_class.mro()[::-1]
    if len(mro) <= 2 or mro[0] != object or mro[1] != Format:
      return False
    for m in mro[2:]:
      if "understand" in m.__dict__ and m.understand(filename) == False:
        return False
    return True

  def _create_record(self, filename, last=None):
    ''' Get any information about the image we can and create a record. '''
    from dxtbx.sweep_filenames import template_regex

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

      # If all models are the same the last models then try to
      # see if we can create a sweep from the scan and template
      if b is last.beam and d is last.detector and g is last.goniometer:
        try:
          if last.template == template and last.index + 1 == index:
            last.scan += s
            s = last.scan
        except Exception:
          pass

    # Create the record and return
    return ImageRecord(
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

  def iter_groups(self):
    ''' Iterate between groups of images. '''
    from itertools import groupby
    groups = groupby(self._images.iteritems(), key=lambda r: r[1].scan)
    for scan, group in groups:
      group = list(group)
      filenames, records = zip(*group)
      assert(len(records) > 0)
      if scan is not None and records[0].template is not None:
        templates = [r.template for r in records]
        indices = [r.index for r in records]
        assert(len(indices) > 0)
        assert(templates.count(templates[0]) == len(templates))
        assert(all(j == i+1 for i, j in zip(indices[:-1], indices[1:])))
        assert(scan.get_image_range() == (indices[0], indices[-1]))
        yield (group, True)
      else:
        yield (group, False)

  def iter_sweep_groups(self):
    ''' Iterate over sweep groups. '''
    for group, sflag in self.iter_groups():
      if sflag == True:
        yield group

  def iter_still_groups(self):
    ''' Iterate over still groups. '''
    for group, sflag in self.iter_groups():
      if sflag == False:
        yield group

  def to_dict(self):
    ''' Convert the datablock to a dictionary '''
    from collections import OrderedDict

    # Get the lists of models
    b, d, g, s = zip(*[(i.beam, i.detector, i.goniometer, i.scan)
      for i in self._images.itervalues()])

    # Get the set of unique models
    b = OrderedDict([(bb, None) for bb in b if bb is not None])
    d = OrderedDict([(dd, None) for dd in d if dd is not None])
    g = OrderedDict([(gg, None) for gg in g if gg is not None])
    s = OrderedDict([(ss, None) for ss in s if ss is not None])

    # Create the data block dictionary
    result = OrderedDict()
    result['__id__'] = 'DataBlock'
    result['imagesets'] = []

    # Convert the sweeps and imagesets to dictionaries
    for group, sflag in self.iter_groups():
      if sflag:
        filename, record = group[0]
        result['imagesets'].append(OrderedDict([
          ('__id__', 'ImageSweep'),
          ('template',   record.template),
          ('beam',       b.keys().index(record.beam)),
          ('detector',   d.keys().index(record.detector)),
          ('goniometer', g.keys().index(record.goniometer)),
          ('scan',       s.keys().index(record.scan))
        ]))
      else:
        image_list = []
        for filename, record in group:
          image_dict = OrderedDict()
          image_dict['filename'] = filename
          if record.beam:
            image_dict['beam'] = b.keys().index(record.beam)
          if record.detector:
            image_dict['detector'] = d.keys().index(record.detector)
          if record.goniometer:
            image_dict['goniometer'] = g.keys().index(record.goniometer)
          if record.scan:
            image_dict['scan'] = s.keys().index(record.scan)
          image_list.append(image_dict)
        result['imagesets'].append(OrderedDict([
          ('__id__', 'ImageSet'), ('images', image_list)]))

    # Add the models to the dictionary
    result['beam'] = [bb.to_dict() for bb in b]
    result['detector'] = [dd.to_dict() for dd in d]
    result['goniometer'] = [gg.to_dict() for gg in g]
    result['scan'] = [ss.to_dict() for ss in s]

    # Return the data block as a dictionary
    return result


class DataBlockFactory(object):
  ''' Class for creating DataBlock instances'''

  @staticmethod
  def from_args(args, verbose=False, unhandled=None):
    ''' Try to load datablocks from any recognised format. '''
    if unhandled is None:
      unhandled = []
    unhandled1 = []
    unhandled2 = []

    # First try as image files
    datablocks = DataBlockFactory.from_filenames(args, verbose, unhandled1)

    # Now try as JSON files
    if len(unhandled1) > 0:
      for filename in unhandled1:
        try:
          datablocks.extend(DataBlockFactory.from_json_file(filename))
          if verbose: print 'Loaded datablock(s) from %s' % filename
        except Exception:
          unhandled2.append(filename)

    # Now try as pickle files
    if len(unhandled2) > 0:
      for filename in unhandled2:
        try:
          datablocks.extend(DataBlockFactory.from_pickle_file(filename))
          if verbose: print 'Loaded datablock(s) from %s' % filename
        except Exception:
          unhandled.append(filename)

    # Return the data blocks
    return datablocks

  @staticmethod
  def from_filenames(filenames, verbose=False, unhandled=None):
    ''' Create a list of data blocks from a list of filenames. '''
    return DataBlockFactory.create_list(filenames, verbose, unhandled)

  @staticmethod
  def create_list(filenames, verbose=False, unhandled=None):
    ''' Create a list of data blocks from a list of filenames. '''
    datablock_list = []
    for f in filenames:
      try:
        datablock_list[-1].append(f)
      except Exception:
        try:
          datablock_list.append(DataBlockFactory.create_single([f]))
        except Exception:
          if unhandled is not None:
            unhandled.append(f)
          continue
        if verbose: print 'Starting datablock %d' % len(datablock_list)
      if verbose: print 'Loaded file: %s' % f
    return datablock_list

  @staticmethod
  def create_single(filenames, verbose=False):
    ''' Create a single data blocks from a list of filenames. '''

    # Ensure we have a list of images
    if len(filenames) < 1:
      raise RuntimeError('Need at least 1 image to create a data block')

    # Create the datablock
    return DataBlock(filenames)

  @staticmethod
  def from_dict(d):
    ''' Create the datablock from a dictionary. '''
    from collections import OrderedDict
    from dxtbx.format.Registry import Registry
    from dxtbx.model import Beam, Detector, Goniometer, Scan

    # If we have a list, extract for each dictionary in the list
    if isinstance(d, list):
      return [DataBlockFactory.from_dict(dd) for dd in d]
    elif not isinstance(d, dict):
      raise RuntimeError('unknown datablock dictionary type')
    assert(d['__id__'] == 'DataBlock')

    # Get the list of models
    blist = d.get('beam', [])
    dlist = d.get('detector', [])
    glist = d.get('goniometer', [])
    slist = d.get('scan', [])

    # Create the dictionary of records
    records = OrderedDict()
    for imageset in d['imagesets']:

      # Add a sweep
      if imageset['__id__'] == 'ImageSweep':
        beam = Beam.from_dict(blist[imageset['beam']])
        detector = Detector.from_dict(dlist[imageset['detector']])
        goniometer = Goniometer.from_dict(glist[imageset['goniometer']])
        scan = Scan.from_dict(slist[imageset['scan']])
        template = imageset['template']
        pfx = template.split('#')[0]
        sfx = template.split('#')[-1]
        template_format = '%s%%0%dd%s' % (pfx, template.count('#'), sfx)
        i0, i1 = scan.get_image_range()
        for i in range(i0, i1 + 1):
          record = ImageRecord()
          record.beam = beam
          record.detector = detector
          record.goniometer = goniometer
          record.scan = scan
          record.template = template
          record.index = i
          records[template_format % i] = record

      # Add an imageset
      elif imageset['__id__'] == 'ImageSet':
        for image in imageset['images']:
          record = ImageRecord()
          if image.get('beam', None):
            record.beam = Beam.from_dict(blist[image['beam']])
          if image.get('detector', None):
            record.detector = Detector.from_dict(dlist[image['detector']])
          if image.get('goniometer', None):
            record.goniometer = Goniometer.from_dict(glist[image['goniometer']])
          if image.get('scan', None):
            record.scan = Scan.from_dict(slist[image['scan']])
          records[image['filename']] = record
      else:
        raise RuntimeError('unknown imageset id %s' % imageset['__id__'])

    # Get the format class from the first image
    format_class = Registry.find(records.keys()[0])

    # Return the datablock
    datablock = DataBlock()
    datablock._images = records
    datablock._format_class = format_class
    return datablock

  @staticmethod
  def from_json(string):
    ''' Decode a datablock from JSON string. '''
    from dxtbx.serialize.load import _decode_dict
    import json
    return DataBlockFactory.from_dict(json.loads(
      string, object_hook=_decode_dict))

  @staticmethod
  def from_json_file(filename):
    ''' Decode a datablock from a JSON file. '''
    with open(filename, 'r') as infile:
      return DataBlockFactory.from_json(infile.read())

  @staticmethod
  def from_pickle_file(filename):
    ''' Decode a datablock from a pickle file. '''
    import cPickle as pickle
    with open(filename, 'rb') as infile:
      obj = pickle.load(infile)
      if isinstance(obj, list):
        assert(all(isinstance(db, DataBlock) for db in obj))
      else:
        assert(isinstance(obj, DataBlock))
      return obj
