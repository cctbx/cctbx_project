#!/usr/bin/env python
# Format.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A top-level class to represent image formats which does little else but
# (i) establish an abstract class for what needs to be implemented and
# (ii) include the format registration code for any image formats which
# inherit from this. This will also contain links to the static methods
# from the X(component)Factories which will allow construction of e.g.
# goniometers etc. from the headers and hence a format specific factory.

from __future__ import absolute_import, division

import sys

try:
  import bz2
except: # intentional
  bz2 = None

try:
  import gzip
except: # intentional
  gzip = None

import dxtbx.filecache_controller
import exceptions

# import access to all of the factories that we will be needing

from dxtbx.model.goniometer import Goniometer, GoniometerFactory
from dxtbx.model.detector import Detector, DetectorFactory
from dxtbx.model.beam import Beam, BeamFactory
from dxtbx.model.scan import Scan, ScanFactory

class _MetaFormat(type):
  '''A metaclass for the Format base class (and hence all format classes)
  to allow autoregistration of the class implementations.'''

  def __init__(self, name, bases, attributes):
    super(_MetaFormat, self).__init__(name, bases, attributes)

    # Do-nothing until the Format module, defining the base class,
    # has been loaded.
    try:
      sys.modules[Format.__module__]
    except NameError:
      return

    # Add the class to the registry if it is directly derived from
    # Format.
    self._children = []
    if Format in bases:
      from dxtbx.format.Registry import Registry
      Registry.add(self)
      return

    # Add the class to the list of children of its superclasses.
    for base in bases:
      base._children.append(self)
    return

  _cache_controller = dxtbx.filecache_controller.simple_controller()
  @classmethod
  def get_cache_controller(cls):
    return cls._cache_controller

class Format(object):
  '''A base class for the representation and interrogation of diffraction
  image formats, from which all classes for reading the header should be
  inherited. This includes: autoregistration of implementation classes,
  stubs which need to be overridden and links to static factory methods
  which will prove to be useful in other implementations.'''

  __metaclass__ = _MetaFormat

  @staticmethod
  def understand(image_file):
    '''Overload this to publish whether this class instance
    understands a given file.  N.B. to say that we understand it,
    return True.  If a subclass also understands the image
    (because, for example, its detector serial number takes a
    certain value), it will by definition understand it better
    than its superclass.  Thus, the preferred class will be the
    deepest subclass in the inheritance hierarchy.  Finally, if
    you are writing this subclass for a given instrument and you
    are given a different example return False.

    Implementing understand() in a subclass, one can safely assume
    that the superclasss understand() function returned True.
    The understand() function of two different classes directly
    derived from the same base should never both return True for
    the same input image.'''

    return False

  @classmethod
  def ignore(cls):
    return False

  def __init__(self, image_file, **kwargs):
    '''Initialize a class instance from an image file.'''

    self._image_file = image_file

    self._goniometer_instance = None
    self._detector_instance = None
    self._beam_instance = None
    self._scan_instance = None

    self._goniometer_factory = GoniometerFactory
    self._detector_factory = DetectorFactory
    self._beam_factory = BeamFactory
    self._scan_factory = ScanFactory

    self.setup()

  def setup(self):
    '''Read the image file, construct the information which we will be
    wanting about the experiment from this. N.B. in your implementation
    of this you will probably want to make use of the static methods
    below and probably add some format parsing code too. Please also keep
    in mind that your implementation may be further subclassed by
    someone else.'''

    self._start()

    try:
      goniometer_instance = self._goniometer()
      #assert(isinstance(goniometer_instance, Goniometer))
      self._goniometer_instance = goniometer_instance

      detector_instance = self._detector()
      #assert(isinstance(detector_instance, Detector))
      self._detector_instance = detector_instance

      beam_instance = self._beam()
      #assert(isinstance(beam_instance, Beam))
      self._beam_instance = beam_instance

      scan_instance = self._scan()
      #assert(isinstance(scan_instance, Scan) or isinstance(scan_instance, list))
      self._scan_instance = scan_instance

    except exceptions.Exception, e:
      # FIXME ideally should not squash the errors here...
      import traceback
      traceback.print_exc()
      pass
    finally:
      self._end()

  def get_goniometer(self):
    '''Get the standard goniometer instance which was derived from the
    image headers.'''

    return self._goniometer_instance

  def get_detector(self):
    '''Get the standard detector instance which was derived from the
    image headers.'''

    return self._detector_instance

  def get_beam(self):
    '''Get the standard beam instance which was derived from the image
    headers.'''

    return self._beam_instance

  def get_scan(self):
    '''Get the standard scan instance which was derived from the image
    headers.'''

    return self._scan_instance

  def get_raw_data(self):
    '''Get the pixel intensities (i.e. read the image and return as a
    flex array.'''
    try:
      image = self.detectorbase
      image.read()
      raw_data = image.get_raw_data()

      return raw_data
    except Exception:
      return None

  def get_vendortype(self):
    return "no dxtbx Format vendortype"

  def detectorbase_start(self):
    print "Overload detectorbase_start"
    raise RuntimeError('Overload!')

  def get_detectorbase(self):
    '''Return the instance of detector base.'''
    self.detectorbase_start()

    # XXX Temporary proxy to aid in the transition to dxtbx.  Remove
    # once completed.
    class _detectorbase_proxy(object):
      def __init__(self, format_instance):
        self._fi = format_instance
        if not hasattr(self, "vendortype"):
          self.vendortype = "generic"

      def __getattribute__(self, name):
        if name == '__class__':
          return self._fi.__class__
        return object.__getattribute__(self, name)

      def __getattr__(self, name):
        try:
          return self._fi.__getattribute__(name)
        except AttributeError:
          #print >> sys.stderr, \
          #    "requesting iotbx.detectors.detectorbase.%s" % name
          try:
            return self._fi.detectorbase.__getattribute__(name)
          except AttributeError:
            return self._fi.detectorbase.__getattr__(name)

    return _detectorbase_proxy(self)

  @classmethod
  def get_instance(Class, filename):
    if not hasattr(Class, "_current_instance_") or Class._current_filename_ != filename:
      Class._current_instance_ = Class(filename)
      Class._current_filename_ = filename
    return Class._current_instance_

  @classmethod
  def get_detectorbase_factory(Class):
    '''
    Return a factory object to create detector base instances

    '''

    class DetectorBaseFactory(object):

      def __init__(self, filenames):
        self._filenames = filenames

      def get(self, index):
        format_instance = Class.get_instance(self._filenames[index])
        return format_instance.get_detectorbase()

      def __len__(self):
        return len(self._filenames)

      def copy(self, filenames):
        return DetectorBaseFactory(filenames)

    return DetectorBaseFactory

  @classmethod
  def get_reader(Class):
    '''
    Return a reader class

    '''

    class Reader(object):

      def __init__(self, filenames):
        self._filenames = filenames

      def read(self, index):
        format_instance = Class.get_instance(self._filenames[index])
        return format_instance.get_raw_data()

      def paths(self):
        return self._filenames

      def __len__(self):
        return len(self._filenames)

      def copy(self, filenames):
        return Reader(filenames)

    return Reader

  @classmethod
  def get_masker(Class):
    '''
    Return a masker class

    '''
    class Masker(object):

      def __init__(self, filenames):
        self._filenames = filenames

      def get(self, index):
        format_instance = Class.get_instance(self._filenames[index])
        return format_instance.get_mask()

      def paths(self):
        return self._filenames

      def __len__(self):
        return len(self._filenames)

      def copy(self, filenames):
        return Masker(filenames)

      def identifiers(self):
        return self.paths()

    return Masker

  @classmethod
  def get_imageset(Class,
                   filenames,
                   beam=None,
                   detector=None,
                   goniometer=None,
                   scan=None,
                   sweep_as_imageset=False):
    '''
    Factory method to create an imageset

    '''
    from dxtbx.format.FormatMultiImage import FormatMultiImage
    from dxtbx.imageset import ImageSet
    from dxtbx.imageset import ImageSweep

    # Get some information from the format class
    reader = Class.get_reader()(filenames)
    masker = Class.get_masker()(filenames)
    dbfact = Class.get_detectorbase_factory()(filenames)

    # Get the format instance
    format_instance = Class(filenames[0])
    if isinstance(format_instance, FormatMultiImage):
      assert len(filenames) == 1
      is_multi_image_format = True
    else:
      is_multi_image_format = False

    # Read the vendor type
    vendor = format_instance.get_vendortype()

    # Check if we have a sweep
    scan = format_instance.get_scan()
    if scan is not None and scan.get_oscillation()[1] != 0:
      is_sweep = True
    else:
      is_sweep = False

    # Create an imageset or sweep
    if not is_sweep or sweep_as_imageset == True:

      # Create the imageset
      iset = ImageSet(
        reader = reader,
        masker = masker,
        properties = {
          "vendor" : vendor
        },
        detectorbase_factory = dbfact)

      # If any are None then read from format
      if [beam, detector, goniometer, scan].count(None) != 0:

        # Get list of models
        beam = []
        detector = []
        goniometer = []
        scan = []
        if is_multi_image_format:
          for i in range(format_instance.get_num_images()):
            beam.append(format_instance.get_beam(i))
            detector.append(format_instance.get_detector(i))
            goniometer.append(format_instance.get_goniometer(i))
            scan.append(format_instance.get_scan(i))
        else:
          for f in filenames:
            format_instance = Class(f)
            beam.append(format_instance.get_beam())
            detector.append(format_instance.get_detector())
            goniometer.append(format_instance.get_goniometer())
            scan.append(format_instance.get_scan())

      # Set the list of models
      for i in range(len(filenames)):
        format_instance = Class(filenames[i])
        iset.set_beam(beam[i], i)
        iset.set_detector(detector[i], i)
        iset.set_goniometer(goniometer[i], i)
        iset.set_scan(scan[i], i)

    else:

      # If any are None then read from format
      if [beam, detector, goniometer, scan].count(None) != 0:
        beam       = format_instance.get_beam()
        detector   = format_instance.get_detector()
        goniometer = format_instance.get_goniometer()
        scan       = format_instance.get_scan()

        # Get the scan model
        if is_multi_image_format:
          for i in range(format_instance.get_num_images()):
            scan += format_instance.get_scan(i)
        else:
          for f in filenames[1:]:
            format_instance = Class(f)
            scan += format_instance.get_scan()

      # Create the sweep
      iset = ImageSweep(
        reader     = reader,
        masker     = masker,
        beam       = beam,
        detector   = detector,
        goniometer = goniometer,
        scan       = scan,
        properties = {
          "vendor" : vendor
        },
        detectorbase_factory = dbfact)

    # Return the imageset
    return iset

  def get_image_file(self):
    '''Get the image file provided to the constructor.'''

    return self._image_file

  # methods which must be overloaded in order to produce a useful Format
  # class implementation

  def _start(self):
    '''Start code for handling this image file, which may open a link
    to it once, say, and pass this around within the implementation.
    How you use this is up to you, though you probably want to overload
    it...'''

    return

  def _end(self):
    '''Clean up things - keeping in mind that this should run even in the
    case of an exception being raised.'''

    return

  def _goniometer(self):
    '''Overload this method to read the image file however you like so
    long as the result is an goniometer.'''
    return None

  def _detector(self):
    '''Overload this method to read the image file however you like so
    long as the result is an detector.'''
    return None

  def _beam(self):
    '''Overload this method to read the image file however you like so
    long as the result is an beam.'''
    return None

  def _scan(self):
    '''Overload this method to read the image file however you like so
    long as the result is an scan.'''
    return None

  def get_mask(self, index=None, goniometer=None):
    '''Overload this method to provide dynamic masks to be used during
    spotfinding or integration.'''

    return None

  def get_goniometer_shadow_masker(self, goniometer=None):
    '''Overload this method to allow generation of dynamic goniometer shadow
    masks to be used during spotfinding or integration.'''

    return None

  ####################################################################
  #                                                                  #
  # Helper functions for dealing with compressed images.             #
  #                                                                  #
  ####################################################################

  @staticmethod
  def is_url(path):
    '''See if the file is a URL.'''

    from urlparse import urlparse

    # Windows file paths can get caught up in this - check that the
    # first letter is one character (which I think should be safe: all
    # URL types are longer than this right?)

    scheme = urlparse(path).scheme
    if scheme and len(scheme) != 1:
      return True

    return False

  @staticmethod
  def is_bz2(filename):
    '''Check if a file pointed at by filename is bzip2 format.'''

    if not '.bz2' in filename[-4:]:
      return False

    return 'BZh' in open(filename, 'rb').read(3)

  @staticmethod
  def is_gzip(filename):
    '''Check if a file pointed at by filename is gzip compressed.'''

    if not '.gz' in filename[-3:]:
      return False

    magic = open(filename, 'rb').read(2)

    return ord(magic[0]) == 0x1f and ord(magic[1]) == 0x8b

  @classmethod
  def open_file(cls, filename, mode='rb', url=False):
    '''Open file for reading, decompressing silently if necessary,
       caching transparently if possible.'''

    if url and Format.is_url(filename):
      import urllib2
      fh_func = lambda: urllib2.urlopen(filename)

    elif Format.is_bz2(filename):
      if bz2 is None:
        raise RuntimeError, 'bz2 file provided without bz2 module'
      fh_func = lambda: bz2.BZ2File(filename, mode)

    elif Format.is_gzip(filename):
      if gzip is None:
        raise RuntimeError, 'gz file provided without gzip module'
      fh_func = lambda: gzip.GzipFile(filename, mode)

    else:
      fh_func = lambda: open(filename, mode)

##  To disable caching logic:
    #return fh_func()
    return cls.__metaclass__.get_cache_controller() \
      .check(filename, fh_func)
