#!/usr/bin/env python
# Registry.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A registry class to handle Format classes and provide lists of them when
# this is useful for i.e. identifying the best tool to read a given range
# of image formats.

from __future__ import absolute_import, division

from dxtbx.format.RegistryHelpers import LoadFormatClasses
from libtbx.utils import Sorry
from .Format import Format

class SorryIOError(IOError, Sorry):
  """Fusion of Sorry and IO Errors.

  This is so we can send a better user-error, but still allow the simplicity
  of just catching an IOError in any parent scopes.
  """
  def __init__(self, ioerror):
    self.errno = ioerror.errno
    self.strerror = ioerror.strerror
    self.filename = ioerror.filename
    super(SorryIOError, self).__init__("{}: {}".format(self.strerror, self.filename))

  def __str__(self):
    return self.message

class _Registry:
  '''A class to handle all of the recognised image formats within xia2
  working towards the generalization project in #1555 and specifically
  to address the requirements in #1573.'''

  def __init__(self):
    '''Set myself up - N.B. this could trigger a directory search for
    different Format classes to auto-register them.'''

    self._formats = []

    self._setup = False

  def setup(self):
    '''Look to import format defining modules from around the place -
    this will look in dxtbx/format and $HOME/.xia2/ for files starting
    with Format and ending in .py to allow user extensibility.'''

    if self._setup:
      return

    LoadFormatClasses()

    self._setup = True

  def add(self, format):
    '''Register a new image format with the registry. N.B. to work
    this object must be inherited from the base Format class.'''

    from dxtbx.format.Format import Format

    assert(issubclass(format, Format))
    if not format in self._formats:
      self._formats.append(format)

  def get(self):
    '''Get a list of image formats registered here.'''

    self.setup()

    return tuple(self._formats)

  def find(self, image_file):
    '''More useful - find the best format handler in the registry for your
    image file. N.B. this is in principle a factory function.'''

    self.setup()

    # Recursively check whether any of the children understand
    # image_file, in which case they are preferred over the parent
    # format.
    def recurse(format, image_file):
      for child in format._children:
        try:
          if child.understand(image_file):
            return recurse(child, image_file)
        except Exception:
          pass
      return format

    for format in self._formats:
      try:
        if format.understand(image_file):
          return recurse(format, image_file)
      except Exception:
        pass

    # Try opening the file; this could be an easy reason for failure
    if not Format.is_url(image_file):
      try:
        with open(image_file) as f:
          pass
      except IOError as e:
        # Assume that the OS file error makes sense to the user
        raise SorryIOError(e)

    raise IOError('no format support found for %s' % image_file)

class Registry:
  '''A class to turn the registry above into a singleton, so that we
  can work around some of the more knotty dependency loops which come
  out of things like this. This is something of a boiler plate.'''

  __instance = None

  def __init__(self):
    if Registry.__instance is None:
      Registry.__instance = _Registry()

  def __call__(self):
    return self

  def __getattr__(self, attr):
    return getattr(self.__instance, attr)

  def __setattr__(self, attr, value):
    return setattr(self.__instance, attr, value)

Registry = Registry()
