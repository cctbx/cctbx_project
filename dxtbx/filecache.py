#!/usr/bin/env python
#
#   Copyright (C) 2015 Diamond Light Source, Markus Gerstel
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A shared caching layer for file-like objects.
# pseudo_file objects can be used as drop-in replacements for actual file
# handles to provide a transparent caching layer to avoid reading multiple
# times from disk or network.
#
# To create a pseudo_file instance, encapsulate a 'real' file handler
# inside a lazy cache object:
#   from dxtbx.filecache import lazy_file_cache
#   cache = lazy_file_cache(open(filename, 'rb'))
#
# Finally use a reference to the cache object to create one or many pseudo_file
# instances:
#   fh1 = cache.open()
#   from dxtbx.filecache import pseudo_file
#   fh2 = pseudo_file(cache) # equivalent
#   fh3 = pseudo_file(cache)
#   ...
#
# Each pseudo_file instance can then be treated as a proper read-only
# file handle, but will benefit from a shared cache:
#   with cache.open() as fh:
#     fh.read(100)
#     fh.readline()
#     fh.seek(500)
#     fh.read()
#     fh.readlines()
#     fh.close()
#
# To flush the cache and free the memory you can use
#     cache.close()
# This will drop the cache when all associated file handles are closed.
# To instantly drop the cache you can use
#     cache.force_close()
# Any further access attempts will then result in an exception.

from __future__ import absolute_import, division, print_function
from cStringIO import StringIO
import os
from threading import Lock

class lazy_file_cache():
  '''An object providing shared cached access to files'''

  def __init__(self, file_object):
    '''Create a shared cache based on a single file handle.'''
    # Uncomment this line to print debug information
    # self._debug = self._debug_enable

    # Size of a block to read. This should not be smaller than 4k, which is the
    # default block size on many systems.
    self._page_size = 4096

    # Reference to the underlying file object. When no further information can
    # be gained from the file (ie. it has been read once completely), it may
    # be closed.
    self._all_cached = False
    self._file = file_object
    self._file_lock = Lock()

    # StringIO object containing cached information
    self._cache_object = StringIO()
    self._cache_size = 0
    self._cache_limit = 4 * 1024 * 1024
    self._cache_limit_reached = False

    # Current status of lazy cache towards client objects.
    # Opening new file handles is disallowed when the object is closing.
    self._closing = False
    # When the lazy cache object is closed no further access is allowed,
    # and cached information is dropped.
    self._closed = False

    # Number of currently registered client objects
    self._reference_counter = 0
    self._reference_counter_lock = Lock()

    self._debug("Created cache object for %s: %s" % (str(file_object), str(self)))

  def _debug(self, string):
    pass

  def _debug_enable(self, string):
    print("%s: %s" % (format(id(self), '#x'), string))

  def __del__(self):
    '''Close file handles and drop cache on garbage collection.'''
    self.force_close()

  def _cache_up_to(self, position):
    '''Ensure that the file has been read up to "position"'''

    # Limit reads to upper cache size limit
    if position > self._cache_limit:
      position = self._cache_limit

    # Is read actually necessary?
    if self._all_cached or (position < self._cache_size):
      return

    with self._file_lock:
      read_bytes = position - self._cache_size
      # This looks like an unnecessary check, but required for concurrency
      if read_bytes <= 0:
        return

      # Do not read less than a memory page, round up read size to a
      # multiple of page sizes if necessary.
      read_bytes = self._page_size * ((read_bytes + self._page_size - 1) // self._page_size)
      self._debug("Reading %d bytes from file" % read_bytes)

      expected_cache_size = self._cache_size + read_bytes

      data = self._file.read(read_bytes)
      self._cache_object.seek(self._cache_size)
      self._cache_object.write(data)
      self._cache_size = self._cache_object.tell()

      self._debug("Read %d bytes from file, cache size %d" % (len(data), self._cache_size))

      if (expected_cache_size != self._cache_size):
        # must have reached end of file
        self._debug("Lazy cache reached EOF (%d != %d)" % (expected_cache_size, self._cache_size))
        self._all_cached = True
        self._close_file()

      if self._cache_size >= self._cache_limit:
        self._cache_limit_reached = True
        self._debug("Cache limit reached with %d bytes" % self._cache_size)

  def _cache_all(self):
    '''Read entire remaining file into cache.'''

    # Is read actually necessary?
    if self._all_cached:
      return

    # Is read actually possible?
    if self._cache_limit_reached:
      return

    with self._file_lock:
      # Check again with lock held, required for concurrency
      if self._all_cached or self._cache_limit_reached: return

      self._debug("Reading remaining file into cache")

      data = self._file.read(self._cache_limit - self._cache_size)
      self._cache_object.seek(self._cache_size)
      self._cache_object.write(data)
      self._debug("Read %d bytes" % len(data))
      self._cache_size += len(data)

      if self._cache_size >= self._cache_limit:
        # Don't cache more than the set limit. In this case keep file handler
        # open and pass read requests through.
        self._cache_limit_reached = True
        self._debug("Cache size limit reached")
        return

      self._all_cached = True
      self._close_file()

  def _check_not_closed(self):
    if self._closed:
      self._debug("Instance tried to access closed cache")
      raise IOError('Accessing lazy file cache %s after closing is not allowed' % str(self))

  def _close_file(self):
    if self._file is not None:
      self._debug("Closing lazy cache internal file handle (%d bytes read)" % self._cache_size)
      self._file.close()
      self._file = None

  def open(self):
    '''Create and return a new pseudo_file object for this cache.'''
    return pseudo_file(self)

  def close(self):
    if not self._closing:
      self._closing = True
      self._debug("Closing lazy cache %s" % str(self))
      with self._reference_counter_lock:
        if self._reference_counter == 0:
          self.force_close()

  def force_close(self):
    '''Close encapsulated file handle, drop cache and prevent further reads.'''
    if not self._closed:
      self._closing = True
      self._closed = True
      self._close_file()
      if self._reference_counter > 0:
        self._debug("Warning: %d connected instances remain" % self._reference_counter)
      if self._cache_object is not None:
        self._cache_object.close()
        self._cache_object = None

  def register(self):
    '''Register a client object. Reference counting for debug purposes.'''
    with self._reference_counter_lock:
      self._check_not_closed()
      if self._closing:
        self._debug("Instance tried to connect to closing cache")
        raise IOError('Cannot open new file handle: lazy file cache is closing')
      self._reference_counter += 1
    self._debug("Instance connected to lazy cache")

  def unregister(self):
    '''Unregister a client object. Reference counting for debug purposes.'''
    self._debug("Instance disconnected from lazy cache")
    with self._reference_counter_lock:
      self._reference_counter -= 1
      if self._closing and (self._reference_counter == 0):
        self.force_close()

  def pass_read(self, start=0, maxbytes=None):
    '''Read from position start up to maxbytes bytes from file.
       If maxbytes is not set, read the entire file.'''
    self._check_not_closed()

    # Do we need to pass the read request to the underlying file object?
    passthrough = False

    if not self._all_cached:
      # Ensure that relevant data is in cache
      if maxbytes is None:
        self._cache_all()
        if not self._all_cached:
          passthrough = True
      else:
        self._cache_up_to(start + maxbytes)
        if self._cache_limit_reached:
          passthrough = start + maxbytes > self._cache_size

    if passthrough:
      with self._file_lock:
        self._file.seek(start)
        if maxbytes is None:
          self._debug("Passing through read from %d" % start)
          return self._file.read(), self._file.tell()
        else:
          self._debug("Passing through read from %d to %d" % (start, start + maxbytes - 1))
          return self._file.read(maxbytes), self._file.tell()
    else:
      self._cache_object.seek(start)
      if maxbytes is None:
        return self._cache_object.read(), self._cache_object.tell()
      else:
        return self._cache_object.read(maxbytes), self._cache_object.tell()

  def pass_readline(self, start=0, maxbytes=None):
    '''Read a line from file, but no more than maxbytes bytes.'''
    self._check_not_closed()

    if self._all_cached:
      self._cache_object.seek(start)
      if maxbytes is None:
        return self._cache_object.readline(), self._cache_object.tell()
      else:
        return self._cache_object.readline(maxbytes), self._cache_object.tell()

    if self._cache_size <= start:
      self._cache_up_to(start + self._page_size)

    self._cache_object.seek(start)
    if maxbytes is None:
      line_candidate = self._cache_object.readline()
    else:
      line_candidate = self._cache_object.readline(maxbytes)

    end_position = self._cache_object.tell()

    if end_position < self._cache_size or line_candidate.endswith('\n'):
      # Found a complete line within the cache
      return line_candidate, end_position

    if (maxbytes is not None) and (end_position == start + maxbytes):
      # Fulfilled maxbytes condition within the cache
      return line_candidate, end_position

    # Need more data
    while end_position == self._cache_size and not self._cache_limit_reached:
      # Do we have a complete line?
      if line_candidate.endswith('\n'):
        return line_candidate, end_position

      # Ran against cache limit. Extend cache
      self._cache_up_to(self._cache_size + self._page_size)
      self._cache_object.seek(end_position)

      # Continue reading
      if maxbytes is None:
        line_candidate += self._cache_object.readline()
      else:
        foundbytes = end_position - start
        line_candidate += self._cache_object.readline(maxbytes - foundbytes)
      end_position = self._cache_object.tell()

    # Do we have a complete line?
    if line_candidate.endswith('\n') or self._all_cached:
      return line_candidate, end_position

    assert self._cache_limit_reached # Only legitimate way of reaching here

    with self._file_lock:
      self._file.seek(end_position)
      if maxbytes is None:
        return line_candidate + self._file.readline(), self._file.tell()
      else:
        foundbytes = end_position - start
        return line_candidate + self._file.readline(maxbytes - foundbytes), self._file.tell()


class pseudo_file():
  '''A file-like object that serves as frontend to a dxtbx lazy file cache.'''

  def __init__(self, lazy_cache_object):
    self._cache_object = lazy_cache_object
    self._cache_object.register()
    self._closed = False
    self._seek = 0

  def __del__(self):
    self.close()

  def __enter__(self):
    return self

  def __iter__(self):
    return self

  def __exit__(self, exc_type, exc_val, exc_tb):
    self.close()
    return False

  def _check_not_closed(self):
    if self._closed:
      raise IOError('Accessing lazy file cache after closing is not allowed')

  def close(self):
    self._closed = True
    if self._cache_object is not None:
      self._cache_object.unregister()
      self._cache_object = None

  def flush(self):
    self._check_not_closed()

  def next(self):
    data = self.readline()
    if data == "":
      raise StopIteration()
    return data
  __next__ = next

  def read(self, size=-1):
    self._check_not_closed()
    if size > 0:
      data, self._seek = self._cache_object.pass_read(start=self._seek, maxbytes=size)
    elif size == 0:
      data = ''
    else:
      data, self._seek = self._cache_object.pass_read(start=self._seek)
    return data

  def readline(self, size=-1):
    self._check_not_closed()
    if (size > 0):
      data, self._seek = self._cache_object.pass_readline(start=self._seek, maxbytes=size)
    else:
      data, self._seek = self._cache_object.pass_readline(start=self._seek)
    return data

  def readlines(self, sizehint=-1):
    self._check_not_closed()
    if (sizehint > 0):
      data, self._seek = self._cache_object.pass_read(start=self._seek, maxbytes=sizehint)
    else:
      data, self._seek = self._cache_object.pass_read(start=self._seek)
    return data.splitlines(True)

  def seek(self, offset, whence=os.SEEK_SET):
    self._check_not_closed()
    if whence == os.SEEK_SET:
      self._seek = offset
    elif whence == os.SEEK_CUR:
      self._seek += offset
    else:
      raise NotImplementedError('Seeking relative to file length is not supported')

  def tell(self):
    self._check_not_closed()
    return self._seek

  def truncate(self, size=0):
    raise NotImplementedError('Truncating lazy file caches is not allowed')

  def write(self, string):
    raise NotImplementedError('Writing to lazy file caches is not allowed')

  def writelines(self, sequence):
    raise NotImplementedError('Writing to lazy file caches is not allowed')
