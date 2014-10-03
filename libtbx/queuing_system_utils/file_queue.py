from __future__ import division

from Queue import Empty
import pickle
import fcntl
import tempfile
import os


class Instant(object):
  """
  Timeout immediately
  """

  def __call__(self):

    raise Empty, "No data found in queue"


class Timed(object):
  """
  Timeout after given time
  """

  def __init__(self, timeout, waittime):

    self.timeout = timeout
    self.waittime = waittime


  def __call__(self):

    if 0 < self.timeout:
      import time
      waittime = min( self.timeout, self.waittime )
      self.timeout -= self.waittime
      time.sleep( waittime )

    else:
      raise Empty, "No data found in queue within timeout"


class Eternal(object):
  """
  No timeout
  """

  def __init__(self, waittime):

    self.waittime = waittime


  def __call__(self):

    import time
    time.sleep( self.waittime )


def get_timeout_object(block, timeout, waittime):

  if not block:
    return Instant()

  else:
    if timeout is not None:
      return Timed( timeout = timeout, waittime = waittime )

    else:
      return Eternal( waittime = waittime )


def temp_name(prefix, suffix, folder):

  ( fd, fname ) = tempfile.mkstemp( suffix = suffix, prefix = prefix, dir = folder )
  os.close( fd )
  return fname


def _lock(fobj, timeout):

  while True:
    try:
      fcntl.lockf( fobj.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB )
      break

    except IOError:
      timeout()


def _unlock(fobj):

  fcntl.lockf( fobj.fileno(), fcntl.LOCK_UN )


class Queue(object):
  """
  File-based queue. Can be used cross-host via an NFS 3.0 filesystem.
  """

  def __init__(self, prefix = "tmp", folder = ".", waittime = 0.1):

    self.waittime = waittime

    self._put_file = temp_name( prefix = prefix, suffix = ".put", folder = folder )
    self._get_file = temp_name( prefix = prefix, suffix = ".get", folder = folder )
    self._offset_file = temp_name( prefix = prefix, suffix = ".offset", folder = folder )

    # No need for locking, as queue has not been constructed yet
    with open( self._offset_file, "w" ) as foff:
      pickle.dump( 0, foff )


  def put(self, value, block = True, timeout = None):

    tobj = get_timeout_object(
      block = block,
      timeout = timeout,
      waittime = self.waittime,
      )

    with open( self._put_file, "r+b" ) as fput:
      _lock( fput, tobj )

      try:
        fput.seek( 0, os.SEEK_END )
        pickle.dump( value, fput )
        fput.flush()

      finally:
        _unlock( fput )


  def get_next_item(self, offset, timeout):

    with open( self._get_file, "r+b" ) as fget:
      fget.seek( offset, os.SEEK_SET )

      try:
        value = pickle.load( fget )
        offset = fget.tell()

      except EOFError:
        while os.path.getsize( self._put_file ) == 0:
          timeout()

        fget.seek( 0, os.SEEK_SET )
        value = self.read_one_and_transfer_remanining_put_contents(
          fget = fget,
          timeout = timeout,
          )
        fget.truncate()
        offset = 0

    return ( value, offset )


  def read_one_and_transfer_remanining_put_contents(self, fget, timeout):

    with open( self._put_file, "r+b" ) as fput:
      _lock( fput, timeout )

      try:
        value = pickle.load( fput )
        chunksize = 1024

        while True:
          data = fput.read( chunksize )
          fget.write( data )

          if len( data ) < chunksize:
            break

        fput.seek( 0, os.SEEK_SET )
        fput.truncate()
        fput.flush()

      finally:
        _unlock( fput )

    return value


  def get(self, block = True, timeout = None):

    tobj = get_timeout_object(
      block = block,
      timeout = timeout,
      waittime = self.waittime,
      )

    with open( self._offset_file, "r+b" ) as foff:
      _lock( foff, tobj )

      try:
        offset = pickle.load( foff )
        ( value, offset ) = self.get_next_item( offset = offset, timeout = tobj )
        foff.seek( 0 )
        pickle.dump( offset, foff )
        foff.truncate()

      finally:
        _unlock( foff )

    return value


  def close(self):

    for fname in [ self._get_file, self._put_file, self._offset_file ]:
      if os.path.exists( fname ):
        os.remove( fname )


class QFactory(object):
  """
  Creator pattern for file queue, also include destruction
  """

  def __init__(self, prefix = "tmp", folder = ".", waittime = 0.1):

    self.prefix = prefix
    self.folder = folder
    self.waittime = waittime


  def create(self):

    return Queue( prefix = self.prefix, folder = self.folder, waittime = self.waittime )


  @staticmethod
  def destroy(queue):

    queue.close()
