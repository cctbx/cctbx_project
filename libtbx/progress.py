from __future__ import absolute_import, division, print_function


class ProgressError(Exception):
  """
  Module exception
  """


class TimeoutError(ProgressError):
  """
  Timeout exceeded in wait
  """


class streamprint(object):
  """
  Prints a character into a stream to show progress
  """

  def __init__(self, stream, character = "."):

    self.stream = stream
    self.character = character


  def __call__(self):

    self.stream.write( self.character )
    self.stream.flush()


class complete_on_success(object):
  """
  A condition that is complete when no expected errors occur. The result is
  stored to avoid problems with concurrency, i.e. that a subsequent call to the
  function returns no result
  """

  def __init__(self, func, excspec):

    self.func = func
    self.excspec = excspec


  def __call__(self):

    try:
      self.result = self.func()

    except self.excspec:
      return False

    return True


def wait(condition, waittime = 2, timeout = 600, callback = lambda: None):
  """
  Waits for a condition to become True
  """

  import time
  start = time.time()

  while not condition():
    if timeout < time.time() - start:
      raise TimeoutError("Timeout (%s s) exceeded" % timeout)

    time.sleep( 2 )
    callback()
