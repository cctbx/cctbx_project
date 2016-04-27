from __future__ import division

import multiprocessing

class stderr_capturing_process(multiprocessing.Process):
  """
  A Process-specialization that adds the .err attribute
  """

  def __init__(
    self,
    group = None,
    target = None,
    name = None,
    propagate_error_message = False,
    args = (),
    kwargs = {},
    ):

    super( stderr_capturing_process, self ).__init__(group, target, name, args, kwargs)
    self.propagate_error_message = propagate_error_message

    import tempfile
    self.stderr = tempfile.TemporaryFile()


  def run(self):

    import sys
    stderr_fileno = sys.stderr.fileno()
    sys.stderr = self.stderr # adjust Python level
    import os
    os.dup2( self.stderr.fileno(), stderr_fileno )

    super( stderr_capturing_process, self ).run()


  def join(self, timeout = None):

    super( stderr_capturing_process, self ).join( timeout )

    if not self.stderr:
      return

    self.stderr.seek( 0 )
    self.stacktrace = self.stderr.read()
    self.stderr.close()
    self.stderr = None

    if self.propagate_error_message:
      import sys
      sys.stderr.write( self.stacktrace )
      sys.stderr.write( "\n" )


class fifo_qfactory(object):
  """
  Creator pattern for multiprocessing.Queue, also include destruction

  Note this is a singleton object
  """

  @staticmethod
  def create():

    return multiprocessing.Queue()


  @staticmethod
  def destroy(queue):

    pass


class managed_qfactory(object):
  """
  Creator pattern for multiprocessing.Manager.Queue, also include destruction

  Note this is not a singleton object
  """

  def __init__(self):

    self.manager = multiprocessing.Manager()


  def create(self):

    return self.manager.Queue()


  def destroy(self, queue):

    pass
