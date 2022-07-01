from __future__ import absolute_import, division, print_function

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
    import os
    ( fd, self.errfile ) = tempfile.mkstemp()
    os.close( fd )


  def run(self):

    stderr = open( self.errfile, "w" )
    import sys
    #  XXX CATCH CASE WHERE SYS.STDERR HAS BEEN MODIFIED
    stderr_fileno = sys.stderr.fileno() if hasattr(sys.stderr, 'fileno') else 2
    sys.stderr = stderr # adjust Python level
    import os
    os.dup2( stderr.fileno(), stderr_fileno ) # adjust OS level

    super( stderr_capturing_process, self ).run()


  def join(self, timeout = None):

    super( stderr_capturing_process, self ).join( timeout )

    if not self.errfile:
      return

    try:
      stderr = open( self.errfile )
      self.stacktrace = stderr.read()

    except IOError:
      self.stacktrace = ""

    else:
      stderr.close()

      import os
      os.remove( self.errfile )

    finally:
      self.errfile = None

    if self.propagate_error_message and self.stacktrace:
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
