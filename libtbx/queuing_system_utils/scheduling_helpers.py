from __future__ import division

import threading
import multiprocessing

class Thread(threading.Thread):
  """
  A Thread-specialization that adds the exitcode attribute

  http://stackoverflow.com/questions/986616/python-thread-exit-code
  """

  def run(self):

    try:
      super( Thread, self ).run()

    except Exception, e:
      import traceback
      print traceback.format_exc()
      self.exitcode = 1
      self.err = e

    else:
      self.exitcode = 0
      self.err = None


class Process(multiprocessing.Process):
  """
  A Process-specialization that adds the .err attribute
  """

  def __init__(
    self,
    group = None,
    target = None,
    name = None,
    display_stderr = False,
    args = (),
    kwargs = {},
    ):

    super( Process, self ).__init__(group, target, name, args, kwargs)
    self.display_stderr = display_stderr

    import tempfile
    self.stderr = tempfile.TemporaryFile()


  def run(self):

    import sys
    stderr_fileno = sys.stderr.fileno()
    sys.stderr = self.stderr # adjust Python level
    import os
    os.dup2( self.stderr.fileno(), stderr_fileno )

    super( Process, self ).run()


  def join(self, timeout = None):

    super( Process, self ).join( timeout )

    if not self.stderr:
      return

    self.stderr.seek( 0 )
    stderr = self.stderr.read()
    self.stderr = None

    if self.display_stderr:
      import sys
      sys.stderr.write( stderr )
      sys.stderr.write( "\n" )

    if self.exitcode != 0:
      self.err = RuntimeError( stderr )

    else:
      self.err = None

