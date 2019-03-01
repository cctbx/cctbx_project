from __future__ import absolute_import, division, print_function

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

    except Exception as e:
      #import traceback
      #print traceback.format_exc()
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
    self.stderr.close()
    self.stderr = None

    if self.display_stderr:
      import sys
      sys.stderr.write( stderr )
      sys.stderr.write( "\n" )

    if self.exitcode != 0:
      self.err = RuntimeError( stderr )

    else:
      self.err = None


class ThreadFactory(object):
  """
  Creator for threading.Thread or scheduling_helpers.Thread
  """

  def __init__(self, name = None, preserve_exception_message = False, **kwargs):

    if preserve_exception_message:
      self.factory = Thread

    else:
      self.factory = threading.Thread

    self.name = name


  def __call__(self, target, args = (), kwargs = {}):

    return self.factory(
      name = self.name,
      target = target,
      args = args,
      kwargs = kwargs,
      )


class ProcessFactory(object):
  """
  Creator for multiprocessing.Process or scheduling_helpers.Process
  """

  def __init__(self, name = None, preserve_exception_message = False, **kwargs):

    if preserve_exception_message:
      self.factory = Process

    else:
      self.factory = multiprocessing.Process

    self.name = name


  def __call__(self, target, args = (), kwargs = {}):

    return self.factory(
      name = self.name,
      target = target,
      args = args,
      kwargs = kwargs,
      )


class QQFactory(object):
  """
  Creator pattern for Queue.Queue, also include destruction

  Note this is a singleton object
  """

  @staticmethod
  def create():

    from six.moves import queue
    return queue.Queue()


  @staticmethod
  def destroy(queue):

    pass


class MPQFactory(object):
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


class MPManagerQFactory(object):
  """
  Creator pattern for multiprocessing.Manager.Queue, also include destruction

  Note this is a not singleton object
  """

  def __init__(self):

    self.manager = multiprocessing.Manager()


  def create(self):

    return self.manager.Queue()


  def destroy(self, queue):

    pass
