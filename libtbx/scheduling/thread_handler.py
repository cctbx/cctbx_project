from __future__ import absolute_import, division, print_function

import threading
from six.moves import queue


class exception_capturing_thread(threading.Thread):
  """
  A Thread-specialization that adds the exitcode attribute

  http://stackoverflow.com/questions/986616/python-thread-exit-code
  """

  def __init__(
    self,
    group=None,
    target=None,
    name=None,
    propagate_error_message = False,
    args=(),
    kwargs={},
    ):

    super( exception_capturing_thread, self ).__init__(group, target, name, args, kwargs)
    self.propagate_error_message = propagate_error_message


  def run(self):

    try:
      super( exception_capturing_thread, self ).run()

    except Exception as e:
      self.exitcode = 1
      self.err = e

      import traceback
      import sys
      self.stacktrace = "".join( traceback.format_tb( tb = sys.exc_info()[2] ) )

      if self.propagate_error_message:
        traceback.print_exc()

    else:
      self.exitcode = 0
      self.err = None
      self.stacktrace = None


class qfactory(object):
  """
  Creator pattern for queue.Queue, also include destruction

  Note this is a singleton object
  """

  @staticmethod
  def create():

    return queue.Queue()


  @staticmethod
  def destroy(q):

    pass
