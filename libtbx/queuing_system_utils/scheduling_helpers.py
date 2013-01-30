from __future__ import division

import threading

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

