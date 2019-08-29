from __future__ import absolute_import, division, print_function
from six.moves import cStringIO as StringIO
import sys

class LoggingFramework:
  def __init__(self):
    self.k = StringIO()
    self.current_out = sys.stdout
    self.current_err = sys.stderr
    sys.stdout = self.k
    sys.stderr = self.k

  def __del__(self):
    sys.stdout = self.current_out
    sys.stderr = self.current_err
    self.k.flush()
    self.k.close()

  def getvalue(self): return self.k.getvalue()

