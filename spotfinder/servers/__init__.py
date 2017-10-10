from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import object
import io as StringIO,sys

class LoggingFramework(object):
  def __init__(self):
    self.k = StringIO.StringIO()
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

