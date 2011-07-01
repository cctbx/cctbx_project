
from libtbx import str_utils
from libtbx import group_args, adopt_init_args
import os
import sys

class manager (object) :
  def __init__ (self) :
    self.handlers = []
    self._pid = os.getpid()
    self._log = None

  def set_warning_log (self, log) :
    self._log = log

  def add_piped_callback (self, connection) :
    cb = piped_callback(connection)
    self.register_handler(cb)

  def register_handler (self, handler) :
    assert hasattr(handler, "__call__")
    self.handlers.append(handler)

  def __call__ (self, message, data, accumulate=True, cached=True) :
    for handler in self.handlers :
      handler(message=message,
              data=data,
              accumulate=accumulate,
              cached=cached)

  def add_citation (self, citation_info) :
    self.__call__("citation", citation_info, accumulate=True)

  def showwarning (self, message, category, filename, lineno, file=None) :
    self.warn(message)

  def warn (self, message) :
    self.__call__(message="warn", data=message, accumulate=True, cached=True)
    if (self._log is not None) :
      log = self._log
    else :
      log = sys.stdout
    msg = "WARNING: %s\n" % message
    print >> log, ""
    for line in str_utils.line_breaker(msg, 72) :
      print >> log, "  " + line
    #print >> log, str_utils.wordwrap(msg, 79)

class piped_callback (object) :
  def __init__(self, connection) :
    adopt_init_args(self, locals())

  def __call__ (self, message, data, accumulate=False, cached=True) :
    self.connection.send(
      group_args(message=message,
                 data=data,
                 accumulate=accumulate,
                 cached=cached))

import libtbx
if not hasattr(libtbx, "call_back") :
  libtbx.call_back = manager()
  libtbx.warn = libtbx.call_back.warn
