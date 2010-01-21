
import os
from libtbx import group_args

class manager (object) :
  def __init__ (self) :
    self.handlers = []
    self._pid = os.getpid()

  def register_handler (self, handler) :
    assert hasattr(handler, "__call__")
    self.handlers.append(handler)

  def __call__ (self, message, data, accumulate=False, cached=True) :
#    assert os.getpid() == self._pid
    for handler in self.handlers :
      handler(message=message,
              data=data,
              accumulate=accumulate,
              cached=cached)

  def add_citation (self, citation_info) :
    self.__call__("citation", citation_info, accumulate=True)

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
if not hasattr(libtbx, "callback") :
  libtbx.call_back = manager()
