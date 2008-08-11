from libtbx import adopt_init_args # XXX backward compatibility 2005_07_29
from libtbx.utils import plural_s # XXX backward compatibility 2006-10-16
from libtbx.utils import user_plus_sys_time # XXX backward compatibility 2006-10-16
from libtbx.utils import time_log # XXX backward compatibility 2006-10-16
from libtbx.utils import human_readable_time # XXX backward compatibility 2006-10-16
from libtbx.utils import human_readable_time_as_seconds # XXX backward compatibility 2006-10-16
from libtbx.utils import show_total_time # XXX backward compatibility 2006-11-08
from libtbx.utils import input_with_prompt # XXX backward compatibility 2006-11-08
from libtbx.str_utils import line_feeder # XXX backward compatibility 2006-11-08
import sys, os

class store(object):

  def __init__(self, **kw):
    self.__dict__.update(kw)

class sorted_store(object):

  def keys(self):
    raise RuntimeError, "Programming error: derived class must override keys()"

  def __init__(self, *args, **kw):
    assert len(args) + len(kw) == len(self.keys())
    for key,value in zip(self.keys()[:len(args)], args):
      setattr(self, key, value)
    for key,value in kw.items():
      assert key in self.keys()
      assert getattr(self, key, None) is None
      setattr(self, key, value)

  def show(self, f=None, indentation=""):
    if (f is None): f = sys.stdout
    for key in self.keys():
      print >> f, "%s%s:" % (indentation, key), getattr(self, key)

def get_caller_name(n_back=2):
  try: raise Exception
  except:
    t = sys.exc_info()[2]
    f = t.tb_frame
    for i in xrange(n_back):
      if (f.f_back is None): return None
      f = f.f_back
    return f.f_code.co_name
