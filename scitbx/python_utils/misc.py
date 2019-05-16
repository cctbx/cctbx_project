from __future__ import absolute_import, division, print_function
import sys
from six.moves import range
from six.moves import zip

class store(object):

  def __init__(self, **kw):
    self.__dict__.update(kw)

class sorted_store(object):

  def keys(self):
    raise RuntimeError("Programming error: derived class must override keys()")

  def __init__(self, *args, **kw):
    assert len(args) + len(kw) == len(list(self.keys()))
    for key,value in zip(list(self.keys())[:len(args)], args):
      setattr(self, key, value)
    for key,value in kw.items():
      assert key in list(self.keys())
      assert getattr(self, key, None) is None
      setattr(self, key, value)

  def show(self, f=None, indentation=""):
    if (f is None): f = sys.stdout
    for key in self.keys():
      print("%s%s:" % (indentation, key), getattr(self, key), file=f)

def get_caller_name(n_back=2):
  try: raise Exception
  except Exception:
    t = sys.exc_info()[2]
    f = t.tb_frame
    for i in range(n_back):
      if (f.f_back is None): return None
      f = f.f_back
    return f.f_code.co_name
