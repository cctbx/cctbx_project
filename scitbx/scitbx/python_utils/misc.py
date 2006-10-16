from libtbx import adopt_init_args # XXX backward compatibility 2005_07_29
from libtbx.utils import plural_s # XXX backward compatibility 2006-10-16
from libtbx.utils import user_plus_sys_time # XXX backward compatibility 2006-10-16
from libtbx.utils import time_log # XXX backward compatibility 2006-10-16
from libtbx.utils import human_readable_time # XXX backward compatibility 2006-10-16
from libtbx.utils import human_readable_time_as_seconds # XXX backward compatibility 2006-10-16
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

def show_total_time(
      out=None,
      show_micro_seconds_per_bytecode_instruction=True):
  if (out == None): out = sys.stdout
  total_time = user_plus_sys_time().get()
  try: python_ticker = sys.gettickeraccumulation()
  except AttributeError: pass
  else:
    print >> out, "Time per interpreted Python bytecode instruction:",
    print >> out, "%.3f micro seconds" % (total_time / python_ticker * 1.e6)
  print >> out, "Total CPU time: %.2f %s" % human_readable_time(total_time)

class line_feeder(object):

  def __init__(self, f):
    self.f = iter(f)
    self.eof = False

  def __iter__(self):
    return self

  def next(self):
    if (not self.eof):
      try:
        return self.f.next()[:-1]
      except StopIteration:
        self.eof = True
    return ""

  def next_non_empty(self):
    while 1:
      result = self.next()
      if (self.eof or len(result.strip()) != 0):
        return result

class input_with_prompt(object):

  def __init__(self, prompt, tracebacklimit=0):
    try: import readline
    except: pass
    try: self.previous_tracebacklimit = sys.tracebacklimit
    except: self.previous_tracebacklimit = None
    if (tracebacklimit is not None):
      sys.tracebacklimit = tracebacklimit
    self.input = raw_input(prompt)

  def __del__(self):
    if (self.previous_tracebacklimit is None):
      del sys.tracebacklimit
    else:
      sys.tracebacklimit = self.previous_tracebacklimit

def get_caller_name(n_back=2):
  try: raise Exception
  except:
    t = sys.exc_info()[2]
    f = t.tb_frame
    for i in xrange(n_back):
      if (f.f_back is None): return None
      f = f.f_back
    return f.f_code.co_name

def exercise():
  from libtbx.test_utils import approx_equal
  time_in_seconds = 1.1
  for i_trial in xrange(55):
    time_in_seconds = time_in_seconds**1.1
    time_units, time_unit = human_readable_time(
      time_in_seconds=time_in_seconds)
    assert approx_equal(
      human_readable_time_as_seconds(time_units, time_unit), time_in_seconds)
  print "OK"

if (__name__ == "__main__"):
  exercise()
