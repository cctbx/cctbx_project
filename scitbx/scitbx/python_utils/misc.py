from libtbx import adopt_init_args # XXX backward compatibility 2005_07_29
import sys, os

class store:

  def __init__(self, **kw):
    self.__dict__.update(kw)

class sorted_store:

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

class user_plus_sys_time:

  def __init__(self):
    self.t = self.get()

  def get(self):
    t = os.times()
    return t[0] + t[1]

  def elapsed(self):
    t = self.get()
    d = t - self.t
    return d

  def delta(self):
    t = self.get()
    d = t - self.t
    self.t = t
    return d

class time_log:

  def __init__(self, label):
    self.label = label
    self.accumulation = 0
    self.n = 0
    self.delta = 0
    self.timer = None

  def start(self):
    self.timer = user_plus_sys_time()
    return self

  def stop(self):
    self.delta = self.timer.delta()
    self.timer = None
    self.accumulation += self.delta
    self.n += 1

  def average(self):
    return self.accumulation / max(1,self.n)

  def log(self):
    self.stop()
    return self.report()

  def log_elapsed(self, local_label):
    return "time_log: %s: %.2f elapsed %s" % (
      self.label, self.timer.elapsed(), local_label)

  legend = "time_log: label: n accumulation delta average"

  def report(self):
    assert self.timer is None
    return "time_log: %s: %d %.2f %.3g %.3g" % (
      self.label, self.n, self.accumulation,
      self.delta, self.average())

def human_readable_time(time_in_seconds):
  time_units = time_in_seconds
  time_unit = "seconds"
  if (time_units > 120):
    time_units /= 60
    time_unit = "minutes"
    if (time_units > 120):
      time_units /= 60
      time_unit = "hours"
      if (time_units > 48):
        time_units /= 24
        time_unit = "days"
  return time_units, time_unit

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

class line_feeder:

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

class input_with_prompt:

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

def plural_s(n):
  if (n == 1): return n, ""
  return n, "s"

def get_caller_name(n_back=2):
  try: raise Exception
  except:
    t = sys.exc_info()[2]
    f = t.tb_frame
    for i in xrange(n_back):
      if (f.f_back is None): return None
      f = f.f_back
    return f.f_code.co_name
