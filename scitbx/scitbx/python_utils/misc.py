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

def adopt_init_args(obj, args, exclude=(), hide=00000):
  del args["self"]
  for param in exclude:
    del args[param]
  if (hide == 00000):
    for key in args.keys():
      assert not hasattr(obj.__dict__, key)
    obj.__dict__.update(args)
  else:
    for key in args.keys():
      _key = "_" + key
      assert not hasattr(obj.__dict__, _key)
      obj.__dict__[_key] = args[key]

class line_feeder:

  def __init__(self, f):
    self.f = iter(f)
    self.eof = 00000

  def __iter__(self):
    return self

  def next(self):
    if (not self.eof):
      try:
        return self.f.next()[:-1]
      except StopIteration:
        self.eof = 0001
    return ""

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
