import sys, os

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

  def log(self):
    self.stop()
    return self.report()

  legend = "time_log: label: n accumulation delta average"

  def report(self):
    assert self.timer == None
    return "time_log: %s: %d %.2f %.2f %.2f" % (
      self.label, self.n, self.accumulation,
      self.delta, self.accumulation/max(1,self.n))

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

def import_regular_symbols(dict_target, dict_source):
  for key, value in dict_source.items():
    if (not key.startswith("_") and not key in dict_target):
      dict_target[key] = value

def import_ext(name):
  components = name.split(".")
  if (len(components) > 1):
    __import__(".".join(components[:-1]))
  previous_dlopenflags = None
  if (sys.platform == "linux2"):
    previous_dlopenflags = sys.getdlopenflags()
    sys.setdlopenflags(0x100|0x2)
  mod = __import__(name)
  for comp in components[1:]:
    mod = getattr(mod, comp)
  if (previous_dlopenflags != None):
    sys.setdlopenflags(previous_dlopenflags)
  return mod

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
