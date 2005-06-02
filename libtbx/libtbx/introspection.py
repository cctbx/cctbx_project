from __future__ import division
import sys, os

def frame_object(frames_back=0):
  try: raise Exception
  except:
    t = sys.exc_info()[2]
    f = t.tb_frame
    for i in xrange(frames_back+1):
      if (f.f_back is None): break
      f = f.f_back
    return f

def varnames(frames_back=0):
  return frame_object(frames_back=frames_back+1).f_code.co_varnames

def adopt_init_args(exclusions=[], prefix="", frames_back=0):
  frame = frame_object(frames_back=frames_back+1)
  varnames = frame.f_code.co_varnames
  exclusions.append(varnames[0]) # self
  init_locals = frame.f_locals
  self = init_locals[varnames[0]]
  for varname in varnames:
    if (varname not in exclusions):
      setattr(self, prefix+varname, init_locals[varname])
  if ("__init__varnames__" not in exclusions):
    setattr(self, "__init__varnames__", varnames)

class caller_location:

  def __init__(self, frames_back=0):
    f = frame_object(frames_back=frames_back+1)
    self.file_name = f.f_code.co_filename
    self.line_number = f.f_lineno

  def __str__(self):
    return "%s(%d)" % (self.file_name, self.line_number)

def check_point(frames_back=0):
  print caller_location(frames_back=frames_back+1)
  sys.stdout.flush()

try:
  _proc_status = "/proc/%d/status" % os.getpid()
except AttributeError:
  _proc_status = None

class virtual_memory_info:

  def __init__(self):
    try:
      self.proc_status = open(_proc_status).read()
    except IOError:
      self.proc_status = None

  def get_bytes(self, vm_key):
    if (self.proc_status is None):
      return None
    try:
      i = self.proc_status.index(vm_key)
    except ValueError:
      return None
    flds = self.proc_status[i:].split(None, 3)
    unit = flds[2].upper()
    if (len(flds) < 3 or unit not in ["KB", "MB"]):
      return None
    try:
      result = int(flds[1])
    except ValueError:
      return None
    result *= 1024
    if (unit == "MB"):
      result *= 1024
    return result

  def virtual_memory_size(self):
    return self.get_bytes('VmSize:')

  def resident_set_size(self):
    return self.get_bytes('VmRSS:')

  def stack_size(self):
    return self.get_bytes('VmStk:')

  def show(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    def size_as_string(sz):
      if (sz is None): return "unknown"
      result = []
      while (sz > 0):
        if (sz >= 1000):
          result.insert(0, "%03d" % (sz % 1000))
          sz //= 1000
        else:
          result.insert(0, "%d" % sz)
          break
      return ",".join(result)
    vms = size_as_string(self.virtual_memory_size())
    rss = size_as_string(self.resident_set_size())
    sts = size_as_string(self.stack_size())
    fmt = "%%%ds" % max(len(vms), len(rss), len(sts))
    print >> out, prefix+"Virtual memory size:", fmt % vms
    print >> out, prefix+"Resident set size:  ", fmt % rss
    print >> out, prefix+"Stack size:         ", fmt % sts

if (__name__ == "__main__"):
  virtual_memory_info().show()
