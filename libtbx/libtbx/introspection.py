from __future__ import division
from libtbx.str_utils import size_as_string_with_commas
import sys, os

def varnames(frames_back=0):
  f_code = sys._getframe(frames_back+1).f_code
  return f_code.co_varnames[:f_code.co_argcount]

class caller_location(object):

  def __init__(self, frames_back=0):
    f = sys._getframe(frames_back+1)
    self.file_name = f.f_code.co_filename
    self.line_number = f.f_lineno

  def __str__(self):
    return "%s(%d)" % (self.file_name, self.line_number)

def check_point(frames_back=0):
  print caller_location(frames_back=frames_back+1)
  sys.stdout.flush()

def show_stack(max_frames_back=None, out=None):
  if (out is None): out = sys.stdout
  try:
    frames_back = 0
    while True:
      if (max_frames_back is not None and frames_back == max_frames_back):
        break
      f = sys._getframe(frames_back+1)
      print >> out, "show_stack(%d): %s(%d) %s" % (
        frames_back, f.f_code.co_filename, f.f_lineno, f.f_code.co_name)
      frames_back += 1
  except ValueError:
    print >> out

def print_trace(frame, event, arg):
  if (event == "line"):
    sys.stderr.flush()
    print "%s(%d)" % (frame.f_code.co_filename, frame.f_lineno)
    sys.stdout.flush()
  return print_trace

def start_print_trace():
  sys.settrace(print_trace)

def stop_print_trace():
  sys.settrace(None)

class proc_file_reader(object):

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

try:
  _proc_status = "/proc/%d/status" % os.getpid()
except AttributeError:
  _proc_status = None

class virtual_memory_info(proc_file_reader):

  def __init__(self):
    try:
      self.proc_status = open(_proc_status).read()
    except IOError:
      self.proc_status = None

  def virtual_memory_size(self):
    return self.get_bytes('VmSize:')

  def resident_set_size(self):
    return self.get_bytes('VmRSS:')

  def stack_size(self):
    return self.get_bytes('VmStk:')

  def show(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    vms = size_as_string_with_commas(self.virtual_memory_size())
    rss = size_as_string_with_commas(self.resident_set_size())
    sts = size_as_string_with_commas(self.stack_size())
    fmt = "%%%ds" % max(len(vms), len(rss), len(sts))
    print >> out, prefix+"Virtual memory size:", fmt % vms
    print >> out, prefix+"Resident set size:  ", fmt % rss
    print >> out, prefix+"Stack size:         ", fmt % sts

  def show_if_available(self, out=None, prefix=""):
    if (self.proc_status is not None):
      self.show(out=out, prefix=prefix)

try:
  _proc_meminfo = "/proc/meminfo"
except AttributeError:
  _proc_meminfo = None

class machine_memory_info(proc_file_reader):

  def __init__(self):
    try:
      self.proc_status = open(_proc_meminfo).read()
    except IOError:
      self.proc_status = None

  def memory_total(self):
    return self.get_bytes("MemTotal:")

  def memory_free(self):
    return self.get_bytes("MemFree:")

  def show(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    mt = size_as_string_with_commas(self.memory_total())
    mf = size_as_string_with_commas(self.memory_free())
    fmt = "%%%ds" % max(len(mt), len(mf))
    print >> out, prefix+"Memory total: ", fmt % mt
    print >> out, prefix+"Memory free:  ", fmt % mf

if (__name__ == "__main__"):
  def exercise_varnames(a, b, c):
    d = 0
    return varnames()
  assert exercise_varnames(1,2,3) == ("a", "b", "c")
  virtual_memory_info().show()
  machine_memory_info().show()
  assert str(caller_location()).find("introspection") > 0
  print "OK"
