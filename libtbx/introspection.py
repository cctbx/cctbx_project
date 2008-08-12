from __future__ import division
from libtbx import Auto
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

def show_stack(
      max_frames_back=None,
      frames_back=0,
      reverse=False,
      out=None,
      prefix=""):
  if (out is None): out = sys.stdout
  lines = []
  try:
    while True:
      if (max_frames_back is not None and frames_back == max_frames_back):
        break
      f = sys._getframe(frames_back+1)
      lines.append(prefix+"show_stack(%d): %s(%d) %s" % (
        frames_back, f.f_code.co_filename, f.f_lineno, f.f_code.co_name))
      frames_back += 1
  except ValueError:
    pass
  if (reverse): lines.reverse()
  for line in lines:
    print >> out, line

def show_stack_true_stderr():
  sys.__stdout__.flush()
  show_stack(out=sys.__stderr__, frames_back=1)
  sys.__stderr__.flush()

def print_trace(frame, event, arg):
  if (event == "line"):
    sys.stderr.flush()
    print "%s(%d)" % (frame.f_code.co_filename, frame.f_lineno)
    sys.stdout.flush()
  return print_trace

def start_print_trace():
  if ("pydoc" in sys.modules):
    from cStringIO import StringIO
    s = StringIO()
    show_stack(out=s)
    for line in s.getvalue().splitlines():
      if (line.find("pydoc.py") >= 0 and line.endswith(" cli")):
        print "libtbx.introspection.start_print_trace():" \
          " pydoc.cli() detected: no tracing"
        return
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

  have_vmpeak = False
  max_virtual_memory_size = 0
  max_resident_set_size = 0
  max_stack_size = 0

  def __init__(self):
    try:
      self.proc_status = open(_proc_status).read()
    except IOError:
      self.proc_status = None

  def virtual_memory_peak_size(self):
    result = self.get_bytes('VmPeak:')
    if (result is not None):
      virtual_memory_info.have_vmpeak = True
    virtual_memory_info.max_virtual_memory_size = max(
    virtual_memory_info.max_virtual_memory_size, result)
    return result

  def virtual_memory_size(self):
    result = self.get_bytes('VmSize:')
    virtual_memory_info.max_virtual_memory_size = max(
    virtual_memory_info.max_virtual_memory_size, result)
    return result

  def resident_set_size(self):
    result = self.get_bytes('VmRSS:')
    virtual_memory_info.max_resident_set_size = max(
    virtual_memory_info.max_resident_set_size, result)
    return result

  def stack_size(self):
    result = self.get_bytes('VmStk:')
    virtual_memory_info.max_stack_size = max(
    virtual_memory_info.max_stack_size, result)
    return result

  def update_max(self):
    if (self.proc_status is not None):
      self.virtual_memory_peak_size()
      self.virtual_memory_size()
      self.resident_set_size()
      self.stack_size()

  def show(self, out=None, prefix="", show_max=False):
    if (out is None): out = sys.stdout
    from libtbx.str_utils import size_as_string_with_commas
    vms = size_as_string_with_commas(self.virtual_memory_size())
    rss = size_as_string_with_commas(self.resident_set_size())
    sts = size_as_string_with_commas(self.stack_size())
    fmt = "%%%ds" % max(len(vms), len(rss), len(sts))
    lvms = prefix + "Virtual memory size:"
    lrss = prefix + "Resident set size:  "
    lsts = prefix + "Stack size:         "
    if (not show_max):
      print >> out, lvms, fmt % vms
      print >> out, lrss, fmt % rss
      print >> out, lsts, fmt % sts
    else:
      self.virtual_memory_peak_size()
      vmi = virtual_memory_info
      max_vms = size_as_string_with_commas(vmi.max_virtual_memory_size)
      max_rss = size_as_string_with_commas(vmi.max_resident_set_size)
      max_sts = size_as_string_with_commas(vmi.max_stack_size)
      max_fmt = "%%%ds" % max(len(max_vms), len(max_rss), len(max_sts))
      if (vmi.have_vmpeak):
        vms_what_max = "    exact max:"
      else:
        vms_what_max = "  approx. max:"
      print >> out, lvms, fmt % vms, vms_what_max,     max_fmt % max_vms
      print >> out, lrss, fmt % rss, "  approx. max:", max_fmt % max_rss
      print >> out, lsts, fmt % sts, "  approx. max:", max_fmt % max_sts

  def show_if_available(self, out=None, prefix="", show_max=False):
    if (self.proc_status is not None):
      self.show(out=out, prefix=prefix, show_max=show_max)

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
    from libtbx.str_utils import size_as_string_with_commas
    mt = size_as_string_with_commas(self.memory_total())
    mf = size_as_string_with_commas(self.memory_free())
    fmt = "%%%ds" % max(len(mt), len(mf))
    print >> out, prefix+"Memory total: ", fmt % mt
    print >> out, prefix+"Memory free:  ", fmt % mf

_number_of_processors = Auto

def number_of_processors(return_value_if_unknown=None):
  global _number_of_processors
  if (_number_of_processors is Auto):
    _number_of_processors = None
    try: import boost.python
    except ImportError: pass
    else:
      n = boost.python.ext.number_of_processors()
      if (n != 0):
        _number_of_processors = n
    if (_number_of_processors is None):
      cpuinfo = "/proc/cpuinfo" # Linux
      if (os.path.isfile(cpuinfo)):
        n = 0
        for line in open(cpuinfo).read().splitlines():
          if (not line.startswith("processor")): continue
          line = line[9:].replace(" ", "").replace("\t", "")
          if (not line.startswith(":")): continue
          n += 1
        if (n != 0):
          _number_of_processors = n
    if (_number_of_processors is None):
      cmd = "/usr/sbin/system_profiler" # Mac OS X
      if (os.path.isfile(cmd)):
        keys = [
          "Total Number Of Cores: ",
          "Number Of CPUs: ",
          "Number Of Processors: "]
        ns = [None] * len(keys)
        from libtbx import easy_run
        for line in easy_run.fully_buffered(
                      command=cmd+" SPHardwareDataType").stdout_lines:
          line = line.strip()
          for i,key in enumerate(keys):
            if (line.startswith(key)):
              try: n = int(line[len(key):])
              except ValueError: continue
              if (n > 0 and ns[i] is None):
                ns[i] = n
        for n in ns:
          if (n is not None):
            _number_of_processors = n
            break
    if (_number_of_processors is None):
      n = os.environ.get("NUMBER_OF_PROCESSORS") # Windows
      if (n is not None):
        try: n = int(n)
        except ValueError: pass
        else: _number_of_processors = n
    if (_number_of_processors is None):
      cmd = "/sbin/hinv" # IRIX
      if (os.path.isfile(cmd)):
        from libtbx import easy_run
        for line in easy_run.fully_buffered(command=cmd).stdout_lines:
          if (line.endswith(" Processors")):
            try: n = int(line.split(" ", 1)[0])
            except ValueError: continue
            if (n > 0):
              _number_of_processors = n
              break
  if (_number_of_processors is not None):
    return _number_of_processors
  return return_value_if_unknown

if (__name__ == "__main__"):
  def exercise_varnames(a, b, c):
    d = 0
    return varnames()
  assert exercise_varnames(1,2,3) == ("a", "b", "c")
  #
  assert _number_of_processors is Auto
  print "number of processors:", number_of_processors(
    return_value_if_unknown="unknown")
  assert _number_of_processors is not Auto
  assert number_of_processors() is _number_of_processors
  #
  virtual_memory_info().show()
  buffer = [0]*1000000
  virtual_memory_info().update_max()
  del buffer
  virtual_memory_info().show(show_max=True)
  #
  machine_memory_info().show()
  #
  assert str(caller_location()).find("introspection") > 0
  #
  print "OK"
