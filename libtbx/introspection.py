from __future__ import division
from libtbx import Auto
from libtbx import group_args
import sys, os
op = os.path

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
  if (out == "return_lines"):
    return lines
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

kb_exponents = {
  "KB": 1,
  "MB": 2,
  "GB": 3,
  "TB": 4,
  "PB": 5}

class proc_file_reader(object):

  def get_bytes(self, vm_key):
    if (self.proc_status is None):
      return None
    try:
      i = self.proc_status.index(vm_key)
    except ValueError:
      return None
    flds = self.proc_status[i:].split(None, 3)
    if (len(flds) < 3):
      return None
    exponent = kb_exponents.get(flds[2].upper())
    try:
      num = int(flds[1])
    except ValueError:
      return None
    return num * 1024**exponent

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

  def current_max_sizes_legend(self):
    return ("Virtual memory", "Resident set", "Stack")

  def current_max_sizes(self):
    return group_args(
      virtual_memory=self.max_virtual_memory_size,
      resident_set=self.max_resident_set_size,
      stack=self.max_stack_size)

def get_mac_os_x_memory_total():
  cmd = "/usr/sbin/system_profiler"
  if (not op.isfile(cmd)):
    return None
  from libtbx import easy_run
  for line in easy_run.fully_buffered(
                command=cmd+" SPHardwareDataType").stdout_lines:
    line = line.strip()
    if (not line.startswith("Memory:")):
      continue
    flds = line.split()
    if (len(flds) != 3):
      continue
    try:
      num = int(flds[1])
    except ValueError:
      continue
    if (num <= 0):
      continue
    exponent = kb_exponents.get(flds[2].upper())
    if (exponent is None):
      continue
    return num * 1024**exponent
  return None

class machine_memory_info(proc_file_reader):

  def __init__(self):
    self.proc_status = None
    self._memory_total = Auto
    self._memory_free = Auto
    if (op.isfile("/proc/meminfo")):
      try:
        self.proc_status = open("/proc/meminfo").read()
      except IOError:
        pass
    else:
      self._memory_total = get_mac_os_x_memory_total()
      self._memory_free = None

  def memory_total(self):
    result = self._memory_total
    if (result is Auto):
      result = self._memory_total = self.get_bytes("MemTotal:")
    return result

  def memory_free(self):
    result = self._memory_free
    if (result is Auto):
      result = self._memory_free = self.get_bytes("MemFree:")
    return result

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
      if (op.isfile(cpuinfo)):
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
      if (op.isfile(cmd)):
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
      if (op.isfile(cmd)):
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

class method_debug_log (object) :
  """For Python 2.4 or greater.  Use an instance of this class as a
  decorator for class methods, and it will print the call signature and
  call location before the method is executed.

  Example:
  debug = libtbx.introspection.method_debug_log()
  class a (object) :
    @debug
    def foo (self, x) :
      print x

  def main () :
    my_object = a()
    a.foo(1)
  main()

  Running this results in the following output when LIBTBX_DEBUG_LOG is set:
a.foo(1) @ test.py(13) main
1
  """
  def __init__ (self) :
    self.debug = False
    if os.environ.get("LIBTBX_DEBUG_LOG") is not None :
      self.debug = True

  def __call__ (self, f) :
    def log_wrapper (O, *args, **kwds) :
      if self.debug :
        _args = list(args)
        _kwds = dict(kwds)
        str_args = ", ".join([ str(arg) for arg in _args ])
        str_kwds = ", ".join([ "%s=%s" % (kwd, str(val))
                               for kwd, val in _kwds.iteritems() ])
        call_signature = []
        if str_args != "" : call_signature.append(str_args)
        if str_kwds != "" : call_signature.append(str_kwds)
        caller = sys._getframe(1)
        print "%s.%s(%s) @ %s(%d) %s" % (O.__class__.__name__, f.__name__,
          ", ".join(call_signature), caller.f_code.co_filename,
          caller.f_lineno, caller.f_code.co_name)
        sys.stdout.flush()
      return f(O, *args, **kwds)
    return log_wrapper


class current_process_status(object):
  """
  An interface to the *NIX utility 'ps' to get info on the current process
  (only tested on MacOS X till further notice, so beware dragons)

  SYNOPSIS:
    ps = current_process_status() # <1>
    .....
    print ps['RSS'] # resident size at the time of <1>
    .....
    ps.refresh() # <2>
    .....
    print ps['%CPU'] # % CPU at the time of <2>
  """
  conversions = {
    'RSS': int,
    'VSZ': int,
    '%CPU': float,
    }

  def __init__(self):
    self.id = str(os.getpid())
    self.refresh()

  def refresh(self):
    import subprocess
    ps = subprocess.Popen(
      args=['/bin/ps', 'cux'],
      stderr=subprocess.STDOUT,
      stdout=subprocess.PIPE)
    cols = dict(
      [ (field, i) for i, field in enumerate(ps.stdout.readline().split()) ])
    i_pid = cols['PID']
    for li in ps.stdout:
      field = li.split()
      if field[i_pid] == self.id: break
    else:
      return
    self.field = dict(
      [ (name, field[i_col]) for name, i_col in cols.iteritems() ])
    for name, conv in self.conversions.iteritems():
      self.field[name] = conv(self.field[name])

  def __getitem__(self, field_name):
    return self.field[field_name]


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
