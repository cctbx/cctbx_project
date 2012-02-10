from libtbx.str_utils import show_string
from libtbx.math_utils import ifloor
from libtbx import Auto
from cStringIO import StringIO
import traceback
import os
import sys

_problem_cache = Auto

def detect_problem():
  global _problem_cache
  if (_problem_cache is Auto):
    import os
    if (os.name == "nt"):
      _problem_cache = "libtbx.easy_mp: Windows is not a supported platform."
    else:
      import libtbx.utils
      _problem_cache = libtbx.utils.detect_multiprocessing_problem()
  return _problem_cache

def enable_multiprocessing_if_possible (nproc=Auto, log=None) :
  if (nproc == 1) or (nproc == 0) :
    return 1
  if (log is None) :
    from libtbx.utils import null_out
    log = null_out()
  problems = detect_problem()
  if (problems is not None) and (problems is not Auto) :
    if (nproc is Auto) or (nproc is None) :
      return 1
    else :
      from libtbx.utils import Sorry
      raise Sorry("%s.  Please use nproc=1 or nproc=Auto." % str(problems))
  else :
    print >> log, """
 ******************************************************************
 INFO: Some parts of this job will make use of multiple processors:
 ******************************************************************

   nproc = %s

 Please ask your system administrator for advice about this, in particular if
 you run this job through a queuing system.
""" % str(nproc)
    return nproc

def get_processes (processes) :
  if (processes is None) or (processes is Auto) :
    from libtbx import introspection
    auto_adjust = (processes is Auto)
    processes = introspection.number_of_processors()
    if (auto_adjust) :
      processes = max(ifloor(processes - os.getloadavg()[0]), 1)
  else :
    assert (processes > 0)
  return processes

from weakref import WeakValueDictionary as _
fixed_func_registry = _()

class fixed_func_proxy(object):

  def __init__(self, key, func):
    self.key = key
    fixed_func_registry[key] = func

  def __call__(self, arg):
    key = self.key
    func = fixed_func_registry[key]
    assert func is not None
    return func(arg)

from itertools import count as _
fixed_func_registry_key_generator = _()

try: # cannot use detect_problem() here (hangs in pool.map())
  from multiprocessing.pool import Pool as multiprocessing_Pool
except Exception:
  multiprocessing_Pool = object

class Pool(multiprocessing_Pool):

  def __init__(self,
        processes=None,
        initializer=None,
        initargs=(),
        maxtasksperchild=None,
        fixed_func=None):
    if (multiprocessing_Pool is object):
      mp_problem = detect_problem()
      assert mp_problem is not None
      raise RuntimeError(mp_problem)
    self.processes = get_processes(processes)
    if (fixed_func is not None):
      key = fixed_func_registry_key_generator.next()
      self.fixed_func_proxy = fixed_func_proxy(key, fixed_func)
    else:
      self.fixed_func_proxy = None
    init = super(Pool, self).__init__
    if (maxtasksperchild is None):
      init(
        processes=processes,
        initializer=initializer,
        initargs=initargs)
    else:
      init(
        processes=processes,
        initializer=initializer,
        initargs=initargs,
        maxtasksperchild=maxtasksperchild) # New in Python 2.7

  def map_fixed_func(self, iterable, chunksize=None):
    '''
    Uses fixed_func as passed to the constructor.
    Avoids repeated pickling/unpickling of func, which can be rate-limiting
    if func is large and the amount of work per call is relatively small.
    '''
    assert self.fixed_func_proxy is not None
    return self.map(
      func=self.fixed_func_proxy,
      iterable=iterable,
      chunksize=chunksize)

def show_caught_exception(index, arg):
  print "CAUGHT EXCEPTION: (argument #%d)" % index
  try:
    r = repr(arg)
  except: # intentional
    pass
  else:
    if (len(r) > 256):
      r = r[:127] + "..." + r[-126:]
    print "ARGUMENT LEADING TO EXCEPTION:", r
  traceback.print_exc(file=sys.stdout)

class func_wrapper_default(object):

  def __init__(O, func, buffer_stdout_stderr=False):
    O.func = func
    O.buffer_stdout_stderr = buffer_stdout_stderr

  def __call__(O, index_and_arg):
    assert len(index_and_arg) == 2
    index, arg = index_and_arg
    if (O.buffer_stdout_stderr):
      sys.stderr = sys.stdout = sio = StringIO()
    try:
      result = O.func(arg)
    except: # intentional
      result = None
      show_caught_exception(index, arg)
    if (O.buffer_stdout_stderr):
      return (sio.getvalue(), result)
    return result

class func_wrapper_sub_directories(object):

  def __init__(O, func, sub_name_format="mp%03d", makedirs_mode=0777):
    assert isinstance(sub_name_format, str)
    assert len(sub_name_format % 0) != 0
    O.func = func
    O.sub_name_format = sub_name_format
    O.makedirs_mode = makedirs_mode

  def __call__(O, index_and_arg):
    assert len(index_and_arg) == 2
    index, arg = index_and_arg
    sub_name = O.sub_name_format % index
    op = os.path
    if (op.exists(sub_name)):
      print "ERROR: sub-directory exists already: %s" % show_string(sub_name)
      return
    try:
      os.makedirs(sub_name, mode=O.makedirs_mode)
    except: # intentional
      print "ERROR: cannot create sub-directory: %s" % show_string(sub_name)
      return
    if (not op.isdir(sub_name)):
      print "ERROR: failure creating sub-directory: %s" % show_string(sub_name)
      return
    try:
      os.chdir(sub_name)
    except: # intentional
      print "ERROR: cannot chdir to sub-directory: %s" % show_string(sub_name)
      return
    try:
      _ = open("log", "w")
    except: # intentional
      print "ERROR: cannot open file: %s" % show_string(
        op.join(sub_name, "log"))
      return
    sys.stderr = sys.stdout = _
    try:
      result = O.func(arg)
    except: # intentional
      result = None
      show_caught_exception(index, arg)
    return result

def pool_map(
      processes=None,
      initializer=None,
      initargs=(),
      maxtasksperchild=None,
      func=None,
      fixed_func=None,
      iterable=None,
      args=None,
      chunksize=None,
      func_wrapper="default",
      index_args=True,
      log=None,
      buffer_stdout_stderr=False):
  assert [func, fixed_func].count(None) == 1
  assert [iterable, args].count(None) == 1
  if (isinstance(func_wrapper, str)):
    if (func_wrapper == "default"):
      func_wrapper = func_wrapper_default
    elif (func_wrapper == "sub_directories"):
      func_wrapper = func_wrapper_sub_directories
    else:
      raise RuntimeError("Unknown func_wrapper keyword: %s" % func_wrapper)
  if (func_wrapper is not None):
    if (func is not None):
      if (buffer_stdout_stderr):
        func = func_wrapper(func, buffer_stdout_stderr)
      else:
        func = func_wrapper(func)
    else:
      if (buffer_stdout_stderr):
        fixed_func = func_wrapper(fixed_func, buffer_stdout_stderr)
      else:
        fixed_func = func_wrapper(fixed_func)
  processes = get_processes(processes)
  if (args is not None):
    iterable = args
    if (processes is not None):
      processes = min(processes, len(args))
  if (index_args):
    iterable = enumerate(iterable)
  if (log is not None):
    print >> log, "multiprocessing pool size:", processes
    flush = getattr(log, "flush", None)
    if (flush is not None):
      flush()
    import time
    time_start = time.time()
  pool = Pool(
    processes=processes,
    initializer=initializer,
    initargs=initargs,
    maxtasksperchild=maxtasksperchild,
    fixed_func=fixed_func)
  try:
    if (func is not None):
      result = pool.map(func=func, iterable=iterable, chunksize=chunksize)
    else:
      result = pool.map_fixed_func(iterable=iterable, chunksize=chunksize)
  finally:
    pool.close()
    pool.join()
  if (log is not None):
    from libtbx.utils import show_wall_clock_time
    show_wall_clock_time(seconds=time.time()-time_start, out=log)
  return result

del _
