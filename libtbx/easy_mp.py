from __future__ import division
from libtbx.str_utils import show_string
from libtbx.math_utils import ifloor
from libtbx import Auto
from cStringIO import StringIO
import traceback
import os
import sys

_have_maxtasksperchild = (sys.version_info[:2] >= (2,7))

_problem_cache = Auto

def detect_problem():
  """
  Identify situations where multiprocessing will not work as required.
  """
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
  """
  Switch for using multiple CPUs, usually called at the beginning of an app.
  If nproc is Auto or None and we are running Windows, it will be reset to 1.
  """
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

# FIXME should be more flexible on Windows
def get_processes (processes) :
  """
  Determine number of processes dynamically: number of CPUs minus the current
  load average (with a minimum of 1).
  """
  if (processes in [None, Auto]) :
    if (os.name == "nt") or (sys.version_info < (2,6)) :
      return 1
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
  """Implementation detail"""
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
  """Subclass of multiprocessing.Pool, used internally by pool_map."""
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
    if (maxtasksperchild is Auto):
      maxtasksperchild = None
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

class func_wrapper_simple_impl(object):
  """Implementation detail"""

  def __init__(O, options, func):
    O.options = options
    O.func = func

  def __call__(O, index_and_arg):
    assert len(index_and_arg) == 2
    index, arg = index_and_arg
    if (O.options.buffer_stdout_stderr):
      sys.stderr = sys.stdout = sio = StringIO()
    try:
      result = O.func(arg)
    except: # intentional
      result = None
      show_caught_exception(index, arg)
    if (O.options.buffer_stdout_stderr):
      return (sio.getvalue(), result)
    return result

class func_wrapper_simple(object):
  """Implementation detail"""

  def __init__(O, buffer_stdout_stderr=False):
    O.buffer_stdout_stderr = buffer_stdout_stderr

  def wrap(O, func):
    return func_wrapper_simple_impl(options=O, func=func)

class func_wrapper_sub_directories_impl(object):
  """Implementation detail"""

  def __init__(O, options, func):
    O.options = options
    O.func = func

  def __call__(O, index_and_arg):
    assert len(index_and_arg) == 2
    index, arg = index_and_arg
    sub_name = O.options.sub_name_format % index
    op = os.path
    if (op.exists(sub_name)):
      return (
        "sub-directory exists already: %s" % show_string(sub_name),
        None)
    try:
      os.makedirs(sub_name, mode=O.options.makedirs_mode)
    except: # intentional
      return (
        "cannot create sub-directory: %s" % show_string(sub_name),
        None)
    if (not op.isdir(sub_name)):
      return (
        "failure creating sub-directory: %s" % show_string(sub_name),
        None)
    initial_cwd = os.getcwd()
    try:
      try:
        os.chdir(sub_name)
      except: # intentional
        return (
          "cannot chdir to sub-directory: %s" % show_string(sub_name),
          None)
      def sub_log(): return show_string(op.join(sub_name, "log"))
      try:
        log = open("log", "w")
      except: # intentional
        return ("cannot open file: %s" % sub_log(), None)
      initial_out = sys.stdout
      initial_err = sys.stderr
      try:
        sys.stderr = sys.stdout = log
        try:
          result = O.func(arg)
        except: # intentional
          show_caught_exception(index, arg)
          return ("CAUGHT EXCEPTION: %s" % sub_log(), None)
      finally:
        sys.stdout = initial_out
        sys.stderr = initial_err
        log.close()
    finally:
      os.chdir(initial_cwd)
    return (None, result)

class func_wrapper_sub_directories(object):
  """Implementation detail"""

  def __init__(O, sub_name_format="mp%03d", makedirs_mode=0777):
    assert isinstance(sub_name_format, str)
    s = sub_name_format
    if (s.find("%") < 0):
      i = s.find("#")
      if (i >= 0):
        c = s.count("#", i)
        if (i + c == len(s)):
          s = s[:i] + "%0" + str(c) + "d"
        else:
          i = -1
      if (i < 0):
        s += "%03d"
    sub_name_format = s
    assert len(sub_name_format) != 0
    assert len(sub_name_format % 0) != 0
    O.sub_name_format = sub_name_format
    O.makedirs_mode = makedirs_mode

  def wrap(O, func):
    return func_wrapper_sub_directories_impl(options=O, func=func)

def pool_map(
      processes=None,
      initializer=None,
      initargs=(),
      maxtasksperchild=Auto,
      func=None,
      fixed_func=None,
      iterable=None,
      args=None,
      chunksize=Auto,
      func_wrapper="simple",
      index_args=True,
      log=None,
      call_back_for_serial_run=None):
  """
  Parallelized map() using subclassed multiprocessing.Pool.  If func is not
  None, this function essentially calls the Pool's own map method; this means
  that both func and iterable/args must be pickle-able.  If fixed_func is not
  None, it will not be pickled but instead saved as an attribute of the Pool,
  which will be preserved after the fork() call.  Additional features include
  optional redirection of output and automatic process number determination.

  Note that because of the reliance on fork(), this function will run in serial
  on Windows, regardless of how many processors are available.

  :param processes: number of processes to spawn; if None or Auto, the
    get_processes() function will be used.
  :param func: target function (will be pickled)
  :param fixed_func: "fixed" target function, which will be be propagated to
    the child process when forked (instead of pickling)
  :param iterable: argument list
  :param args: same as iterable (alternate keyword)
  :param chunksize: number of arguments to process at once
  """
  assert [func, fixed_func].count(None) == 1
  assert [iterable, args].count(None) == 1
  assert ((call_back_for_serial_run is None) or
          hasattr(call_back_for_serial_run, "__call__"))
  if (isinstance(func_wrapper, str)):
    if (func_wrapper == "simple"):
      func_wrapper = func_wrapper_simple()
    else:
      if (func_wrapper == "buffer_stdout_stderr"):
        func_wrapper = func_wrapper_simple(buffer_stdout_stderr=True)
      elif (func_wrapper == "sub_directories"):
        func_wrapper = func_wrapper_sub_directories()
      elif (func_wrapper.startswith("sub_directories:")):
        func_wrapper = func_wrapper_sub_directories(
          sub_name_format=func_wrapper[16:])
      else:
        raise RuntimeError("Unknown func_wrapper keyword: %s" % func_wrapper)
      if (maxtasksperchild is Auto and _have_maxtasksperchild):
        maxtasksperchild = 1
      if (chunksize is Auto):
        chunksize = 1
  if (func_wrapper is not None):
    wrap = getattr(func_wrapper, "wrap", None)
    if (wrap is None):
      raise RuntimeError("func_wrapper must have a .wrap() method.")
    if (func is not None):
      func = wrap(func)
    else:
      fixed_func = wrap(fixed_func)
  processes = get_processes(processes)
  # XXX since we want to be able to call this function on Windows too, reset
  # processes to 1
  if (os.name == "nt") or (sys.version_info < (2,6)) :
    processes = 1
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
  result = None
  # XXX this allows the function to be used even when parallelization is
  # not enabled or supported, which should keep calling code simpler.
  if (processes == 1) or (os.name == "nt") :
    result = []
    for args in iterable :
      if (func is not None) :
        result.append(func(args))
      else :
        result.append(fixed_func(args))
      if (call_back_for_serial_run is not None) :
        call_back_for_serial_run(result[-1])
  else :
    pool = Pool(
      processes=processes,
      initializer=initializer,
      initargs=initargs,
      maxtasksperchild=maxtasksperchild,
      fixed_func=fixed_func)
    if (chunksize is Auto):
      chunksize = None
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

#-----------------------------------------------------------------------
# application support for parallelization across multiple cores or a
# queuing system (also Unix-only)
parallel_phil_str_base = """
nproc = 1
  .type = int
  .short_caption = Number of processes
  .style = bold renderer:draw_nproc_widget
technology = %s
  .type = choice
  .short_caption = Parallelization method
  .caption = %s
qsub_command = None
  .type = str
  .short_caption = qsub command
"""

parallel_methods = ["*multiprocessing", "sge", "lsf", "pbs", "condor"]
parallel_captions = ["Multiprocessing", "Sun_Grid_Engine", "LSF", "PBS", "Condor"]

parallel_phil_str = parallel_phil_str_base % (
  " ".join(parallel_methods + ["threading"]),
  " ".join(parallel_captions + ["Threading"]))
parallel_phil_str_no_threading = parallel_phil_str_base % (
  " ".join(parallel_methods), " ".join(parallel_captions))

def parallel_map (
    func,
    iterable,
    params=None,
    processes=1,
    method="multiprocessing",
    qsub_command=None,
    asynchronous=False,
    callback=None,
    preserve_order=True) :
  """
  Generic parallel map() implementation for a variety of platforms, including
  the multiprocessing module and supported queuing systems, via the module
  libtbx.queuing_system_utils.scheduling.  This is less flexible than pool_map
  above, since it does not provide a way to use a non-pickleable target
  function, but it provides a consistent API for programs where multiple
  execution methods are desired.  It will also work on Windows (if the method
  is multiprocessing or threading).

  Note that for most applications, the threading method will be constrained
  by the Global Interpreter Lock, therefore multiprocessing is prefered for
  parallelizing across a single multi-core system.

  See Computational Crystallography Newsletter 3:37-42 (2012) for details of
  the underlying method.

  :param func: target function (must be pickleable)
  :param iterable: list of arguments for func
  :param processes: number of processes/threads to start
  :param method: parallelization method (multiprocessing|threading|sge|lsf|pbs)
  :param qsub_command: command to submit queue jobs (optional)
  :param asynchronous: run queue jobs asynchronously
  :returns: a list of result objects
  """
  if (params is not None) :
    method = params.technology
    processes = params.nproc
    qsub_command = params.qsub_command
  from libtbx.queuing_system_utils import scheduling
  from libtbx.queuing_system_utils import processing
  assert ((method in ["multiprocessing", "threading"]) or
          (method in processing.INTERFACE_FOR.keys()))
  assert (callback is None) or (hasattr(callback, "__call__"))
  results = []
  # if we aren't actually going to run multiple processes or use a queuing
  # system, just loop over all arguments and skip the rest of the setup.
  if (processes == 1) and (method in ["threading", "multiprocessing"]) :
    for args in iterable :
      result = func(args)
      if (callback is not None) :
        callback(result)
      results.append(result)
    return results
  units = []
  factory, queue_factory = None, None
  if (method == "multiprocessing") :
    import multiprocessing
    factory = multiprocessing.Process
    # XXX this is essential - using multiprocessing.Queue will not work for
    # large result objects.  explanation here:
    #   http://bugs.python.org/issue8426
    # for reasons which are opaque to me, using the Manager object to create
    # the Queue will circumvent the problem.
    mp_manager = multiprocessing.Manager()
    queue_factory = mp_manager.Queue
  elif (method == "threading") :
    import threading
    import Queue
    factory = threading.Thread
    queue_factory = Queue.Queue
  for i_proc in range(processes) :
    queue = None
    factory_tmp = factory
    if (factory_tmp is None) :
      qhandler_function, evaluator = processing.INTERFACE_FOR[method]
      qhandler = qhandler_function(
        command=qsub_command,
        asynchronous=asynchronous)
      factory_tmp = qhandler.Job
      queue = processing.Queue(identifier="parallel_map")
    else :
      queue = queue_factory()
    assert (not None in [factory_tmp, queue])
    unit = scheduling.ExecutionUnit(
      factory=factory_tmp,
      processor=scheduling.RetrieveProcessor(queue=queue))
    units.append(unit)
  manager = scheduling.Manager(units=units)
  if preserve_order:
    orderer = scheduling.SubmissionOrder
  else:
    orderer = scheduling.FinishingOrder
  parallel_for = scheduling.ParallelForIterator(
    calculations = ( ( func, ( args, ), {} ) for args in iterable ),
    manager = manager,
    )
  for ( params, result_wrapper ) in orderer( parallel_for = parallel_for ) :
    result = result_wrapper()
    results.append(result) # raises exception in case of error
    if (callback is not None) :
      callback(result)
  return results

