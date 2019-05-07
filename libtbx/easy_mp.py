from __future__ import absolute_import, division, print_function
from six.moves import range
from libtbx.str_utils import show_string
from libtbx.math_utils import ifloor
from libtbx import Auto
from six.moves import cStringIO as StringIO
import traceback
import os
import sys

_have_maxtasksperchild = (sys.version_info[:2] >= (2,7))

_problem_cache = Auto

# Patch Python 2.7 multiprocessing module to avoid unnecessary file operations
# on non-Windows systems.
sys.modules['multiprocessing.random'] = None

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

def enable_multiprocessing_if_possible(nproc=Auto, log=None):
  """
  Switch for using multiple CPUs with the pool_map function, usually called at
  the beginning of an app.  If nproc is Auto or None and we are running
  Windows, it will be reset to 1.

  :param nproc: default number of processors to use
  :returns: number of processors to use (None or Auto means automatic)
  """
  if (nproc == 1) or (nproc == 0):
    return 1
  if (log is None):
    from libtbx.utils import null_out
    log = null_out()
  problems = detect_problem()
  if (problems is not None) and (problems is not Auto):
    if (nproc is Auto) or (nproc is None):
      return 1
    else :
      from libtbx.utils import Sorry
      raise Sorry("%s.  Please use nproc=1 or nproc=Auto." % str(problems))
  else :
    print("""
 ******************************************************************
 INFO: Some parts of this job will make use of multiple processors:
 ******************************************************************

   nproc = %s

 Please ask your system administrator for advice about this, in particular if
 you run this job through a queuing system.
""" % str(nproc), file=log)
    return nproc

# FIXME should be more flexible on Windows
def get_processes(processes):
  """
  Determine number of processes dynamically: number of CPUs minus the current
  load average (with a minimum of 1).

  :param processes: default number of processes (may be None or Auto)
  :returns: actual number of processes to use
  """
  if (processes in [None, Auto]):
    if (os.name == "nt") or (sys.version_info < (2,6)):
      return 1
    from libtbx import introspection
    auto_adjust = (processes is Auto)
    processes = introspection.number_of_processors()
    if (auto_adjust):
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
      key = next(fixed_func_registry_key_generator)
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
  print("CAUGHT EXCEPTION: (argument #%d)" % index)
  try:
    r = repr(arg)
  except: # intentional
    pass
  else:
    if (len(r) > 256):
      r = r[:127] + "..." + r[-126:]
    print("ARGUMENT LEADING TO EXCEPTION:", r)
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

  def __init__(O, sub_name_format="mp%03d", makedirs_mode=0o777):
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

  Examples
  --------
  >>> def f(x):
  ...   return some_long_running_method(x)
  ...
  >>> args = range(1000)
  >>> result = easy_mp.pool_map(
  ...   func=f,
  ...   args=args)
  ...
  >>> print len(result)
  ... 1000

  >>> class f_caller(object):
  ...   def __init__(self, non_pickleable_object):
  ...     self._obj = non_pickleable_object
  ...   def __call__(self, x):
  ...     return some_long_running_method(x, self._obj)
  ...
  >>> args = range(1000)
  >>> f = f_caller(processed_pdb_file)
  >>> result = easy_mp.pool_map(
  ...   fixed_func=f,
  ...   args=args)
  ...
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
  if (os.name == "nt") or (sys.version_info < (2,6)):
    processes = 1
  if (args is not None):
    iterable = args
    if (processes is not None):
      processes = min(processes, len(args))
  if (index_args):
    iterable = enumerate(iterable)
  if (log is not None):
    print("multiprocessing pool size:", processes, file=log)
    flush = getattr(log, "flush", None)
    if (flush is not None):
      flush()
    import time
    time_start = time.time()
  result = None
  # XXX this allows the function to be used even when parallelization is
  # not enabled or supported, which should keep calling code simpler.
  if (processes == 1) or (os.name == "nt"):
    result = []
    for args in iterable :
      if (func is not None):
        result.append(func(args))
      else :
        result.append(fixed_func(args))
      if (call_back_for_serial_run is not None):
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

parallel_methods = ["*multiprocessing", "sge", "lsf", "pbs", "condor", "slurm"]
parallel_captions = ["Multiprocessing", "Sun_Grid_Engine", "LSF", "PBS", "Condor", "SLURM"]

parallel_phil_str = parallel_phil_str_base % (
  " ".join(parallel_methods + ["threading"]),
  " ".join(parallel_captions + ["Threading"]))
parallel_phil_str_no_threading = parallel_phil_str_base % (
  " ".join(parallel_methods), " ".join(parallel_captions))


def single_argument(func):

  return func


class posiargs(object):

  def __init__(self, func):

    self.func = func


  def __call__(self, arg):

    return self.func( *arg )


class kwargs(object):

  def __init__(self, func):

    self.func = func


  def __call__(self, arg):

    return self.func( **arg )


class posi_and_kwargs(object):

  def __init__(self, func):

    self.func = func


  def __call__(self, arg):

    ( args, kwargs ) = arg
    return self.func( *args, **kwargs )


def parallel_map(
    func,
    iterable,
    iterable_type = single_argument,
    params=None,
    processes=1,
    method="multiprocessing",
    qsub_command=None,
    asynchronous=True,
    callback=None,
    preserve_order=True,
    preserve_exception_message=False,
    use_manager=False,
    stacktrace_handling = "ignore"):
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
  :param preserve_exception_message: keeps original exception message
  :returns: a list of result objects
  """
  if (params is not None):
    method = params.technology
    processes = params.nproc
    qsub_command = params.qsub_command

  from libtbx.utils import Sorry
  from libtbx.scheduling import SetupError

  if processes == 1 and "LIBTBX_FORCE_PARALLEL" not in os.environ:
    from libtbx.scheduling import mainthread
    creator = mainthread.creator

  else:
    from libtbx.scheduling import philgen
    from libtbx.scheduling import job_scheduler

    if method == "threading":
      technology = philgen.threading(
        capture_exception = preserve_exception_message,
        )
      jfactory = technology.jfactory()
      qfactory = technology.qfactory()[0]
      capacity = job_scheduler.limited(
        njobs = get_processes( processes = processes )
        )

    elif method == "multiprocessing":
      technology = philgen.multiprocessing(
        capture_stderr = preserve_exception_message,
        qtype = philgen.mp_managed_queue if use_manager else philgen.mp_fifo_queue,
        )
      jfactory = technology.jfactory()
      qfactory = technology.qfactory()[0]
      capacity = job_scheduler.limited(
        njobs = get_processes( processes = processes )
        )

    else:
      technology = philgen.cluster(
        asynchronous = asynchronous,
        capture_stderr = preserve_exception_message,
        )

      assert method in technology.platforms # perhaps something less intrusive

      try:
        jfactory = technology.jfactory( platform = method, command = qsub_command )

      except SetupError as e:
        raise Sorry(e)

      from libtbx.scheduling import file_queue
      qfactory = file_queue.qfactory()

      if processes is Auto or processes is None:
        capacity = job_scheduler.unlimited

      else:
        capacity = job_scheduler.limited( njobs = processes )

    creator = job_scheduler.creator(
      job_factory = jfactory,
      queue_factory = qfactory,
      capacity = capacity,
      )

  import libtbx.scheduling

  if stacktrace_handling == "ignore":
    sthandler = libtbx.scheduling.ignore

  elif stacktrace_handling == "excepthook":
    sthandler = libtbx.scheduling.excepthook

  elif stacktrace_handling == "decorate":
    sthandler = libtbx.scheduling.decorate

  else:
    raise Sorry("Unknown stacktrace handling method: %s" % stacktrace_handling)

  from libtbx.scheduling import parallel_for

  if callback is None:
    callback = lambda r: None

  results = []

  with libtbx.scheduling.holder( creator = creator, stacktrace = sthandler ) as manager:
    adfunc = iterable_type( func )

    try:
      pfi = parallel_for.iterator(
        calculations = ( ( adfunc, ( args, ), {} ) for args in iterable ),
        manager = manager,
        keep_input_order = preserve_order,
        )

      for ( calc, res ) in pfi:
        result = res()
        results.append( result )
        callback( result )

    except SetupError as e:
      raise Sorry(e)

    manager.shutdown()
    manager.join()

  return results



def multi_core_run( myfunction, argstuples, nproc ):
  """
  Run myfunction on many cpu cores using multiprocessing.
  A simplified version of parallel_map() above

  myfunction: name of the function to be parallelised,
  argstuples: list of tuples of associated input arguments,
  nproc: number of cores to run on.

  Both myfunction and its input arguments must be pickleable.

  Output is an iterator where each element is a tuple that contains:
  a tuple of arguments for one particular calculation with myfunction,
  the result of this calculation,
  the stacktrace if myfunction crashed

  Example:

  # define RunMyJob() in a file testjob.py
  def RunMyJob( foo, bar):
    import math
    return math.sqrt(foo)/bar

  # then one can start RunMyJob in parallel like:
  >>> import testjob
  >>> from libtbx import easy_mp
  >>>
  >>> argstuples = [( 3, 4), (2, 3) ] # define tuples of arguments
  >>>
  >>> for args, res, errstr in easy_mp.multi_core_run( testjob.RunMyJob, argstuples, 2):
  ...   print "arguments: %s \nresult: %s \nerrorstring: %s\n" %(args, res, errstr)
  ...
  arguments: (2, 3)
  result: 0.471404520791
  errorstring: None

  arguments: (3, 4)
  result: 0.433012701892
  errorstring: None

  >>>

  """
  from libtbx.scheduling import philgen
  from libtbx.scheduling import job_scheduler
  from libtbx.scheduling import parallel_for
  import libtbx.scheduling
  from libtbx.scheduling import stacktrace

  technology = philgen.multiprocessing(
    capture_stderr = True, # catch each individual error message and stack trace
    qtype = philgen.mp_managed_queue,
    )

  jfactory = technology.jfactory()
  qfactory = technology.qfactory()[0]
  capacity = job_scheduler.limited(
    njobs = get_processes( processes = nproc )
    )

  creator = job_scheduler.creator(
    job_factory = jfactory,
    queue_factory = qfactory,
    capacity = capacity,
    )

  manager = creator.create()
  pfi = parallel_for.iterator(
      calculations = ( ( myfunction, args, {} ) for args in argstuples ),
      manager = manager,
      keep_input_order = False,
      )

  for i, ( calc, res ) in enumerate(pfi):
    result = None
    errstr =  None
    try:
      result = res()
    except Exception as e:
      tracestr = ""
      if stacktrace.exc_info()[1]:
        for inf in stacktrace.exc_info()[1]:
          tracestr += inf
      errstr = str(e) + "\n" + tracestr
    #calc[0] is the function name, calc[1] is the tuple of function arguments
    parmres = ( calc[1], result, errstr )

    if i >= len(argstuples)-1:
      manager.shutdown() # clean up once the last calculation has returned
      manager.join()
      creator.destroy( manager = manager )

    yield parmres # spit out results as they emerge



#  -------  SIMPLE INTERFACE TO MULTIPROCESSING -------------

#  For example of use, see phenix_regression/libtbx/tst_easy_mp.py

#  run_parallel

# Simple interface to run any target function in parallel, allowing
# specification of keyword inputs that may be different for each run

#  Class to just run a method. This is part of run_parallel below.

class run_anything(object):
  def __init__(self,kw_list=None,target_function=None):
    self.kw_list=kw_list
    self.target_function=target_function

  def __call__(self, i):
    kw=self.kw_list[i]
    return self.target_function(**kw)

def run_parallel(
   method='multiprocessing',  # multiprocessing, only choice for now
   qsub_command='qsub',       # queue command, not supported yet
   nproc=1,                   # number of processors to use
   target_function=None,      # the method to run
   kw_list=None):             # list of kw dictionaries for target_function

  n=len(kw_list)  # number of jobs to run, one per kw dict

  if nproc==1 or n<=1: # just run it for each case in list, no multiprocessing
    results=[]
    ra=run_anything(kw_list=kw_list,target_function=target_function)
    for i in range(n):
      results.append(ra(i))
  elif 0:  #(method == "multiprocessing") and (sys.platform != "win32"):
    # XXX Can crash 2015-10-13 TT so don't use it
    from libtbx.easy_mp import  pool_map
    results = pool_map(
      func=run_anything(target_function=target_function,kw_list=kw_list),
      iterable=list(range(n)),
      processes=nproc)
  else :
    from libtbx.easy_mp import parallel_map
    results=parallel_map(
      func=run_anything(target_function=target_function,kw_list=kw_list),
      iterable=list(range(n)),
      method=method,
      processes=nproc,
      callback=None,
      preserve_exception_message=True, # 2016-08-17
      qsub_command=qsub_command,
      use_manager=True )#  Always use manager 2015-10-13 TT (sys.platform == "win32"))
  return results

#  -------  END OF SIMPLE INTERFACE TO MULTIPROCESSING -------------
