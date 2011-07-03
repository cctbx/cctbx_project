from libtbx import Auto
from cStringIO import StringIO
import traceback
import sys

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

from multiprocessing.pool import Pool as _

class Pool(_):

  def __init__(self,
        processes=None,
        initializer=None,
        initargs=(),
        maxtasksperchild=None,
        fixed_func=None):
    if (processes is None):
      from libtbx import introspection
      processes = introspection.number_of_processors()
    self.processes = processes
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

class func_wrapper(object):

  def __init__(O, func, buffer_stdout_stderr):
    O.func = func
    O.buffer_stdout_stderr = buffer_stdout_stderr

  def __call__(O, arg):
    if (O.buffer_stdout_stderr):
      sys.stderr = sys.stdout = sio = StringIO()
    try:
      result = O.func(arg)
    except: # intentional
      result = None
      print "CAUGHT EXCEPTION:"
      traceback.print_exc(file=sys.stdout)
    if (O.buffer_stdout_stderr):
      return (sio.getvalue(), result)
    return result

default_func_wrapper = func_wrapper

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
      func_wrapper=Auto,
      log=None,
      buffer_stdout_stderr=False):
  assert [func, fixed_func].count(None) == 1
  assert [iterable, args].count(None) == 1
  if (func_wrapper is Auto):
    func_wrapper = default_func_wrapper
  if (func_wrapper is not None):
    if (func is not None):
      func = func_wrapper(func, buffer_stdout_stderr)
    else:
      fixed_func = func_wrapper(fixed_func, buffer_stdout_stderr)
  if (processes is None):
    from libtbx import introspection
    processes = introspection.number_of_processors()
  if (args is not None):
    iterable = args
    if (processes is not None):
      processes = min(processes, len(args))
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
