import traceback

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
    try:
      return func(arg)
    except: # intentional
      print "CAUGHT EXCEPTION:"
      traceback.print_exc()

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
    if (fixed_func is not None):
      key = fixed_func_registry_key_generator.next()
      self.fixed_func_proxy = fixed_func_proxy(key, fixed_func)
    else:
      self.fixed_func_proxy = None
    super(Pool, self).__init__(
      processes=processes,
      initializer=initializer,
      initargs=initargs,
      maxtasksperchild=maxtasksperchild)

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

def pool_map(
      processes=None,
      initializer=None,
      initargs=(),
      maxtasksperchild=None,
      func=None,
      fixed_func=None,
      iterable=None,
      args=None,
      chunksize=None):
  assert [func, fixed_func].count(None) == 1
  assert [iterable, args].count(None) == 1
  if (processes is None):
    from libtbx import introspection
    processes = introspection.number_of_processors()
  if (args is not None):
    iterable = args
    if (processes is not None):
      processes = min(processes, len(args))
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
  return result

del _
