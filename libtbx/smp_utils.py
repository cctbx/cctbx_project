
import libtbx.phil
from libtbx import thread_utils, easy_pickle
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import time
import sys
import os

multiprocessing_params = libtbx.phil.parse("""
enable_multiprocessing = True
  .type = bool
method = *mp sge
  .type = choice
  .short_caption = Parallelization method
  .caption = Multi-processing Sun_Grid_Engine
  .style = bold
nproc = None
  .type = int
  .short_caption = Number of processes
  .style = bold renderer:draw_nproc_widget
tmp_dir = None
  .type = path
""")

class Pool (object) :
  def __init__ (self, params) : #enable_multiprocessing=True, nproc=None) :
    self.enable_multiprocessing = params.enable_multiprocessing
    self.nproc = params.nproc
    self.smp_method = params.method
    self.tmp_dir = params.tmp_dir
    adopt_init_args(self, locals())
    self._pool = None
    if self.enable_multiprocessing :
      if (params.method == "mp") :
        try :
          import multiprocessing
        except ImportError :
          self.enable_multiprocessing = False
        else :
          if self.nproc is None :
            self.nproc = multiprocessing.cpu_count()
          if self.nproc > 1 :
            self._pool = multiprocessing.Pool(processes=self.nproc)
          else :
            self.enable_multiprocessing = False
      elif (params.method == "sge") :
        from libtbx.queuing_system_utils import sge
        if (not sge.is_enabled) :
          print "Can't use SGE"
          self.enable_multiprocessing = False
        else :
          self._pool = sge.SGEPool(use_grid_engine=True)
    else :
      self.enable_multiprocessing = False

  def get_tmp_dir (self) :
    if (self.tmp_dir is not None) :
      return self.tmp_dir
    else :
      return os.getcwd()

  def join (self) :
    if (self._pool is not None) :
      self._pool.join()

  def show_summary (self, out=sys.stdout) :
    if self.enable_multiprocessing :
      print >> out, "Multiprocessing is ENABLED on %d CPUs" % self.nproc
    else :
      print >> out, "Multiprocessing is DISABLED"

  def map_async (self, func, iterable, chunksize=None, callback=None) :
    if self.enable_multiprocessing :
      self._pool.map_async(_run_wrapper(func), iterable, chunksize, callback)
    else :
      return map(func, iterable)

  def map_with_async_callback (self, func, iterable, callback_async,
      chunksize=None, callback=None) :
    if self.enable_multiprocessing :
      if (self.smp_method == "mp") :
        import multiprocessing
        manager = multiprocessing.Manager()
        q = manager.Queue()
        f = _wrapper_with_queue(func, q)
        t = thread_utils.queue_monitor_thread(q, callback_async)
        t.start()
        result = self._pool.map(f, iterable, chunksize)
        t.exit()
      elif (self.smp_method == "sge") :
        if (chunksize is None) :
          chunksize = 1
        task_ids = make_ids(len(iterable))
        new_iterable = list(zip(task_ids, iterable))
        tmp_dir = os.path.join(self.get_tmp_dir(), ".results")
        if not os.path.isdir(tmp_dir) :
          os.makedirs(tmp_dir)
        f = _wrapper_with_result_file(func, tmp_dir)
        file_names = [ "%s.pkl" % id for id in task_ids ]
        t = thread_utils.file_monitor_thread(
          dir_name=tmp_dir,
          file_names=file_names,
          callback=callback_async)
        t.start()
        try :
          try :
            result = self._pool.map(f, new_iterable, chunksize)
          except KeyboardInterrupt :
            t.exit()
            sys.exit(1)
        finally :
          t.exit()
      else :
        raise RuntimeError("Can't use SMP method %s" % self.smp_method)
      return result
    else :
      results = []
      for item in iterable :
        item_result = func(item)
        callback_async(item_result)
        results.append(item_result)
      return results

  def map (self, func, iterable, chunksize=None) :
    if self.enable_multiprocessing :
      return self._pool.map(func, iterable, chunksize)
    else :
      return map(func, iterable)

  def run_many (self, objects) :
    return self.map(_run_many, objects)

  def run_many_async (self, objects, callback=None) :
    self.map_async(_run_many, objects, callback=callback)

def make_ids (n) :
  import random
  pid = os.getpid()
  pid2 = int(random.random() * 100000)
  ids = [ "%d_%d.%d" % (pid, pid2, x) for x in range(n) ]
  return ids

class _run_wrapper (object) :
  def __init__ (self, f) :
    self.f = f

  def __call__ (self, args) :
    try :
      return self.f(args)
    except Exception, e :
      if hasattr(e, "reset_module") :
        e.reset_module()
      raise

class _wrapper_with_queue (_run_wrapper) :
  def __init__ (self, f, q) :
    self.f = f
    self.q = q

  def __call__ (self, args) :
    result = _run_wrapper.__call__(self, args)
    self.q.put(result)
    return result

class _wrapper_with_result_file (_run_wrapper) :
  def __init__ (self, f, tmp_dir) :
    self.f = f
    self.tmp_dir = tmp_dir

  def __call__ (self, args) :
    task_id = args[0]
    new_args = args[1]
    result = _run_wrapper.__call__(self, new_args)
    pkl_file = os.path.join(self.tmp_dir, "%s.pkl" % task_id)
    easy_pickle.dump(pkl_file, result)
    return result

def _run_many (run_object) :
  return run_object.run()

########################################################################
# tests
#
def exercise (timeout=0.1, verbose=False) :
  import random
  data = [ random.random() for n in range(20) ]
  params = multiprocessing_params.fetch().extract()
  params.enable_multiprocessing = False
  p = Pool(params)
  r1 = p.map(_test_f, data)
  data_objects = [ _exercise_run(x) for x in data ]
  r2 = p.run_many(data_objects)
  t1 = time.time()
  r3 = p.map(_test_function(timeout), data)
  t2 = time.time()
  cb = _exercise_callback()
  r4 = p.map_with_async_callback(_test_f, data, cb)
  assert (cb.c == len(data) == len(r4))
  cb.reset()
  assert (r1 == r2 == r3 == r4)
  params.enable_multiprocessing = True
  params.method = "mp"
  p = Pool(params)
  r5 = p.map(_test_f, data)
  r6 = p.run_many(data_objects)
  t3 = time.time()
  r7 = p.map(_test_function(timeout), data)
  t4 = time.time()
  r8 = p.map_with_async_callback(_test_f, data, cb)
  # XXX this does not work!
  #assert (cb.c == len(data) == len(r8))
  assert (r8 == r7 == r6 == r5 == r4)
  if verbose :
    print "CPU count: %d  Speedup: %.1fx" % (p.nproc, (t2 - t1) / (t4 - t3))
  from libtbx.queuing_system_utils import sge_utils
  try :
    q = sge_utils.qstat_parse()
  except RuntimeError :
    print "Skipping SGE-specific tests in smp_utils"
    return
  params.method = "sge"
  p = Pool(params)
  r9 = p.map(_test_f, data)
  r10 = p.run_many(data_objects)
  t5 = time.time()
  r11 = p.map(_test_function(timeout), data)
  t6 = time.time()
  r12 = p.map_with_async_callback(_test_f, data, cb)
  assert (r12 == r11 == r10 == r9 == r8)

class _exercise_callback (object) :
  def __init__ (self) :
    self.c = 0

  def __call__ (self, result) :
    self.c += 1

  def reset (self) :
    self.c = 0

class _test_function (object) :
  def __init__ (self, timeout=1) :
    self.timeout = timeout

  def __call__ (self, n) :
    return _test_f(n, timeout=self.timeout)

def _test_f (n, timeout=0.1) :
  time.sleep(timeout)
  #print "%.5f" % n
  return n**2

class _exercise_run (object) :
  def __init__ (self, n) :
    self.n = n

  def run (self) :
    return _test_f(self.n)

if __name__ == "__main__" :
  if (len(sys.argv) > 1) :
    exercise(timeout=float(sys.argv[1]))
  else :
    exercise()
  print "OK"
