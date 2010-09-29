
import libtbx.phil
from libtbx import thread_utils, easy_pickle
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import time
import sys
import os

#--- backwards-compatability check
if sys.version_info[0] > 2 or sys.version_info[1] >= 6 :
  import multiprocessing
  cpu_count = multiprocessing.cpu_count

  def smp_map (func, iterable, chunksize=None, callback=None, nproc=None) :
    pool = multiprocessing.Pool(processes=nproc)
    result = pool.map_async(func, iterable, chunksize, callback)
    return result.get()

#--- pre-2.6 equivalents
else :
  cpu_count = lambda: 1
  def smp_map (func, iterable, chunksize=None, callback=None, nproc=None) :
    return map(func, iterable)

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

#--- test functions
def exercise () :
  print "Running smp_utils tests on %d processors" % cpu_count()
  i = [ (0.1 * x, 100000) for x in xrange(0, 16) ]
  (m_single, t_single) = benchmark_function(map, (t, i), "map")
  (m_smp, t_smp) = benchmark_function(smp_map, (t, i), "smp_map")
  assert m_single == m_smp
  print "Speedup on %d cpus: %.1fx" % (cpu_count(), t_single / t_smp)

def benchmark_function (func, args, func_name) :
  t1 = time.time()
  r = func(*args)
  t2 = time.time()
  print "%s(): %.3f seconds" % (func_name, (t2 - t1))
  return (r, t2-t1)

# see http://dan.corlan.net/bench.html#PYTHON
# "The program we benchmarked computes the same 100-term polynomial 500,000
# times, using exactly the same algorithm. In all programs the indices of the
# polynomial are kept in a local float vector. In this, the program only tests
# the quality of code which accesses local vectors and performs simple
# arithmetics in loops, and is free from differences in the standard library,
# operating system calls and, indeed, the presence of any advanced language
# features."
# TODO: benchmark using a suitable cctbx function instead
def t (args) :
  (x, n) = args
  mu = 10.0
  pu = 0.0
  pol =[0] * 100
  r = range(0,100)

  for i in range(0,n):
    for j in r:
      pol[j] = mu = (mu + 2.0) / 2.0
    su = 0.0
    for j in r:
      su = x * su + pol[j]
    pu = pu + su
  return pu
