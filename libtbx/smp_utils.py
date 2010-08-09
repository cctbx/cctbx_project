import sys, time
from libtbx.utils import Sorry
from libtbx import adopt_init_args

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

class manager (object) :
  def __init__ (self, enable_multiprocessing=True, nproc=None) :
    adopt_init_args(self, locals())
    self.pool = None
    if enable_multiprocessing and (self.nproc is None or self.nproc > 1) :
      try :
        import multiprocessing
      except ImportError :
        self.enable_multiprocessing = False
      else :
        if self.nproc is None :
          self.nproc = multiprocessing.cpu_count()
        if self.nproc > 1 :
          self.pool = multiprocessing.Pool(processes=self.nproc)
        else :
          self.enable_multiprocessing = False
    else :
      self.enable_multiprocessing = False

  def show_summary (self, out=sys.stdout) :
    if self.enable_multiprocessing :
      print >> out, "Multiprocessing is ENABLED on %d CPUs" % self.nproc
    else :
      print >> out, "Multiprocessing is DISABLED"

  def map_async (self, func, iterable, chunksize=None, callback=None) :
    if self.enable_multiprocessing :
      result = self.pool.map_async(func, iterable, chunksize, callback)
      return result.get()
    else :
      return map(func, iterable)

  def map (self, func, iterable, chunksize=None) :
    if self.enable_multiprocessing :
      return self.pool.map(func, iterable, chunksize)
    else :
      return map(func, iterable)

  def run_many (self, objects) :
    return self.map(_run_many, objects)

  def run_many_async (self, objects, callback=None) :
    self.map_async(_run_many, objects, callback=callback)

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
