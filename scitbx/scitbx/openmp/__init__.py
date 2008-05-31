import boost.python
try:
  ext = boost.python.import_ext("scitbx_openmp_ext")
  from scitbx_openmp_ext import *
  available = True
except ImportError:
  available = False

if (available):

  class environment(object):

    def num_threads(self):
      return omp_get_max_threads()
    def set_num_threads(self, n):
      omp_set_num_threads(n)
      return self
    num_threads = property(
      num_threads, set_num_threads,
      doc="Number of threads to distribute the work over")

    def dynamic(self):
      return omp_get_dynamic()
    def set_dynamic(self, flag):
      omp_set_dynamic(flag)
      return self
    dynamic = property(
      dynamic, set_dynamic,
      doc="Whether the number of threads is dynamically allocated")

    def nested(self):
      return omp_get_nested()
    def set_nested(self, flag):
      omp_set_nested(flag)
      return self
    nested = property(
      nested, set_nested,
      doc="Whether nested parallelism is enabled")

    def num_procs(self):
      return omp_get_num_procs()
    num_procs = property(
      num_procs, doc="Number of available processors")

  environment = environment()
