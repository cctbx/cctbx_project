import boost.python
ext = boost.python.import_ext("omptbx_ext")
from omptbx_ext import *

import libtbx.introspection

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
    omp_set_dynamic(int(flag))
    return self
  dynamic = property(
    dynamic, set_dynamic,
    doc="Whether the number of threads is dynamically allocated")

  def nested(self):
    return omp_get_nested()
  def set_nested(self, flag):
    omp_set_nested(int(flag))
    return self
  nested = property(
    nested, set_nested,
    doc="Whether nested parallelism is enabled")

  def is_nested_available(self):
    try:
      saved = self.nested
      self.nested = True
      if self.nested: result = True
      else: result = False
      return result
    finally:
      self.nested = saved
  is_nested_available = property(
    is_nested_available,
    doc="Whether nested parallelism is available at all")

  def num_procs(self):
    return omp_get_num_procs()
  num_procs = property(
    num_procs, doc="Number of available processors")

env = environment()
env.dynamic = False
env.num_threads = libtbx.introspection.number_of_processors()
