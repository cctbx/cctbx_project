import os
import time
import boost.optional
import boost.python
boost.python.import_ext("scitbx_random_ext")
import scitbx_random_ext

builtin_int = __builtins__["int"]
builtin_long = __builtins__["long"]

def get_random_seed():
  try:
    result = builtin_long(os.getpid() * (2**16)) \
           + builtin_long(time.time() * (2**8))
  except KeyboardInterrupt: raise
  except:
    result = time.time()
  return builtin_int(result % (2**31-1))

def set_random_seed(value):
  mt19937.seed(value)

from scitbx_random_ext import *
mt19937 = mersenne_twister_19937(value=get_random_seed())


class variate_factory(object):
  """
  The corner stone of this package.

  Synopsis

    g = scitbx.random.variate(scitbx.random.normal_distribution(mean, sigma))

  How to use it from other modules?

    scitbx.random.variate.register(xxxx_ext)

  where xxxx_ext is a Boost.Python module featuring overloaded 'variate'
  functions, variate(engine, distribution), where engine is the like
  of mt19937 and distribution is the like of scitbx.random.normal_distribution.
  """

  def __init__(self):
    self.modules = []

  def register_module(self, module):
    self.modules.append(module)

  def variate_functions(self):
    yield scitbx_random_ext.variate
    for module in self.modules: yield module.variate

  def __call__(self, distribution, engine=mt19937):
    exceptions = []
    for variate in self.variate_functions():
      try:
        return variate(engine, distribution)
      except Exception, e:
        if str(e.__class__) == "<class 'Boost.Python.ArgumentError'>":
          exceptions.append(e)
          continue
        else:
          raise
    else:
      raise RuntimeError('\n'.join([ str(e) for e in exceptions ]))

# Instantiate the one single variate factory doing it all
variate = variate_factory()
