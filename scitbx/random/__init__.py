import os
import time
from scitbx.random import ext

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

mt19937 = ext.mt19937(value=get_random_seed())


class normal_variate_generator(object):
  def __init__(self, mean=0, sigma=1):
    from scitbx.array_family import flex
    self._distribution = ext.normal_distribution(mean=mean, sigma=sigma)
    self._generator = ext.normal_variate_generator(mt19937, self._distribution)

  def __call__(self, size_t=None):
    if size_t is None:
      return self._generator()
    else:
      return self._generator(size_t)
