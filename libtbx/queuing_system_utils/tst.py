from __future__ import absolute_import, division, print_function
from six.moves import range

# XXX this is intended to be a simple template for debugging queueing system
# support issues, not a full regression test.

class target(object):
  def __init__(self, x):
    self.x = x

  def __call__(self):
    import math
    results = []
    for n in range(x):
      nn = math.sqrt(n**3)
      print(nn)
      results.append(nn)

def exercise():
  from libtbx.queuing_system_utils import generic as queuing
  t = target(1000000)
  job = queuing.qsub(
    target=t,
    platform="sge")
  job.start()
  assert (isinstance(job.jobid, int))
  while job.is_alive():
    pass
  print("done")

if (__name__ == "__main__"):
  exercise()
