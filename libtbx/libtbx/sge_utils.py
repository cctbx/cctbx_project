"Sun Grid Engine utilities"

from __future__ import generators
import sys, os

def int_or_none(v):
  if (v is None or v == "undefined"): return None
  return int(v)

def job_id():
  return int_or_none(os.environ.get("JOB_ID"))

def arch():
  return os.environ.get("SGE_ARCH")

class task_info(object):

  def __init__(self):
    self.first = int_or_none(os.environ.get("SGE_TASK_FIRST"))
    self.last = int_or_none(os.environ.get("SGE_TASK_LAST"))
    self.id = int_or_none(os.environ.get("SGE_TASK_ID"))
    assert [self.first, self.last, self.id].count(None) in [0, 3]

  def show(self, out=None, prefix="", even_if_none=False):
    if (out is None): out = sys.stdout
    if (self.first is not None or even_if_none):
      print >> out, prefix+"SGE_TASK_FIRST =", self.first
    if (self.last is not None or even_if_none):
      print >> out, prefix+"SGE_TASK_LAST =", self.last
    if (self.id is not None or even_if_none):
      print >> out, prefix+"SGE_TASK_ID =", self.id
    return self

  def as_n_i_pair(self):
    if (self.first is None): return 1, 0
    assert self.first == 1
    return self.last, self.id-1

  def skip_loop_iteration(self):
    return skip_loop_iteration(*self.as_n_i_pair())

def skip_loop_iteration(n, i):
  j = 0
  while True:
    yield (j % n != i)
    j += 1

class info(task_info):

  def __init__(self):
    task_info.__init__(self)
    self.job_id = job_id()
    self.arch = arch()

  def show(self, out=None, prefix="", even_if_none=False):
    if (out is None): out = sys.stdout
    if (self.job_id is not None or even_if_none):
      print >> out, prefix+"JOB_ID =", self.job_id
    if (self.arch is not None or even_if_none):
      print >> out, prefix+"SGE_ARCH =", self.arch
    task_info.show(self, out=out, prefix=prefix, even_if_none=even_if_none)

if (__name__ == "__main__"):
  info().show(prefix="*** ", even_if_none=True)
  print "OK"
