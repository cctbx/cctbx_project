"Sun Grid Engine utilities"

import sys, os

def int_or_none(v):
  if (v is None or v == "undefined"): return None
  return int(v)

def arch():
  return os.environ.get("SGE_ARCH")

def job_id():
  return int_or_none(os.environ.get("JOB_ID"))

class task_info:

  def __init__(self):
    self.first = int_or_none(os.environ.get("SGE_TASK_FIRST"))
    self.last = int_or_none(os.environ.get("SGE_TASK_LAST"))
    self.id = int_or_none(os.environ.get("SGE_TASK_ID"))
    assert [self.first, self.last, self.id].count(None) in [0, 3]

  def show(self, out=None):
    if (out is None): out = sys.stdout
    print >> out, "SGE_TASK_FIRST:", self.first
    print >> out, "SGE_TASK_LAST: ", self.last
    print >> out, "SGE_TASK_ID:   ", self.id
    return self

  def as_n_i_pair(self):
    if (self.first is None): return 1, 0
    assert self.first == 1
    return self.last, self.id-1
