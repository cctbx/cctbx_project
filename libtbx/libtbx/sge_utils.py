"Sun Grid Engine utilities"

import sys, os

def int_or_none(v):
  if (v is None): return None
  return int(v)

def job_id():
  return int_or_none(os.environ.get("JOB_ID", None))

class task_info:

  def __init__(self):
    self.first = int_or_none(os.environ.get("SGE_TASK_FIRST", None))
    self.last = int_or_none(os.environ.get("SGE_TASK_LAST", None))
    self.id = int_or_none(os.environ.get("SGE_TASK_ID", None))
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
