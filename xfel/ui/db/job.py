from __future__ import division
from xfel.ui.db import db_proxy

class Job(db_proxy):
  def __init__(self, dbobj, job_id = None, **kwargs):
    db_proxy.__init__(self, dbobj, id = job_id, **kwargs)

    if job_id is None:
      job_id = 4 # create a new job

    self.job_id = job_id

    # dummy values:
    self.status = "Unqueued"
    self.trial_id = 0
    self.rungroup_id = 0
    self.run_id = 0
