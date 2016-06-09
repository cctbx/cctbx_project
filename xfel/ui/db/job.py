from __future__ import division

class job(object):
  def __init__(self, job_id = None, **kwargs):
    if job_id is None:
      job_id = 4 # create a new job

    self.job_id = job_id

    # dummy values:
    self.status = "Unqueued"
    self.trial_id = 0
    self.rungroup_id = 0
    self.run_id = 0
