from __future__ import division
from xfel.ui.db import db_proxy

def submit_all_jobs(app):
  jobs = app.get_all_jobs()
  trials = app.get_all_trials(active = True)
  rungroups = app.get_all_rungroups(active = True)

  print "Jobs and trials", jobs, trials, rungroups

class Job(db_proxy):
  def __init__(self, app, job_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_job" % app.params.experiment_tag, id = job_id, **kwargs)
    self.job_id = self.id

