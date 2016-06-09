from __future__ import division

from xfel.ui.db.trial import Trial
from xfel.ui.db.run import Run
from xfel.ui.db.rungroup import Rungroup
from xfel.ui.db.tag import Tag
from xfel.ui.db.job import Job
from xfel.ui.db.stats import Stats

class xfel_db_application(object):
  def __init__(self, dbobj):
    self.dbobj = dbobj
    self.cursor = dbobj.cursor()

  def create_trial(self, **kwargs):
    return Trial(**kwargs)

  def get_trial(self, trial_id):
    return Trial(trial_id)

  def get_all_trials(self):
    return [Trial(i) for i in xrange(3)]

  def create_run(self, **kwargs):
    return Run(**kwargs)

  def get_run(self, run_id):
    return Run(run_id)

  def get_all_runs(self):
    return [Run(i) for i in xrange(3)]

  def create_rungroup(self, **kwargs):
    return Rungroup(**kwargs)

  def get_rungroup(self, rungroup_id):
    return Rungroup(rungroup_id)

  def get_all_rungroups(self):
    return [Rungroup(i) for i in xrange(3)]

  def create_tag(self, **kwargs):
    return Tag(**kwargs)

  def get_tag(self, tag_id):
    return Tag(tag_id)

  def get_all_tags(self):
    return [Tag(i) for i in xrange(3)]

  def get_run_tags(self, run = None, run_id = None):
    if run_id is None:
      run_id = run.run_id

    return [Tag(i) for i in xrange(3)]

  def delete_tag(self, tag = None, tag_id = None):
    if tag_id is None:
      tag_id = tag.tag_id

  def create_job(self, **kwargs):
    return Job(**kwargs)

  def get_job(self, job_id):
    return Job(job_id)

  def get_all_jobs(self):
    return [Job(i) for i in xrange(3)]

  def delete_tag(self, job = None, job_id = None):
    if job_id is None:
      job_id = job.job_id

  def get_stats(self, **kwargs):
    return Stats(**kwargs)
