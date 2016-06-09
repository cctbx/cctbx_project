from __future__ import division

from xfel.ui.db.trial import trial
from xfel.ui.db.run import run
from xfel.ui.db.rungroup import rungroup
from xfel.ui.db.tag import tag
from xfel.ui.db.job import job
from xfel.ui.db.stats import stats

class xfel_db_application(object):
  def __init__(self, dbobj):
    self.dbobj = dbobj
    self.cursor = dbobj.cursor()

  def create_trial(self, **kwargs):
    return trial(**kwargs)

  def get_trial(self, trial_id):
    return trial(trial_id)

  def get_all_trials(self):
    return [trial(i) for i in xrange(3)]

  def create_run(self, **kwargs):
    return run(**kwargs)

  def get_run(self, run_id):
    return run(run_id)

  def get_all_runs(self):
    return [run(i) for i in xrange(3)]

  def create_rungroup(self, **kwargs):
    return rungroup(**kwargs)

  def get_rungroup(self, rungroup_id):
    return rungroup(rungroup_id)

  def get_all_rungroups(self):
    return [rungroup(i) for i in xrange(3)]

  def create_tag(self, **kwargs):
    return tag(**kwargs)

  def get_tag(self, tag_id):
    return tag(tag_id)

  def get_all_tags(self):
    return [tag(i) for i in xrange(3)]

  def delete_tag(self, tag = None, tag_id = None):
    if tag_id is None:
      tag_id = tag.tag_id

  def create_job(self, **kwargs):
    return job(**kwargs)

  def get_job(self, job_id):
    return job(job_id)

  def get_all_jobs(self):
    return [job(i) for i in xrange(3)]

  def delete_tag(self, job = None, job_id = None):
    if job_id is None:
      job_id = job.job_id

  def get_stats(self, **kwargs):
    return stats(**kwargs)
