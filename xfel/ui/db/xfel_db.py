from __future__ import division

import os
import libtbx.load_env

from xfel.ui.db.trial import Trial
from xfel.ui.db.run import Run
from xfel.ui.db.rungroup import Rungroup
from xfel.ui.db.tag import Tag
from xfel.ui.db.job import Job
from xfel.ui.db.stats import Stats

from xfel.command_line.experiment_manager import initialize as initialize_base
class initialize(initialize_base):
  expected_tables = ["run", "job", "rungroup", "trial", "tag", "run_tag", "event", "trial_rungroup",
                     "imageset", "imageset_frame", "beam", "detector", "experiment",
                     "crystal", "cell", "cell_bin", "bin"]

  def create_tables(self, sql_path = None):
    if sql_path is None:
      sql_path = os.path.join(libtbx.env.find_in_repositories("xfel/ui/db"), "schema.sql")

    return initialize_base.create_tables(self, sql_path)

class xfel_db_application(object):
  def __init__(self, dbobj, params):
    self.params = params
    self._dbobj = dbobj
    self.cursor = dbobj.cursor()

    self.init_tables = initialize(params, dbobj)

  def __getattr__(self, name):
    if name == "dbobj":
      while True:
        try:
          self._dbobj.ping()
        except Exception, e:
          print "Lost connection to db, reconnecting..."
          del self._dbobj
          from xfel.ui.db import get_db_connection
          self._dbobj = get_db_connection(self.params)
          continue
        return self._dbobj
    else:
      raise AttributeError()

  def list_lcls_runs(self):
    from xfel.xpp.simulate import file_table
    query = "https://pswww.slac.stanford.edu/ws-auth/dataexport/placed?exp_name=%s" % (self.params.experiment)
    FT = file_table(self.params, query)
    return FT.get_runs()

  def verify_tables(self):
    return self.init_tables.verify_tables()

  def create_tables(self):
    return self.init_tables.create_tables()

  def drop_tables(self):
    return self.init_tables.drop_tables()

  def create_trial(self, **kwargs):
    return Trial(self, **kwargs)

  def get_all_x(self, cls, name):
    query = "SELECT id from `%s_%s`" % (self.params.experiment_tag, name)
    cursor = self.dbobj.cursor()
    cursor.execute(query)
    return [cls(self, i[0]) for i in cursor.fetchall()]

  def get_trial(self, trial_id):
    return Trial(self, trial_id)

  def get_trial_rungroups(self, rungroup_id, only_active = False):
    query = "SELECT rungroup_id from `%s_trial_rungroup` WHERE `%s_trial_rungroup`.trial_id = %d" % \
            (self.params.experiment_tag, self.params.experiment_tag, rungroup_id)
    cursor = self.dbobj.cursor()
    cursor.execute(query)
    rungroups = [Rungroup(self, i[0]) for i in cursor.fetchall()]
    if only_active:
      return [rg for rg in rungroups if rg.active]
    else:
      return rungroups

  def get_all_trials(self, only_active = False):
    if only_active:
      return [t for t in self.get_all_x(Trial, "trial") if t.active]
    else:
      return self.get_all_x(Trial, "trial")

  def create_run(self, **kwargs):
    return Run(self, **kwargs)

  def get_run(self, run_id):
    return Run(self, run_id)

  def get_all_runs(self):
    return self.get_all_x(Run, "run")

  def create_rungroup(self, **kwargs):
    return Rungroup(self, **kwargs)

  def get_rungroup(self, rungroup_id):
    return Rungroup(self, rungroup_id)

  def get_all_rungroups(self, only_active = False):
    if only_active:
      return [rg for rg in self.get_all_x(Rungroup, "rungroup") if rg.active]
    else:
      return self.get_all_x(Rungroup, "rungroup")

  def create_tag(self, **kwargs):
    return Tag(self, **kwargs)

  def get_tag(self, tag_id):
    return Tag(self, tag_id)

  def get_run_tags(self, run_id):
    query = "SELECT tag_id from `%s_run_tag` WHERE `%s_run_tag`.run_id = %d" % \
            (self.params.experiment_tag, self.params.experiment_tag, run_id)
    cursor = self.dbobj.cursor()
    cursor.execute(query)
    return [Tag(self, i[0]) for i in cursor.fetchall()]

  def get_all_tags(self):
    return self.get_all_x(Tag, "tag")

  def delete_tag(self, tag = None, tag_id = None):
    if tag_id is None:
      tag_id = tag.tag_id

  def create_job(self, **kwargs):
    return Job(self, **kwargs)

  def get_job(self, job_id):
    return Job(self, job_id)

  def get_all_jobs(self):
    return self.get_all_x(Job, "job")

  def delete_job(self, job = None, job_id = None):
    if job_id is None:
      job_id = job.job_id

  def get_stats(self, **kwargs):
    return Stats(**kwargs)
