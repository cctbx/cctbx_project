from __future__ import division

import os
import libtbx.load_env

from xfel.ui.db.trial import Trial
from xfel.ui.db.run import Run
from xfel.ui.db.rungroup import Rungroup
from xfel.ui.db.tag import Tag
from xfel.ui.db.job import Job
from xfel.ui.db.stats import Stats
from xfel.ui.db.experiment import Cell #Cell_Bin

from xfel.ui.db import get_db_connection

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
  def __init__(self, params):
    self.params = params
    dbobj = get_db_connection(params)
    self.init_tables = initialize(params, dbobj) # only place where a connection is held

  def execute_query(self, query, commit = False):
    try:
      dbobj = get_db_connection(self.params)
      cursor = dbobj.cursor()
      cursor.execute(query)
      if commit:
        dbobj.commit()
      return cursor
    except Exception, e:
      print "Couldn't execute MYSQL query.  Query:"
      print query
      print "Exception:"
      print str(e)
      raise e

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
    trial = Trial(self, **kwargs)

    from xfel.command_line.xtc_process import phil_scope
    from iotbx.phil import parse
    trial_params = phil_scope.fetch(parse(trial.target_phil_str)).extract()
    if len(trial_params.indexing.stills.isoforms) > 0:
      for isoform in trial_params.indexing.stills.isoforms:
        cell = self.get_cell(name = isoform.name)
        if cell is None:
          print "Creating isoform", isoform.name
          a, b, c, alpha, beta, gamma = isoform.cell.parameters()
          cell = self.create_cell(name = isoform.name,
                                  cell_a = a,
                                  cell_b = b,
                                  cell_c = c,
                                  cell_alpha = alpha,
                                  cell_beta = beta,
                                  cell_gamma = gamma,
                                  lookup_symbol = isoform.lookup_symbol)
          #TODO add isoform bins
    return trial

  def create_cell(self, **kwargs):
    return Cell(self, **kwargs)

  def get_cell(self, cell_id = None, name = None):
    assert [cell_id, name].count(None) == 1
    if name is not None:
      query = "SELECT id FROM `%s_cell` WHERE name = '%s'"%(self.params.experiment_tag, name)
      cursor = self.execute_query(query)
      results = cursor.fetchall()
      assert len(results) in [0,1]
      if len(results) == 0:
        return None
      cell_id = int(results[0][0])

    return Cell(self, cell_id=cell_id)

  def get_all_x(self, cls, name):
    query = "SELECT id FROM `%s_%s`" % (self.params.experiment_tag, name)
    cursor = self.execute_query(query)
    return [cls(self, i[0]) for i in cursor.fetchall()]

  def get_trial(self, trial_id):
    return Trial(self, trial_id)

  def get_trial_rungroups(self, rungroup_id, only_active = False):
    query = "SELECT rungroup_id FROM `%s_trial_rungroup` WHERE `%s_trial_rungroup`.trial_id = %d" % \
            (self.params.experiment_tag, self.params.experiment_tag, rungroup_id)
    cursor = self.execute_query(query)
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

  def get_run(self, run_id = None, run_number = None):
    assert [run_id, run_number].count(None) == 1
    if run_id is not None:
      return Run(self, run_id)
    runs = [r for r in self.get_all_runs() if r.run == run_number]
    assert len(runs) == 1
    return runs[0]

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
    cursor = self.execute_query(query)
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
