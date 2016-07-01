from __future__ import division

import os, time
import libtbx.load_env

from xfel.ui.db.trial import Trial
from xfel.ui.db.run import Run
from xfel.ui.db.rungroup import Rungroup
from xfel.ui.db.tag import Tag
from xfel.ui.db.job import Job
from xfel.ui.db.stats import Stats
from xfel.ui.db.experiment import Cell, Bin, Isoform

from xfel.ui.db import get_db_connection

try:
  from MySQLdb import OperationalError
except ImportError:
  from libtbx.utils import Sorry
  raise Sorry('Mysql not available')

from xfel.command_line.experiment_manager import initialize as initialize_base

class initialize(initialize_base):
  expected_tables = ["run", "job", "rungroup", "trial", "tag", "run_tag", "event", "trial_rungroup",
                     "imageset", "imageset_event", "beam", "detector", "experiment",
                     "crystal", "cell", "cell_bin", "bin", "isoform"]

  def __init__(self, params, dbobj):
    initialize_base.__init__(self, params, dbobj, interactive = False, drop_tables = None)

  def create_tables(self, sql_path = None):
    if sql_path is None:
      sql_path = os.path.join(libtbx.env.find_in_repositories("xfel/ui/db"), "schema.sql")

    return initialize_base.create_tables(self, sql_path)

  def set_up_columns_dict(self, app):
    columns_dict = {}
    for table in self.expected_tables:
      table_name = "%s_%s" % (self.params.experiment_tag, table)
      query = "SHOW COLUMNS FROM `%s`" % (table_name)
      cursor = app.execute_query(query)
      columns_dict[table_name] = [c[0] for c in cursor.fetchall() if c[0] != 'id']
    return columns_dict

class xfel_db_application(object):
  def __init__(self, params):
    self.params = params
    dbobj = get_db_connection(params)
    self.init_tables = initialize(params, dbobj) # only place where a connection is held

    self.columns_dict = self.init_tables.set_up_columns_dict(self)

  def execute_query(self, query, commit = False, verbose = False):
    if verbose:
      from time import time
      st = time()

    retry_count = 0
    retry_max = 10
    sleep_time = 0.1
    while retry_count < retry_max:
      try:
        dbobj = get_db_connection(self.params)
        cursor = dbobj.cursor()
        cursor.execute(query)
        if commit:
          dbobj.commit()

        if verbose:
          print 'SQLTime Taken = % 10.6f seconds' % (time() - st), query[:min(len(query),160)]
        return cursor
      except OperationalError, e:
        if "Can't connect to MySQL server" not in str(e):
          raise e
        retry_count += 1
        print "Couldn't connect to MYSQL, retry", retry_count
        time.sleep(sleep_time)
        sleep_time *= 2
      except Exception, e:
        print "Couldn't execute MYSQL query.  Query:"
        print query
        print "Exception:"
        print str(e)
        raise e
    raise Sorry("Couldn't execute MYSQL query. Too many reconnects. Query: %s"%query)

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

  def create_trial(self, d_min = 1.5, n_bins = 10, **kwargs):
    # d_min and n_bins only used if isoforms are in this trial

    trial = Trial(self, **kwargs)
    if trial.target_phil_str is not None:
      from xfel.command_line.xtc_process import phil_scope
      from iotbx.phil import parse
      trial_params = phil_scope.fetch(parse(trial.target_phil_str)).extract()
      if len(trial_params.indexing.stills.isoforms) > 0:
        for isoform in trial_params.indexing.stills.isoforms:
          print "Creating isoform", isoform.name
          db_isoform = Isoform(self,
                               name = isoform.name,
                               trial_id = trial.id)
          a, b, c, alpha, beta, gamma = isoform.cell.parameters()
          cell = self.create_cell(cell_a = a,
                                  cell_b = b,
                                  cell_c = c,
                                  cell_alpha = alpha,
                                  cell_beta = beta,
                                  cell_gamma = gamma,
                                  lookup_symbol = isoform.lookup_symbol,
                                  isoform_id = db_isoform.id)
          from cctbx.crystal import symmetry

          cs = symmetry(unit_cell = isoform.cell,space_group_symbol=str(isoform.lookup_symbol))
          mset = cs.build_miller_set(anomalous_flag=False, d_min=d_min)
          binner = mset.setup_binner(n_bins=n_bins)
          for i in binner.range_used():
            d_max, d_min = binner.bin_d_range(i)
            Bin(self,
                number = i,
                d_min = d_min,
                d_max = d_max,
                total_hkl = binner.counts_complete()[i],
                cell_id = cell.id)
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

  def get_cell_bins(self, cell_id):
    query = "SELECT id FROM `%s_bin` WHERE cell_id = %d" % \
            (self.params.experiment_tag, cell_id)
    cursor = self.execute_query(query)
    return [Bin(self, bin_id = i[0]) for i in cursor.fetchall()]

  def get_all_x(self, cls, name):
    query = "SELECT id FROM `%s_%s`" % (self.params.experiment_tag, name)
    cursor = self.execute_query(query)
    return [cls(self, i[0]) for i in cursor.fetchall()]

  def get_trial(self, trial_id = None, trial_number = None):
    assert [trial_id, trial_number].count(None) == 1
    if trial_id is None:
      trials = [t for t in self.get_all_trials() if t.trial == trial_number]
      assert len(trials) == 1
      return trials[0]
    else:
      return Trial(self, trial_id)

  def get_trial_rungroups(self, trial_id, only_active = False):
    query = "SELECT rungroup_id FROM `%s_trial_rungroup` WHERE `%s_trial_rungroup`.trial_id = %d" % \
            (self.params.experiment_tag, self.params.experiment_tag, trial_id)
    cursor = self.execute_query(query)
    rungroups = [Rungroup(self, i[0]) for i in cursor.fetchall()]
    if only_active:
      return [rg for rg in rungroups if rg.active]
    else:
      return rungroups

  def get_trial_runs(self, trial_id):
    rungroups = self.get_trial_rungroups(trial_id, only_active=True)
    runs = []
    run_ids = []
    for rungroup in rungroups:
      for run in rungroup.runs:
        if run.id not in run_ids:
          run_ids.append(run.id)
          runs.append(run)
    return runs

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

  def delete_x(self, item, item_id):
    query = "DELETE FROM `%s` WHERE id = %d"%(item.table_name, item_id)
    self.execute_query(query, commit = True)

  def delete_tag(self, tag = None, tag_id = None):
    assert [tag, tag_id].count(None) == 1

    if tag_id is None:
      tag_id = tag.tag_id
    else:
      tag = self.get_tag(tag_id)

    query = "DELETE FROM `%s_run_tag` WHERE tag_id = %d" % (self.params.experiment_tag, tag_id)
    self.execute_query(query, commit=True)

    self.delete_x(tag, tag_id)

  def create_job(self, **kwargs):
    return Job(self, **kwargs)

  def get_job(self, job_id):
    return Job(self, job_id)

  def get_all_jobs(self):
    return self.get_all_x(Job, "job")

  def delete_job(self, job = None, job_id = None):
    assert [job, job_id].count(None) == 1
    if job_id is None:
      job_id = job.job_id
    else:
      job = self.get_job(job_id)

    self.delete_x(job, job_id)

  def get_stats(self, **kwargs):
    return Stats(self, **kwargs)
