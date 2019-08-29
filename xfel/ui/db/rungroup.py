from __future__ import absolute_import, division, print_function
from xfel.ui.db import db_proxy

class Rungroup(db_proxy):
  def __init__(self, app, rungroup_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_rungroup" % app.params.experiment_tag, id = rungroup_id, **kwargs)
    self.rungroup_id = self.id

  def __getattr__(self, name):
    # Called only if the property cannot be found
    if name == "runs":
      return self.app.get_rungroup_runs(self.id)
    else:
      return super(Rungroup, self).__getattr__(name)

  def __setattr__(self, name, value):
    assert name != "runs"
    super(Rungroup, self).__setattr__(name, value)

  def get_first_and_last_runs(self):
    runs = self.runs
    if len(runs) == 0:
      return (None, None)
    run_ids = [r.id for r in runs]
    first = runs[run_ids.index(min(run_ids))]
    if self.open:
      last = None
    else:
      last = runs[run_ids.index(max(run_ids))]
    return first, last

  def sync_runs(self, first_run, last_run, use_ids = True):
    all_runs = self.app.get_all_runs()
    runs = self.runs
    run_ids = [r.id for r in runs]
    if self.open:
      if use_ids:
        tester = lambda x: x.id >= first_run
      else:
        tester = lambda x: int(x.run) >= first_run
    else:
      if use_ids:
        tester = lambda x: x.id >= first_run and x.id <= last_run
      else:
        tester = lambda x: int(x.run) >= first_run and int(x.run) <= last_run

    for run in all_runs:
      if tester(run):
        if not run.id in run_ids:
          self.add_run(run.id)
      else:
        if run.id in run_ids:
          self.remove_run(run.id)

  def add_run(self, run_id):
    query = "INSERT INTO `%s_rungroup_run` (rungroup_id, run_id) VALUES (%d, %d)" % ( \
      self.app.params.experiment_tag, self.id, run_id)
    self.app.execute_query(query, commit=True)

  def remove_run(self, run_id):
    query = "DELETE FROM `%s_rungroup_run` WHERE rungroup_id = %d AND run_id = %d" % ( \
      self.app.params.experiment_tag, self.id, run_id)
    self.app.execute_query(query, commit=True)
