from __future__ import absolute_import, division, print_function
from xfel.ui.db import db_proxy

def runs_in_range(all_runs, first_run, last_run, use_ids = True):
  """ The runs a run group covers: those from first_run up to last_run, or
      every run from first_run on when last_run is None (an open group).
      Compared by run id, or by run number when use_ids is False. """
  key = (lambda r: r.id) if use_ids else (lambda r: int(r.run))
  if last_run is None:
    return [r for r in all_runs if key(r) >= first_run]
  return [r for r in all_runs if first_run <= key(r) <= last_run]

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
    use_ids = self.app.params.facility.name not in ['lcls']
    runs = self.runs
    if len(runs) == 0:
      return (None, None)
    if use_ids:
      run_ids = [r.id for r in runs]
      first = runs[run_ids.index(min(run_ids))]
      if self.open:
        last = None
      else:
        last = runs[run_ids.index(max(run_ids))]
    else:
      run_numbers = [int(r.run) for r in runs]
      first = runs[run_numbers.index(min(run_numbers))]
      if self.open:
        last = None
      else:
        last = runs[run_numbers.index(max(run_numbers))]
    return first, last

  def sync_runs(self, first_run, last_run, use_ids = True):
    all_runs = self.app.get_all_runs()
    run_ids = [r.id for r in self.runs]
    wanted = [r.id for r in runs_in_range(all_runs, first_run, None if self.open else last_run, use_ids)]
    for run in all_runs:
      if run.id in wanted:
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
