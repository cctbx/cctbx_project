from __future__ import absolute_import, division, print_function
from xfel.ui.db import db_proxy

class Trial(db_proxy):
  def __init__(self, app, trial_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_trial" % app.params.experiment_tag, id = trial_id, **kwargs)
    self.trial_id = self.id

  def __getattr__(self, name):
    # Called only if the property cannot be found
    if name == "rungroups":
      return self.app.get_trial_rungroups(self.id, only_active=True)
    elif name == "runs":
      return self.app.get_trial_runs(self.id)
    elif name == "isoforms":
      return self.app.get_trial_isoforms(self.id)
    elif name == "tags":
      return self.app.get_trial_tags(self.id)
    elif name == "cell":
      return self.app.get_trial_cell(self.id)
    else:
      return super(Trial, self).__getattr__(name)

  def __setattr__(self, name, value):
    assert name not in ["rungroups", "runs", "isoforms"]
    super(Trial, self).__setattr__(name, value)

  def add_rungroup(self, rungroup):
    query = "INSERT INTO `%s_trial_rungroup` (trial_id, rungroup_id) VALUES (%d, %d)" % (
      self.app.params.experiment_tag, self.id, rungroup.id)
    self.app.execute_query(query, commit=True)

  def remove_rungroup(self, rungroup):
    query = "DELETE FROM `%s_trial_rungroup` WHERE trial_id = %d AND rungroup_id = %s" % (
      self.app.params.experiment_tag, self.id, rungroup.id)
    self.app.execute_query(query, commit=True)
