from __future__ import division
from xfel.ui.db import db_proxy

class Trial(db_proxy):
  def __init__(self, app, trial_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_trial" % app.params.experiment_tag, id = trial_id, **kwargs)
    self.trial_id = self.id
    self._rungroups = None  # delay querying the db for the tags until needed

  def __getattr__(self, name):
    # Called only if the property cannot be found
    if name == "rungroups":
      if self._rungroups is None:
        self._rungroups = self.app.get_trial_rungroups(self.id)
      return self._rungroups
    else:
      return super(Trial, self).__getattr__(name)

  def __setattr__(self, name, value):
    assert name != "rungroups"
    super(Trial, self).__setattr__(name, value)

  def add_rungroup(self, rungroup):
    query = "INSERT INTO `%s_trial_rungroup` (trial_id, rungroup_id) VALUES (%d, %d)" % (
      self.app.params.experiment_tag, self.id, rungroup.id)
    self.app.execute_query(query, commit=True)
    self._rungroups = self.app.get_trial_rungroups(self.id)

  def remove_tag(self, rungroup):
    query = "DELETE FROM `%s_trial_rungroup` WHERE trial_id = %d AND rungroup_id = %s" % (
      self.app.params.experiment_tag, self.id, rungroup.id)
    self.app.execute_query(query, commit=True)
    self._rungroups = self.app.get_trial_rungroups(self.id)
