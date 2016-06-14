from __future__ import division
from xfel.ui.db import db_proxy

class Trial(db_proxy):
  def __init__(self, app, trial_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_trial" % app.params.experiment_tag, id = trial_id, **kwargs)
    self.trial_id = self.id

