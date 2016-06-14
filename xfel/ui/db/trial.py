from __future__ import division
from xfel.ui.db import db_proxy

class Trial(db_proxy):
  def __init__(self, dbobj, trial_id = None, **kwargs):
    db_proxy.__init__(self, dbobj, id = trial_id, **kwargs)

    if trial_id is None:
      trial_id = 4 # create a new trial

    self.trial_id = trial_id

    # dummy values:
    self.active = True
    self.phil = ""
    self.comment = ""
