from __future__ import division

class trial(object):
  def __init__(self, trial_id = None, **kwargs):
    if trial_id is None:
      trial_id = 4 # create a new trial

    self.trial_id = trial_id

    # dummy values:
    self.active = True
    self.phil = ""
    self.comment = ""
