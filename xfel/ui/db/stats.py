from __future__ import division

class Stats(object):
  def __init__(self, tags = None,
                     trial = None,
                     d_min = None):
    self.tags = tags
    self.trial = trial
    self.d_min = d_min

    # dummy values:
    self.multiplicity = 10
