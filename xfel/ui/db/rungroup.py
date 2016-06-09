from __future__ import division

class rungroup(object):
  def __init__(self, rungroup_id = None, **kwargs):
    if rungroup_id is None:
      rungroup_id = 4 # create a new rungroup

    self.rungroup_id = rungroup_id

    # dummy values:
    startrun = 1
    endrun = 5
    detz_parameter = 580
    binning = None
    comment = ""
