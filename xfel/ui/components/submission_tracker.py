from __future__ import division

class QueueInterrogator(object):
  def __init__(self, params):
    pass

  def query(self, submission_id):
    pass

class SubmissionTracker(object):
  def __init__(self, params):
    self.params = params
    self.interrogator = QueueInterrogator(self.params)
  def track(self, submission_id):
    pass
