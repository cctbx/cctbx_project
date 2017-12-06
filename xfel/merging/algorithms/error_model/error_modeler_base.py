from __future__ import division

"""
Base class for error modeling. Used after merging to adjust the errors.
"""
class error_modeler_base(object):
  def __init__(self, scaler):
    self.scaler = scaler
    self.log = scaler.log

  def adjust_errors(self):
    raise NotImplementedError("Override!")
