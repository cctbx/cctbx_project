from __future__ import division
from xfel.merging.application.scale.frame_scaler import frame_scaler
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for scaling frames of data. """
  @staticmethod
  def from_parameters(params):
    """ """
    return [frame_scaler(params)]
