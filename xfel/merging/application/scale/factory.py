from __future__ import absolute_import, division, print_function
from xfel.merging.application.scale.experiment_scaler import experiment_scaler
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for scaling frames of data. """
  @staticmethod
  def from_parameters(params, additional_info=None):
    """ """
    return [experiment_scaler(params)]
