from __future__ import division
from xfel.merging.application.filter.experiment_filter import experiment_filter
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for filtering experiments. """
  @staticmethod
  def from_parameters(params):
    """ """
    return [experiment_filter(params)]
