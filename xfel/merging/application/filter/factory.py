from __future__ import absolute_import, division, print_function
from xfel.merging.application.filter.experiment_filter import experiment_filter
from xfel.merging.application.filter.reflection_filter import reflection_filter
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for filtering experiments. """
  @staticmethod
  def from_parameters(params, additional_info=None):
    """ """
    return [experiment_filter(params), reflection_filter(params)]
