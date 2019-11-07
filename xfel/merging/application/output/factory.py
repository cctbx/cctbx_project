from __future__ import division
from xfel.merging.application.output.output import output
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """Factory class for outputting merged multiple measurements of symmetry-reduced HKLs"""
  @staticmethod
  def from_parameters(params, additional_info=None):
    return [output(params)]
