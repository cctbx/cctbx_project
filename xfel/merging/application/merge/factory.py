from __future__ import absolute_import, division, print_function
from xfel.merging.application.merge.merger import merger
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """Factory class for calculating averaged intensities of symmetry-reduced HKLs."""
  @staticmethod
  def from_parameters(params, additional_info=None):
    return [merger(params)]
