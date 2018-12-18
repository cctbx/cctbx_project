from __future__ import division
from xfel.merging.application.merge.merger import merger
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for averaging intensities of symmetry-reduced hkl's. """
  @staticmethod
  def from_parameters(params):
    return [merger(params)]
