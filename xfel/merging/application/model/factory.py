from __future__ import absolute_import, division, print_function
from xfel.merging.application.model.crystal_model import crystal_model
from xfel.merging.application.model.resolution_binner import resolution_binner
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  '''Factory class for modifying errors of measured intensities'''
  @staticmethod
  def from_parameters(params, additional_info=None):
    return [crystal_model(params), resolution_binner(params)]
