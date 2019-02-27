from __future__ import division
from xfel.merging.application.error_model.error_model import error_model
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  '''Factory class for modifying errors of measured intensities'''
  @staticmethod
  def from_parameters(params):
    return [error_model(params)]
