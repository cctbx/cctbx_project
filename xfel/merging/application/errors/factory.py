from __future__ import absolute_import, division, print_function
from xfel.merging.application.errors.error_modifier import error_modifier
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  '''Factory class for modifying errors of measured intensities'''
  @staticmethod
  def from_parameters(params, additional_info=None):
    return [error_modifier(params)]
