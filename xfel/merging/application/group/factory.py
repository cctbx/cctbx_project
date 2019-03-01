from __future__ import division
from xfel.merging.application.group.group_reflections import hkl_group
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  '''Factory class for grouping all measurements of an asu hkl from all ranks at a single rank'''
  @staticmethod
  def from_parameters(params):
    """ """
    return [hkl_group(params)]
