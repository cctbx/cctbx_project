from __future__ import division
from xfel.merging.application.statistics.unit_cell_statistics import unit_cell_statistics
from xfel.merging.application.statistics.intensity_statistics import intensity_statistics
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for calculating statistics of merged measurements. """
  @staticmethod
  def from_parameters(params):
    """ """
    return [unit_cell_statistics(params), intensity_statistics(params)]
