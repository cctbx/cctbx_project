from __future__ import absolute_import, division, print_function
from xfel.merging.application.statistics.unit_cell_statistics import unit_cell_statistics
from xfel.merging.application.statistics.intensity_resolution_statistics import intensity_resolution_statistics
from xfel.merging.application.statistics.experiment_resolution_statistics import experiment_resolution_statistics
from xfel.merging.application.statistics.intensity_histogram import intensity_histogram
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for calculating statistics of merged measurements. """
  @staticmethod
  def from_parameters(params, additional_info=None):
    """ """
    if additional_info == 'experiment':
      return [unit_cell_statistics(params), experiment_resolution_statistics(params)]
    elif additional_info == 'unit_cell':
      return [unit_cell_statistics(params)]
    elif additional_info == 'intensity':
      return [intensity_resolution_statistics(params)]
    elif additional_info == 'intensity_histogram':
      return [intensity_histogram(params)]
