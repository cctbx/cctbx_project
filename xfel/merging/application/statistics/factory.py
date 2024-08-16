from __future__ import absolute_import, division, print_function
from xfel.merging.application.statistics.unit_cell_statistics import unit_cell_statistics
from xfel.merging.application.statistics.beam_statistics import beam_statistics
from xfel.merging.application.statistics.intensity_resolution_statistics import intensity_resolution_statistics
from xfel.merging.application.statistics.intensity_resolution_statistics_cxi import intensity_resolution_statistics_cxi
from xfel.merging.application.statistics.deltaccint import deltaccint
from xfel.merging.application.statistics.experiment_resolution_statistics import experiment_resolution_statistics
from xfel.merging.application.statistics.intensity_histogram import intensity_histogram
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for calculating statistics of merged measurements. """
  @staticmethod
  def from_parameters(params, additional_info=[], mpi_helper=None, mpi_logger=None):
    """ """
    info_count = len(additional_info)
    assert info_count > 0
    if additional_info[0] == 'unitcell':
      return [unit_cell_statistics(params, mpi_helper, mpi_logger)]
    elif additional_info[0] == 'beam':
      return [beam_statistics(params, mpi_helper, mpi_logger)]
    elif additional_info[0] == 'resolution':
      return [experiment_resolution_statistics(params, mpi_helper, mpi_logger)]
    elif additional_info[0] == 'intensity':
      if info_count == 1:
        return [intensity_resolution_statistics(params, mpi_helper, mpi_logger)]
      elif info_count > 1 and additional_info[1] == 'cxi':
        return [intensity_resolution_statistics_cxi(params, mpi_helper, mpi_logger)]
      elif info_count > 1 and additional_info[1] == 'histogram':
        return [intensity_histogram(params, mpi_helper, mpi_logger)]
    elif additional_info[0] == 'deltaccint':
      return [deltaccint(params, mpi_helper, mpi_logger)]

