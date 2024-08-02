from __future__ import absolute_import, division, print_function
from xfel.merging.application.filter.experiment_filter import experiment_filter
from xfel.merging.application.filter.reflection_filter import reflection_filter
from xfel.merging.application.filter.global_filter import GlobalFilter
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for filtering experiments. """
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    """ """
    workers = []
    if additional_info and additional_info[0] == 'global':
      workers.append(GlobalFilter(params, mpi_helper, mpi_logger))
    elif additional_info:
      raise KeyError('Unknown worker: filter_' + '_'.join(additional_info))
    else:  # if not additional_info
      if params.filter.algorithm is not None:
        workers.append(experiment_filter(params, mpi_helper, mpi_logger))
      if params.select.algorithm is not None:
        workers.append(reflection_filter(params, mpi_helper, mpi_logger))
    return workers
