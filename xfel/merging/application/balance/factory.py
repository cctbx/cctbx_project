from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import factory as factory_base
from xfel.merging.application.balance.load_balancer import load_balancer

class factory(factory_base):
  """Factory class for balancing input data load."""
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    if params.input.parallel_file_load.balance != None:
      return [load_balancer(params, mpi_helper, mpi_logger)]
    else:
      return []
