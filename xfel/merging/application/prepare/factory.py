from __future__ import absolute_import, division, print_function
from xfel.merging.application.prepare.prepare_spread import prepare_spread
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """Factory class for preparing data for additional analysis."""
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    assert additional_info[0] in ['spread',]
    return [prepare_spread(params, mpi_helper, mpi_logger)]
