from __future__ import division
from xfel.merging.application.publish.publisher import publisher
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    return [publisher(params, mpi_helper, mpi_logger)]
