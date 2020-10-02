from __future__ import absolute_import, division, print_function
from xfel.merging.application.integrate.integrate import integrate
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for modification of intensites. """
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    return [integrate(params, mpi_helper, mpi_logger)]
