from __future__ import absolute_import, division, print_function

from xfel.merging.application.truncate.truncate import Truncate
from xfel.merging.application.worker import factory as factory_base


class factory(factory_base):
  """Factory class for filtering expts/refls based on their aggregated stats"""
  @staticmethod
  def from_parameters(params, additional_info=(),
                      mpi_helper=None, mpi_logger=None):
    """Initiate a new instance of the truncate worker"""
    return [Truncate(params, mpi_helper, mpi_logger)]
