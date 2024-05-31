from __future__ import absolute_import, division, print_function

from xfel.merging.application.monitor.monitor import MonitorWorker
from xfel.merging.application.worker import factory as factory_base


class factory(factory_base):
  """Factory class for monitoring CPU and GPU resources in the background"""
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    return [MonitorWorker(params, mpi_helper, mpi_logger)]
