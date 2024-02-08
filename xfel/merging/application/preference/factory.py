from __future__ import absolute_import, division, print_function
from xfel.merging.application.preference.preference import PreferenceWorker
from xfel.merging.application.worker import factory as factory_base


class factory(factory_base):
  """Factory for PreferenceWorker which evaluates preferential orientation"""
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    return [PreferenceWorker(params, mpi_helper, mpi_logger)]
