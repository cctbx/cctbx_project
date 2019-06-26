from __future__ import absolute_import, division, print_function
from xfel.merging.application.model.crystal_model import crystal_model
from xfel.merging.application.model.resolution_binner import resolution_binner
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  '''Factory class for creating crystal model and resolution binner'''
  @staticmethod
  def from_parameters(params, additional_info=[], mpi_helper=None, mpi_logger=None):
    assert len(additional_info) > 0
    return [crystal_model(params, additional_info[0], mpi_helper, mpi_logger), resolution_binner(params, mpi_helper, mpi_logger)]
