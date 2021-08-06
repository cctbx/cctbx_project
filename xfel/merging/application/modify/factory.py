from __future__ import absolute_import, division, print_function
from xfel.merging.application.modify.polarization import polarization
from xfel.merging.application.modify.reindex_to_reference import reindex_to_reference
from xfel.merging.application.modify.cosym import cosym
from xfel.merging.application.modify.reindex_to_abc import reindex_to_abc
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for modification of intensites. """
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    if additional_info == "reindex_to_reference".split("_"):
      return [reindex_to_reference(params, mpi_helper, mpi_logger)]
    elif additional_info == "cosym".split("_"):
      return [cosym(params, mpi_helper, mpi_logger)]
    elif additional_info == "reindex_to_abc".split("_"):
      return [reindex_to_abc(params, mpi_helper, mpi_logger)]
    else:
      assert not additional_info, additional_info
      return [polarization(params, mpi_helper, mpi_logger)]
