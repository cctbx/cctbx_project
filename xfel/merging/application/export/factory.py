from __future__ import absolute_import, division, print_function
from xfel.merging.application.export.careless import export_careless_mtz
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for export of batch MTZ file. """
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    if params.export.output_label is not None:
      return [export_careless_mtz(params, mpi_helper, mpi_logger)]
    else:
      raise("Must specify output label if using export worker")
