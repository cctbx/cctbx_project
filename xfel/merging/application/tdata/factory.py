from __future__ import absolute_import, division, print_function
from xfel.merging.application.tdata.cell_listing import simple_cell_listing
from xfel.merging.application.worker import factory as factory_base

""" Factory class for tdata cell listing """

class factory(factory_base):
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
     if params.tdata.output_path is not None:
      return [simple_cell_listing(params, mpi_helper, mpi_logger)]
