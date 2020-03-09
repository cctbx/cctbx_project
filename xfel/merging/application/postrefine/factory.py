from __future__ import absolute_import, division, print_function
from xfel.merging.application.postrefine.postrefinement_rs import postrefinement_rs
from xfel.merging.application.postrefine.postrefinement_rs2 import postrefinement_rs2
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for post-refining experiments. """
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    """ """
    if params.postrefinement.algorithm in ['rs','eta_deff']:
      return [postrefinement_rs(params, mpi_helper, mpi_logger)]
    elif params.postrefinement.algorithm in ['rs2']:
      return [postrefinement_rs2(params, mpi_helper, mpi_logger)]

    '''
    FROM CXI-MERGE
    if params.postrefinement.enable:
      if params.postrefinement.algorithm in ['rs','eta_deff']:
        from xfel.cxi.postrefinement_legacy_rs import legacy_rs as postrefinement
      elif params.postrefinement.algorithm in ['rs2']:
        from xfel.cxi.postrefinement_updated_rs import updated_rs as postrefinement_algorithm
      elif self.params.postrefinement.algorithm in ['rs_hybrid']:
        from xfel.cxi.postrefinement_hybrid_rs import rs_hybrid as postrefinement_algorithm
      return [postrefinement_algorithm]
    return None
    '''
