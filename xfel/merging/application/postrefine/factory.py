from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import factory as factory_base

class factory(factory_base):
  """ Factory class for post-refining experiments. """
  @staticmethod
  def from_parameters(params, additional_info=None, mpi_helper=None, mpi_logger=None):
    """ """
    if params.postrefinement.algorithm == 'ratchet1':
      from xfel.merging.application.postrefine.ratchet1_postrefinement import ratchet1_postrefinement as postrefinement
    else:
      from xfel.merging.application.postrefine.postrefinement import postrefinement
    return [postrefinement(params, mpi_helper, mpi_logger)]

    '''
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
