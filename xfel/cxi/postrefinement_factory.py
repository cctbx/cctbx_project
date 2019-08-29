from __future__ import absolute_import, division, print_function

class factory(object):
  def __init__(self, params):
    self.params = params

  def postrefinement_algorithm(self,):
    if self.params.postrefinement.enable:
      if self.params.postrefinement.algorithm in ['rs','eta_deff']:
        from xfel.cxi.postrefinement_legacy_rs import legacy_rs as postrefinement_algorithm
      elif self.params.postrefinement.algorithm in ['rs2']:
        from xfel.cxi.postrefinement_updated_rs import updated_rs as postrefinement_algorithm
      elif self.params.postrefinement.algorithm in ['rs_hybrid']:
        from xfel.cxi.postrefinement_hybrid_rs import rs_hybrid as postrefinement_algorithm
      return postrefinement_algorithm
    return None

  @staticmethod
  def insert_frame_call(kwargs):
    # legacy behavior
    if kwargs["self"].params.postrefinement.enable==False or \
       kwargs["self"].params.postrefinement.algorithm in ["rs","eta_deff"]:

       return kwargs["db_mgr"].insert_frame_legacy(
         **dict(
           [ (k,kwargs[k]) for k in ["result","wavelength","corr","slope","offset","data"] ]
             )
         )
    # new behavior for rs2 and rs_hybrid
    elif kwargs["self"].params.postrefinement.algorithm in ["rs2", "rs_hybrid"]:
       return kwargs["db_mgr"].insert_frame_updated(
         **dict(
           [ (k,kwargs[k]) for k in ["result","wavelength","data","postx"] ]
             )
         )
