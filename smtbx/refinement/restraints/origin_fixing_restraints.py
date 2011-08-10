from scitbx.array_family import flex
from cctbx.eltbx import tiny_pse
from smtbx.refinement.restraints import origin_fixing as base


class homogeneous_weighting(base):

  def weights(self,
              normal_eqns,
              jacobian_transpose_matching_grad_fc,
              params):
    z_max = max([
      tiny_pse.table(p.scatterer.element_symbol()).atomic_number()
      for p in params ])
    return flex.double(params.size(), z_max**2)


class atomic_number_weighting(base):

  def weights(self,
              normal_eqns,
              jacobian_transpose_matching_grad_fc,
              params):
    w = flex.double([
      tiny_pse.table(p.scatterer.element_symbol()).atomic_number()
      for p in params ])**2
    return w


class flack_schwarzenbach_weighting(base):

  def weights(self,
            normal_eqns,
            jacobian_transpose_matching_grad_fc,
            params):
    raise NotImplementedError()
