from cctbx.array_family import flex
from scitbx import lbfgs
from libtbx.test_utils import approx_equal
from libtbx import adopt_init_args

class minimizer(object):

  def __init__(self,
        fmodel,
        groups,
        call_back_after_minimizer_cycle=None,
        number_of_minimizer_cycles=3,
        lbfgs_max_iterations=20,
        number_of_finite_difference_tests=0):
    adopt_init_args(self, locals())
    self.alpha_w, self.beta_w = self.fmodel.alpha_beta_w(
      only_if_required_by_target=True)
    self.x = flex.double()
    for group in groups:
      if (not group.fix_f_prime): self.x.append(group.f_prime)
      if (not group.fix_f_double_prime): self.x.append(group.f_double_prime)
    fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    for group in groups:
      if (not group.fix_f_prime):
        fmodel.xray_structure.scatterers().flags_set_grad_fp(
          iselection=group.iselection)
      if (not group.fix_f_double_prime):
        fmodel.xray_structure.scatterers().flags_set_grad_fdp(
          iselection=group.iselection)
    for self.i_cycle in xrange(number_of_minimizer_cycles):
      self.lbfgs = lbfgs.run(
        target_evaluator=self,
        termination_params=lbfgs.termination_parameters(
          max_iterations=lbfgs_max_iterations),
        exception_handling_params=lbfgs.exception_handling_parameters(
          ignore_line_search_failed_step_at_lower_bound = True))
      if (call_back_after_minimizer_cycle is not None):
        self.unpack()
        if (not call_back_after_minimizer_cycle(minimizer=self)):
          break
    if (call_back_after_minimizer_cycle is None):
      self.unpack()
    del self.i_cycle
    del self.lbfgs
    del self.x
    del self.alpha_w
    del self.beta_w
    del self.fmodel
    del self.groups

  def unpack(self):
    xi = iter(self.x)
    for group in self.groups:
      if (not group.fix_f_prime): group.f_prime = xi.next()
      if (not group.fix_f_double_prime): group.f_double_prime = xi.next()
    for group in self.groups:
      group.copy_to_scatterers_in_place(
        scatterers=self.fmodel.xray_structure.scatterers())
    self.fmodel.update_xray_structure(update_f_calc=True)

  def compute_functional_and_gradients(self):
    self.unpack()
    fmodel = self.fmodel
    t_w = fmodel.xray_target_functor_result(
      compute_gradients=True,
      alpha=self.alpha_w,
      beta=self.beta_w,
      flag="work")
    f = t_w.target()
    sfg = fmodel.structure_factor_gradients_w(
      u_iso_refinable_params=None,
      d_target_d_f_calc=t_w.derivatives() * fmodel.core.fb_cart_w,
      xray_structure=fmodel.xray_structure,
      n_parameters=0,
      miller_set=fmodel.f_obs_w,
      algorithm=fmodel.sf_algorithm)
    d_t_d_fp = sfg.d_target_d_fp()
    d_t_d_fdp = sfg.d_target_d_fdp()
    del sfg
    g = flex.double()
    for group in self.groups:
      if (not group.fix_f_prime):
        g.append(flex.sum(d_t_d_fp.select(group.iselection)))
      if (not group.fix_f_double_prime):
        g.append(flex.sum(d_t_d_fdp.select(group.iselection)))
    if (self.number_of_finite_difference_tests != 0):
      self.number_of_finite_difference_tests -= 1
      g_fin = []
      eps = 1.e-5
      x = self.x
      for i in xrange(x.size()):
        fs = []
        xi0 = x[i]
        for signed_eps in [eps,-eps]:
          x[i] = xi0 + signed_eps
          self.unpack()
          x[i] = xi0
          t_w = fmodel.xray_target_functor_result(
            compute_gradients=False,
            alpha=self.alpha_w,
            beta=self.beta_w,
            flag="work")
          fs.append(t_w.target())
        g_fin.append((fs[0]-fs[1])/(2*eps))
      self.unpack()
      assert approx_equal(g_fin, g)
    return f, g
