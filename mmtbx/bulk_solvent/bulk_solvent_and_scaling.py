from cctbx.array_family import flex
import math, time
from cctbx import miller
from cctbx import crystal
from cctbx import adptbx
from scitbx import lbfgs
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
from mmtbx import bulk_solvent
from cctbx import xray
from mmtbx import max_lik
from mmtbx.max_lik import maxlik
from mmtbx.refinement import print_statistics
from scitbx import matrix
from mmtbx.max_lik import max_like_non_uniform
import iotbx.phil
from libtbx.utils import Sorry

import boost.python
ext = boost.python.import_ext("mmtbx_f_model_ext")

master_params = iotbx.phil.parse("""\
  bulk_solvent = True
    .type = bool
  anisotropic_scaling = True
    .type = bool
  k_sol_b_sol_grid_search = True
    .type = bool
    .expert_level=2
  minimization_k_sol_b_sol = True
    .type = bool
    .expert_level=2
  minimization_b_cart = True
    .type = bool
    .expert_level=2
  target = ls_wunit_k1 *ml
    .type = choice
  symmetry_constraints_on_b_cart = True
    .type = bool
    .expert_level=2
  k_sol_max = 0.6
    .type = float
    .expert_level=2
  k_sol_min = 0.0
    .type = float
    .expert_level=2
  b_sol_max = 300.0
    .type = float
    .expert_level=2
  b_sol_min = 0.0
    .type = float
    .expert_level=2
  k_sol_grid_search_max = 0.6
    .type = float
    .expert_level=2
  k_sol_grid_search_min = 0.0
    .type = float
    .expert_level=2
  b_sol_grid_search_max = 80.0
    .type = float
    .expert_level=2
  b_sol_grid_search_min = 20.0
    .type = float
    .expert_level=2
  k_sol_step = 0.2
    .type = float
    .expert_level=2
  b_sol_step = 20.0
    .type = float
    .expert_level=2
  number_of_macro_cycles = 2
    .type = int
    .expert_level=2
  max_iterations = 25
    .type = int
    .expert_level=3
  min_iterations = 25
    .type = int
    .expert_level=3
  fix_k_sol = None
    .type = float
    .expert_level=2
  fix_b_sol = None
    .type = float
    .expert_level=2
  fix_b_cart
    .expert_level=2
    .style = box
  {
    b11 = None
      .type = float
    b22 = None
      .type = float
    b33 = None
      .type = float
    b12 = None
      .type = float
    b13 = None
      .type = float
    b23 = None
      .type = float
  }
  apply_back_trace_of_b_cart = True
    .type = bool
    .expert_level=2
  verbose = -1
    .type = int
    .expert_level=3
  ignore_bulk_solvent_and_scale_failure = False
    .type = bool
    .expert_level = 3
""")

class kbu_minimizer(object):
  def __init__(self,
               fmodel_core_data,
               f_obs,
               k_initial,
               b_initial,
               u_initial,
               refine_k,
               refine_b,
               refine_u,
               min_iterations,
               max_iterations,
               fmodel_core_data1 = None,
               twin_fraction = None,
               symmetry_constraints_on_b_cart = True,
               u_min_max = 500.,
               u_min_min =-500.,
               k_sol_max = 10.,
               k_sol_min =-10.,
               b_sol_max = 500.,
               b_sol_min =-500.):
    adopt_init_args(self, locals())
    assert [fmodel_core_data1,twin_fraction].count(None) in [0,2]
    self.k_min = self.k_initial
    self.b_min = self.b_initial
    self.u_min = self.u_initial
    self.u_factor = self.fmodel_core_data.uc.volume()**(2/3.)
    if(self.symmetry_constraints_on_b_cart):
      self.adp_constraints = self.f_obs.space_group().adp_constraints()
      u_star = self.f_obs.space_group().average_u_star(u_star = self.u_initial)
      self.dim_u = self.adp_constraints.n_independent_params()
      assert self.dim_u <= 6
      independent_params = self.adp_constraints.independent_params(u_star)
      self.x = self.pack(
        u=independent_params,
        k=self.k_min,
        b=self.b_min,
        u_factor=self.u_factor)
    else:
      self.dim_u = len(self.u_initial)
      assert self.dim_u == 6
      self.x = self.pack(
        u=flex.double(self.u_min),
        k=self.k_min,
        b=self.b_min,
        u_factor=self.u_factor)
    lbfgs_exception_handling_params = lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True,
      ignore_line_search_failed_step_at_upper_bound = True,
      ignore_line_search_failed_maxfev              = True)
    self.minimizer = lbfgs.run(
      target_evaluator = self,
      core_params = lbfgs.core_parameters(),
      termination_params = lbfgs.termination_parameters(
        min_iterations            = min_iterations,
        max_iterations            = max_iterations),
        exception_handling_params = lbfgs_exception_handling_params)
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, u, k, b, u_factor):
    v = []
    if (self.refine_u): v += [ui*u_factor for ui in u]
    if (self.refine_k): v.append(k)
    if (self.refine_b): v.append(b)
    return flex.double(v)

  def unpack_x(self):
    i = 0
    if(self.refine_u):
      if(self.symmetry_constraints_on_b_cart):
        self.u_min = list(self.adp_constraints.all_params(
          iter(self.x[i:self.dim_u]/self.u_factor)))
      else:
        self.u_min = list(iter(self.x[i:self.dim_u]/self.u_factor))
      i = self.dim_u
    if(self.refine_k):
      self.k_min = self.x[i]
      i += 1
    if(self.refine_b):
      self.b_min = self.x[i]

  def compute_functional_and_gradients(self):
    self.unpack_x()
    self.u_min = tuple([max(self.u_min_min, min(self.u_min_max, v))
      for v in self.u_min])
    if(self.b_min > self.b_sol_max): self.b_min = self.b_sol_max
    if(self.b_min < self.b_sol_min): self.b_min = self.b_sol_min
    if(self.k_min > self.k_sol_max): self.k_min = self.k_sol_max
    if(self.k_min < self.k_sol_min): self.k_min = self.k_sol_min
    if(self.twin_fraction is None):
      self.fmodel_core_data.update(k_sol = self.k_min, b_sol = self.b_min,
        u_star = self.u_min)
      tg = bulk_solvent.bulk_solvent_and_aniso_scale_target_and_grads_ls(
        fm                  = self.fmodel_core_data.data,
        fo                  = self.f_obs.data(),
        compute_k_sol_grad  = self.refine_k,
        compute_b_sol_grad  = self.refine_b,
        compute_u_star_grad = self.refine_u)
    else:
      self.fmodel_core_data.update(k_sol = self.k_min, b_sol = self.b_min,
        u_star = self.u_min)
      self.fmodel_core_data1.update(k_sol = self.k_min, b_sol = self.b_min,
        u_star = self.u_min)
      tg = bulk_solvent.bulk_solvent_and_aniso_scale_target_and_grads_ls(
        fm1                 = self.fmodel_core_data.data,
        fm2                 = self.fmodel_core_data1.data,
        twin_fraction       = self.twin_fraction,
        fo                  = self.f_obs.data(),
        compute_k_sol_grad  = self.refine_k,
        compute_b_sol_grad  = self.refine_b,
        compute_u_star_grad = self.refine_u)
    self.f = tg.target()
    gk=0
    gb=0
    gu=[0,0,0,0,0,0]
    if(self.refine_k or self.refine_b):
      gk = tg.grad_k_sol()
      gb = tg.grad_b_sol()
    if(self.refine_u): gu = list(tg.grad_u_star())
    if(self.symmetry_constraints_on_b_cart and self.refine_u):
      independent_params = flex.double(
        self.adp_constraints.independent_gradients(all_gradients=gu))
      self.g = self.pack(
        u=independent_params,
        k=gk,
        b=gb,
        u_factor=1/self.u_factor)
    else:
      self.g = self.pack(u=gu, k=gk, b=gb, u_factor=1/self.u_factor)
    return self.f, self.g

def k_sol_b_sol_b_cart_minimizer(
      fmodel,
      params        = None,
      refine_k_sol  = False,
      refine_b_sol  = False,
      refine_u_star = False):
  if(params is None): params = master_params.extract()
  fmodel_core_data_work = fmodel.core_data_work()
  return kbu_minimizer(
    fmodel_core_data = fmodel_core_data_work,
    f_obs            = fmodel.f_obs_w,
    k_initial        = fmodel_core_data_work.k_sol,
    b_initial        = fmodel_core_data_work.b_sol,
    u_initial        = fmodel_core_data_work.u_star,
    refine_k         = refine_k_sol,
    refine_b         = refine_b_sol,
    refine_u         = refine_u_star,
    max_iterations   = params.max_iterations,
    min_iterations   = params.min_iterations,
    symmetry_constraints_on_b_cart = params.symmetry_constraints_on_b_cart)

def _approx_le(x,y):
  x = float("%.4f"%x)
  y = float("%.4f"%y)
  assert x <= y
  return x <= y

def _approx_lt(x,y):
  x = float("%.4f"%x)
  y = float("%.4f"%y)
  assert x < y
  return x < y

def _extract_fix_b_cart(fix_b_cart_scope):
  fbs = fix_b_cart_scope
  b_cart = [fbs.b11,fbs.b22,fbs.b33,fbs.b12,fbs.b13,fbs.b23]
  if(b_cart.count(None) > 0): return None
  else: return b_cart


class bulk_solvent_and_scales(object):

  def __init__(self, fmodel, params = None, log = None):
    start_target = fmodel.r_work()
    self.params = params
    self.log = log
    if(params is None): params = master_params.extract()
    if([params.bulk_solvent, params.anisotropic_scaling].count(True) > 0):
       if(not params.bulk_solvent):
         params.k_sol_b_sol_grid_search = False
         params.minimization_k_sol_b_sol = False
       if(not params.anisotropic_scaling):
         params.minimization_b_cart = False
       params_target = params.target
       fmodel_target = fmodel.target_name
       if(fmodel.alpha_beta_params is not None):
          save_interpolation_flag = fmodel.alpha_beta_params.interpolation
          fmodel.alpha_beta_params.interpolation = False
       m = "macro_cycle= "
       if(params.bulk_solvent):
         assert not fmodel.check_f_mask_all_zero()
       macro_cycles = range(1, params.number_of_macro_cycles+1)
       self.show(fmodel = fmodel, message = m+str(0)+" (start)")
       mask_ok = not fmodel.check_f_mask_all_zero()
       if(params.fix_k_sol is not None and mask_ok):
         print params.bulk_solvent
         assert params.bulk_solvent
         assert not params.k_sol_b_sol_grid_search
         assert not params.minimization_k_sol_b_sol
         fmodel.update(k_sol = params.fix_k_sol)
       if(params.fix_b_sol is not None and mask_ok):
         assert params.bulk_solvent
         assert not params.k_sol_b_sol_grid_search
         assert not params.minimization_k_sol_b_sol
         fmodel.update(b_sol = params.fix_b_sol)
       fix_b_cart = _extract_fix_b_cart(fix_b_cart_scope = params.fix_b_cart)
       if(fix_b_cart is not None):
         assert params.anisotropic_scaling
         assert not params.minimization_b_cart
         fmodel.update(b_cart = fix_b_cart)
       target = fmodel.r_work()
       for mc in macro_cycles:
         if(params.k_sol_b_sol_grid_search and mc == macro_cycles[0] and
            not self._is_within_grid_search(fmodel = fmodel)):
           ksol,bsol,b_cart,target = self._ksol_bsol_grid_search(fmodel=fmodel)
           fmodel.update(k_sol = ksol, b_sol = bsol, b_cart = b_cart)
           self.ERROR_MESSAGE(status=approx_equal(target, fmodel.r_work()))
           if(not params.apply_back_trace_of_b_cart):
             self.ERROR_MESSAGE(status=_approx_le(target, start_target))
           self.show(fmodel = fmodel, message=m+str(mc)+": k & b: grid search")
         if(params.minimization_k_sol_b_sol):
           ksol, bsol, target = self._ksol_bsol_cart_minimizer(fmodel = fmodel)
           fmodel.update(k_sol = ksol, b_sol = bsol)
           self.ERROR_MESSAGE(status=approx_equal(target, fmodel.r_work()))
           if(not params.apply_back_trace_of_b_cart):
             self.ERROR_MESSAGE(status=_approx_le(target, start_target))
           self.show(fmodel = fmodel,message=m+str(mc)+": k & b: minimization")
         if(params.minimization_b_cart):
           b_cart,target = self._b_cart_minimizer(fmodel = fmodel)
           fmodel.update(b_cart = b_cart)
           if(not params.apply_back_trace_of_b_cart):
             self.ERROR_MESSAGE(status=_approx_le(target, start_target))
           self.ERROR_MESSAGE(status=approx_equal(target, fmodel.r_work()))
           self.show(fmodel = fmodel, message =m+str(mc)+": anisotropic scale")
         if(params.apply_back_trace_of_b_cart and abs(fmodel.b_iso()) > 0.0):
            fmodel.apply_back_b_iso()
            self.show(fmodel = fmodel,
              message = m+str(mc)+": apply back trace(b_cart)")
       if(params.apply_back_trace_of_b_cart and abs(fmodel.b_iso()) > 0.0):
         fmodel.apply_back_b_iso()
         self.show(fmodel = fmodel,
           message = m+str(mc)+": apply back trace(b_cart)")
       fmodel.update(target_name = fmodel_target)
       if(abs(fmodel.k_sol()) < 0.01):
         fmodel.update(k_sol = 0.0, b_sol = 0.0)

  def show(self, fmodel, message):
    if(self.params.verbose > 0):
      fmodel.info().show_rfactors_targets_scales_overall(
        header = message, out = self.log)

  def _ksol_bsol_grid_search(self, fmodel):
    start_r_work = fmodel.r_work()
    final_ksol = fmodel.k_sol()
    final_bsol = fmodel.b_sol()
    final_b_cart = fmodel.b_cart()
    final_r_work = start_r_work
    k_sols = kb_range(self.params.k_sol_grid_search_max,
                      self.params.k_sol_grid_search_min,
                      self.params.k_sol_step)
    b_sols = kb_range(self.params.b_sol_grid_search_max,
                      self.params.b_sol_grid_search_min,
                      self.params.b_sol_step)
    for ksol in k_sols:
      for bsol in b_sols:
        for bc in fmodel.b_cart():
          if(abs(bc) > 300.):
            fmodel.update(b_cart = [0,0,0,0,0,0])
            break
        fmodel.update(k_sol = ksol, b_sol = bsol)
        if(self.params.minimization_k_sol_b_sol):
          ksol_, bsol_, dummy = self._ksol_bsol_cart_minimizer(fmodel = fmodel)
          fmodel.update(k_sol = ksol_, b_sol = bsol_)
        if(self.params.minimization_b_cart):
          b_cart, dummy = self._b_cart_minimizer(fmodel = fmodel)
          fmodel.update(b_cart = b_cart)
        r_work = fmodel.r_work()
        if(r_work < final_r_work):
          final_r_work = r_work
          final_ksol = fmodel.k_sol()
          final_bsol = fmodel.b_sol()
          final_b_cart = fmodel.b_cart()
          fmodel.update(k_sol=final_ksol,b_sol=final_bsol,b_cart=final_b_cart)
    fmodel.update(k_sol  = final_ksol,
                  b_sol  = final_bsol,
                  b_cart = final_b_cart)
    self.ERROR_MESSAGE(status=approx_equal(fmodel.r_work(), final_r_work))
    self.ERROR_MESSAGE(status=_approx_le(fmodel.r_work(), start_r_work))
    return final_ksol, final_bsol, final_b_cart, final_r_work

  def _ksol_bsol_cart_minimizer(self, fmodel):
    start_r_work = fmodel.r_work()
    final_ksol = fmodel.k_sol()
    final_bsol = fmodel.b_sol()
    final_r_work = fmodel.r_work()
    ksol, bsol = self._k_sol_b_sol_minimization_helper(fmodel = fmodel)
    fmodel.update(k_sol = ksol, b_sol = bsol)
    r_work = fmodel.r_work()
    if(r_work < final_r_work):
      final_ksol = ksol
      final_bsol = bsol
      final_r_work = r_work
    self.ERROR_MESSAGE(status=_approx_le(final_r_work, start_r_work))
    return final_ksol, final_bsol, final_r_work

  def _b_cart_minimizer(self, fmodel):
    start_r_work = fmodel.r_work()
    final_b_cart = fmodel.b_cart()
    final_r_work = fmodel.r_work()
    b_cart = self._b_cart_minimizer_helper(fmodel = fmodel)
    fmodel.update(b_cart = b_cart)
    r_work = fmodel.r_work()
    if(r_work < final_r_work):
      final_b_cart = b_cart
      final_r_work = r_work
    self.ERROR_MESSAGE(status=_approx_le(final_r_work, start_r_work))
    return final_b_cart, final_r_work

  def _b_cart_minimizer_helper(self, fmodel, n_macro_cycles = 2):
    r_start = fmodel.r_work()
    b_start = fmodel.b_cart()
    for u_cycle in xrange(n_macro_cycles):
      u_min = k_sol_b_sol_b_cart_minimizer(fmodel = fmodel,
        params = self.params, refine_u_star = True).u_min
      b_cart = adptbx.u_as_b(
        adptbx.u_star_as_u_cart(fmodel.f_obs_w.unit_cell(),u_min))
      fmodel.update(b_cart = b_cart)
    r_final = fmodel.r_work()
    if(r_final >= r_start):
       fmodel.update(b_cart = b_start)
       return b_start
    else: return b_cart

  def _k_sol_b_sol_minimization_helper(self, fmodel):
    ksol_orig, bsol_orig = fmodel.k_sol_b_sol()
    r_start = fmodel.r_work()
    minimizer_obj = k_sol_b_sol_b_cart_minimizer(fmodel = fmodel,
      params = self.params, refine_k_sol = True, refine_b_sol = True)
    ksol, bsol = minimizer_obj.k_min, minimizer_obj.b_min
    if(ksol < self.params.k_sol_min): ksol = self.params.k_sol_min
    if(ksol > self.params.k_sol_max): ksol = self.params.k_sol_max
    if(bsol < self.params.b_sol_min): bsol = self.params.b_sol_min
    if(bsol > self.params.b_sol_max): bsol = self.params.b_sol_max
    fmodel.update(k_sol = ksol, b_sol = bsol)
    r_end = fmodel.r_work()
    if(r_end >= r_start):
       fmodel.update(k_sol = ksol_orig, b_sol = bsol_orig)
       return ksol_orig, bsol_orig
    else: return ksol, bsol

  def _is_within_grid_search(self, fmodel):
    result = False
    ksol, bsol = fmodel.k_sol_b_sol()
    keps = 0.05
    beps = 5.0
    if(ksol >= self.params.k_sol_grid_search_min+keps and
       ksol <= self.params.k_sol_grid_search_max-keps and
       bsol >= self.params.b_sol_grid_search_min+beps and
       bsol <= self.params.b_sol_grid_search_max-beps): result = True
    return result

  def ERROR_MESSAGE(self, status):
    if(not status):
      if(not self.params.ignore_bulk_solvent_and_scale_failure):
        raise Sorry("""

   Internal error in bulk solvent and scaling. To ignore this problem please
   run again with the following keyword:

      ignore_bulk_solvent_and_scale_failure=True

   Please report this to PAfonine@lbl.gov
""")
      else:
        print >> self.log, \
          "WARNING -- Internal error in bulk solvent and scaling."

def kb_range(x_max, x_min, step):
  sc = 1000.
  return [i/sc for i in range(int(x_min*sc), int(x_max*sc)+1, int(step*sc))]
