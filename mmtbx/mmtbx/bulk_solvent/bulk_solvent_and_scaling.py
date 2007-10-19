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
  b_sol_max = 150.0
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
  k_sol_step = 0.3
    .type = float
    .expert_level=2
  b_sol_step = 30.0
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
  apply_back_trace_of_b_cart = False
    .type = bool
    .expert_level=2
  verbose = -1
    .type = int
    .expert_level=3
  ignore_bulk_solvent_and_scale_failure = False
    .type = bool
    .expert_level = 3
""")

def k_sol_b_sol_b_cart_minimizer(fmodel,
                                 params,
                                 refine_k_sol = False,
                                 refine_b_sol = False,
                                 refine_b_cart = False):
  if(fmodel.target_name == "ml"):
    alpha, beta = fmodel.alpha_beta_w()
    alpha_data, beta_data = alpha.data(), beta.data()
  elif(fmodel.target_name == "ls_wunit_k1"):
    alpha_data, beta_data = None, None
  else:
    raise RuntimeError("requested target for aniso scaling is not available")
  if(refine_b_cart): symm_constr = params.symmetry_constraints_on_b_cart
  else: symm_constr = False
  return uaniso_ksol_bsol_scaling_minimizer(
         fc            = fmodel.f_calc_w(),
         fo            = fmodel.f_obs_w,
         fm            = fmodel.f_mask_w(),
         k_initial     = fmodel.k_sol(),
         b_initial     = fmodel.b_sol(),
         u_initial     = fmodel.b_cart(),
         scale_initial = 1.0,
         refine_k      = refine_k_sol,
         refine_b      = refine_b_sol,
         refine_u      = refine_b_cart,
         refine_scale  = False,
         alpha         = alpha_data,
         beta          = beta_data,
         max_iterations= params.max_iterations,
         min_iterations= params.min_iterations,
         lbfgs_exception_handling_params = lbfgs.exception_handling_parameters(
                         ignore_line_search_failed_step_at_lower_bound = True,
                         ignore_line_search_failed_step_at_upper_bound = True,
                         ignore_line_search_failed_maxfev              = True),
         symmetry_constraints_on_b_cart = symm_constr)

class uaniso_ksol_bsol_scaling_minimizer(object):
  def __init__(self,
               fc,
               fo,
               fm,
               k_initial,
               b_initial,
               u_initial,
               scale_initial,
               refine_k,
               refine_b,
               refine_u,
               refine_scale,
               min_iterations,
               max_iterations,
               alpha = None,
               beta = None,
               lbfgs_exception_handling_params = None,
               symmetry_constraints_on_b_cart = False,
               u_min_max = 500.,
               u_min_min =-500.,
               k_sol_max = 10.,
               k_sol_min =-10.,
               b_sol_max = 500.,
               b_sol_min =-500.):
    adopt_init_args(self, locals())
    self.xxx = 0
    assert self.fc.indices().all_eq(self.fm.indices()) == 1
    assert self.fc.indices().all_eq(self.fo.indices()) == 1
    self.gradient_flags = [self.refine_k,self.refine_b,self.refine_u]
    self.sg = fc.space_group()
    self.fc = fc.data()
    self.fo = fo.data()
    self.fm = fm.data()
    self.uc = fo.unit_cell()
    self.hkl = fo.indices()
    self.k_min = self.k_initial
    self.b_min = self.b_initial
    self.u_min = self.u_initial
    self.scale_min = self.scale_initial
    self.flag = (self.alpha is not None and self.beta is not None)
    if(self.flag):
      self.eps = fc.epsilons().data()
      self.cs  = fc.centric_flags().data()
    ################################
    if(self.symmetry_constraints_on_b_cart == True):
       self.adp_constraints = self.sg.adp_constraints()
       u_star = adptbx.u_cart_as_u_star(
         self.uc, self.sg.average_u_star(u_star= self.u_min))
       self.dim_u = self.adp_constraints.n_independent_params()
       assert self.dim_u <= 6
       independent_params = self.adp_constraints.independent_params(u_star)
       self.u_factor = self.uc.volume()**(2/3.)
       self.x = self.pack(
         u=independent_params,
         k=self.k_min,
         b=self.b_min,
         scale=self.scale_min,
         u_factor=self.u_factor)
    else:
       self.u_factor = 1.0
       self.dim_u = len(self.u_initial)
       assert self.dim_u == 6
       self.x = self.pack(
         u=flex.double(self.u_min),
         k=self.k_min,
         b=self.b_min,
         scale=self.scale_min,
         u_factor=self.u_factor)
    ################################
    self.minimizer = lbfgs.run(
                             target_evaluator = self,
                             core_params = lbfgs.core_parameters(),
                             termination_params = lbfgs.termination_parameters(
                                  min_iterations = min_iterations,
                                  max_iterations = max_iterations),
                                  exception_handling_params =
                                   self.lbfgs_exception_handling_params
                              )
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, u, k, b, scale, u_factor):
    v = []
    if (self.refine_u): v += [ui*u_factor for ui in u]
    if (self.refine_k): v.append(k)
    if (self.refine_b): v.append(b)
    if (self.refine_scale): v.append(scale)
    return flex.double(v)

  def unpack_x(self):
    i = 0
    if (self.refine_u):
      if(self.symmetry_constraints_on_b_cart == True):
         self.u_min = adptbx.u_star_as_u_cart(self.uc,
           list(self.adp_constraints.all_params(
             iter(self.x[i:self.dim_u]/self.u_factor))))
      else:
         self.u_min = tuple(self.x)[i:self.dim_u]
      i = self.dim_u
    if (self.refine_k):
      self.k_min = self.x[i]
      i += 1
    if (self.refine_b):
      self.b_min = self.x[i]
      i += 1
    if (self.refine_scale):
      self.scale_min = self.x[i]

  def compute_functional_and_gradients(self):
    self.unpack_x()
    self.u_min = tuple([max(self.u_min_min, min(self.u_min_max, v))
      for v in self.u_min])
    if(self.b_min > self.b_sol_max): self.b_min = self.b_sol_max
    if(self.b_min < self.b_sol_min): self.b_min = self.b_sol_min
    if(self.k_min > self.k_sol_max): self.b_min = self.k_sol_max
    if(self.k_min < self.k_sol_min): self.b_min = self.k_sol_min
    if(self.flag):
      manager = bulk_solvent.target_gradients_aniso_ml(
               self.fo,
               self.fc,
               self.fm,
               self.u_min,
               self.k_min,
               self.b_min,
               self.hkl,
               self.uc,
               self.sg,
               flex.bool(self.gradient_flags),
               self.alpha,
               self.beta,
               self.scale_min)
    else:
      manager = bulk_solvent.target_gradients_aniso(
                                       self.fo,
                                       self.fc,
                                       self.fm,
                                       self.u_min,
                                       self.k_min,
                                       self.b_min,
                                       self.hkl,
                                       self.uc,
                                       self.refine_u,
                                       self.refine_k,
                                       self.refine_b)
    self.f = manager.target()
    gk = manager.grad_ksol()
    gb = manager.grad_bsol()
    try: gscale = manager.grad_k()
    except KeyboardInterrupt: pass
    except: gscale = 0.0
    if(self.symmetry_constraints_on_b_cart == True):
       gu = adptbx.grad_u_cart_as_u_star(self.uc, list(manager.grad_b_cart()))
       independent_params = flex.double(
         self.adp_constraints.independent_gradients(all_gradients=gu))
       self.g = self.pack(
         u=independent_params,
         k=gk,
         b=gb,
         scale=gscale,
         u_factor=1/self.u_factor)
    else:
       gu = manager.grad_b_cart()
       self.g = self.pack(u=gu, k=gk, b=gb, scale=gscale, u_factor=1.0)
    return self.f, self.g

################################################################################

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
         assert abs(flex.max(flex.abs(fmodel.f_mask().data()))) > 1.e-3
       macro_cycles = range(1, params.number_of_macro_cycles+1)
       self.show(fmodel = fmodel, message = m+str(0)+" (start)")
       if(params.fix_k_sol is not None):
         assert params.bulk_solvent
         assert not params.k_sol_b_sol_grid_search
         assert not params.minimization_k_sol_b_sol
         fmodel.update(k_sol = params.fix_k_sol)
       if(params.fix_b_sol is not None):
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
       if(not params.apply_back_trace_of_b_cart):
         self.ERROR_MESSAGE(status=_approx_le(fmodel.r_work(), start_target))
       if(abs(fmodel.k_sol()) < 0.01):
         fmodel.update(k_sol = 0.0, b_sol = 0.0)

  def show(self, fmodel, message):
    if(self.params.verbose > 0):
      fmodel.info().show_rfactors_targets_scales_overall(
        header = message, out = self.log)

  def _ksol_bsol_grid_search(self, fmodel):
    save_target = None
    if(fmodel.target_name != "ls_wunit_k1"):
      save_target = fmodel.target_name
      fmodel.update(target_name = "ls_wunit_k1")
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
        #fmodel.update(b_cart = [0,0,0,0,0,0])
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
    if(save_target is not None): fmodel.update(target_name = save_target)
    return final_ksol, final_bsol, final_b_cart, final_r_work

  def _ksol_bsol_cart_minimizer(self, fmodel):
    original_target = fmodel.target_name
    start_r_work = fmodel.r_work()
    save_target = None
    final_ksol = fmodel.k_sol()
    final_bsol = fmodel.b_sol()
    final_r_work = fmodel.r_work()
    if(fmodel.target_name != "ls_wunit_k1"):
      save_target = "ml"
      fmodel.update(target_name = "ls_wunit_k1")
    ksol, bsol = self._k_sol_b_sol_minimization_helper(fmodel = fmodel)
    fmodel.update(k_sol = ksol, b_sol = bsol)
    r_work = fmodel.r_work()
    if(r_work < final_r_work):
      final_ksol = ksol
      final_bsol = bsol
      final_r_work = r_work
    if(save_target is not None):
      fmodel.update(target_name = save_target)
      ksol, bsol = self._k_sol_b_sol_minimization_helper(fmodel = fmodel)
      fmodel.update(k_sol = ksol, b_sol = bsol)
      r_work = fmodel.r_work()
      if(r_work < final_r_work):
        final_ksol = ksol
        final_bsol = bsol
        final_r_work = r_work
    self.ERROR_MESSAGE(status=_approx_le(final_r_work, start_r_work))
    fmodel.update(target_name = original_target)
    return final_ksol, final_bsol, final_r_work

  def _b_cart_minimizer(self, fmodel):
    original_target = fmodel.target_name
    start_r_work = fmodel.r_work()
    save_target = None
    final_b_cart = fmodel.b_cart()
    final_r_work = fmodel.r_work()
    if(fmodel.target_name != "ls_wunit_k1"):
      save_target = "ml"
      fmodel.update(target_name = "ls_wunit_k1")
    b_cart = self._b_cart_minimizer_helper(fmodel = fmodel)
    fmodel.update(b_cart = b_cart)
    r_work = fmodel.r_work()
    if(r_work < final_r_work):
      final_b_cart = b_cart
      final_r_work = r_work
    if(save_target is not None):
      fmodel.update(target_name = save_target)
      b_cart = self._b_cart_minimizer_helper(fmodel = fmodel)
      fmodel.update(b_cart = b_cart)
      r_work = fmodel.r_work()
      if(r_work < final_r_work):
        final_b_cart = b_cart
        final_r_work = r_work
    self.ERROR_MESSAGE(status=_approx_le(final_r_work, start_r_work))
    fmodel.update(target_name = original_target)
    return final_b_cart, final_r_work

  def _b_cart_minimizer_helper(self, fmodel, n_macro_cycles = 2):
    r_start = fmodel.r_work()
    b_start = fmodel.b_cart()
    for u_cycle in xrange(n_macro_cycles):
      b_cart = k_sol_b_sol_b_cart_minimizer(fmodel = fmodel,
        params = self.params, refine_b_cart = True).u_min
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
