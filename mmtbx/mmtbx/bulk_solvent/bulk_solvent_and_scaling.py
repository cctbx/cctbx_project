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


master_params = iotbx.phil.parse("""\
  bulk_solvent = True
    .type = bool
  anisotropic_scaling = True
    .type = bool
  statistical_solvent = False
    .type = bool
  k_sol_b_sol_grid_search = True
    .type = bool
  minimization_k_sol_b_sol = True
    .type = bool
  minimization_b_cart = True
    .type = bool
  target = ls_wunit_k1 *ml
    .type = choice
  symmetry_constraints_on_b_cart = True
    .type = bool
  k_sol_max = 0.6
    .type = float
  k_sol_min = 0.0
    .type = float
  b_sol_max = 80.0
    .type = float
  b_sol_min = 10.0
    .type = float
  k_sol_grid_search_max = 0.6
    .type = float
  k_sol_grid_search_min = 0.0
    .type = float
  b_sol_grid_search_max = 80.0
    .type = float
  b_sol_grid_search_min = 10.0
    .type = float
  k_sol_step = 0.05
    .type = float
  b_sol_step = 5.0
    .type = float
  number_of_macro_cycles = 2
    .type = int
  number_of_minimization_macro_cycles = 3
    .type = int
  number_of_cycles_for_anisotropic_scaling = 3
    .type = int
  fix_k_sol = None
    .type = float
  fix_b_sol = None
    .type = float
  fix_b_cart {
    u11 = None
      .type = float
    u22 = None
      .type = float
    u33 = None
      .type = float
    u12 = None
      .type = float
    u13 = None
      .type = float
    u23 = None
      .type = float
  }
  apply_back_trace_of_b_cart = False
    .type = bool
  start_minimization_from_k_sol = 0.35
    .type = float
  start_minimization_from_b_sol = 46.0
    .type = float
  start_minimization_from_b_cart {
    u11 = 0.0
      .type = float
    u22 = 0.0
      .type = float
    u33 = 0.0
      .type = float
    u12 = 0.0
      .type = float
    u13 = 0.0
      .type = float
    u23 = 0.0
      .type = float
  }
  nu_fix_n_atoms = None
    .type = float
  nu_fix_b_atoms = None
    .type = float
  verbose = 1
    .type = int
""")


class solvent_and_scale_params(object):
  def __init__(self, params = None,
                     verbose = -1):
    if (params is None): params = master_params.extract()
    else: o=params
    o=params
    self.__dict__.update(params.__dict__)
    if(self.fix_b_cart is not None):
       try: self.fix_b_cart = [self.fix_b_cart.u11,
                                self.fix_b_cart.u22,
                                self.fix_b_cart.u33,
                                self.fix_b_cart.u12,
                                self.fix_b_cart.u13,
                                self.fix_b_cart.u23]
       except: self.fix_b_cart = self.fix_b_cart
    prefix = "solvent_and_scale_params: "
    if(1):#(self.overwrite is not None):
       #o = self.overwrite
       self.bulk_solvent                    = o.bulk_solvent
       self.anisotropic_scaling             = o.anisotropic_scaling
       self.statistical_solvent             = o.statistical_solvent
       self.k_sol_b_sol_grid_search         = o.k_sol_b_sol_grid_search
       self.minimization_k_sol_b_sol        = o.minimization_k_sol_b_sol
       self.minimization_b_cart             = o.minimization_b_cart
       self.target                          = o.target
       self.symmetry_constraints_on_b_cart  = o.symmetry_constraints_on_b_cart
       self.k_sol_grid_search_max           = o.k_sol_grid_search_max
       self.k_sol_grid_search_min           = o.k_sol_grid_search_min
       self.b_sol_grid_search_max           = o.b_sol_grid_search_max
       self.b_sol_grid_search_min           = o.b_sol_grid_search_min
       self.k_sol_max                       = o.k_sol_max
       self.k_sol_min                       = o.k_sol_min
       self.b_sol_max                       = o.b_sol_max
       self.b_sol_min                       = o.b_sol_min
       self.k_sol_step                      = o.k_sol_step
       self.b_sol_step                      = o.b_sol_step
       self.number_of_macro_cycles          = o.number_of_macro_cycles
       self.number_of_minimization_macro_cycles      = \
                                          o.number_of_minimization_macro_cycles
       self.number_of_cycles_for_anisotropic_scaling = \
                                     o.number_of_cycles_for_anisotropic_scaling
       self.fix_k_sol                         = o.fix_k_sol
       self.fix_b_sol                         = o.fix_b_sol
       self.apply_back_trace_of_b_cart       = o.apply_back_trace_of_b_cart
       self.start_minimization_from_k_sol     = o.start_minimization_from_k_sol
       self.start_minimization_from_b_sol     = o.start_minimization_from_b_sol
       self.verbose                           = o.verbose
       self.nu_fix_n_atoms                    = o.nu_fix_n_atoms
       self.nu_fix_b_atoms                    = o.nu_fix_b_atoms
       uasc = o.fix_b_cart
       if(uasc is not None):
          try:
              if([uasc.u11, uasc.u22, uasc.u33, uasc.u12, uasc.u13,
                                                    uasc.u23].count(None) > 0):
                 self.fix_b_cart = None
              else:
                 self.fix_b_cart = [uasc.u11, uasc.u22, uasc.u33, uasc.u12,
                                                            uasc.u13, uasc.u23]
          except:
              self.fix_b_cart = uasc
       else:
          self.fix_b_cart = None
       uasc = o.start_minimization_from_b_cart
       try: self.start_minimization_from_b_cart = \
                   [uasc.u11, uasc.u22, uasc.u33, uasc.u12, uasc.u13, uasc.u23]
       except: self.start_minimization_from_b_cart = uasc

    if(self.target not in ("ls_wunit_k1","ml")):
       raise RuntimeError(prefix+"target name unavailable: "+self.target)

    if(self.bulk_solvent == True):
       if((self.fix_k_sol, self.fix_b_sol) != (None, None)):
          if(self.fix_k_sol is None or self.fix_b_sol is None):
             raise RuntimeError(prefix+"one of k_sol or b_sol is not fixed")
          if(self.minimization_k_sol_b_sol == True or \
                                         self.k_sol_b_sol_grid_search == True):
             raise RuntimeError(prefix+"ambiguous parameter combination")
          self.k_sol_max  = None
          self.k_sol_min  = None
          self.b_sol_max  = None
          self.b_sol_min  = None
          self.k_sol_step = None
          self.b_sol_step = None
          self.k_sol_grid_search_max = None
          self.k_sol_grid_search_min = None
          self.b_sol_grid_search_max = None
          self.b_sol_grid_search_min = None
          self.start_minimization_from_k_sol = None
          self.start_minimization_from_b_sol = None
       elif(self.k_sol_b_sol_grid_search == True):
          assert self.k_sol_grid_search_max > self.k_sol_grid_search_min
          assert self.k_sol_step > 0.0
          assert self.b_sol_grid_search_max > self.b_sol_grid_search_min
          assert self.b_sol_step > 0.0
          self.start_minimization_from_k_sol = None
          self.start_minimization_from_b_sol = None
       elif(self.minimization_k_sol_b_sol == True):
          if(self.number_of_minimization_macro_cycles < 1):
             raise RuntimeError(prefix+\
                                     "number_of_minimization_macro_cycles < 1")
          assert self.start_minimization_from_k_sol is not None
          assert self.start_minimization_from_b_sol is not None
       else: raise RuntimeError(prefix+"ambiguous parameter combination")
    if(self.bulk_solvent == False):
       self.k_sol_b_sol_grid_search                  = False
       self.minimization_k_sol_b_sol                 = False
       self.k_sol_max                                = None
       self.k_sol_min                                = None
       self.b_sol_max                                = None
       self.b_sol_min                                = None
       self.k_sol_step                               = None
       self.b_sol_step                               = None
       self.k_sol_grid_search_max = None
       self.k_sol_grid_search_min = None
       self.b_sol_grid_search_max = None
       self.b_sol_grid_search_min = None
       self.fix_k_sol                                = None
       self.fix_b_sol                                = None
       self.start_minimization_from_k_sol            = None
       self.start_minimization_from_b_sol            = None
    if(self.anisotropic_scaling == False):
       self.minimization_b_cart                     = False
       self.symmetry_constraints_on_b_cart          = False
       self.number_of_cycles_for_anisotropic_scaling = None
       self.fix_b_cart                              = None
       self.apply_back_trace_of_b_cart              = None
       self.start_minimization_from_b_cart          = None
    if(self.anisotropic_scaling == True):
       if(self.minimization_b_cart == True):
          self.fix_b_cart = None
          if(self.number_of_cycles_for_anisotropic_scaling < 1):
             raise RuntimeError(prefix+\
                                "number_of_cycles_for_anisotropic_scaling < 1")
          if(self.number_of_minimization_macro_cycles > 0):
             if(self.minimization_k_sol_b_sol == False):
                self.number_of_minimization_macro_cycles = 0
          if(self.fix_b_cart is not None):
             raise RuntimeError(prefix+"ambiguous parameter combination")
       else:
          if(self.fix_b_cart is None):
             raise RuntimeError(prefix+"no option for anisotropic scaling")
          try: dim = len(self.fix_b_cart)
          except:
               try: dim = self.fix_b_cart.size()
               except:
                    dim = 0
          if(dim != 6): raise RuntimeError(prefix+"wrong fix_b_cart size")
          self.number_of_cycles_for_anisotropic_scaling = None
          self.start_minimization_from_b_cart = None
          self.apply_back_trace_of_b_cart = None
          self.symmetry_constraints_on_b_cart = None
    if(self.statistical_solvent):
       fix_flag = [self.nu_fix_n_atoms,self.nu_fix_b_atoms].count(None)
       assert fix_flag == 2 or fix_flag == 0

def aniso_scale_minimizer(fmodel, symm_constr, alpha=None, beta=None):
  if(fmodel.target_name == "ml"):
     if([alpha,beta].count(None) == 2):
        alpha, beta = fmodel.alpha_beta_w()
        alpha_data, beta_data = alpha.data(), beta.data()
     elif([alpha,beta].count(None) == 0):
        alpha_data, beta_data = alpha.data(), beta.data()
     else:
        raise RuntimeError("aniso_scale_minimizer: alpha & beta are not available")
  elif(fmodel.target_name == "ls_wunit_k1"):
     alpha_data, beta_data = None, None
  else:
     raise RuntimeError("requested target for aniso scaling is not navailable")
  f_c_d = fmodel.f_calc_w().data()
  f_o_d = fmodel.f_ordered_solvent_w().data()
  f_calc_new = miller.array(miller_set = fmodel.f_calc_w(), data = f_c_d + f_o_d)
  return uaniso_ksol_bsol_scaling_minimizer(
         fc            = f_calc_new,
         fo            = fmodel.f_obs_w,
         fm            = fmodel.f_mask_w(),
         k_initial     = fmodel.k_sol_b_sol()[0],
         b_initial     = fmodel.k_sol_b_sol()[1],
         u_initial     = fmodel.b_cart(),
         scale_initial = 1.0,
         refine_k      = False,
         refine_b      = False,
         refine_u      = True,
         refine_scale  = False,
         alpha         = alpha_data,
         beta          = beta_data,
         lbfgs_exception_handling_params = lbfgs.exception_handling_parameters(
                         ignore_line_search_failed_step_at_lower_bound = True,
                         ignore_line_search_failed_step_at_upper_bound = True,
                         ignore_line_search_failed_maxfev              = True),
         symmetry_constraints_on_b_cart = symm_constr).u_min


def k_sol_b_sol_minimizer(fmodel):
  if(fmodel.target_name == "ml"):
     alpha, beta = fmodel.alpha_beta_w()
     alpha_data, beta_data = alpha.data(), beta.data()
  elif(fmodel.target_name == "ls_wunit_k1"):
     alpha_data, beta_data = None, None
  else:
     print fmodel.target_name
     raise RuntimeError("requested target for aniso scaling is not navailable")
  f_c_d = fmodel.f_calc_w().data()
  f_o_d = fmodel.f_ordered_solvent_w().data()
  f_calc_new = miller.array(miller_set = fmodel.f_calc_w(), data = f_c_d + f_o_d)
  manager = uaniso_ksol_bsol_scaling_minimizer(
         fc            = f_calc_new,
         fo            = fmodel.f_obs_w,
         fm            = fmodel.f_mask_w(),
         k_initial     = fmodel.k_sol_b_sol()[0],
         b_initial     = fmodel.k_sol_b_sol()[1],
         u_initial     = fmodel.b_cart(),
         scale_initial = 1.0,
         refine_k      = True,
         refine_b      = True,
         refine_u      = False,
         refine_scale  = False,
         alpha         = alpha_data,
         beta          = beta_data,
         lbfgs_exception_handling_params = lbfgs.exception_handling_parameters(
                         ignore_line_search_failed_step_at_lower_bound = True,
                         ignore_line_search_failed_step_at_upper_bound = True,
                         ignore_line_search_failed_maxfev              = True),
         symmetry_constraints_on_b_cart = False)
  return manager.k_min, manager.b_min


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
               alpha = None,
               beta = None,
               min_iterations=50,
               max_iterations=50,
               lbfgs_exception_handling_params = None,
               symmetry_constraints_on_b_cart = False,
               u_min_max = 1000.,
               u_min_min =-1000.,
               k_sol_max = 10.,
               k_sol_min =-10.,
               b_sol_max = 500.,
               b_sol_min =-500.):
    adopt_init_args(self, locals())
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
                             core_params = lbfgs.core_parameters(
                                               maxfev = 20, #Rfactor strongly depends on this number
                                               stpmin = 1.e-10,
                                               stpmax = 1.e10),
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


class bulk_solvent_and_scales(object):

  def __init__(self, fmodel, params = None, log = None):
    if(params is None): params = solvent_and_scale_params()
    else:               params = solvent_and_scale_params(params = params)
    to_do = [params.bulk_solvent,
             params.anisotropic_scaling,
             params.statistical_solvent]
    if(to_do.count(False) != 3):
       params_target = params.target
       fmodel_target = fmodel.target_name
       if(fmodel.alpha_beta_params is not None):
          save_interpolation_flag = fmodel.alpha_beta_params.interpolation
          fmodel.alpha_beta_params.interpolation = False
       if(params_target != "ls_wunit_k1"): params.target = "ls_wunit_k1"
       if(to_do == [False,False,True]): params.target = params_target
       fmodel.update(target_name = params.target)
       m = "macro_cycle= "
       minimization_macro_cycles = \
                         range(1, params.number_of_minimization_macro_cycles+1)
       if(params.bulk_solvent):
          assert abs(flex.max(flex.abs(fmodel.f_mask().data()))) > 1.e-3
       macro_cycles = range(1, params.number_of_macro_cycles+1)
       if(params.k_sol_b_sol_grid_search):
          k_sols = kb_range(params.k_sol_grid_search_max,
                            params.k_sol_grid_search_min,
                            params.k_sol_step)
          b_sols = kb_range(params.b_sol_grid_search_max,
                            params.b_sol_grid_search_min,
                            params.b_sol_step)
       if(params.verbose > 0):
          fmodel.show_k_sol_b_sol_b_cart_target(header = m + str(0)+\
                             " (start) target= "+fmodel.target_name, out = log)
       if(params.fix_k_sol is not None):
          fmodel.update(k_sol = params.fix_k_sol, b_sol = params.fix_b_sol)
       if(params.fix_b_cart is not None):
          fmodel.update(b_cart = params.fix_b_cart)
       if(to_do.count(False) == 2): macro_cycles = range(1,2)
       target = fmodel.target_w()
       ksol   = fmodel.k_sol()
       bsol   = fmodel.b_sol()
       grid_search_done = False
       for mc in macro_cycles:
           do_grid_search = (ksol <= params.k_sol_grid_search_min or
                             ksol >= params.k_sol_grid_search_max or
                             bsol <= params.b_sol_grid_search_min or
                             ksol >= params.b_sol_grid_search_max)
           if(params.k_sol_b_sol_grid_search and do_grid_search and
                                                         not grid_search_done):
              grid_search_done = True
              for ksol_ in k_sols:
                  for bsol_ in b_sols:
                      fmodel.update(k_sol = ksol_, b_sol = bsol_)
                      target_ = fmodel.target_w()
                      if(target_ < target):
                         target = target_
                         ksol = ksol_
                         bsol = bsol_
              fmodel.update(k_sol = ksol, b_sol = bsol)
              if(params.verbose > 0):
                 h=m+str(mc)+": k & b: grid search; T= "+fmodel.target_name
                 fmodel.show_k_sol_b_sol_b_cart_target(header = h, out = log)
           if(params.minimization_k_sol_b_sol):
              self._k_sol_b_sol_minimization_helper(params, fmodel)
              if(params.verbose > 0):
                 h=m+str(mc)+": k & b: minimization; T= "+fmodel.target_name
                 fmodel.show_k_sol_b_sol_b_cart_target(header = h, out = log)
           if(params.minimization_b_cart):
              self._b_cart_minimizer_helper(params, fmodel)
              if(params.verbose > 0):
                 h=m+str(mc)+": anisotropic scale; T= "+fmodel.target_name
                 fmodel.show_k_sol_b_sol_b_cart_target(header = h, out = log)
           if(params.statistical_solvent):
              pass
              #self._set_f_ordered_solvent(params = params)
              #target = self.target_w()
              #if(target > target_start):
              #   print "ordered solvent: T start=, end= ",target_start,target
              #target_start = target
              #if(params.verbose > 0):
              #   h=m+str(mc)+": (ordered solvent) T= "+self.target_name
              #   self.show_k_sol_b_sol_b_cart_target(header = h, out = log)
       ### start ml optimization
       if(abs(fmodel.k_sol()) < 0.01 or abs(fmodel.b_sol()) < 1.0):
          fmodel.update(k_sol = 0, b_sol = 0)

       ksol = fmodel.k_sol()
       bsol = fmodel.b_sol()
       bcart= fmodel.b_cart()
       r_work = fmodel.r_work()
       if(params_target == "ml"):
          params.target = params_target
          fmodel.update(target_name = params_target)
          if(params.minimization_k_sol_b_sol):
             for mc in minimization_macro_cycles:
                 self._k_sol_b_sol_minimization_helper(params, fmodel)
                 if(params.verbose > 0):
                    h=m+str(mc)+": k_sol & b_sol min.; T= "+fmodel.target_name
                    fmodel.show_k_sol_b_sol_b_cart_target(header=h, out = log)
          if(fmodel.alpha_beta_params is not None):
             fmodel.alpha_beta_params.interpolation = save_interpolation_flag
       if(params.apply_back_trace_of_b_cart and abs(fmodel.u_iso()) > 0.0):
          fmodel.apply_back_b_iso()
          if(params.verbose > 0):
             h=m+str(mc)+": apply back trace(b_cart): T= "+fmodel.target_name
             fmodel.show_k_sol_b_sol_b_cart_target(header = h, out = log)
       fmodel.update(target_name = fmodel_target)
       if(abs(fmodel.k_sol()) < 0.01 or abs(fmodel.b_sol()) < 1.0):
          fmodel.update(k_sol = 0, b_sol = 0)

  def _set_f_ordered_solvent(self):
    pass

  def _b_cart_minimizer_helper(self, params, fmodel):
    symm_constr = params.symmetry_constraints_on_b_cart
    u_cycles = range(1, params.number_of_cycles_for_anisotropic_scaling+1)
    r_start = fmodel.r_work()
    u_start = fmodel.b_cart()
    for u_cycle in u_cycles:
        b_cart = aniso_scale_minimizer(fmodel      = fmodel,
                                        symm_constr = symm_constr)
        fmodel.update(b_cart = b_cart)
    r_final = fmodel.r_work()
    if(r_final - r_start > 0.01):
       # XXX reproduce and fix
       #print "Warning: r went up after anisotropic scaling:"
       #print "   r_start = ", r_start
       #print "   r_final = ", r_final
       #print "   u_start = ", ["%9.3f"%u for u in u_start]
       #print "   u_final = ", ["%9.3f"%u for u in fmodel.b_cart()]
       fmodel.update(b_cart = u_start)

  def _k_sol_b_sol_minimization_helper(self, params, fmodel):
    ksol_orig, bsol_orig = fmodel.k_sol_b_sol()
    r_start = fmodel.r_work()
    ksol, bsol = k_sol_b_sol_minimizer(fmodel = fmodel)
    if(ksol <= params.k_sol_min or ksol >= params.k_sol_max):
       k1 = abs(abs(ksol) - abs(params.k_sol_min))
       k2 = abs(abs(ksol) - abs(params.k_sol_max))
       if(k1 >= k2): ksol = params.k_sol_max
       if(k1 <= k2): ksol = params.k_sol_min
    if(bsol <= params.b_sol_min or bsol >= params.b_sol_max):
       b1 = abs(abs(bsol) - abs(params.b_sol_min))
       b2 = abs(abs(bsol) - abs(params.b_sol_max))
       if(b1 >= b2): bsol = params.b_sol_max
       if(b1 <= b2): bsol = params.b_sol_min
    fmodel.update(k_sol = ksol, b_sol = bsol)
    r_end = fmodel.r_work()
    if(r_end - r_start > 0.01):
       fmodel.update(k_sol = ksol_orig, b_sol = bsol_orig)

def kb_range(x_max, x_min, step):
  x_range = []
  x = x_min
  while x <= x_max + 0.0001:
    x_range.append(x)
    x += step
  return x_range
