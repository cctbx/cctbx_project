from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx import lbfgs
import copy, math
from cctbx import adptbx


class manager(object):
  def __init__(self, fmodel,
                     tan_b_iso_max            = 0,
                     selections               = None,
                     u_initial                = None,
                     max_number_of_iterations = 50,
                     number_of_macro_cycles   = 5,
                     convergence_test         = True,
                     convergence_delta        = 0.00001,
                     log                      = None):
    if(log is None): log = sys.stdout
    if(selections is None):
       selections = []
       selections.append(flex.bool(fmodel.xray_structure.scatterers().size(),
                                                                         True))
    else: assert len(selections) > 0
    if(u_initial is None):
       u_initial = []
       eps = math.pi**2*8
       for sel in selections:
           u = fmodel.xray_structure.select(sel).extract_u_iso_or_u_equiv()
           u_initial.append(flex.mean(u))
    fmodel_copy = fmodel.deep_copy()
    rworks = flex.double()
    minimized = None
    for macro_cycle in xrange(1,number_of_macro_cycles+1,1):
        if(minimized is not None):
           u = flex.double(minimized.u_min)
           un_sel = u <= 0.0
           up_sel = u >  0.0
           up = u.select(up_sel)
           u.set_selected(un_sel, flex.mean(up))
           assert (u <= 0).count(True) == 0
           u_initial = list(u)
        minimized = group_u_iso_minimizer(
                           fmodel                   = fmodel_copy,
                           selections               = selections,
                           u_initial                = u_initial,
                           max_number_of_iterations = max_number_of_iterations,
                           tan_b_iso_max            = tan_b_iso_max)
        new_xrs = apply_transformation(
                              xray_structure = minimized.fmodel.xray_structure,
                              u              = u_initial,
                              selections     = selections)
        fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                          update_f_calc  = True,
                                          out            = log)
        rwork = minimized.fmodel.r_work()
        rfree = minimized.fmodel.r_free()
        self.show(f     = fmodel_copy.f_obs_w(),
                  rw    = rwork,
                  rf    = rfree,
                  tw    = minimized.fmodel.target_w(),
                  mc    = macro_cycle,
                  it    = minimized.counter,
                  ct    = convergence_test,
                  out   = log)
        if(convergence_test):
           rworks.append(rwork)
           if(rworks.size() > 1):
              size = rworks.size() - 1
              if(abs(rworks[size]-rworks[size-1])<convergence_delta):
                 break
    fmodel.update_xray_structure(xray_structure = fmodel_copy.xray_structure,
                                 update_f_calc  = True)
    self.fmodel = fmodel

  def rotation(self):
    return self.total_rotation

  def translation(self):
    return self.total_translation

  def show(self, f,
                 rw,
                 rf,
                 tw,
                 mc,
                 it,
                 ct,
                 out = None):
    if(out is None): out = sys.stdout
    d_max, d_min = f.d_max_min()
    nref = f.data().size()
    mc = str(mc)
    it = str(it)
    part1 = "|-group b-factor refinement (macro cycle = "
    part2 = "; iterations = "
    n = 77 - len(part1 + part2 + mc + it)
    part3 = ")"+"-"*n+"|"
    print >> out, part1 + mc + part2 + it + part3
    part1 = "| "
    if(ct): ct = "on"
    else:   ct = "off"
    part4 = " convergence test = "+str("%s"%ct)
    rw = "| r_work = "+str("%.6f"%rw)
    rf = " r_free = "+str("%.6f"%rf)
    tw = " target = "+str("%.6f"%tw)
    n = 78 - len(rw+rf+tw+part4)
    end = " "*n+"|"
    print >> out, rw+rf+tw+part4+end
    print >> out, "|" +"-"*77+"|"

class group_u_iso_minimizer(object):
  def __init__(self,
               fmodel,
               selections,
               u_initial,
               max_number_of_iterations,
               tan_b_iso_max):
    adopt_init_args(self, locals())
    if(self.fmodel.target_name in ["ml","lsm"]):
       # XXX looks like alpha & beta must be recalcuclated in line search
       #self.alpha, self.beta = self.fmodel.alpha_beta_w()
       self.alpha, self.beta = None, None
    else:
       self.alpha, self.beta = None, None
    self.fmodel_copy = self.fmodel.deep_copy()
    self.counter=0
    assert len(self.selections) == len(self.u_initial)
    self.u_min = copy.deepcopy(self.u_initial)
    self.x = self.pack(self.u_min)
    self.n = self.x.size()
    self.minimizer = lbfgs.run(
               target_evaluator = self,
               termination_params = lbfgs.termination_parameters(
                    max_iterations = max_number_of_iterations),
               exception_handling_params = lbfgs.exception_handling_parameters(
                    ignore_line_search_failed_step_at_lower_bound = True)
                              )
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, u):
    if(self.tan_b_iso_max > 0):
       for i, ui in enumerate(u):
           u[i] =math.tan(math.pi*(ui/adptbx.b_as_u(self.tan_b_iso_max)-1./2.))
    return flex.double(tuple(u))

  def unpack_x(self):
    if(self.tan_b_iso_max > 0):
       for i, ui in enumerate(self.x):
           self.u_min[i] = adptbx.b_as_u(self.tan_b_iso_max)*(math.atan(ui)+\
                                                            math.pi/2.)/math.pi
    else:
       self.u_min = tuple(self.x)

  def compute_functional_and_gradients(self):
    self.unpack_x()
    self.counter += 1
    new_xrs = apply_transformation(xray_structure = self.fmodel.xray_structure,
                                   u              = self.u_min,
                                   selections     = self.selections)
    self.fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                           update_f_calc  = True)
    tg_obj = target_and_grads(fmodel        = self.fmodel_copy,
                              alpha         = self.alpha,
                              beta          = self.beta,
                              selections    = self.selections,
                              tan_b_iso_max = 0)
    self.f = tg_obj.target()
    ##########################################################################
    #eps = 0.001
    #new_xrs = apply_transformation(xray_structure = self.fmodel.xray_structure,
    #                               u              = [self.u_min[0]+eps, self.u_min[1],self.u_min[2]],
    #                               selections     = self.selections)
    #self.fmodel_copy.update_xray_structure(xray_structure = new_xrs,
    #                                       update_f_calc  = True)
    #tg_obj = target_and_grads(fmodel        = self.fmodel_copy,
    #                          alpha         = self.alpha,
    #                          beta          = self.beta,
    #                          selections    = self.selections,
    #                          tan_b_iso_max = self.tan_b_iso_max)
    #t1 = tg_obj.target()
    ##################
    #new_xrs = apply_transformation(xray_structure = self.fmodel.xray_structure,
    #                               u              = [self.u_min[0]-eps, self.u_min[1],self.u_min[2]],
    #                               selections     = self.selections)
    #self.fmodel_copy.update_xray_structure(xray_structure = new_xrs,
    #                                       update_f_calc  = True)
    #tg_obj = target_and_grads(fmodel        = self.fmodel_copy,
    #                          alpha         = self.alpha,
    #                          beta          = self.beta,
    #                          selections    = self.selections,
    #                          tan_b_iso_max = self.tan_b_iso_max)
    #t2 = tg_obj.target()
    ##########################################################################
    grads = tg_obj.gradients_wrt_u()

    #print grads[0], (t1-t2)/(eps*2)
    self.g = flex.double(tuple(grads))
    return self.f, self.g

def apply_transformation(xray_structure,
                         u,
                         selections):
  assert len(selections) == len(u)
  new_us = flex.double()
  for sel, ui in zip(selections, u):
      xrs = xray_structure.select(sel)
      new_us.extend( flex.double(xrs.scatterers().size(), ui) )
  return xray_structure.set_u_iso(values = new_us)

class target_and_grads(object):
  def __init__(self, fmodel,
                     alpha,
                     beta,
                     selections,
                     tan_b_iso_max):
    self.grads_wrt_u = []
    target_grads_wrt_adp = fmodel.gradient_wrt_atomic_parameters(
                                                 sites         = False,
                                                 u_iso         = True,
                                                 alpha         = alpha,
                                                 beta          = beta,
                                                 tan_b_iso_max = tan_b_iso_max)
    target_grads_wrt_adp = target_grads_wrt_adp.packed()
    self.f = fmodel.target_w(alpha = alpha, beta = beta)
    for sel in selections:
        target_grads_wrt_adp_sel = target_grads_wrt_adp.select(sel)
        self.grads_wrt_u.append(flex.sum(target_grads_wrt_adp_sel))

  def target(self):
    return self.f

  def gradients_wrt_u(self):
    return self.grads_wrt_u
