from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx import lbfgs
import copy
from cctbx import adptbx


class manager(object):
  def __init__(self, fmodel,
                     tan_b_iso_max,
                     selections        = None,
                     u_initial         = None,
                     max_iterations    = 20,
                     convergence_test  = True,
                     convergence_delta = 0.00001,
                     log               = None):
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
           u_initial.append(flex.mean(u)/2)
       print >> log, "Start B-values for group B-factor refinement:"
       if(len(u_initial) < 10):
          print >> log, [float("%10.3f"%adptbx.u_as_b(u)) for u in u_initial]
       else:
          print >> log, "first 10: ", \
                     [float("%10.3f"%adptbx.u_as_b(u)) for u in u_initial[:10]]
    fmodel_copy = fmodel.deep_copy()
    print >> log
    rworks = flex.double()
    minimized = None
    for macro_cycle in xrange(5):
        if(minimized is not None):
           u = flex.double(minimized.u_min)
           un_sel = u <= 0.0
           up_sel = u >  0.0
           up = u.select(up_sel)
           u.set_selected(un_sel, flex.mean(up))
           assert (u <= 0).count(True) == 0
           u_initial = list(u)
        minimized = group_u_iso_minimizer(fmodel         = fmodel_copy,
                                          selections     = selections,
                                          u_initial      = u_initial,
                                          max_iterations = max_iterations,
                                          tan_b_iso_max  = tan_b_iso_max)
        new_xrs = apply_transformation(
                              xray_structure = minimized.fmodel.xray_structure,
                              u              = u_initial,
                              selections     = selections)
        fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                          update_f_calc  = True,
                                          out            = log)
        rwork = minimized.fmodel.r_work()
        rfree = minimized.fmodel.r_free()
        assert approx_equal(rwork, fmodel_copy.r_work())
        print "*"*30
        print "rwork, rfree = ",rwork, rfree
        print adptbx.u_as_b(minimized.u_min[0]), adptbx.u_as_b(minimized.u_min[1]), adptbx.u_as_b(minimized.u_min[2])
        print "*"*30
        #self.show(f     = fmodel_copy.f_obs_w(),
        #          r_mat = self.total_rotation,
        #          t_vec = self.total_translation,
        #          rw    = rwork,
        #          rf    = rfree,
        #          tw    = minimized.fmodel.target_w(),
        #          mc    = macro_cycle,
        #          it    = minimized.counter,
        #          ct    = convergence_test,
        #          out   = log)
        #if(convergence_test):
        #   rworks.append(rwork)
        #   if(rworks.size() > 1):
        #      size = rworks.size() - 1
        #      if(abs(rworks[size]-rworks[size-1])<convergence_delta):
        #         break
    fmodel.update_xray_structure(xray_structure = fmodel_copy.xray_structure,
                                 update_f_calc  = True)
    self.fmodel = fmodel

  def rotation(self):
    return self.total_rotation

  def translation(self):
    return self.total_translation

  def show(self, f,
                 r_mat,
                 t_vec,
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
    part1 = "|-rigid body refinement (macro cycle = "
    part2 = "; iterations = "
    n = 77 - len(part1 + part2 + mc + it)
    part3 = ")"+"-"*n+"|"
    print >> out, part1 + mc + part2 + it + part3
    part1 = "| resolution range: "
    d_max = str("%.3f"%d_max)
    part2 = " - "
    d_min = str("%.3f"%d_min)
    part3 = " ("
    nref = str("%d"%nref)
    if(ct): ct = "on"
    else:   ct = "off"
    part4 = " reflections) convergence test = "+str("%s"%ct)
    n = 78 - len(part1+d_max+part2+d_min+part3+nref+part4)
    part5 = " "*n+"|"
    print >> out, part1+d_max+part2+d_min+part3+nref+part4+part5
    rw = "| r_work = "+str("%.6f"%rw)
    rf = " r_free = "+str("%.6f"%rf)
    tw = " target = "+str("%.6f"%tw)
    n = 78 - len(rw+rf+tw)
    end = " "*n+"|"
    print >> out, rw+rf+tw+end
    print >> out, "|                         rotation (deg.)             "\
                  "   translation (A)      |"
    i = 1
    for r,t in zip(r_mat,t_vec):
        part1 = "| group"+str("%5d:  "%i)
        part2 = str("%8.4f"%r[0])+" "+str("%8.4f"%r[1])+" "+str("%8.4f"%r[2])
        part3 = "     "
        part4 = str("%8.4f"%t[0])+" "+str("%8.4f"%t[1])+" "+str("%8.4f"%t[2])
        n = 78 - len(part1 + part2 + part3 + part4)
        part5 = " "*n+"|"
        print >> out, part1 + part2 + part3 + part4 + part5
        i += 1
    print >> out, "|" +"-"*77+"|"

class group_u_iso_minimizer(object):
  def __init__(self,
               fmodel,
               selections,
               u_initial,
               max_iterations,
               tan_b_iso_max):
    adopt_init_args(self, locals())
    # XXX ???
    if(self.fmodel.target_name in ["ml","lsm"]):
       self.alpha, self.beta = None, None#self.fmodel.alpha_beta_w()
    else:
       self.alpha, self.beta = None, None
    self.fmodel_copy = self.fmodel.deep_copy()
    self.n_groups = len(self.selections)
    assert self.n_groups > 0
    self.counter=0
    assert len(self.selections) == len(self.u_initial)
    self.u_min = copy.deepcopy(self.u_initial)
    print self.u_min
    for i in xrange(len(self.u_min)):
        self.u_min[i] = self.u_min[i]
    self.x = self.pack(self.u_min)
    self.n = self.x.size()
    self.minimizer = lbfgs.run(
               target_evaluator = self,
               termination_params = lbfgs.termination_parameters(
                    max_iterations = max_iterations),
               exception_handling_params = lbfgs.exception_handling_parameters(
                    ignore_line_search_failed_step_at_lower_bound = True)
                              )
    self.compute_functional_and_gradients()
    del self.x

  def pack(self, u):
    return flex.double(tuple(u))

  def unpack_x(self):
    self.u_min = tuple(self.x)

  def compute_functional_and_gradients(self):
    self.unpack_x()
    self.counter += 1
    new_xrs = apply_transformation(xray_structure = self.fmodel.xray_structure,
                                   u              = self.u_min,
                                   selections     = self.selections)
    self.fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                           update_f_calc  = True)
    tg_obj = target_and_grads(fmodel     = self.fmodel_copy,
                              alpha      = self.alpha,
                              beta       = self.beta,
                              selections = self.selections,
                              tan_b_iso_max = self.tan_b_iso_max)
    self.f = tg_obj.target()
    self.g = self.pack( tg_obj.gradients_wrt_u() )
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
        self.grads_wrt_u.append(flex.sum(target_grads_wrt_adp.select(sel)))

  def target(self):
    return self.f

  def gradients_wrt_u(self):
    return self.grads_wrt_u
