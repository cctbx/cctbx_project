from cctbx.array_family import flex
from libtbx import adopt_init_args
import math, sys
from libtbx.test_utils import approx_equal
from scitbx import matrix
from scitbx import lbfgs
import copy, math
from cctbx import adptbx
from cctbx import xray
from libtbx.utils import user_plus_sys_time

time_group_py  = 0.0

def show_times(out = None):
  if(out is None): out = sys.stdout
  total = time_group_py
  if(total > 0.01):
     print >> out, "Group ADP refinement:"
     print >> out, "  time_group_py                          = %-7.2f" % time_group_py
  return total


class manager(object):
  def __init__(self, fmodel,
                     selections                  = None,
                     max_number_of_iterations    = 50,
                     number_of_macro_cycles      = 5,
                     convergence_test            = True,
                     convergence_delta           = 0.00001,
                     run_finite_differences_test = False,
                     refine_adp                  = False,
                     refine_occ                  = False,
                     log                         = None):
    global time_group_py
    timer = user_plus_sys_time()
    if(log is None): log = sys.stdout
    assert [refine_adp, refine_occ].count(True) == 1
    xray.set_scatterer_grad_flags(
                               scatterers = fmodel.xray_structure.scatterers(),
                               u_iso      = refine_adp,
                               occupancy  = refine_occ)
    if(selections is None):
       selections = []
       selections.append(
                    flex.bool(fmodel.xray_structure.scatterers().size(), True))
    else: assert len(selections) > 0
    par_initial = []
    selections_ = []
    for sel in selections:
        if(refine_adp): par_initial.append(adptbx.b_as_u(0.0))
        if(refine_occ): par_initial.append(0.0)
        if(str(type(sel).__name__) == "bool"):
           selections_.append(sel.iselection())
        else:
           selections_.append(sel)
    selections = selections_
    fmodel_copy = fmodel.deep_copy()
    rworks = flex.double()
    sc_start = fmodel.xray_structure.scatterers().deep_copy()
    minimized = None
    for macro_cycle in xrange(1,number_of_macro_cycles+1,1):
        if(minimized is not None): par_initial = minimized.par_min
        minimized = group_minimizer(
                     fmodel                      = fmodel_copy,
                     sc_start                    = sc_start,
                     selections                  = selections,
                     par_initial                 = par_initial,
                     refine_adp                  = refine_adp,
                     refine_occ                  = refine_occ,
                     max_number_of_iterations    = max_number_of_iterations,
                     run_finite_differences_test = run_finite_differences_test)
        if(minimized is not None): par_initial = minimized.par_min
        new_xrs = apply_transformation(xray_structure = fmodel.xray_structure,
                                       par            = par_initial,
                                       sc_start       = sc_start,
                                       selections     = selections,
                                       refine_adp     = refine_adp,
                                       refine_occ     = refine_occ)
        fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                          update_f_calc  = True,
                                          out            = log)
        rwork = minimized.fmodel_copy.r_work()
        rfree = minimized.fmodel_copy.r_free()
        assert approx_equal(rwork, fmodel_copy.r_work())
        self.show(f     = fmodel_copy.f_obs_w,
                  rw    = rwork,
                  rf    = rfree,
                  tw    = minimized.fmodel_copy.target_w(),
                  mc    = macro_cycle,
                  it    = minimized.counter,
                  ct    = convergence_test,
                  refine_adp = refine_adp,
                  refine_occ = refine_occ,
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
    time_group_py += timer.elapsed()

  def show(self, f,
                 rw,
                 rf,
                 tw,
                 mc,
                 it,
                 ct,
                 refine_adp,
                 refine_occ,
                 out = None):
    if(out is None): out = sys.stdout
    d_max, d_min = f.d_max_min()
    nref = f.data().size()
    mc = str(mc)
    it = str(it)
    if(refine_adp): part1 = "|-group b-factor refinement (macro cycle = "
    if(refine_occ): part1 = "|-group occupancy refinement (macro cycle = "
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

class group_minimizer(object):
  def __init__(self,
               fmodel,
               sc_start,
               selections,
               par_initial,
               refine_adp,
               refine_occ,
               max_number_of_iterations,
               run_finite_differences_test = False):
    adopt_init_args(self, locals())
    if(self.fmodel.target_name in ["ml","lsm", "mlhl"]):
       self.alpha, self.beta = self.fmodel.alpha_beta_w()
    else:
       self.alpha, self.beta = None, None
    self.fmodel_copy = self.fmodel.deep_copy()
    self.counter=0
    assert len(self.selections) == len(self.par_initial)
    self.par_min = copy.deepcopy(self.par_initial)
    self.x = self.pack(self.par_min)
    self.n = self.x.size()
    if (run_finite_differences_test):
      self.buffer_ana = []
      self.buffer_fin = []
    self.minimizer = lbfgs.run(
               target_evaluator = self,
               termination_params = lbfgs.termination_parameters(
                    max_iterations = max_number_of_iterations),
               exception_handling_params = lbfgs.exception_handling_parameters(
                    ignore_line_search_failed_step_at_lower_bound = True)
                              )
    self.compute_functional_and_gradients()
    del self.x
    if (run_finite_differences_test):
      print "analytical gradients:", self.buffer_ana
      print "finite differences:  ", self.buffer_fin
      if (len(self.buffer_ana) >= 3):
        corr = flex.linear_correlation(
          flex.double(self.buffer_ana),
          flex.double(self.buffer_fin))
        assert corr.is_well_defined()
        assert corr.coefficient() > 1-1.e-5

  def pack(self, par):
    return flex.double(tuple(par))

  def unpack_x(self):
    self.par_min = tuple(self.x)

  def compute_functional_and_gradients(self):
    self.unpack_x()
    self.counter += 1
    new_xrs = apply_transformation(xray_structure = self.fmodel.xray_structure,
                                   par            = self.par_min,
                                   sc_start       = self.sc_start,
                                   selections     = self.selections,
                                   refine_adp     = self.refine_adp,
                                   refine_occ     = self.refine_occ)
    self.fmodel_copy.update_xray_structure(xray_structure = new_xrs,
                                           update_f_calc  = True)
    tg_obj = target_and_grads(fmodel        = self.fmodel_copy,
                              alpha         = self.alpha,
                              beta          = self.beta,
                              selections    = self.selections,
                              refine_adp    = self.refine_adp,
                              refine_occ    = self.refine_occ)
    self.f = tg_obj.target()
    ##########################################################################
    if(self.run_finite_differences_test):
       eps = 1.e-5
       fmodel = self.fmodel_copy.deep_copy()
       par = []
       for i_seq, pari in enumerate(self.par_min):
         if(i_seq == 0): pari += eps
         par.append(pari)
       new_xrs = apply_transformation(xray_structure = fmodel.xray_structure,
                                      par            = par,
                                      sc_start       = self.sc_start,
                                      selections     = self.selections,
                                      refine_adp     = self.refine_adp,
                                      refine_occ     = self.refine_occ)
       fmodel.update_xray_structure(xray_structure = new_xrs,
                                    update_f_calc  = True)
       tg_obj = target_and_grads(fmodel        = fmodel,
                                 alpha         = self.alpha,
                                 beta          = self.beta,
                                 selections    = self.selections,
                                 refine_adp    = self.refine_adp,
                                 refine_occ    = self.refine_occ)
       t1 = tg_obj.target()
       #################
       par = []
       for i_seq, pari in enumerate(self.par_min):
         if(i_seq == 0): pari -= eps
         par.append(pari)
       new_xrs = apply_transformation(xray_structure = fmodel.xray_structure,
                                      par            = par,
                                      sc_start       = self.sc_start,
                                      selections     = self.selections,
                                      refine_adp     = self.refine_adp,
                                      refine_occ     = self.refine_occ)
       fmodel.update_xray_structure(xray_structure = new_xrs,
                                    update_f_calc  = True)
       tg_obj = target_and_grads(fmodel        = fmodel,
                                 alpha         = self.alpha,
                                 beta          = self.beta,
                                 selections    = self.selections,
                                 refine_adp    = self.refine_adp,
                                 refine_occ    = self.refine_occ)
       t2 = tg_obj.target()
    ##########################################################################
    grads = tg_obj.gradients_wrt_par()

    if(self.run_finite_differences_test):
       compare = True
       for pari in self.par_min:
         if(self.refine_adp): pari = abs(adptbx.u_as_b(pari))
         if(self.refine_adp and (pari < 0.5 or pari > 70.0)): compare = False
         if(self.refine_occ and (pari < 0.0 or pari > 10.0)): compare = False
       if(compare):
          self.buffer_ana.append(grads[0])
          self.buffer_fin.append((t1-t2)/(eps*2))
    self.g = flex.double(grads)
    return self.f, self.g


def apply_transformation(xray_structure,
                         par,
                         selections,
                         sc_start,
                         refine_adp,
                         refine_occ):
  assert len(selections) == len(par)
  assert [refine_adp, refine_occ].count(True) == 1
  new_sc = sc_start.deep_copy()
  for sel, pari in zip(selections, par):
      if(refine_adp):
         xray.shift_us(scatterers = new_sc,
                       unit_cell  = xray_structure.unit_cell(),
                       u_shift    = pari,
                       selection  = sel)
      if(refine_occ):
         xray.shift_occupancies(scatterers = new_sc,
                                q_shift    = pari,
                                selection  = sel)
  xray_structure.replace_scatterers(new_sc)
  return xray_structure

class target_and_grads(object):
  def __init__(self, fmodel,
                     alpha,
                     beta,
                     selections,
                     refine_adp,
                     refine_occ):
    self.grads_wrt_par = []
    assert [refine_adp, refine_occ].count(True) == 1
    target_grads_wrt_par = fmodel.gradient_wrt_atomic_parameters(
                                                    u_iso         = refine_adp,
                                                    occupancy     = refine_occ,
                                                    alpha         = alpha,
                                                    beta          = beta,
                                                    tan_b_iso_max = 0.0)
    target_grads_wrt_par = target_grads_wrt_par.packed()
    self.f = fmodel.target_w(alpha = alpha, beta = beta)
    for sel in selections:
        target_grads_wrt_par_sel = target_grads_wrt_par.select(sel)
        self.grads_wrt_par.append(flex.sum(target_grads_wrt_par_sel))

  def target(self):
    return self.f

  def gradients_wrt_par(self):
    return self.grads_wrt_par
