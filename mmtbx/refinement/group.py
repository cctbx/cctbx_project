from cctbx.array_family import flex
from libtbx import adopt_init_args
import sys
from libtbx.test_utils import approx_equal
from scitbx import lbfgs
import copy
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
                     log                         = None,
                     occupancy_max               = None,
                     occupancy_min               = None):
    global time_group_py
    timer = user_plus_sys_time()
    self.show(rw         = fmodel.r_work(),
              rf         = fmodel.r_free(),
              tw         = fmodel.target_w(),
              mc         = 0,
              it         = 0,
              ct         = convergence_test,
              refine_adp = refine_adp,
              refine_occ = refine_occ,
              out        = log)
    if(log is None): log = sys.stdout
    assert [refine_adp, refine_occ].count(True) == 1
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
    scatterers = fmodel.xray_structure.scatterers()
    scatterers.flags_set_grads(state=False)
    # XXX very inefficient: same code is in driver.py file. fix asap. Pavel.
    save_use_u_iso = fmodel.xray_structure.use_u_iso()
    save_use_u_aniso = fmodel.xray_structure.use_u_aniso()
    for sel in selections:
      if(refine_adp):
         for s in sel:
           sc = scatterers[s]
           if(not sc.flags.use_u_iso()):
             sc.flags.set_use_u_iso(True)
             if(sc.u_iso == -1): sc.u_iso = 0
         scatterers.flags_set_grad_u_iso(iselection = sel)
      if(refine_occ):
         scatterers.flags_set_grad_occupancy(iselection = sel)
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
        apply_transformation(
          xray_structure = fmodel.xray_structure,
          par            = par_initial,
          sc_start       = sc_start,
          selections     = selections,
          refine_adp     = refine_adp,
          refine_occ     = refine_occ)
        fmodel_copy.update_xray_structure(
          xray_structure = fmodel.xray_structure,
          update_f_calc  = True,
          out            = log)
        rwork = minimized.fmodel.r_work()
        rfree = minimized.fmodel.r_free()
        assert approx_equal(rwork, fmodel_copy.r_work())
        self.show(rw    = rwork,
                  rf    = rfree,
                  tw    = minimized.fmodel.target_w(),
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
    # XXX absence of tidy_us will lead to crash in very rare cases; fixing requires too
    # XXX much of time and brain investment.
    #fmodel_copy.xray_structure.tidy_us()
    fmodel.update_xray_structure(xray_structure = fmodel_copy.xray_structure,
                                 update_f_calc  = True)
    if(refine_occ):
      i_selection = flex.size_t()
      for sel in selections:
        i_selection.extend(sel)
      fmodel.xray_structure.adjust_occupancy(occ_max   = occupancy_max,
                                             occ_min   = occupancy_min,
                                             selection = i_selection)
    self.fmodel = fmodel
    time_group_py += timer.elapsed()

  def show(self, rw,
                 rf,
                 tw,
                 mc,
                 it,
                 ct,
                 refine_adp,
                 refine_occ,
                 out = None):
    if(out is None): out = sys.stdout
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
    rw = "| r_work = "+str("%.4f"%rw)
    rf = " r_free = "+str("%.4f"%rf)
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
    self.target_functor = fmodel.target_functor()
    self.target_functor.prepare_for_minimization()
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
                    ignore_line_search_failed_step_at_lower_bound = True,
                    ignore_line_search_failed_step_at_upper_bound = True,)
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
    apply_transformation(
      xray_structure = self.fmodel.xray_structure,
      par            = self.par_min,
      sc_start       = self.sc_start,
      selections     = self.selections,
      refine_adp     = self.refine_adp,
      refine_occ     = self.refine_occ)
    self.fmodel.update_xray_structure(update_f_calc=True)
    tg_obj = target_and_grads(
      target_functor=self.target_functor,
      selections=self.selections,
      refine_adp=self.refine_adp,
      refine_occ=self.refine_occ)
    self.f = tg_obj.target()
    self.g = flex.double(tg_obj.gradients_wrt_par())
    compare = False
    if (self.run_finite_differences_test):
      for pari in self.par_min:
        if (self.refine_adp):
          pari = adptbx.u_as_b(pari)
          if (pari < 0 or pari > 100):
            break
        if (self.refine_occ and (pari < 0 or pari > 10)):
          break
      else:
        compare = True
    if (compare):
       i_g_max = flex.max_index(flex.abs(self.g))
       eps = 1.e-5
       par_eps = list(self.par_min)
       par_eps[i_g_max] = self.par_min[i_g_max] + eps
       apply_transformation(
         xray_structure = self.fmodel.xray_structure,
         par            = par_eps,
         sc_start       = self.sc_start,
         selections     = self.selections,
         refine_adp     = self.refine_adp,
         refine_occ     = self.refine_occ)
       self.fmodel.update_xray_structure(update_f_calc=True)
       t1 = target_and_grads(
         target_functor=self.target_functor,
         selections=self.selections,
         refine_adp=self.refine_adp,
         refine_occ=self.refine_occ,
         compute_gradients=False).target()
       par_eps[i_g_max] = self.par_min[i_g_max] - eps
       apply_transformation(
         xray_structure = self.fmodel.xray_structure,
         par            = par_eps,
         sc_start       = self.sc_start,
         selections     = self.selections,
         refine_adp     = self.refine_adp,
         refine_occ     = self.refine_occ)
       del par_eps
       self.fmodel.update_xray_structure(update_f_calc=True)
       t2 = target_and_grads(
         target_functor=self.target_functor,
         selections=self.selections,
         refine_adp=self.refine_adp,
         refine_occ=self.refine_occ,
         compute_gradients=False).target()
       apply_transformation(
         xray_structure = self.fmodel.xray_structure,
         par            = self.par_min,
         sc_start       = self.sc_start,
         selections     = self.selections,
         refine_adp     = self.refine_adp,
         refine_occ     = self.refine_occ)
       self.fmodel.update_xray_structure(update_f_calc=True)
       self.buffer_ana.append(self.g[i_g_max])
       self.buffer_fin.append((t1-t2)/(eps*2))
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

class target_and_grads(object):
  def __init__(self, target_functor,
                     selections,
                     refine_adp,
                     refine_occ,
                     compute_gradients=True):
    assert [refine_adp, refine_occ].count(True) == 1
    t_r = target_functor(compute_gradients=compute_gradients)
    self.f = t_r.target_work()
    if (compute_gradients):
      target_grads_wrt_par = t_r.gradients_wrt_atomic_parameters(
        u_iso=refine_adp,
        occupancy=refine_occ)
      self.grads_wrt_par = []
      for sel in selections:
        target_grads_wrt_par_sel = target_grads_wrt_par.select(sel)
        self.grads_wrt_par.append(flex.sum(target_grads_wrt_par_sel))
    else:
      self.grads_wrt_par = None

  def target(self):
    return self.f

  def gradients_wrt_par(self):
    return self.grads_wrt_par
