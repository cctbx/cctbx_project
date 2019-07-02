from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx import adopt_init_args
import sys
from libtbx.test_utils import approx_equal
from scitbx import lbfgs
import copy
from cctbx import adptbx
from cctbx import xray
from libtbx.utils import user_plus_sys_time
from cctbx import crystal
import random
from six.moves import zip
from six.moves import range

time_group_py  = 0.0

def show_times(out = None):
  if(out is None): out = sys.stdout
  total = time_group_py
  if(total > 0.01):
     print("Group ADP refinement:", file=out)
     print("  time_group_py                          = %-7.2f" % time_group_py, file=out)
  return total

class sphere_similarity_restraints(object):
  def __init__(
        self,
        xray_structure,
        selection,
        refine_adp,
        refine_occ,
        sphere_radius = 5.0):
    self.selection = selection
    assert str(type(selection).__name__) == "bool"
    self.sphere_radius = sphere_radius
    self.refine_adp = refine_adp
    self.refine_occ = refine_occ
    pair_asu_table = xray_structure.pair_asu_table(
      distance_cutoff = self.sphere_radius)
    self.pair_sym_table = pair_asu_table.extract_pair_sym_table()
    assert [self.refine_adp, self.refine_occ].count(True) == 1
    self.sites_frac = xray_structure.sites_frac()
    self.orthogonalization_matrix = \
      xray_structure.unit_cell().orthogonalization_matrix()

  def target_and_gradients(self, xray_structure, to_compute_weight=False):
    if(to_compute_weight):
      xrs = xray_structure.deep_copy_scatterers()
      # This may be useful to explore:
      #xrs.shake_adp_if_all_equal(b_iso_tolerance = 1.e-3)
      #xrs.shake_adp(spread=10, keep_anisotropic= False)
    else:
      xrs = xray_structure
    if(self.refine_adp): params = xrs.extract_u_iso_or_u_equiv()
    if(self.refine_occ): params = xrs.scatterers().extract_occupancies()
    if(to_compute_weight):
      pmin = flex.min(params)
      pmax = flex.max(params)
      if(abs(pmin+pmax)!=0. and abs(pmin-pmax)/abs(pmin+pmax)*2*100<1.e-3):
        pmean = flex.mean(params)
        n_par = params.size()
        params = flex.double()
        for i in range(n_par):
          params.append(pmean + 0.1 * pmean * random.choice([-1,0,1]))
    return crystal.adp_iso_local_sphere_restraints_energies(
      pair_sym_table           = self.pair_sym_table,
      orthogonalization_matrix = self.orthogonalization_matrix,
      sites_frac               = self.sites_frac,
      u_isos                   = params,
      selection                = self.selection,
      use_u_iso                = self.selection,
      grad_u_iso               = self.selection,
      sphere_radius            = self.sphere_radius,
      distance_power           = 2,
      average_power            = 1,
      min_u_sum                = 1.e-6,
      compute_gradients        = True,
      collect                  = False)

class manager(object):
  def __init__(
        self,
        fmodel,
        selections                  = None,
        max_number_of_iterations    = 50,
        number_of_macro_cycles      = 5,
        use_restraints              = False,
        restraints_weight           = None,
        convergence_test            = True,
        convergence_delta           = 0.00001,
        run_finite_differences_test = False,
        refine_adp                  = False,
        refine_occ                  = False,
        log                         = None,
        occupancy_max               = None,
        occupancy_min               = None):
    global time_group_py
    #
    tmp = flex.size_t()
    for s in selections: tmp.extend(s)
    selections_as_bool = flex.bool(fmodel.xray_structure.scatterers().size(),
      tmp)
    #
    if(log is None): log = sys.stdout
    timer = user_plus_sys_time()
    self.show(
      rw         = fmodel.r_work(),
      rf         = fmodel.r_free(),
      tw         = fmodel.target_w(),
      mc         = 0,
      it         = 0,
      refine_adp = refine_adp,
      refine_occ = refine_occ,
      weight     = restraints_weight,
      out        = log)
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
    restraints_manager = None
    if(use_restraints):
      restraints_manager = sphere_similarity_restraints(
        xray_structure = fmodel.xray_structure,
        selection      = selections_as_bool,
        refine_adp     = refine_adp,
        refine_occ     = refine_occ)
    fmodel_copy = fmodel.deep_copy()
    rworks = flex.double()
    sc_start = fmodel.xray_structure.scatterers().deep_copy()
    minimized = None
    self.tested = 0
    for macro_cycle in range(1,number_of_macro_cycles+1,1):
      if(minimized is not None): par_initial = minimized.par_min
      minimized = group_minimizer(
        fmodel                      = fmodel_copy,
        sc_start                    = sc_start,
        selections                  = selections,
        par_initial                 = par_initial,
        refine_adp                  = refine_adp,
        refine_occ                  = refine_occ,
        max_number_of_iterations    = max_number_of_iterations,
        run_finite_differences_test = run_finite_differences_test,
        restraints_manager          = restraints_manager,
        restraints_weight           = restraints_weight)
      if(minimized is not None):
        par_initial = minimized.par_min
        self.tested += minimized.tested
      apply_transformation(
        xray_structure = fmodel.xray_structure,
        par            = par_initial,
        sc_start       = sc_start,
        selections     = selections,
        refine_adp     = refine_adp,
        refine_occ     = refine_occ)
      fmodel_copy.update_xray_structure(
        xray_structure = fmodel.xray_structure,
        update_f_calc  = True)
      rwork = minimized.fmodel.r_work()
      rfree = minimized.fmodel.r_free()
      assert approx_equal(rwork, fmodel_copy.r_work())
      self.show(
        rw         = rwork,
        rf         = rfree,
        tw         = minimized.fmodel.target_w(),
        mc         = macro_cycle,
        it         = minimized.counter,
        refine_adp = refine_adp,
        refine_occ = refine_occ,
        weight     = minimized.weight,
        out        = log)
      if(convergence_test):
        rworks.append(rwork)
        if(rworks.size() > 1):
          size = rworks.size() - 1
          if(abs(rworks[size]-rworks[size-1])<convergence_delta):
             break
    fmodel_copy.xray_structure.tidy_us()
    fmodel.update_xray_structure(
      xray_structure = fmodel_copy.xray_structure, update_f_calc  = True)
    if(refine_occ):
      i_selection = flex.size_t()
      for sel in selections:
        i_selection.extend(sel)
      fmodel.xray_structure.adjust_occupancy(
        occ_max   = occupancy_max,
        occ_min   = occupancy_min,
        selection = i_selection)
    self.fmodel = fmodel
    time_group_py += timer.elapsed()

  def show(
        self,
        rw,
        rf,
        tw,
        mc,
        it,
        refine_adp,
        refine_occ,
        weight,
        out):
    if(out is None): return
    mc = str(mc)
    it = str(it)
    if(refine_adp): part1 = "|-group b-factor refinement (macro cycle = "
    if(refine_occ): part1 = "|-group occupancy refinement (macro cycle = "
    part2 = "; iterations = "
    n = 77 - len(part1 + part2 + mc + it)
    part3 = ")"+"-"*n+"|"
    print(part1 + mc + part2 + it + part3, file=out)
    part1 = "| "
    if(weight is None):
      part4 = " restraints weight = "+str(weight)
    else:
      part4 = " restraints weight = "+str("%10.3f"%weight).strip()
    rw = "| r_work = "+str("%.4f"%rw)
    rf = " r_free = "+str("%.4f"%rf)
    tw = " target = "+str("%.6f"%tw)
    n = 78 - len(rw+rf+tw+part4)
    end = " "*n+"|"
    print(rw+rf+tw+part4+end, file=out)
    print("|" +"-"*77+"|", file=out)

class group_minimizer(object):
  def __init__(
        self,
        fmodel,
        sc_start,
        selections,
        par_initial,
        refine_adp,
        refine_occ,
        max_number_of_iterations,
        run_finite_differences_test = False,
        restraints_weight = None,
        restraints_manager = None):
    adopt_init_args(self, locals())
    self.target_functor = fmodel.target_functor()
    self.target_functor.prepare_for_minimization()
    self.counter=0
    assert len(self.selections) == len(self.par_initial)
    self.par_min = copy.deepcopy(self.par_initial)
    self.x = self.pack(self.par_min)
    self.n = self.x.size()
    self.weight = restraints_weight
    if(self.restraints_manager is not None and self.weight is None):
      gx = self.target_functor(
        compute_gradients=True).gradients_wrt_atomic_parameters(
          u_iso     = refine_adp,
          occupancy = refine_occ)
      rtg = self.restraints_manager.target_and_gradients(
        xray_structure    = self.fmodel.xray_structure,
        to_compute_weight = True)
      gx_norm = gx.norm()
      if(gx_norm != 0):
        self.weight = rtg.gradients.norm()/gx_norm
      else: self.weight = 1.0
    if(self.weight is not None):
      assert self.restraints_manager is not None
    if(run_finite_differences_test):
      self.buffer_ana = flex.double()
      self.buffer_fin = flex.double()
    self.minimizer = lbfgs.run(
      target_evaluator = self,
      termination_params = lbfgs.termination_parameters(
        max_iterations = max_number_of_iterations),
      exception_handling_params = lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = True,
        ignore_line_search_failed_step_at_upper_bound = True))
    self.compute_functional_and_gradients()
    del self.x
    self.tested = 0
    if(run_finite_differences_test):
      # For debugging only
      #for a,f in zip(self.buffer_ana, self.buffer_fin):
      #  print a, f
      diff = flex.abs(self.buffer_ana - self.buffer_fin)
      s = diff < 1.e-3
      if(s.size()>0 and s.count(True)*100./s.size()>50):
        self.tested += 1

  def get_restraints_tg(self):
    result = None
    if(self.restraints_manager is not None):
      result = self.restraints_manager.target_and_gradients(
        xray_structure = self.fmodel.xray_structure)
    return result

  def get_tg(self, compute_gradients):
    return target_and_grads(
      target_functor    = self.target_functor,
      selections        = self.selections,
      refine_adp        = self.refine_adp,
      refine_occ        = self.refine_occ,
      compute_gradients = compute_gradients,
      rtg               = self.get_restraints_tg(),
      weight            = self.weight)

  def apply_shifts(self, par):
    apply_transformation(
      xray_structure = self.fmodel.xray_structure,
      par            = par,
      sc_start       = self.sc_start,
      selections     = self.selections,
      refine_adp     = self.refine_adp,
      refine_occ     = self.refine_occ)

  def pack(self, par):
    return flex.double(tuple(par))

  def unpack_x(self):
    self.par_min = tuple(self.x)

  def finite_difference_test(self):
    if(self.fmodel.r_work()>1.e-3):
      i_g_max = flex.max_index(flex.abs(self.g))
      eps = 1.e-5
      par_eps = list(self.par_min)
      par_eps[i_g_max] = self.par_min[i_g_max] + eps
      self.apply_shifts(par = par_eps)
      self.fmodel.update_xray_structure(update_f_calc=True)
      t1 = self.get_tg(compute_gradients=False).target()
      par_eps[i_g_max] = self.par_min[i_g_max] - eps
      self.apply_shifts(par = par_eps)
      del par_eps
      self.fmodel.update_xray_structure(update_f_calc=True)
      t2 = self.get_tg(compute_gradients=False).target()
      self.apply_shifts(par = self.par_min)
      self.fmodel.update_xray_structure(update_f_calc=True)
      self.buffer_ana.append(self.g[i_g_max])
      self.buffer_fin.append((t1-t2)/(eps*2))

  def compute_functional_and_gradients(self):
    self.unpack_x()
    self.counter += 1
    self.apply_shifts(par = self.par_min)
    self.fmodel.update_xray_structure(update_f_calc=True)
    tg_obj = self.get_tg(compute_gradients=True)
    self.f = tg_obj.target()
    self.g = flex.double(tg_obj.gradients_wrt_par())
    if(self.run_finite_differences_test):
      self.finite_difference_test()
    return self.f, self.g

def apply_transformation(
      xray_structure,
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
      xray.shift_us(
        scatterers = new_sc,
        unit_cell  = xray_structure.unit_cell(),
        u_shift    = pari,
        selection  = sel)
    if(refine_occ):
      xray.shift_occupancies(
        scatterers = new_sc,
        q_shift    = pari,
        selection  = sel)
  xray_structure.replace_scatterers(new_sc)

class target_and_grads(object):
  def __init__(
        self,
        target_functor,
        selections,
        refine_adp,
        refine_occ,
        compute_gradients=True,
        rtg=None,
        weight=None):
    assert [refine_adp, refine_occ].count(True) == 1
    t_r = target_functor(compute_gradients=compute_gradients)
    self.f = t_r.target_work()
    if(rtg is not None): self.f = self.f*weight+rtg.residual_sum
    if(compute_gradients):
      target_grads_wrt_par = t_r.gradients_wrt_atomic_parameters(
        u_iso     = refine_adp,
        occupancy = refine_occ)
      if(rtg is not None):
        target_grads_wrt_par = target_grads_wrt_par*weight+rtg.gradients
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
