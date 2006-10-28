from cctbx import xray
from cctbx import crystal
from cctbx.xray.structure import structure as cctbx_xray_structure
from cctbx.array_family import flex
import scitbx.lbfgs
from libtbx import adopt_init_args
from stdlib import math
import sys, time
from libtbx.test_utils import approx_equal
from libtbx.utils import user_plus_sys_time

time_site_individual = 0.0

class lbfgs(object):

  def __init__(self, restraints_manager,
                     fmodel,
                     model,
                     wx,
                     wc                       = None,
                     wu                       = None,
                     wilson_b                 = None,
                     tan_b_iso_max            = None,
                     refine_xyz               = False,
                     refine_adp               = False,
                     refine_occ               = False,
                     lbfgs_termination_params = None,
                     f_a_to_be_r_i            = None,
                     use_fortran              = False,
                     verbose                  = 0,
                     iso_restraints           = None,
                     alpha_w                  = None,
                     beta_w                   = None,
                     fmodel_neutron           = None,
                     wn                       = None,
                     neutron_scattering_dict  = None,
                     xray_scattering_dict     = None,
                     wxnc_scale               = None,
                     wxnu_scale               = None,
                     selection                = None):
    global time_site_individual
    timer = user_plus_sys_time()
    adopt_init_args(self, locals())
    assert [refine_xyz, refine_adp, refine_occ].count(False) == 2
    assert [refine_xyz, refine_adp, refine_occ].count(True)  == 1
    self.xray_structure = self.fmodel.xray_structure
    if(refine_xyz): self.wr = wc
    if(refine_adp): self.wr = wu
    if(refine_occ): self.wr = None
    del self.wc, self.wu
    if(f_a_to_be_r_i is None):
       f_a_to_be_r_i = False
    else:
       # XXX temorary: while we do not have aniso restraints
       f_a_to_be_r_i = not f_a_to_be_r_i
    xray.set_scatterer_grad_flags(scatterers = self.xray_structure.scatterers(),
                                  site       = refine_xyz,
                                  u_iso      = refine_adp,
                                  u_aniso    = f_a_to_be_r_i,
                                  occupancy  = refine_occ)
    self.neutron_refinement = (self.fmodel_neutron is not None and
                               self.wn is not None)
    if(self.neutron_refinement):
       assert self.neutron_scattering_dict is not None and \
              self.xray_scattering_dict is not None
       if(refine_xyz):
          assert wxnc_scale is not None
          self.wn *= wxnc_scale
       if(refine_adp):
          assert wxnu_scale is not None
          self.wn *= wxnu_scale
       self.collector = minimization_history(
                                         wx= self.wx, wn= self.wn, wr= self.wr)
       self.collector.collect(rnw = self.fmodel_neutron.r_work(),
                              rnf = self.fmodel_neutron.r_free(),
                              rxw = self.fmodel.r_work(),
                              rxf = self.fmodel.r_free())
    else:
       self.collector = minimization_history(wx = self.wx, wr = self.wr)
       self.collector.collect(rxw = self.fmodel.r_work(),
                              rxf = self.fmodel.r_free())
    # XXX
    self.d_selection = self.model.atoms_selection(scattering_type = "D")
    if(self.fmodel.alpha_beta_params.method == "calc"):
       if(self.fmodel.alpha_beta_params.fix_scale_for_calc_option == None):
          self.scale_ml = self.fmodel.scale_ml()
       else:
          self.scale_ml = self.fmodel.alpha_beta_params.fix_scale_for_calc_option
    if(self.fmodel.alpha_beta_params.method == "est"):
       self.scale_ml = 1.0
    self.f_obs_w = self.fmodel.f_obs_w
    self.xray_structure.tidy_us(u_min = 1.e-2)
    self.target_name = self.fmodel.target_name
    assert self.target_name in ("ml","mlhl") or self.target_name.count("ls") == 1
    if(self.target_name in ("ml","mlhl", "lsm")):
       if(self.alpha_w is None or self.beta_w is None):
          self.alpha_w, self.beta_w = self.fmodel.alpha_beta_w()
          if(self.neutron_refinement):
             self.alpha_w_neutron, self.beta_w_neutron = \
                                             self.fmodel_neutron.alpha_beta_w()
       else:
          assert self.alpha_w.data().size() == self.f_obs_w.data().size()
          assert self.beta_w.data().size() == self.f_obs_w.data().size()
    self.x = flex.double(self.xray_structure.n_parameters_XXX(), 0)
    self._scatterers_start = self.xray_structure.scatterers()
    self._lock_for_line_search = False
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator          = self,
      termination_params        = lbfgs_termination_params,
      use_fortran               = use_fortran,
      exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
                         ignore_line_search_failed_step_at_lower_bound = True))
    self.apply_shifts()
    del self._scatterers_start
    self._lock_for_line_search = False
    self.compute_target(compute_gradients = False,u_iso_reinable_params = None)
    del self._lock_for_line_search
    self.collector.collect(et   = self.f,
                           iter = self.minimizer.iter(),
                           nfun = self.minimizer.nfun(),
                           rxw  = self.fmodel.r_work(),
                           rxf  = self.fmodel.r_free())
    if(self.neutron_refinement):
       self.collector.collect(rnw = self.fmodel_neutron.r_work(),
                              rnf = self.fmodel_neutron.r_free())
    time_site_individual += timer.elapsed()

  def apply_shifts(self):
    selection = None
    if(self.refine_occ):
       if(self.model.refinement_flags.occupancies_individual is not None):
          selection = self.model.refinement_flags.occupancies_individual[0]
    if(self.refine_xyz):
      if(self.model.refinement_flags.sites_individual is not None):
         selection = self.model.refinement_flags.sites_individual[0]
    if(self.refine_adp):
       if(self.model.refinement_flags.adp_individual is not None):
         selection = self.model.refinement_flags.adp_individual[0]
    apply_shifts_result = xray.ext.minimization_apply_shifts(
                              unit_cell      = self.xray_structure.unit_cell(),
                              scatterers     = self._scatterers_start,
                              shifts         = self.x,
                              refine_xyz     = self.refine_xyz,
                              refine_adp     = self.refine_adp,
                              refine_occ     = self.refine_occ,
                              selection      = selection)
    scatterers_shifted = apply_shifts_result.shifted_scatterers
    if(self.refine_xyz):
       site_symmetry_table = self.xray_structure.site_symmetry_table()
       for i_seq in site_symmetry_table.special_position_indices():
         scatterers_shifted[i_seq].site = crystal.correct_special_position(
                crystal_symmetry = self.xray_structure,
                special_op       = site_symmetry_table.get(i_seq).special_op(),
                site_frac        = scatterers_shifted[i_seq].site)
    self.xray_structure.replace_scatterers(scatterers = scatterers_shifted)
    if(self.refine_adp):
       return apply_shifts_result.u_iso_reinable_params
    else:
       return None

  def compute_target(self, compute_gradients, u_iso_reinable_params):
    selection = None
    if(self.refine_occ):
       if(self.model.refinement_flags.occupancies_individual is not None):
          selection = self.model.refinement_flags.occupancies_individual[0]
    if(self.refine_xyz):
      if(self.model.refinement_flags.sites_individual is not None):
         selection = self.model.refinement_flags.sites_individual[0]
    if(self.refine_adp):
       self.xray_structure.adjust_u_iso()
       if(self.model.refinement_flags.adp_individual is not None):
         selection = self.model.refinement_flags.adp_individual[0]
    #if(self.refine_occ): self.xray_structure.adjust_occupancy()
    self.stereochemistry_residuals = None
    if(self.neutron_refinement):
       self.xray_structure.scattering_type_registry(
                                       custom_dict = self.xray_scattering_dict)
    self.fmodel.update_xray_structure(xray_structure = self.xray_structure,
                                      update_f_calc  = True)
    ex = self.fmodel.target_w(alpha    = self.alpha_w,
                              beta     = self.beta_w,
                              scale_ml = self.scale_ml)
    self.collector.collect(ex = ex)
    self.f = ex * self.wx
    if(compute_gradients):
       sf = self.fmodel.gradient_wrt_atomic_parameters(
                        alpha                 = self.alpha_w,
                        beta                  = self.beta_w,
                        u_iso_reinable_params = u_iso_reinable_params).packed()
       if(self.selection is not None):
          sf_v3d = flex.vec3_double(sf)
          sf_v3d_sel = sf_v3d.set_selected(self.selection, [0,0,0])
          sf = sf_v3d_sel.as_double()
       self.g = sf * self.wx

    if(self.neutron_refinement):
       self.xray_structure.scattering_type_registry(
                                    custom_dict = self.neutron_scattering_dict)
       self.fmodel_neutron.update_xray_structure(
                                         xray_structure = self.xray_structure,
                                         update_f_calc  = True)
       en = self.fmodel_neutron.target_w(alpha    = self.alpha_w_neutron,
                                         beta     = self.beta_w_neutron,
                                         scale_ml = self.scale_ml)
       self.collector.collect(en = en)
       self.f += en * self.wn
       if(compute_gradients):
          sf = self.fmodel_neutron.gradient_wrt_atomic_parameters(
                        alpha                 = self.alpha_w_neutron,
                        beta                  = self.beta_w_neutron,
                        u_iso_reinable_params = u_iso_reinable_params).packed()
          if(self.selection is not None):
             sf = sf.set_selected(self.selection, 0.0)
          self.g = self.g + sf * self.wn

          ii = 0
          if(self.refine_xyz):
             tt = flex.vec3_double(self.g)
             rr = flex.vec3_double(sf)
             for i, j, fd in zip(tt, rr, self.d_selection):
               angle = flex.double(i).angle(flex.double(j), deg = True)
               if(angle >= 90.0):
                  if(fd): tt[ii] = [0,0,0]
                  else:   rr[ii] = [0,0,0]
                  #tt[ii] = [0,0,0]
                  #rr[ii] = [0,0,0]
                  #if(fd): tt[ii] = [tt[ii][0]/angle,tt[ii][1]/angle,tt[ii][2]/angle]
                  #else:   rr[ii] = [rr[ii][0]/angle,rr[ii][1]/angle,rr[ii][2]/angle]
                  #tt[ii] = [tt[ii][0]/angle,tt[ii][1]/angle,tt[ii][2]/angle]
                  #rr[ii] = [rr[ii][0]/angle,rr[ii][1]/angle,rr[ii][2]/angle]
               ii += 1
             self.g = tt.as_double() + rr.as_double() * self.wn
          if(self.refine_adp):
             tt = self.g
             rr = sf
             for i, j, fd in zip(tt, rr, self.d_selection):
               if(i*j < 0):
                  if(fd): tt[ii] = 0.
                  else:   rr[ii] = 0.
                  #tt[ii] = 0.0
                  #rr[ii] = 0.0
                  #if(fd): tt[ii] = tt[ii]/3.
                  #else:   rr[ii] = rr[ii]/3.
                  #tt[ii] = tt[ii]/3.
                  #rr[ii] = rr[ii]/3.
               ii += 1
             self.g = tt + rr * self.wn

    if(self.refine_xyz and self.restraints_manager is not None and self.wr > 0.0):
       self.stereochemistry_residuals = self.restraints_manager.energies_sites(
                       sites_cart        = self.xray_structure.sites_cart(),
                       compute_gradients = compute_gradients,
                       lock_for_line_search = self._lock_for_line_search)
       self._lock_for_line_search = True
       er = self.stereochemistry_residuals.target
       self.collector.collect(er = er)
       self.f += er * self.wr
       if(compute_gradients):
          sgc = self.stereochemistry_residuals.gradients
          if(self.selection is not None):
             sgc = sgc.set_selected(self.selection, [0,0,0])
          xray.minimization.add_gradients(
                             scatterers     = self.xray_structure.scatterers(),
                             xray_gradients = self.g,
                             site_gradients = sgc*self.wr,
                             refine_xyz     = self.refine_xyz,
                             refine_adp     = self.refine_adp,
                             refine_occ     = self.refine_occ,
                             selection      = selection)

    if(self.refine_adp and self.restraints_manager.geometry is not None
                        and self.wr > 0.0 and self.iso_restraints is not None):
       energies_adp_iso = self.restraints_manager.energies_adp_iso(
                    xray_structure    = self.xray_structure,
                    parameters        = self.iso_restraints,
                    wilson_b          = self.wilson_b,
                    tan_b_iso_max     = self.tan_b_iso_max,
                    use_u_local_only  = self.iso_restraints.use_u_local_only,
                    compute_gradients = compute_gradients)
       er = energies_adp_iso.target
       self.collector.collect(er = er)
       self.f += er * self.wr
       if(compute_gradients):
          sgu = energies_adp_iso.gradients
          if(self.selection is not None):
             sgu = sgu.set_selected(self.selection, 0.0)
          xray.minimization.add_gradients(
                        scatterers      = self.xray_structure.scatterers(),
                        xray_gradients  = self.g,
                        u_iso_gradients = sgu * self.wr,
                        refine_xyz     = self.refine_xyz,
                        refine_adp     = self.refine_adp,
                        refine_occ     = self.refine_occ,
                        selection      = selection)

  def callback_after_step(self, minimizer):
    self._lock_for_line_search = False
    if (self.verbose > 0):
      print "refinement.minimization step: f,iter,nfun:",
      print self.f,minimizer.iter(),minimizer.nfun()

  def compute_functional_and_gradients(self):
    u_iso_reinable_params = self.apply_shifts()
    self.compute_target(compute_gradients     = True,
                        u_iso_reinable_params = u_iso_reinable_params)
    self.collector.collect(et = self.f)
    if (self.verbose > 1):
      print "xray.minimization line search: f,rms(g):",
      print self.f, math.sqrt(flex.mean_sq(self.g))
    return self.f, self.g

class minimization_history(object):
  def __init__(self, wx = None,
                     wn = None,
                     wr = None,
                     iter = None,
                     nfun = None):
    adopt_init_args(self, locals())
    self.ex = []
    self.en = []
    self.er = []
    self.et = []
    self.rxw = []
    self.rxf = []
    self.rnw = []
    self.rnf = []

  def collect(self, ex=None,en=None,er=None,rxw=None,rxf=None,rnw=None,
              rnf=None,iter=None,nfun=None,et=None):
    if(ex is not None): self.ex.append(ex)
    if(en is not None): self.en.append(en)
    if(er is not None): self.er.append(er)
    if(et is not None): self.et.append(et)
    if(rxw is not None): self.rxw.append(rxw)
    if(rxf is not None): self.rxf.append(rxf)
    if(rnw is not None): self.rnw.append(rnw)
    if(rnf is not None): self.rnf.append(rnf)
    if(iter is not None): self.iter = iter
    if(nfun is not None): self.nfun = nfun

  def first_last(self, a):
    if(len(a) > 1):
       first = a[0]
       last  = a[len(a)-1]
    elif(len(a) == 1):
       first = a[0]
       last  = a[0]
    else:
       first = None
       last  = None
    return first, last

  def check_size(self, a, b):
    if(len(a) > 0 and len(b) > 0): assert len(a) == len(b)

  def check_consistency(self):
    assert len(self.rxw) == len(self.rxf)
    assert len(self.rnw) == len(self.rnf)
    self.check_size(a = self.rxw, b = self.rnw)
    self.check_size(a = self.rxf, b = self.rnf)
    self.check_size(a = self.ex, b = self.er)
    self.check_size(a = self.en, b = self.er)
    self.check_size(a = self.ex, b = self.en)
    self.check_size(a = self.ex, b = self.et)
    self.check_size(a = self.en, b = self.et)
    self.ex_first,  self.ex_last  = self.first_last(a = self.ex)
    self.en_first,  self.en_last  = self.first_last(a = self.en)
    self.er_first,  self.er_last  = self.first_last(a = self.er)
    self.et_first,  self.et_last  = self.first_last(a = self.et)
    self.rxw_first, self.rxw_last = self.first_last(a = self.rxw)
    self.rxf_first, self.rxf_last = self.first_last(a = self.rxf)
    self.rnw_first, self.rnw_last = self.first_last(a = self.rnw)
    self.rnf_first, self.rnf_last = self.first_last(a = self.rnf)
    neutron_refinement = len(self.en) > 0
    if(not neutron_refinement):
       for et, ex, er in zip(self.et, self.ex, self.er):
           et_sum = ex * self.wx + er * self.wr
           assert approx_equal(et, et_sum)
    else:
       for et, ex, en, er in zip(self.et, self.ex, self.en, self.er):
           et_sum = ex * self.wx + en * self.wn + er * self.wr
           assert approx_equal(et, et_sum)

  def show(self, text=None, out=None):
    if (out is None): out = sys.stdout
    self.check_consistency()
    neutron_refinement = len(self.en) > 0
    line_len = len("| "+text+"|")
    fill_len = 80 - line_len-1
    print >> out, "| "+text+"-"*(fill_len)+"|"
    if(neutron_refinement):
       print >> out, "|                                xray data:                                   |"
    print >> out, "| start r-factor (work) = %6.4f      final r-factor "\
     "(work) = %6.4f          |"%(self.rxw_first, self.rxw_last)
    print >> out, "| start r-factor (free) = %6.4f      final r-factor "\
     "(free) = %6.4f          |"%(self.rxf_first, self.rxf_last)
    if(neutron_refinement):
       print >> out, "|                              neutron data:                                  |"
       print >> out, "| start r-factor (work) = %6.4f      final r-factor "\
       "(work) = %6.4f          |"%(self.rnw_first, self.rnw_last)
       print >> out, "| start r-factor (free) = %6.4f      final r-factor "\
       "(free) = %6.4f          |"%(self.rnf_first, self.rnf_last)
    print >> out, "|"+"-"*77+"|"
    et = self.nn(self.et_first)
    wx = self.nn(self.wx)
    wn = self.nn(self.wn)
    ex = self.nn(self.ex_first)
    en = self.nn(self.en_first)
    wc = self.nn(self.wr)
    ec = self.nn(self.er_first)
    if(neutron_refinement):
       nr = "+ wnc"+" "*abs(len("+ wnc")-len(wn)-3)+\
                            "* Eneutron "+" "*abs(len("* Eneutron ")-len(en)-3)
    else: nr = ""
    h = "| T_start "+" "*abs(len("T_start")-len(et))+\
        "= wxc"+" "*abs(len("= wxc")-len(wx)-3)+\
        "* Exray "+" "*abs(len("* Exray ")-len(ex)-3)+\
        nr+\
        "+ wc "+" "*abs(len("+ wc ")-len(wc)-3)+\
        "* Echem "+" "*abs(len("* Echem ")-len(ec)-3)
    print >> out, h + " "*(78-len(h))+"|"
    if(not neutron_refinement):
       format = "| %s = %s * %s + %s * %s"
       st = format % (et,wx,ex,wc,ec)
    else:
       format = "| %s = %s * %s + %s * %s + %s * %s"
       st = format % (et,wx,ex,wn,en,wc,ec)
    print >> out, st + " "*(78-len(st))+"|"

    et = self.nn(self.et_last)
    ex = self.nn(self.ex_last)
    en = self.nn(self.en_last)
    ec = self.nn(self.er_last)
    if(neutron_refinement):
       nr = "+ wnc"+" "*abs(len("+ wnc")-len(wn)-3)+\
                            "* Eneutron "+" "*abs(len("* Eneutron ")-len(en)-3)
    else: nr = ""
    print >> out, "|"+"-"*77+"|"
    h = "| T_final "+" "*abs(len("T_final")-len(et))+\
        "= wxc"+" "*abs(len("= wxc")-len(wx)-3)+\
        "* Exray "+" "*abs(len("* Exray ")-len(ex)-3)+\
        nr+\
        "+ wc "+" "*abs(len("+ wc ")-len(wc)-3)+\
        "* Echem "+" "*abs(len("* Echem ")-len(ec)-3)
    print >> out, h + " "*(78-len(h))+"|"
    if(not neutron_refinement):
       format = "| %s = %s * %s + %s * %s"
       st = format % (et,wx,ex,wc,ec)
    else:
       format = "| %s = %s * %s + %s * %s + %s * %s"
       st = format % (et,wx,ex,wn,en,wc,ec)
    print >> out, st + " "*(78-len(st))+"|"
    print >> out, "|"+"-"*77+"|"
    print >> out, "| number of iterations = %4d    |    number of function "\
                  "evaluations = %4d   |"%(self.iter, self.nfun)
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def nn(self, x):
    x = str(x)
    try:
      result = x[:len(x[:x.index(".")])+5]
      while len(result) < 8: result += "0"
    except: result = "undef"
    return result
