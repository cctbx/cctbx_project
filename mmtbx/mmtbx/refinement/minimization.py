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
import cctbx.adp_restraints
from cctbx import adptbx

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
                     refine_dbe               = False,
                     lbfgs_termination_params = None,
                     use_fortran              = False,
                     verbose                  = 0,
                     iso_restraints           = None,
                     fmodel_neutron           = None,
                     wn                       = None,
                     neutron_scattering_dict  = None,
                     xray_scattering_dict     = None,
                     wxnc_scale               = None,
                     wxnu_scale               = None,
                     occupancy_max            = None,
                     occupancy_min            = None,
                     h_params                 = None,
                     use_xn_grads_filtering   = None,
                     u_min                    = adptbx.b_as_u(-100.0),
                     u_max                    = adptbx.b_as_u(1000.0)):
    global time_site_individual
    timer = user_plus_sys_time()
    adopt_init_args(self, locals())
    assert [refine_xyz, refine_adp, refine_occ].count(False) == 2
    assert [refine_xyz, refine_adp, refine_occ].count(True)  == 1
    self.xray_structure = self.fmodel.xray_structure
    self.xray_structure.tidy_us()
    if(refine_xyz): self.wr = wc
    if(refine_adp): self.wr = wu
    if(refine_occ): self.wr = None
    self.wxc_dbe = None
    del self.wc, self.wu
    self.hd_selection = self.xray_structure.hd_selection()
    self.hd_flag = self.hd_selection.count(True) > 0
    if(self.hd_selection.count(True) > 0):
       self.xh_connectivity_table = xh_connectivity_table(
                                    geometry       = restraints_manager,
                                    xray_structure = self.xray_structure).table
    self.xray_structure.scatterers().flags_set_grads(state=False)
    if (refine_xyz):
      sel = flex.bool(self.model.refinement_flags.sites_individual[0].size(), False)
      for m in self.model.refinement_flags.sites_individual:
         sel = sel | m
      self.hd_selection = self.hd_selection.select(sel)
      #if (self.h_params.mode == "riding"):
      #  sel.set_selected(self.hd_selection, False)
      self.xray_structure.scatterers().flags_set_grad_site(
        iselection=sel.iselection())
      del sel
    if (refine_occ):
      sel = flex.bool(self.model.refinement_flags.occupancies_individual[0].size(), False)
      for m in self.model.refinement_flags.occupancies_individual:
         sel = sel | m
      self.xray_structure.scatterers().flags_set_grad_occupancy(
        iselection=sel.iselection())
      del sel
    if (refine_adp):
      sel = self.model.refinement_flags.adp_individual_iso[0]
      if (self.h_params.mode == "riding"):
        sel.set_selected(self.hd_selection, False)
      self.xray_structure.scatterers().flags_set_grad_u_iso(
        iselection=sel.iselection())
      #
      sel = self.model.refinement_flags.adp_individual_aniso[0]
      if (self.h_params.mode == "riding"):
        sel.set_selected(self.hd_selection, False)
      self.xray_structure.scatterers().flags_set_grad_u_aniso(
        iselection=sel.iselection())
      del sel
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
       self.collector = minimization_history(wx = self.wx,
                                             wn = self.wn,
                                             wr = self.wr,
                                             refine_xyz = self.refine_xyz,
                                             refine_adp = self.refine_adp,
                                             refine_occ = self.refine_occ)
       self.collector.collect(rnw = self.fmodel_neutron.r_work(),
                              rnf = self.fmodel_neutron.r_free(),
                              rxw = self.fmodel.r_work(),
                              rxf = self.fmodel.r_free())
    else:
       self.collector = minimization_history(wx = self.wx,
                                             wr = self.wr,
                                             refine_xyz = self.refine_xyz,
                                             refine_adp = self.refine_adp,
                                             refine_occ = self.refine_occ)
       self.collector.collect(rxw = self.fmodel.r_work(),
                              rxf = self.fmodel.r_free())
    self.target_functor_xray = self.fmodel.target_functor()
    if (self.neutron_refinement):
      self.target_functor_neutron = self.fmodel_neutron.target_functor()
    self.x = flex.double(self.xray_structure.n_parameters_XXX(), 0)
    self._scatterers_start = self.xray_structure.scatterers()
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator          = self,
      termination_params        = lbfgs_termination_params,
      use_fortran               = use_fortran,
      exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
                         ignore_line_search_failed_step_at_lower_bound = True))
    self.apply_shifts()
    del self._scatterers_start
    self.compute_target(compute_gradients = False,u_iso_refinable_params = None)
    self.xray_structure.tidy_us()
    if(refine_occ):
       self.xray_structure.adjust_occupancy(occ_max = occupancy_max,
                                            occ_min = occupancy_min)
    if(self.neutron_refinement):
       self.xray_structure.scattering_type_registry(
                                       custom_dict = self.xray_scattering_dict)
    self.fmodel.update_xray_structure(update_f_calc=True)
    self.collector.collect(et   = self.f,
                           iter = self.minimizer.iter(),
                           nfun = self.minimizer.nfun(),
                           rxw  = self.fmodel.r_work(),
                           rxf  = self.fmodel.r_free())
    if(self.neutron_refinement):
       # XXX UNSAFE: invalidates self.fmodel
       self.xray_structure.scattering_type_registry(
                                    custom_dict = self.neutron_scattering_dict)
       self.fmodel_neutron.update_xray_structure(
         xray_structure = self.xray_structure,
         update_f_calc  = True)
       self.collector.collect(rnw = self.fmodel_neutron.r_work(),
                              rnf = self.fmodel_neutron.r_free())
       # XXX UNSAFE: invalidates self.fmodel_neutron (but fixes self.fmodel)
       self.xray_structure.scattering_type_registry(
                                       custom_dict = self.xray_scattering_dict)
    del self.target_functor_xray
    if (self.neutron_refinement):
      del self.target_functor_neutron
    time_site_individual += timer.elapsed()

  def apply_shifts(self):
    # XXX inefficient ?
    if(self.refine_adp):
       sel = self.x < self.u_min
       if(sel.count(True) > 0): self.x.set_selected(sel, self.u_min)
       sel = self.x > self.u_max
       if(sel.count(True) > 0): self.x.set_selected(sel, self.u_max)
    # XXX inefficient ?
    # XXX Fix for normal cases at normal resolutions
    if(self.refine_xyz and self.h_params.fix_xh_distances and self.hd_flag and
       self.neutron_scattering_dict is None):
    # THIS LOOKS AS desired to be but does not work!
    #if(self.refine_xyz and self.hd_flag and self.neutron_scattering_dict is None):
       v3d_x = flex.vec3_double(self.x)
       for bond in self.xh_connectivity_table:
           xsh = v3d_x[bond[0]]
           v3d_x[bond[1]] = xsh
       sel = flex.bool(self.x.size(), True)
       self.x.set_selected(sel, v3d_x.as_double())
    apply_shifts_result = xray.ext.minimization_apply_shifts(
                              unit_cell      = self.xray_structure.unit_cell(),
                              scatterers     = self._scatterers_start,
                              shifts         = self.x)
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
       return apply_shifts_result.u_iso_refinable_params
    else:
       return None

  def compute_target(self, compute_gradients, u_iso_refinable_params):
    self.stereochemistry_residuals = None
    if(self.neutron_refinement):
       self.xray_structure.scattering_type_registry(
                                       custom_dict = self.xray_scattering_dict)
    self.fmodel.update_xray_structure(xray_structure = self.xray_structure,
                                      update_f_calc  = True)
    t_r = self.target_functor_xray(compute_gradients=compute_gradients)
    ex = t_r.target_work()
    self.collector.collect(ex = ex)
    self.f = ex * self.wx
    if(compute_gradients):
       sf = t_r.gradients_wrt_atomic_parameters(
         u_iso_refinable_params=u_iso_refinable_params).packed()
       #XXX do not count grads for H or D:
       #if(self.hd_selection.count(True) > 0 and not self.model.use_dbe):
       if(self.hd_flag and self.h_params.mode != "full"):
          if(self.refine_xyz):
             sf_v3d = flex.vec3_double(sf)
             sf_v3d_sel = sf_v3d.set_selected(self.hd_selection, [0,0,0])
             sf = sf_v3d_sel.as_double()
          #if(self.refine_adp):
          #   sf = sf.set_selected(self.hd_selection, 0.0)
       self.g = sf * self.wx

    if(self.neutron_refinement):
       self.xray_structure.scattering_type_registry(
                                    custom_dict = self.neutron_scattering_dict)
       self.fmodel_neutron.update_xray_structure(
                                         xray_structure = self.xray_structure,
                                         update_f_calc  = True)
       t_r = self.target_functor_neutron(compute_gradients=compute_gradients)
       en = t_r.target_work()
       self.collector.collect(en = en)
       self.f += en * self.wn
       if(compute_gradients):
          sf = t_r.gradients_wrt_atomic_parameters(
            u_iso_refinable_params=u_iso_refinable_params).packed()
          self.g = self.g + sf * self.wn

          ii = 0
          if(self.refine_xyz):
             tt = flex.vec3_double(self.g)
             rr = flex.vec3_double(sf)
             if(self.use_xn_grads_filtering):
               for i, j, fd in zip(tt, rr, self.hd_selection):
                 angle = flex.double(i).angle(flex.double(j), deg = True)
                 if(angle is not None and angle >= 90.0):
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
             if(self.use_xn_grads_filtering):
               for i, j, fd in zip(tt, rr, self.hd_selection):
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
       self.stereochemistry_residuals = self.model.restraints_manager_energies_sites(
                       sites_cart           = self.xray_structure.sites_cart(),
                       compute_gradients    = compute_gradients)
       er = self.stereochemistry_residuals.target
       self.collector.collect(er = er)
       self.f += er * self.wr
       if(compute_gradients):
          sgc = self.stereochemistry_residuals.gradients
          xray.minimization.add_gradients(
                             scatterers     = self.xray_structure.scatterers(),
                             xray_gradients = self.g,
                             site_gradients = sgc*self.wr)

    if(self.refine_adp and self.restraints_manager.geometry is not None
                        and self.wr > 0.0 and self.iso_restraints is not None):
       if(self.model.refinement_flags.adp_individual_aniso[0].count(True) == 0):
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
       else:
          energies_adp_aniso = self.restraints_manager.energies_adp_aniso(
                       xray_structure    = self.xray_structure,
                       compute_gradients = compute_gradients)
          er = energies_adp_aniso.target
          self.collector.collect(er = er)
          self.f += er * self.wr

       if(compute_gradients):
          if(self.model.refinement_flags.adp_individual_aniso[0].count(True) == 0):
             sgu = energies_adp_iso.gradients
             #if(self.hd_selection.count(True) > 0):
             #   sgu = sgu.set_selected(self.hd_selection, 0.0)
             xray.minimization.add_gradients(
                           scatterers      = self.xray_structure.scatterers(),
                           xray_gradients  = self.g,
                           u_iso_gradients = sgu * self.wr)
          else:
             ####################################################################
             energies_adp_aniso.gradients_aniso_star *= self.wr
             if(energies_adp_aniso.gradients_iso is not None):
                energies_adp_aniso.gradients_iso *= self.wr
             xray.minimization.add_gradients(
                           scatterers      = self.xray_structure.scatterers(),
                           xray_gradients  = self.g,
                           u_aniso_gradients = energies_adp_aniso.gradients_aniso_star,
                           u_iso_gradients   = energies_adp_aniso.gradients_iso)
             energies_adp_aniso.gradients_aniso_star = None # just for safety
             energies_adp_aniso.gradients_iso = None
             ####################################################################

  def callback_after_step(self, minimizer):
    if (self.verbose > 0):
      print "refinement.minimization step: f,iter,nfun:",
      print self.f,minimizer.iter(),minimizer.nfun()

  def compute_functional_and_gradients(self):
    u_iso_refinable_params = self.apply_shifts()
    self.compute_target(compute_gradients     = True,
                        u_iso_refinable_params = u_iso_refinable_params)
    self.collector.collect(et = self.f)
    if (self.verbose > 1):
      print "xray.minimization line search: f,rms(g):",
      print self.f, math.sqrt(flex.mean_sq(self.g))
    return self.f, self.g

class minimization_history(object):
  def __init__(self, wx        = None,
                     wn        = None,
                     wr        = None,
                     iter      = None,
                     nfun      = None,
                     refine_xyz=False,
                     refine_adp=False,
                     refine_occ=False):
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
    if(self.refine_xyz):
       wrestrx = "= wxc"
       wrestr  = "+ wc "
       energy  = "* Echem "
    elif(self.refine_adp):
       wrestrx = "= wxu"
       wrestr  = "+ wu "
       energy  = "* Eadp "
    else:
       wrestrx = "=    "
       wrestr  = "     "
       energy  = "       "
    h = "| T_start "+" "*abs(len("T_start")-len(et))+\
        wrestrx+" "*abs(len(wrestrx)-len(wx)-3)+\
        "* Exray "+" "*abs(len("* Exray ")-len(ex)-3)+\
        nr+\
        wrestr+" "*abs(len(wrestr)-len(wc)-3)+\
        energy+" "*abs(len(energy)-len(ec)-3)
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
        wrestrx+" "*abs(len(wrestrx)-len(wx)-3)+\
        "* Exray "+" "*abs(len("* Exray ")-len(ex)-3)+\
        nr+\
        wrestr+" "*abs(len(wrestr)-len(wc)-3)+\
        energy+" "*abs(len(energy)-len(ec)-3)
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
    if(x is not None): x = "%15.7f"%x
    x = str(x).strip()
    try:
      if(x is not None): x = "%15.6f"%float(x)
      x = str(x).strip()
      result = x[:len(x[:x.index(".")])+7]
      while len(result) < 8: result += "0"
    except: result = "undef"
    return result

class xh_connectivity_table(object):
  def __init__(self, geometry, xray_structure):
    bond_proxies_simple = geometry.geometry.pair_proxies().bond_proxies.simple
    self.table = []
    scatterers = xray_structure.scatterers()
    for proxy in bond_proxies_simple:
        i_seq, j_seq = proxy.i_seqs
        i_x, i_h = None, None
        if(scatterers[i_seq].element_symbol() == "H"):
           i_h = i_seq
           i_x = j_seq
           self.table.append([i_x, i_h])
        if(scatterers[j_seq].element_symbol() == "H"):
           i_h = j_seq
           i_x = i_seq
           self.table.append([i_x, i_h])
