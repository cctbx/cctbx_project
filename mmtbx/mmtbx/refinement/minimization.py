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
                     fmodels,
                     model,
                     target_weights           = None,
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
                     occupancy_max            = None,
                     occupancy_min            = None,
                     h_params                 = None,
                     u_min                    = adptbx.b_as_u(-100.0),
                     u_max                    = adptbx.b_as_u(1000.0)):
    global time_site_individual
    timer = user_plus_sys_time()
    adopt_init_args(self, locals())
    fmodels.create_target_functors()
    assert [refine_xyz, refine_adp, refine_occ].count(False) == 2
    assert [refine_xyz, refine_adp, refine_occ].count(True)  == 1
    self.xray_structure = self.fmodels.fmodel_xray().xray_structure
    self.xray_structure.tidy_us()
    self.weights = None
    if(refine_xyz):   self.weights = target_weights.xyz_weights_result
    elif(refine_adp): self.weights = target_weights.adp_weights_result
    else:
      from phenix.refinement import weight_xray_chem
      self.weights = weight_xray_chem.weights(wx       = 1,
                                              wx_scale = 1,
                                              angle_x  = None,
                                              wn       = 1,
                                              wn_scale = 1,
                                              angle_n  = None,
                                              w        = 0,
                                              wxn      = 1)
    self.wxc_dbe = None
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
    self.neutron_refinement = (self.fmodels.fmodel_n is not None)
    if(self.neutron_refinement):
       self.collector = minimization_history(wx = self.weights.wx * self.weights.wx_scale,
                                             wn = self.weights.wn * self.weights.wn_scale,
                                             wr = self.weights.w,
                                             refine_xyz = self.refine_xyz,
                                             refine_adp = self.refine_adp,
                                             refine_occ = self.refine_occ)
       self.collector.collect(rnw = self.fmodels.fmodel_neutron().r_work(),
                              rnf = self.fmodels.fmodel_neutron().r_free(),
                              rxw = self.fmodels.fmodel_xray().r_work(),
                              rxf = self.fmodels.fmodel_xray().r_free())
    else:
       self.collector = minimization_history(wx = self.weights.wx*self.weights.wx_scale,
                                             wr = self.weights.w,
                                             refine_xyz = self.refine_xyz,
                                             refine_adp = self.refine_adp,
                                             refine_occ = self.refine_occ)
       self.collector.collect(rxw = self.fmodels.fmodel_xray().r_work(),
                              rxf = self.fmodels.fmodel_xray().r_free())
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
    self.fmodels.update_xray_structure(
      update_f_calc  = True,
      xray_structure = self.xray_structure)
    self.collector.collect(et   = self.f,
                           iter = self.minimizer.iter(),
                           nfun = self.minimizer.nfun(),
                           rxw  = self.fmodels.fmodel_xray().r_work(),
                           rxf  = self.fmodels.fmodel_xray().r_free())
    if(self.neutron_refinement):
       self.collector.collect(rnw = self.fmodels.fmodel_neutron().r_work(),
                              rnf = self.fmodels.fmodel_neutron().r_free())
    time_site_individual += timer.elapsed()

  def apply_shifts(self):
    # XXX inefficient
    if(self.refine_adp):
       sel = self.x < self.u_min
       if(sel.count(True) > 0): self.x.set_selected(sel, self.u_min)
       sel = self.x > self.u_max
       if(sel.count(True) > 0): self.x.set_selected(sel, self.u_max)
    # XXX inefficient
    # XXX Fix for normal cases at normal resolutions
    if(self.refine_xyz and self.h_params.fix_xh_distances and self.hd_flag):
    # THIS LOOKS AS desired to be but does not work!
    #if(self.refine_xyz and self.hd_flag):
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
    h_flag = self.hd_flag and self.h_params.mode != "full" and self.refine_xyz
    self.stereochemistry_residuals = None
    self.fmodels.update_xray_structure(xray_structure = self.xray_structure,
                                       update_f_calc  = True)
    fmodels_target_and_gradients = self.fmodels.target_and_gradients(
      weights                = self.weights,
      compute_gradients      = compute_gradients,
      hd_selection           = self.hd_selection,
      h_flag                 = h_flag,
      u_iso_refinable_params = u_iso_refinable_params)
    self.f = fmodels_target_and_gradients.target()
    self.g = fmodels_target_and_gradients.gradients()
    self.collector.collect(ex = fmodels_target_and_gradients.target_work_xray)
    if(self.neutron_refinement):
      self.collector.collect(en = fmodels_target_and_gradients.target_work_neutron)
    if(self.refine_xyz and self.restraints_manager is not None and
       self.weights.w > 0.0):
      self.stereochemistry_residuals = \
        self.model.restraints_manager_energies_sites(
          compute_gradients = compute_gradients)
      er = self.stereochemistry_residuals.target
      self.collector.collect(er = er)
      self.f += er * self.weights.w
      if(compute_gradients):
        sgc = self.stereochemistry_residuals.gradients
        xray.minimization.add_gradients(
          scatterers     = self.xray_structure.scatterers(),
          xray_gradients = self.g,
          site_gradients = sgc*self.weights.w)
    if(self.refine_adp and self.restraints_manager.geometry is not None
       and self.weights.w > 0.0 and self.iso_restraints is not None):
      energies_adp = self.model.energies_adp(
        iso_restraints    = self.iso_restraints,
        compute_gradients = compute_gradients)
      self.f += energies_adp.target * self.weights.w
      if(compute_gradients):
        if(energies_adp.u_aniso_gradients is None):
          xray.minimization.add_gradients(
            scatterers      = self.xray_structure.scatterers(),
            xray_gradients  = self.g,
            u_iso_gradients = energies_adp.u_iso_gradients * self.weights.w)
        else:
          energies_adp.u_aniso_gradients *= self.weights.w
          if(energies_adp.u_iso_gradients is not None):
            energies_adp.u_iso_gradients *= self.weights.w
          xray.minimization.add_gradients(
            scatterers        = self.xray_structure.scatterers(),
            xray_gradients    = self.g,
            u_aniso_gradients = energies_adp.u_aniso_gradients,
            u_iso_gradients   = energies_adp.u_iso_gradients)
          energies_adp.u_aniso_gradients = None # just for safety
          energies_adp.u_iso_gradients = None

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


#class monitor(object):
#  def __init__(self, weights,
#                     refine_xyz = False,
#                     refine_adp = False,
#                     refine_occ = False):

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
           #assert approx_equal(et, et_sum)
    else:
       for et, ex, en, er in zip(self.et, self.ex, self.en, self.er):
           et_sum = ex * self.wx + en * self.wn + er * self.wr
           #assert approx_equal(et, et_sum)

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
