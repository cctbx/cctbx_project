from __future__ import division, print_function
from mmtbx.refinement import print_statistics
from cctbx.array_family import flex
from libtbx import adopt_init_args
import math
import sys
from libtbx.utils import Sorry, user_plus_sys_time
from libtbx.str_utils import format_value
import iotbx.phil


time_weights_xray_chem_py  = 0.0

core_params_str = """
  optimize_xyz_weight = False
    .alias = optimise_xyz_weight
    .type = bool
    .short_caption = Optimize X-ray/stereochemistry weight
    .style = bold noauto
  optimize_adp_weight = False
    .alias = optimise_adp_weight
    .type = bool
    .short_caption = Optimize X-ray/ADP weight
    .style = bold noauto
  r_free_only = False
    .type = bool
    .help = Use only R-free to choose the best weight
  wxc_scale = 0.5
    .type = float
    .optional = False
    .short_caption = Scale factor for X-ray/stereochemistry weight (wxc_scale)
  wxu_scale = 1.0
    .type = float
    .optional = False
    .short_caption = Scale factor for X-ray/ADP weight (wxu_scale)
  wc = 1.0
    .type = float
    .optional = False
    .short_caption = Stereochemistry weight scale (wc)
  wu = 1.0
    .type = float
    .optional = False
    .short_caption = ADP weight scale (wu)
  fix_wxc = None
    .type = float
    .short_caption = Fix X-ray/stereochemistry weight (wxc)
  fix_wxu = None
    .type = float
    .short_caption = Fix X-ray/ADP weight (wxu)
  shake_sites = True
    .type = bool
    .expert_level=3
  shake_adp = 10.0
    .type = float
    .short_caption = Shake ADPs
    .expert_level=3
  regularize_ncycles = 50
    .type = int
    .expert_level=3
    .short_caption = Number of regularization cycles
  verbose = 1
    .type = int
  wnc_scale = 0.5
    .type = float
    .optional = False
    .short_caption = Scale factor for neutron/stereochemistry weight (wnc_scale)
  wnu_scale = 1.0
    .type = float
    .optional = False
    .short_caption = Scale factor for neutron/ADP weight (wnu_scale)
  rmsd_cutoff_for_gradient_filtering = 3.0
    .type = float
    .expert_level=3
    .short_caption = RMSD cutoff for gradient filtering
  force_optimize_weights = False
    .type = bool
    .expert_level = 3
    .style = hidden
"""

master_params_str = """\
%s
  weight_selection_criteria
    .style = box
  {
    bonds_rmsd = None
      .type=float
      .short_caption = RMS(bonds)
    angles_rmsd = None
      .type=float
      .short_caption = RMS(angles)
    r_free_minus_r_work = None
      .type=float
      .short_caption = R_free - R_work
    r_free_range_width = None
      .type = float
      .short_caption = R-free range width
    mean_diff_b_iso_bonded_fraction = None
      .type = float
    min_diff_b_iso_bonded = None
      .type = float
  }
"""%core_params_str

master_params = iotbx.phil.parse(master_params_str)

class weights(object):
  def __init__(self,
               wx        = None,
               wx_scale  = None,
               angle_x   = None,
               wn        = None,
               wn_scale  = None,
               angle_n   = None,
               w         = None,
               wxn       = None,
               angle_xn  = None,
               angle_xnr = None):
    adopt_init_args(self, locals())

class show(object):
  def __init__(self, adp = None, xyz = None, log = None):
    if(log is None): log = sys.stdout
    if([xyz, adp].count(None)==0):
      print("|"+"-"*77+"|", file=log)
      if(xyz is not None):
        if(xyz.wn is not None):
          print("| XYZ refinement: T = (Exray*wxc_scale + Eneutron*wnc_scale)*wxnc + Echem*wc  |", file=log)
          print("| wxnc = %s   wxc scale = %s   wnc scale = %s   wc = %s |" % self.format_4_values(obj = xyz), file=log)
          print("|-----------------------------------------------------------------------------|", file=log)
        else:
          print("| XYZ refinement: T = Eexperimental * wxc * wxc_scale + Echem * wc"+" "*12+"|", file=log)
          print("| wxc = %s       wxc_scale = %s       wc = %s    |" % self.format_3_values(obj = xyz), file=log)
      if(adp is not None):
        if(adp.wn is not None):
          print("| ADP refinement: T = (Exray*wxu_scale + Eneutron*wnu_scale)*wxnu + Eadp *wu  |", file=log)
          print("| wxnu = %s   wxu scale = %s   wnu scale = %s   wu = %s |" % self.format_4_values(obj = adp), file=log)
        else:
          if(xyz is not None): print("|"+" "*77+"|", file=log)
          print("| ADP refinement: T = Eexperimental * wxu * wxu_scale + Eadp  * wu"+" "*12+"|", file=log)
          print("| wxu = %s       wxu_scale = %s       wu = %s    |" % self.format_3_values(obj = adp), file=log)
      print("|"+"-"*77+"|", file=log)
      print(file=log)

  def format_3_values(self, obj):
    return (format_value("%-15.6f",obj.wx),
            format_value("%-10.3f",obj.wx_scale),
            format_value("%-10.3f",obj.w))

  def format_4_values(self, obj):
    return (format_value("%-8.6f",obj.wxn),
            format_value("%-8.3f",obj.wx_scale),
            format_value("%-8.3f",obj.wn_scale),
            format_value("%-6.3f",obj.w))

class adp_gradients(object):
  def __init__(self, fmodel,
                     model,
                     iso_restraints,
                     shake):
    fmodel_dc = fmodel.deep_copy()
    xray_structure = fmodel_dc.xray_structure
    sel_i = model.refinement_flags.adp_individual_iso
    sel_a = model.refinement_flags.adp_individual_aniso
    if model.ias_manager is not None:
      ias_selection = model.ias_manager.get_ias_selection()
      if ias_selection is not None:
        xray_structure = xray_structure.select(~ias_selection)
        sel_i = sel_i.select(~ias_selection)
        sel_a = sel_a.select(~ias_selection)
    restraints_manager = model.restraints_manager
    hd_sel = xray_structure.hd_selection()
    if(hd_sel.count(True) > 0 and hd_sel.count(True) != hd_sel.size()):
      if(sel_i is not None):
        sel_i = sel_i.select(~hd_sel)
      if(sel_a is not None):
        sel_a = sel_a.select(~hd_sel)
      xray_structure = xray_structure.select(~hd_sel)
      restraints_manager = model.restraints_manager.select(~hd_sel)
      pp = restraints_manager.geometry.pair_proxies(sites_cart =
        xray_structure.sites_cart())
    xray_structure.shake_adp_if_all_equal(b_iso_tolerance = 1.e-3)
    if(shake):
      xray_structure.shake_adp(spread=shake, keep_anisotropic= False)
    scatterers = xray_structure.scatterers()
    scatterers.flags_set_grads(state=False)
    # heavy atoms may overwhelm gradients
    scat_types = xray_structure.scatterers().extract_scattering_types()
    sel_use = (scat_types == "C") | (scat_types == "N") | (scat_types == "O") |\
              (scat_types == "S") | (scat_types == "P")
    if(sel_i is not None):
      if((sel_i & sel_use).count(True)>0):
        sel_i = sel_i & sel_use
      scatterers.flags_set_grad_u_iso(iselection = sel_i.iselection())
    if(sel_a is not None):
      if((sel_a & sel_use).count(True)>0):
        sel_a = sel_a & sel_use
      scatterers.flags_set_grad_u_aniso(iselection = sel_a.iselection())
    fmodel_dc.update_xray_structure(xray_structure = xray_structure,
      update_f_calc = True)
    gxu = fmodel_dc.one_time_gradients_wrt_atomic_parameters().packed()
    gu = None
    if(sel_a is None or sel_a.count(True) == 0):
      restraints_manager.geometry.pair_proxies(sites_cart =
        xray_structure.sites_cart())
      energies_adp = restraints_manager.energies_adp_iso(
        xray_structure    = xray_structure,
        parameters        = iso_restraints,
        use_u_local_only  = iso_restraints.use_u_local_only,
        use_hd            = model.is_neutron(),
        compute_gradients = True)
      if(sel_i is not None):
        gu = energies_adp.gradients.as_double().select(sel_i)
    else:
      restraints_manager.geometry.pair_proxies(sites_cart =
        xray_structure.sites_cart())
      energies_adp = restraints_manager.energies_adp_aniso(
        xray_structure    = xray_structure,
        compute_gradients = True)
      gu_i = None
      if(sel_a is not None):
        gu = energies_adp.gradients_aniso_star.select(sel_a).as_double()
        gu_i = energies_adp.gradients_iso
      if(gu_i is not None):
        if(sel_i is not None):
          gu_i = gu_i.select(sel_i).as_double()
          gu.extend(gu_i)
    self.gu = gu
    self.gxu = gxu
    if([self.gu,self.gxu].count(None)==0):
      assert self.gu.size() == self.gxu.size()

class site_gradients(object):
  def __init__(self, fmodel,
                     model,
                     correct_special_position_tolerance,
                     cartesian_dynamics_parameters,
                     shake,
                     regularize_ncycles,
                     rmsd_cutoff_for_gradient_filtering,
                     gradient_filtering = False,
                     log=None):
    fmodel_dc = fmodel.deep_copy()
    xray_structure = fmodel_dc.xray_structure
    sel_si  = model.refinement_flags.sites_individual
    sel_sta = model.refinement_flags.sites_torsion_angles
    sel_count = [sel_si, sel_sta].count(None)
    assert sel_count != 2
    if (sel_count == 0):
      if (not sel_si.all_eq(sel_sta)):
        raise Sorry(
          "Not implemented: support for different sites.individual"
          " and sites.torsion_angles selections.")
      sel = sel_si
    elif (sel_si is not None):
      sel = sel_si
    else:
      sel = sel_sta
    if model.ias_manager is not None:
      ias_selection = model.ias_manager.get_ias_selection()
      if ias_selection is not None:
        xray_structure = xray_structure.select(~ias_selection)
        sel = sel.select(~ias_selection)
    assert sel.count(True) > 0
    restraints_manager = model.restraints_manager
    if(shake):
      from phenix.refinement import memory_eraser
      geometry = getattr(restraints_manager, "geometry", None)
      xray_structure = memory_eraser.shake_sites(
        cartesian_dynamics_parameters = cartesian_dynamics_parameters,
        restraints_manager            = geometry,
        xray_structure                = xray_structure,
        max_iterations                = regularize_ncycles,
        neutron                       = model.is_neutron(),
        log=log,
        verbose=0)
    scatterers = xray_structure.scatterers()
    scatterers.flags_set_grads(state=False)
    scatterers.flags_set_grad_site(iselection = sel.iselection())
    fmodel_dc.update_xray_structure(xray_structure = xray_structure,
      update_f_calc = True, update_f_mask = True)
    gxc = flex.vec3_double(
      fmodel_dc.one_time_gradients_wrt_atomic_parameters(site = True).packed())
    gco = restraints_manager.energies_sites(
      sites_cart        = xray_structure.sites_cart(),
      compute_gradients = True,
      hd_selection      = xray_structure.hd_selection(), # need for afitt
      )
    gc=gco.gradients
    self.gc = gc.select(sel)
    self.gxc = gxc
    assert self.gc.size() == self.gxc.size()
    assert self.gc.size() == self.gxc.size()
    self.gc_filtered = None
    self.gxc_filtered = None
    if(gradient_filtering):
      self.gc_filtered = flex.sqrt(self.gc.dot())
      gxc_f = flex.sqrt(self.gxc.dot())
      gxc_f_mean = flex.mean(gxc_f)
      #gxc_rms = math.sqrt(flex.mean(self.gxc.dot()))
      if(rmsd_cutoff_for_gradient_filtering is None):
        raise Sorry("rmsd_cutoff_for_gradient_filtering is set to None but a number is expected.")
      gxc_f_selection = (gxc_f < gxc_f_mean*rmsd_cutoff_for_gradient_filtering)
      self.gxc_filtered = gxc_f.select(gxc_f_selection)

class weight:
   def __init__(self, fmodel,
                      model,
                      correct_special_position_tolerance,
                      target_weights_params,
                      macro_cycle,
                      cartesian_dynamics_parameters = None,
                      iso_restraints                = None,
                      amber_params                  = None,
                      log                           = None,
                      show_summary                  = True):
     global time_weights_xray_chem_py
     timer = user_plus_sys_time()
     adopt_init_args(self, locals())
     self.compute_wxc = self.model.refinement_flags.individual_sites \
                     or self.model.refinement_flags.torsion_angles
     self.compute_wxu = self.model.refinement_flags.individual_adp
     self.twp = self.target_weights_params

     if self.iso_restraints is None:
       import mmtbx.refinement.adp_refinement
       self.iso_restraints = mmtbx.refinement.adp_refinement.\
         adp_restraints_master_params.extract().iso

     self.special_case = False
     d_min = self.fmodel.f_obs().d_min()
     if(self.macro_cycle>1 and self.twp.fix_wxc is None and
        d_min>=3.5 and d_min<=4.5):
       self.special_case = True


     if(self.log is None): self.log = sys.stdout
     gxc = None
     gxu = None
     gc = None
     gu = None
     self.adp_weights_result = self.adp_weights()
     wxu       = self.adp_weights_result.wx
     wu        = self.adp_weights_result.w
     self.xyz_weights_result = self.xyz_weights(log=log)
     wxc       = self.xyz_weights_result.wx
     wc        = self.xyz_weights_result.w
     if(show_summary):
       show(xyz = self.xyz_weights_result, adp = self.adp_weights_result,
         log = self.log)
     time_weights_xray_chem_py += timer.elapsed()

   def xyz_weights(self, log=None):
     wxc       = None
     wxc_scale = self.twp.wxc_scale#None
     angle_xc  = None
     wc        = None
     if(self.compute_wxc and self.model.restraints_manager is not None):
        wc = self.twp.wc
        if(self.twp.wxc_scale != 0.0 and
           self.twp.fix_wxc is None and
           self.twp.wc != 0.0
          ):
          if self.special_case:
            result_a = site_gradients(
              fmodel                        = self.fmodel,
              model                         = self.model,
              correct_special_position_tolerance=self.correct_special_position_tolerance,
              cartesian_dynamics_parameters = self.cartesian_dynamics_parameters,
              shake                         = False,
              regularize_ncycles            = self.twp.regularize_ncycles,
              rmsd_cutoff_for_gradient_filtering = self.twp.rmsd_cutoff_for_gradient_filtering,
              gradient_filtering            = False,
              log=self.log)
            angle_xc = self.angle(result_a.gxc, result_a.gc)
          result = site_gradients(
            fmodel                        = self.fmodel,
            model                         = self.model,
            correct_special_position_tolerance=self.correct_special_position_tolerance,
            cartesian_dynamics_parameters = self.cartesian_dynamics_parameters,
            regularize_ncycles            = self.twp.regularize_ncycles,
            shake                         = self.twp.shake_sites,
            rmsd_cutoff_for_gradient_filtering = self.twp.rmsd_cutoff_for_gradient_filtering,
            gradient_filtering            = True,
            log=self.log)
          gc = result.gc_filtered
          gxc = result.gxc_filtered
          gc_norm  = gc.norm()
          gxc_norm = gxc.norm()
          if math.isnan(gc_norm): raise Sorry("Norm of gc gradients is %s" % gc_norm)
          if math.isnan(gxc_norm): raise Sorry("Norm of gxc gradients is %s" % gxc_norm)
          if(gxc_norm != 0.0):
            wxc = gc_norm / gxc_norm
          else:
            wxc = 1.0
          if self.amber_params and self.amber_params.use_amber:
            if wxc > 1000: raise Sorry("wxc is too high: wxc=%f" %wxc)
            print("  Setting wxc_scale using amber.wxc_factor : %0.3f * %0.3f = %0.3f" % (
              wxc_scale,
              self.amber_params.wxc_factor,
              wxc_scale*self.amber_params.wxc_factor,
            ), file=log)
            wxc_scale *= self.amber_params.wxc_factor
        elif(self.twp.fix_wxc is not None):
           wxc = self.twp.fix_wxc
           wxc_scale = 1.0
        elif(self.twp.wxc_scale == 0.0):
           wxc = 0.0
        elif(self.twp.wc == 0.0):
           wxc = 1.0
           wc  = 0.0
           wnc = 1.0
           wxnc = 1.0
           wxc_scale = 1.0
        else:
           raise Sorry("Wrong parameter type, value or combination for "+
            "stereochemistry/X-ray weighting.")
     else:
        wxc = 1.0
        wc  = 0.0
     # case-specific: resolutions 3.5-4.5A
     if(self.special_case and angle_xc is not None):
       wxc_scale_ = math.cos((180-angle_xc)*math.pi/180)
       if(wxc_scale_<=0): wxc_scale_=wxc_scale
       wxc_scale_ = 0.0204 * math.exp(3.7444*wxc_scale_)
       if(wxc_scale_>0): wxc_scale = wxc_scale_
     #
     return weights(wx       = wxc,
                    wx_scale = wxc_scale,
                    w        = wc)


   def adp_weights(self):
     wxu       = None
     wxu_scale = None
     wu        = None
     if(self.compute_wxu):
       wu = self.twp.wu
       if(self.twp.wxu_scale != 0.0 and self.twp.fix_wxu is None and
          self.twp.wu != 0.0):
         result = adp_gradients(
           fmodel         = self.fmodel,
           model          = self.model,
           iso_restraints = self.iso_restraints,
           shake          = self.twp.shake_adp)
         gu = result.gu
         gxu = result.gxu
         gu_norm  = gu.norm()
         gxu_norm = gxu.norm()
         if(gxu_norm != 0.0):
           if(gu_norm != 0):
             wxu = gu_norm / gxu_norm
           else:
             wxu = 1.0
         else:
           wxu = 1.0
       elif(self.twp.fix_wxu is not None):
         wxu = self.twp.fix_wxu
       elif(self.twp.wxu_scale == 0.0):
         wxu = 0.0
       elif(self.twp.wu == 0.0):
         wxu = 1.0
         wu  = 0.0
       else:
         raise Sorry("Wrong parameter type, value or combination for "+
          "ADP/X-ray weighting.")
     else:
       wxu = 1.0
       wu  = 0.0
     return weights(wx       = wxu,
                    wx_scale = self.twp.wxu_scale,
                    w        = wu)

   def angle(self, a, b):
     if([a,b].count(None) > 0): return None
     result = a.as_double().angle(b.as_double())
     if(result is None): return None
     return result * 180 / math.pi

def run(
      params,
      fmodels,
      model,
      macro_cycle,
      prefix,
      log):
  print_statistics.make_header(prefix, out = log)
  result = weight(
    fmodel                        = fmodels.fmodel_xray(),
    model                         = model,
    correct_special_position_tolerance =
      params.main.correct_special_position_tolerance,
    cartesian_dynamics_parameters = params.cartesian_dynamics,
    target_weights_params         = params.target_weights,
    iso_restraints                = params.adp_restraints.iso,
    macro_cycle                   = macro_cycle,
    amber_params                  = getattr(params, "amber", None),
    log                           = log)
  return result
