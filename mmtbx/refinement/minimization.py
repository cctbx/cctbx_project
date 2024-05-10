from __future__ import absolute_import, division, print_function
from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex
import scitbx.lbfgs
from libtbx import adopt_init_args
import math
import sys
from libtbx.utils import user_plus_sys_time
from cctbx import adptbx
from libtbx.str_utils import format_value
from libtbx.utils import Sorry
import mmtbx.refinement.minimization_ncs_constraints

class lbfgs(object):

  def __init__(self, fmodels,
                     restraints_manager       = None,
                     model                    = None,
                     is_neutron_scat_table    = None,
                     target_weights           = None,
                     refine_xyz               = False,
                     refine_adp               = False,
                     lbfgs_termination_params = None,
                     use_fortran              = False,
                     verbose                  = 0,
                     correct_special_position_tolerance = 1.0,
                     iso_restraints           = None,
                     h_params                 = None,
                     qblib_params             = None,
                     macro_cycle              = None,
                     u_min                    = adptbx.b_as_u(-5.0),
                     u_max                    = adptbx.b_as_u(999.99),
                     collect_monitor          = True,
                     log                      = None):
    timer = user_plus_sys_time()
    adopt_init_args(self, locals())
    self.f=None
    self.xray_structure = self.fmodels.fmodel_xray().xray_structure
    self.fmodels.create_target_functors()
    self.fmodels.prepare_target_functors_for_minimization()
    if(self.refine_adp and fmodels.fmodel_neutron() is None):
      self.xray_structure.tidy_us()
      self.fmodels.update_xray_structure(
        xray_structure = self.xray_structure,
        update_f_calc  = True)
    self.weights = None
# QBLIB INSERT
    self.qblib_params = qblib_params
    if(self.qblib_params is not None and self.qblib_params.qblib):
        self.macro = macro_cycle
        self.qblib_cycle_count = 0
        self.tmp_XYZ = None
        self.XYZ_diff_curr=None
# QBLIB END
    self.correct_special_position_tolerance = correct_special_position_tolerance
    if(refine_xyz and target_weights is not None):
      self.weights = target_weights.xyz_weights_result
    elif(refine_adp and target_weights is not None):
      self.weights = target_weights.adp_weights_result
    else:
      from mmtbx.refinement import weights
      self.weights = weights.weights(wx = 1, wx_scale = 1, w = 0)
    if(self.collect_monitor):
      self.monitor = monitor(
        weights        = self.weights,
        fmodels        = fmodels,
        model          = model,
        iso_restraints = iso_restraints,
        refine_xyz     = refine_xyz,
        refine_adp     = refine_adp,
        refine_occ     = False)
    if(self.collect_monitor): self.monitor.collect()
    self.neutron_refinement = (self.fmodels.fmodel_n is not None)
    self.x = flex.double(self.xray_structure.n_parameters(), 0)
    self.riding_h_manager = self.model.riding_h_manager
    self._scatterers_start = self.xray_structure.scatterers()
    #lbfgs_core_params = scitbx.lbfgs.core_parameters(
    #  stpmin = 1.e-9,
    #  stpmax = adptbx.b_as_u(10))
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator          = self,
      termination_params        = lbfgs_termination_params,
      use_fortran               = use_fortran,
    #  core_params = lbfgs_core_params,
    #  gradient_only=True,
      exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
                         ignore_line_search_failed_step_at_lower_bound = True))
    self.apply_shifts()
    del self._scatterers_start
    self.compute_target(compute_gradients = False,u_iso_refinable_params = None)
    if(self.refine_adp and self.fmodels.fmodel_neutron() is None):
      self.xray_structure.tidy_us()
      self.fmodels.update_xray_structure(
        xray_structure = self.xray_structure,
        update_f_calc  = True)
    if(self.collect_monitor):
      self.monitor.collect(iter = self.minimizer.iter(),
                           nfun = self.minimizer.nfun())
    self.fmodels.create_target_functors()
# QBLIB INSERT
    if(self.qblib_params is not None and self.qblib_params.qblib):
       print('{:-^80}'.format(""), file=self.qblib_params.qblib_log)
       print(file=self.qblib_params.qblib_log)
# QBLIB END

  def apply_shifts(self):
    if(self.refine_adp):
      xray.ext.truncate_shifts(
        shifts    = self.x,
        min_value = self.u_min,
        max_value = self.u_max)
    if(self.riding_h_manager is not None and self.refine_xyz):
      # temporary array with zero shifts for all atoms
      tmp = flex.vec3_double(self._scatterers_start.size(), [0,0,0])
      # Set shifts according to self.x
      tmp = tmp.set_selected(self.model.refinement_flags.sites_individual,
        flex.vec3_double(self.x))
      # Set shifts of H atoms to zero (they will be idealized later anyway)
      shifts_h = flex.vec3_double(self.riding_h_manager.hd_selection.count(True),
        [0,0,0])
      tmp = tmp.set_selected(self.riding_h_manager.hd_selection, shifts_h)
      # Select entries for shifted atoms only (e.g. if there is xyz selection)
      tmp_x = tmp.select(self.model.refinement_flags.sites_individual)
      # orig
      #tmp = tmp.set_selected(self.riding_h_manager.not_hd_selection,
      #  flex.vec3_double(self.x))
      # Set entries in self.x corresponding to H atoms to 0
      self.x.set_selected(flex.bool(self.x.size()), tmp_x.as_double())
      # orig
      #self.x.set_selected(flex.bool(self.x.size()), tmp.as_double())

    apply_shifts_result = xray.ext.minimization_apply_shifts(
      unit_cell      = self.xray_structure.unit_cell(),
      scatterers     = self._scatterers_start,
      shifts         = self.x)
    scatterers_shifted = apply_shifts_result.shifted_scatterers
    if(self.refine_xyz):
      site_symmetry_table = self.xray_structure.site_symmetry_table()
      for i_seq in site_symmetry_table.special_position_indices():
        try:
          scatterers_shifted[i_seq].site = crystal.correct_special_position(
            crystal_symmetry = self.xray_structure,
            special_op       = site_symmetry_table.get(i_seq).special_op(),
            site_frac        = scatterers_shifted[i_seq].site,
            site_label       = scatterers_shifted[i_seq].label,
            tolerance        = self.correct_special_position_tolerance)
        except Exception as e:
          print(str(e), file=self.log)
    self.xray_structure.replace_scatterers(scatterers = scatterers_shifted)
    #
    if(self.riding_h_manager is not None and self.refine_xyz):
      sites_cart = self.xray_structure.sites_cart()
      self.riding_h_manager.idealize_riding_h_positions(
        sites_cart = sites_cart,
        selection_bool = self.model.refinement_flags.sites_individual)
      self.xray_structure.set_sites_cart(sites_cart=sites_cart)
    #
    if(self.refine_adp):
      return apply_shifts_result.u_iso_refinable_params
    else:
      return None

  def compute_target(self, compute_gradients, u_iso_refinable_params):
    self.stereochemistry_residuals = None
    self.fmodels.update_xray_structure(
      xray_structure = self.xray_structure,
      update_f_calc  = True)
    fmodels_target_and_gradients = self.fmodels.target_and_gradients(
      weights                = self.weights,
      compute_gradients      = compute_gradients,
      u_iso_refinable_params = u_iso_refinable_params)
    self.f = fmodels_target_and_gradients.target()
    self.g = fmodels_target_and_gradients.gradients()
    if(self.refine_xyz and self.restraints_manager is not None and
       self.weights.w > 0.0):
      self.stereochemistry_residuals = \
        self.model.restraints_manager_energies_sites(
          compute_gradients = compute_gradients)
# QBLIB INSERT
      if(self.qblib_params is not None and self.qblib_params.qblib):
        from qbpy import qb_refinement
        self.qblib_cycle_count +=1
        if(self.tmp_XYZ is not None):
          diff,max_diff,max_elem = qb_refinement.array_difference(
            self.tmp_XYZ,
            self.model.get_sites_cart(),
            )
          if(self.XYZ_diff_curr is not None):
            self.XYZ_diff_curr=max_elem
          else:
            self.XYZ_diff_curr=max_elem
          if(max_elem>=self.qblib_params.skip_divcon_threshold):
               self.tmp_XYZ = self.model.get_sites_cart()
        else:
          self.tmp_XYZ = self.model.get_sites_cart()
        if (self.macro_cycle != self.qblib_params.macro_cycle_to_skip):
          qblib_call = qb_refinement.QBblib_call_manager(
            hierarchy = self.model.get_hierarchy(),
            xray_structure=self.model.get_xray_structure(),
            geometry_residuals = self.stereochemistry_residuals,
            qblib_params=self.qblib_params,
            diff_curr=self.XYZ_diff_curr,
            macro=self.macro_cycle,
            micro=self.qblib_cycle_count,
            )
          try:
            qblib_call.run()
          except Exception as e:
            if e.__class__.__name__=="QB_none_fatal_error":
              raise Sorry(str(e))
            else:
              raise e

          if(qblib_call.QBStatus):
            qblib_g=qblib_call.result_QBlib.gradients
            qblib_f=qblib_call.result_QBlib.target
            self.stereochemistry_residuals.gradients=qblib_g
            self.stereochemistry_residuals.target=qblib_f
# QBLIB END
      er = self.stereochemistry_residuals.target
      self.f += er * self.weights.w
      if(compute_gradients):
        sgc = self.stereochemistry_residuals.gradients
        # ias do not participate in geometry restraints
        if self.model is not None and self.model.ias_manager is not None:
          ias_selection = self.model.ias_manager.get_ias_selection()
          if ias_selection is not None and ias_selection.count(True) > 0:
            sgc.extend(flex.vec3_double(ias_selection.count(True),[0,0,0]))
        xray.minimization.add_gradients(
          scatterers     = self.xray_structure.scatterers(),
          xray_gradients = self.g,
          site_gradients = sgc*self.weights.w)
    if(self.refine_adp and self.restraints_manager is not None and
       self.restraints_manager.geometry is not None
       and self.weights.w > 0.0 and self.iso_restraints is not None):
      use_hd = False
      if(self.fmodels.fmodel_n is not None or
         (self.is_neutron_scat_table is not None and
          self.is_neutron_scat_table == "neutron") or
         self.h_params.refine == "individual"):
        use_hd = True
      energies_adp = self.model.energies_adp(
        iso_restraints    = self.iso_restraints,
        use_hd            = use_hd,
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
      print("refinement.minimization step: f,iter,nfun:", end=' ')
      print(self.f,minimizer.iter(),minimizer.nfun())

  def compute_functional_and_gradients(self):
    u_iso_refinable_params = self.apply_shifts()
    self.compute_target(compute_gradients     = True,
                        u_iso_refinable_params = u_iso_refinable_params)
    if (self.verbose > 1):
      print("xray.minimization line search: f,rms(g):", end=' ')
      print(self.f, math.sqrt(flex.mean_sq(self.g)))
    return self.f, self.g

class monitor(object):
  def __init__(self, weights,
                     fmodels,
                     model,
                     iso_restraints,
                     refine_xyz = False,
                     refine_adp = False,
                     refine_occ = False):
    adopt_init_args(self, locals())
    self.ex = []
    self.en = []
    self.er = []
    self.et = []
    self.rxw = []
    self.rxf = []
    self.rnw = []
    self.rnf = []
    self.iter = None
    self.nfun = None

  def collect(self, iter = None, nfun = None):
    if(iter is not None): self.iter = format_value("%-4d", iter)
    if(nfun is not None): self.nfun = format_value("%-4d", nfun)
    self.fmodels.fmodel_xray().xray_structure == self.model.get_xray_structure()
    fmodels_tg = self.fmodels.target_and_gradients(
      weights           = self.weights,
      compute_gradients = False)
    self.rxw.append(format_value("%6.4f", self.fmodels.fmodel_xray().r_work()))
    self.rxf.append(format_value("%6.4f", self.fmodels.fmodel_xray().r_free()))
    self.ex.append(format_value("%10.4f", fmodels_tg.target_work_xray))
    if(self.fmodels.fmodel_neutron() is not None):
      self.rnw.append(format_value("%6.4f", self.fmodels.fmodel_neutron().r_work()))
      self.rnf.append(format_value("%6.4f", self.fmodels.fmodel_neutron().r_free()))
      self.en.append(format_value("%10.4f", fmodels_tg.target_work_neutron))
    if(self.refine_xyz and self.weights.w > 0):
      er = self.model.restraints_manager_energies_sites(
        compute_gradients = False).target
    elif(self.refine_adp and self.weights.w > 0):
      use_hd = False
      if(self.fmodels.fmodel_n is not None): use_hd = True
      er = self.model.energies_adp(
        iso_restraints = self.iso_restraints,
        use_hd = use_hd,
        compute_gradients = False).target
    elif(self.refine_occ):
      er = 0
    else: er = 0
    self.er.append(format_value("%10.4f", er))
    self.et.append(format_value(
      "%10.4f", fmodels_tg.target()+er*self.weights.w))

  def show(self, message = "", log = None):
    if(log is None): log = sys.stdout
    print("|-"+message+"-"*(79-len("| "+message+"|"))+"|", file=log)
    if(self.fmodels.fmodel_neutron() is not None):
      print("|"+" "*33+"x-ray data:"+" "*33+"|", file=log)
    print("| start r-factor (work) = %s      final r-factor "\
     "(work) = %s          |"%(self.rxw[0], self.rxw[1]), file=log)
    print("| start r-factor (free) = %s      final r-factor "\
     "(free) = %s          |"%(self.rxf[0], self.rxf[1]), file=log)
    if(self.fmodels.fmodel_neutron() is not None):
      print("|"+" "*32+"neutron data:"+" "*32+"|", file=log)
      print("| start r-factor (work) = %s      final r-factor "\
      "(work) = %s          |"%(self.rnw[0], self.rnw[1]), file=log)
      print("| start r-factor (free) = %s      final r-factor "\
      "(free) = %s          |"%(self.rnf[0], self.rnf[1]), file=log)
    print("|"+"-"*77+"|", file=log)
    #
    if(self.fmodels.fmodel_neutron() is None):
      for i_seq, stage in enumerate(["start", "final"]):
        if(self.refine_xyz):
          print("| T_%s = wxc * wxc_scale * Exray + wc * Echem"%stage+" "*30+"|", file=log)
          print(self.xray_line(i_seq), file=log)
        elif(self.refine_adp):
          print("| T_%s = wxu * wxu_scale * Exray + wu * Eadp"%stage+" "*31+"|", file=log)
          print(self.xray_line(i_seq), file=log)
        elif(self.refine_occ):
          print("| T_%s = 1.0 * 1.0 * Exray"%stage+" "*43+"|", file=log)
          print(self.xray_line(i_seq), file=log)
        if(i_seq == 0): print("|"+" "*77+"|", file=log)
    #
    if(self.fmodels.fmodel_neutron() is not None):
      for i_seq, stage in enumerate(["start", "final"]):
        if(self.refine_xyz):
          print("| T_%s = wnxc * (wxc_scale * Exray + wnc_scale * Eneutron) + wc * Echem"%stage+" "*4+"|", file=log)
          print(self.neutron_line(i_seq), file=log)
        elif(self.refine_adp):
          print("| T_%s = wnxu * (wxu_scale * Exray + wnu_scale * Eneutron) + wu * Eadp"%stage+" "*5+"|", file=log)
          print(self.neutron_line(i_seq), file=log)
        elif(self.refine_occ):
          print("| T_%s = 1.0 * (1.0 * Exray + 1.0 * Eneutron)"%stage+" "*30+"|", file=log)
          print(self.neutron_line(i_seq), file=log)
        if(i_seq == 0): print("|"+" "*77+"|", file=log)
    #
    print("|"+"-"*77+"|", file=log)
    print("| number of iterations = %s    |    number of function "\
                  "evaluations = %s   |"%(self.iter, self.nfun), file=log)
    print("|"+"-"*77+"|", file=log)
    log.flush()

  def neutron_line(self, i_seq):
    line = "| %s = %s * (%s * %s + %s * %s) + %s * %s"%(
      self.et[i_seq].strip(),
      format_value("%6.2f",self.weights.wxn).strip(),
      format_value("%6.2f",self.weights.wx_scale).strip(),
      self.ex[i_seq].strip(),
      format_value("%6.2f",self.weights.wn_scale).strip(),
      self.en[i_seq].strip(),
      format_value("%6.2f",self.weights.w).strip(),
      self.er[i_seq].strip())
    line = line + " "*(78 - len(line))+"|"
    return line

  def xray_line(self, i_seq):
    line = "| %s = %s * %s * %s + %s * %s"%(
      self.et[i_seq].strip(),
      format_value("%6.2f",self.weights.wx).strip(),
      format_value("%6.2f",self.weights.wx_scale).strip(),
      self.ex[i_seq].strip(),
      format_value("%6.2f",self.weights.w).strip(),
      self.er[i_seq].strip())
    line = line + " "*(78 - len(line))+"|"
    return line

class run_constrained(object):
  def __init__(self,
               model,
               fmodel,
               target_weight,
               log,
               params,
               prefix,
               refine_sites           = False,
               refine_u_iso           = False,
               refine_transformations = False):
    #
    # XXX Same needs to be done for refine_u_iso
    #
    refine_selection = None
    if(refine_sites):
      refine_selection = model.refinement_flags.sites_individual.iselection()
    #
    tg_object = mmtbx.refinement.minimization_ncs_constraints.\
      target_function_and_grads_reciprocal_space(
        fmodel                    = fmodel,
        ncs_restraints_group_list = model.get_ncs_groups(),
        refine_selection          = refine_selection, # XXX
        restraints_manager        = model.restraints_manager,
        iso_restraints            = params.adp_restraints.iso,
        data_weight               = target_weight,
        refine_sites              = refine_sites,
        refine_u_iso              = refine_u_iso,
        refine_transformations    = refine_transformations)
    self.minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
       target_and_grads_object      = tg_object,
       xray_structure               = fmodel.xray_structure,
       ncs_restraints_group_list    = model.get_ncs_groups(),
       refine_selection             = refine_selection, # XXX
       finite_grad_differences_test = False,
       max_iterations               = params.main.max_number_of_iterations,
       refine_sites                 = refine_sites,
        refine_u_iso                = refine_u_iso,
        refine_transformations      = refine_transformations)
    fmodel.update_xray_structure(update_f_calc=True)
    model.set_xray_structure(fmodel.xray_structure)
