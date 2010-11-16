from cctbx.array_family import flex
import iotbx.phil
from mmtbx.refinement import minimization
import mmtbx.refinement.group
from mmtbx.tls import tools
from mmtbx.refinement import print_statistics
import scitbx.lbfgs
from libtbx.test_utils import approx_equal
from libtbx.utils import user_plus_sys_time

time_adp_refinement_py = 0.0

def show_times(out = None):
  if(out is None): out = sys.stdout
  total = time_adp_refinement_py
  if(total > 0.01):
     print >> out, "ADP refinement:"
     print >> out, "  time spent in adp_refinement.py          = %-7.2f" % time_adp_refinement_py
  return total

group_adp_master_params = iotbx.phil.parse("""\
  number_of_macro_cycles   = 3
    .type = int
  max_number_of_iterations = 25
    .type = int
  convergence_test         = False
    .type = bool
  run_finite_differences_test = False
    .type = bool
""")

tls_master_params = iotbx.phil.parse("""\
  one_residue_one_group       = None
    .type = bool
    .style = tribool
  refine_T                    = True
    .type = bool
  refine_L                    = True
    .type = bool
  refine_S                    = True
    .type = bool
  number_of_macro_cycles      = 2
    .type = int
  max_number_of_iterations    = 25
    .type = int
  start_tls_value             = None
    .type = float
  run_finite_differences_test = False
    .type = bool
  eps                         = 1.e-6
    .type = float
""")

individual_adp_master_params = iotbx.phil.parse("""\
  iso {
    max_number_of_iterations = 25
      .type = int
    automatic_randomization_if_all_equal = True
      .type = bool
    scaling {
      scale_max       = 3.0
        .type = float
      scale_min       = 10.0
        .type = float
    }
  }
""")

adp_restraints_master_params = iotbx.phil.parse("""\
  iso {
    use_u_local_only = False
      .type = bool
    sphere_radius = 5.0
      .type = float
    distance_power = 1.69
      .type = float
    average_power = 1.03
      .type = float
    wilson_b_weight_auto = False
      .type = bool
    wilson_b_weight = None
      .type = float
    plain_pairs_radius = 5.0
      .type = float
    refine_ap_and_dp = False
      .type = bool
  }
""")

class manager(object):
  def __init__(
            self,
            fmodels,
            model,
            all_params,
            group_adp_selections   = None,
            group_adp_selections_h = None,
            group_adp_params       = group_adp_master_params.extract(),
            tls_selections         = None,
            tls_params             = tls_master_params.extract(),
            individual_adp_params  = individual_adp_master_params.extract(),
            adp_restraints_params  = adp_restraints_master_params.extract(),
            refine_adp_individual  = None,
            refine_adp_group       = None,
            refine_tls             = None,
            tan_b_iso_max          = None,
            restraints_manager     = None,
            target_weights         = None,
            macro_cycle            = None,
            log                    = None,
            h_params               = None):
    global time_adp_refinement_py
    scatterers = fmodels.fmodel_xray().xray_structure.scatterers()
    timer = user_plus_sys_time()
    if(log is None): log = sys.stdout
    tan_u_iso = False
    param = 0
    if(tan_b_iso_max > 0.0):
       tan_u_iso = True
       param = int(tan_b_iso_max)
    if(macro_cycle == 1):
       offset = True
    else:
       offset = False

    if(refine_tls):
       print_statistics.make_sub_header(text = "TLS refinement",
                                        out  = log)
       tls_sel_st = flex.size_t()
       for ts in tls_selections:
         tls_sel_st.extend(ts)
       tls_sel_bool = flex.bool(scatterers.size(), flex.size_t(tls_sel_st))
       ### totally ad hoc fix
       tmp_site_t = flex.size_t()
       for gs in group_adp_selections:
         for gs_ in gs:
           tmp_site_t.append(gs_)
       ###
       if(macro_cycle == 1 or tmp_site_t.size() != scatterers.size()):
          gbr_selections = []
          for s in tls_selections:
            gbr_selections.append(s)
       else:
          gbr_selections = []
          for gs in group_adp_selections:
            gbr_selection = flex.size_t()
            for gs_ in gs:
              if(tls_sel_bool[gs_]):
                gbr_selection.append(gs_)
            if(gbr_selection.size() > 0):
              gbr_selections.append(gbr_selection)
       gbr_selections_one_arr = flex.size_t()
       for gbs in gbr_selections:
         gbr_selections_one_arr.extend(gbs)
       scatterers = fmodels.fmodel_xray().xray_structure.scatterers()
       for gbr_selection in gbr_selections_one_arr:
         scatterers[gbr_selection].flags.set_use_u_iso(True)
       group_b_manager = mmtbx.refinement.group.manager(
          fmodel                   = fmodels.fmodel_xray(),
          selections               = gbr_selections,
          convergence_test         = group_adp_params.convergence_test,
          max_number_of_iterations = 50,
          number_of_macro_cycles   = 1,
          refine_adp               = True,
          log                      = log)
       scatterers = fmodels.fmodel_xray().xray_structure.scatterers()
       for tls_selection_ in tls_selections:
         for tls_selection__ in tls_selection_:
           scatterers[tls_selection__].flags.set_use_u_aniso(True)
       model.show_groups(tls = True, out = log)
       current_target_name = fmodels.fmodel_xray().target_name
       fmodels.fmodel_xray().update(target_name = "ls_wunit_k1")
       tools.split_u(fmodels.fmodel_xray().xray_structure, tls_selections, offset)
       self.tls_refinement_manager = tools.tls_refinement(
          fmodel                      = fmodels.fmodel_xray(),
          model                       = model,
          selections                  = tls_selections,
          selections_1d               = tls_sel_st,
          refine_T                    = tls_params.refine_T,
          refine_L                    = tls_params.refine_L,
          refine_S                    = tls_params.refine_S,
          number_of_macro_cycles      = tls_params.number_of_macro_cycles,
          max_number_of_iterations    = tls_params.max_number_of_iterations,
          start_tls_value             = tls_params.start_tls_value,
          run_finite_differences_test = tls_params.run_finite_differences_test,
          eps                         = tls_params.eps,
          out                         = log,
          macro_cycle = macro_cycle)
       fmodels.fmodel_xray().update(target_name = current_target_name)
       fmodels.update_xray_structure(
            xray_structure = self.tls_refinement_manager.fmodel.xray_structure,
            update_f_calc  = True)
       model.xray_structure = fmodels.fmodel_xray().xray_structure

    if(refine_adp_individual):
       refine_adp(model, fmodels, target_weights, individual_adp_params, adp_restraints_params, h_params, log,
       all_params = all_params)

    if(refine_adp_group):
       print_statistics.make_sub_header(text= "group isotropic ADP refinement",
                                        out = log)
       group_b_manager = mmtbx.refinement.group.manager(
          fmodel                   = fmodels.fmodel_xray(),
          selections               = group_adp_selections,
          convergence_test         = group_adp_params.convergence_test,
          max_number_of_iterations = group_adp_params.max_number_of_iterations,
          number_of_macro_cycles   = group_adp_params.number_of_macro_cycles,
          run_finite_differences_test = group_adp_params.run_finite_differences_test,
          refine_adp               = True,
          log                      = log)
    time_adp_refinement_py += timer.elapsed()

def refine_adp(model,
               fmodels,
               target_weights,
               individual_adp_params,
               adp_restraints_params,
               h_params,
               log,
               all_params):
  model.set_refine_individual_adp()
  print_statistics.make_sub_header(text="Individual ADP refinement", out = log)
  lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
    max_iterations = individual_adp_params.iso.max_number_of_iterations)
  assert fmodels.fmodel_xray().xray_structure is model.xray_structure
  if(target_weights is not None and target_weights.twp.optimize_wxu):
    save_scatterers = fmodels.fmodel_xray().xray_structure.\
      deep_copy_scatterers().scatterers()
    print >> log, "Start r_free = %6.4f"%fmodels.fmodel_xray().r_free()
    new_scatterers = None
    scaler_values = [1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,6.,7.,8.,9.,10.]
  else: scaler_values = [1.,]
  r_free = 1.e+9
  cntr = 0
  if(target_weights is not None):
    w = target_weights.adp_weights_result.w
  for scaler in scaler_values:
    if(target_weights is not None and target_weights.twp.optimize_wxu):
      fmodels.fmodel_xray().xray_structure.replace_scatterers(
        save_scatterers.deep_copy())
      fmodels.update_xray_structure(update_f_calc = True)
      target_weights.adp_weights_result.w = w / scaler
    is_neutron_scat_table = False
    if(all_params.main.scattering_table == "neutron"):
      is_neutron_scat_table = True
    minimized = minimization.lbfgs(
      restraints_manager       = model.restraints_manager,
      fmodels                  = fmodels,
      model                    = model,
      refine_adp               = True,
      is_neutron_scat_table    = is_neutron_scat_table,
      lbfgs_termination_params = lbfgs_termination_params,
      iso_restraints           = adp_restraints_params.iso,
      verbose                  = 0,
      target_weights           = target_weights,
      h_params                 = h_params)
    r_free_ = fmodels.fmodel_xray().r_free()
    if(target_weights is not None and target_weights.twp.optimize_wxu):
      print >> log, "scale= %8.4f total_weight= %8.4f r_free= %6.4f"%(scaler,
        target_weights.adp_weights_result.wx, r_free_), cntr
    if(r_free_ < r_free):
      r_free = r_free_
      new_scatterers = fmodels.fmodel_xray().xray_structure.scatterers(
        ).deep_copy()
      cntr = 0
    elif(r_free != 1.e+9):
      cntr += 1
    if(cntr == 3): break
    assert fmodels.fmodel_xray().xray_structure is model.xray_structure
    assert minimized.xray_structure is model.xray_structure
  #########
  if(target_weights is not None and target_weights.twp.optimize_wxu):
    scaler_values = [1./1.5,1./2.,1./2.5,1./3.,1./3.5,1./4.,1./4.5,1./5.]
    cntr = 0
    for scaler in scaler_values:
      if(target_weights is not None and target_weights.twp.optimize_wxu):
        fmodels.fmodel_xray().xray_structure.replace_scatterers(
          save_scatterers.deep_copy())
        fmodels.update_xray_structure(update_f_calc = True)
        target_weights.adp_weights_result.w = w / scaler
      is_neutron_scat_table = False
      if(all_params.main.scattering_table == "neutron"):
        is_neutron_scat_table = True
      minimized = minimization.lbfgs(
        restraints_manager       = model.restraints_manager,
        fmodels                  = fmodels,
        model                    = model,
        refine_adp               = True,
        lbfgs_termination_params = lbfgs_termination_params,
        iso_restraints           = adp_restraints_params.iso,
        verbose                  = 0,
        target_weights           = target_weights,
        is_neutron_scat_table    = is_neutron_scat_table,
        h_params                 = h_params)
      r_free_ = fmodels.fmodel_xray().r_free()
      if(target_weights.twp.optimize_wxu):
        print >> log, "scale= %8.4f total_weight= %8.4f r_free= %6.4f"%(scaler,
          target_weights.adp_weights_result.wx, r_free_), cntr
      if(r_free_ < r_free):
        r_free = r_free_
        new_scatterers = fmodels.fmodel_xray().xray_structure.scatterers(
          ).deep_copy()
        cntr = 0
      elif(r_free != 1.e+9):
        cntr += 1
      if(cntr == 3): break
      assert fmodels.fmodel_xray().xray_structure is model.xray_structure
      assert minimized.xray_structure is model.xray_structure
  #########
  if(target_weights is not None and target_weights.twp.optimize_wxu):
    print >> log, "Best r_free =  %6.4f"%r_free
    fmodels.fmodel_xray().xray_structure.replace_scatterers(new_scatterers)
    fmodels.update_xray_structure(update_f_calc = True)
    assert approx_equal(fmodels.fmodel_xray().r_free(), r_free)
  else: fmodels.update_xray_structure(update_f_calc = True)
  if(target_weights is not None and not target_weights.twp.optimize_wxu):
    minimized.monitor.show(message = "LBFGS minimization", log  = log)
    #monitors.collect(step = str(macro_cycle)+"_adp:") XXX
