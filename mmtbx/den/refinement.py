from mmtbx import utils
from mmtbx.dynamics import simulated_annealing
from mmtbx.refinement import tardy
from libtbx import easy_mp
from mmtbx.refinement import print_statistics
from mmtbx.refinement import adp_refinement
import sys

class manager(object):
  def __init__(
            self,
            fmodels,
            model,
            params,
            target_weights,
            monitors,
            macro_cycle,
            log=None):
    if log is None:
      log = sys.stdout
    self.log = log
    self.fmodels = fmodels
    self.model = model
    self.params = params
    self.target_weights = target_weights
    self.monitors = monitors
    self.macro_cycle = macro_cycle
    self.tan_b_iso_max = 0
    self.save_scatterers_local = fmodels.fmodel_xray().\
        xray_structure.deep_copy_scatterers().scatterers()
    den_manager = model.restraints_manager. \
      geometry.generic_restraints_manager.den_manager
    print_statistics.make_header("DEN refinement", out=self.log)
    if params.den.optimize:
      grid = den_manager.get_optimization_grid()
      print >> log, \
        "Running DEN torsion optimization on %d processors..." % \
        params.main.nproc
    else:
      grid = [(params.den.gamma, params.den.weight)]

    if "torsion" in params.den.annealing_type:
      print >> self.log, "Running torsion simulated annealing"
      stdout_and_results = easy_mp.pool_map(
        processes=params.main.nproc,
        fixed_func=self.try_den_weight_torsion,
        args=grid,
        buffer_stdout_stderr=True)
      if params.den.optimize:
        self.show_den_opt_summary_torsion(stdout_and_results)
    elif "cartesian" in params.den.annealing_type:
      print >> self.log, "Running Cartesian simulated annealing"
      stdout_and_results = easy_mp.pool_map(
        processes=params.main.nproc,
        fixed_func=self.try_den_weight_cartesian,
        args=grid,
        buffer_stdout_stderr=True)
      if params.den.optimize:
        self.show_den_opt_summary_cartesian(stdout_and_results)
    else:
      raise "error in DEN annealing type"
    low_r_free = 1.0
    best_xray_structure = None
    best_gamma =  None
    best_weight = None
    best_out = None
    for out, result in stdout_and_results:
      cur_r_free = result[2]
      if cur_r_free < low_r_free:
        low_r_free = cur_r_free
        best_gamma = result[0]
        best_weight = result[1]
        best_xray_structure = result[3]
        best_out = out
    assert best_xray_structure is not None
    if params.den.optimize:
      print >> self.log, "\nbest gamma: %.1f" % best_gamma
      print >> self.log, "best weight: %.1f" % best_weight
      if params.den.verbose:
        print >> self.log, "\nBest annealing results:\n"
        print >> self.log, best_out
    fmodels.fmodel_xray().xray_structure.replace_scatterers(
      best_xray_structure.deep_copy())
    fmodels.update_xray_structure(
      xray_structure = fmodels.fmodel_xray().xray_structure,
      update_f_calc  = True)
    utils.assert_xray_structures_equal(
      x1 = fmodels.fmodel_xray().xray_structure,
      x2 = model.xray_structure)

  def try_den_weight_torsion(self, grid_pair):
    self.fmodels.fmodel_xray().xray_structure.replace_scatterers(
      self.save_scatterers_local.deep_copy())
    self.fmodels.update_xray_structure(
      xray_structure = self.fmodels.fmodel_xray().xray_structure,
      update_f_calc  = True)
    utils.assert_xray_structures_equal(
      x1 = self.fmodels.fmodel_xray().xray_structure,
      x2 = self.model.xray_structure)
    gamma_local = grid_pair[0]
    weight_local = grid_pair[1]
    self.model.restraints_manager.geometry.\
      generic_restraints_manager.den_manager.gamma = \
        gamma_local
    self.model.restraints_manager.geometry.\
      generic_restraints_manager.den_manager.weight = \
        weight_local
    cycle = 0
    num_den_cycles = self.model.restraints_manager.geometry.\
      generic_restraints_manager.den_manager.num_cycles
    if self.params.den.optimize:
      local_log = sys.stdout
    else:
      local_log = self.log
    print >> self.log, "  ...trying gamma %.1f, weight %.1f" % (
      gamma_local, weight_local)
    while cycle < num_den_cycles:
      print >> local_log, "DEN cycle %s" % (cycle+1)
      r_free = self.fmodels.fmodel_xray().r_free()
      print >> local_log, "rfree at start of SA cycle: %f" % r_free
      tardy.run(
        fmodels=self.fmodels,
        model=self.model,
        target_weights=self.target_weights,
        params=self.params.tardy,
        log=local_log,
        format_for_phenix_refine=True,
        monitor=self.monitors.monitor_xray)
      if self.params.den.refine_adp:
        self.adp_refinement(log=local_log)
      cycle += 1
      r_free = self.fmodels.fmodel_xray().r_free()
      print >> local_log, "rfree at end of SA cycle: %f" % r_free
    r_free = self.fmodels.fmodel_xray().r_free()
    step_xray_structure = self.fmodels.fmodel_xray().\
      xray_structure.deep_copy_scatterers().scatterers()
    return (gamma_local,
            weight_local,
            r_free,
            step_xray_structure)

  def try_den_weight_cartesian(self, grid_pair):
    self.fmodels.fmodel_xray().xray_structure.replace_scatterers(
      self.save_scatterers_local.deep_copy())
    self.fmodels.update_xray_structure(
      xray_structure = self.fmodels.fmodel_xray().xray_structure,
      update_f_calc  = True)
    utils.assert_xray_structures_equal(
      x1 = self.fmodels.fmodel_xray().xray_structure,
      x2 = self.model.xray_structure)
    gamma_local = grid_pair[0]
    weight_local = grid_pair[1]
    self.model.restraints_manager.geometry.\
      generic_restraints_manager.den_manager.gamma = \
        gamma_local
    self.model.restraints_manager.geometry.\
      generic_restraints_manager.den_manager.weight = \
        weight_local
    cycle = 0
    num_den_cycles = self.model.restraints_manager.geometry.\
      generic_restraints_manager.den_manager.num_cycles
    if self.params.den.optimize:
      local_log = sys.stdout
    else:
      local_log = self.log
    print >> self.log, "  ...trying gamma %f, weight %f" % (
      gamma_local, weight_local)
    while cycle < num_den_cycles:
      print >> local_log, "DEN cycle %s" % (cycle+1)
      r_free = self.fmodels.fmodel_xray().r_free()
      print >> local_log, "rfree at start of SA cycle: %f" % r_free
      simulated_annealing.manager(
        simulated_annealing_params = self.params.simulated_annealing,
        bulk_solvent_parameters    = self.params.bulk_solvent_and_scale,
        refinement_parameters      = self.params.main,
        alpha_beta_parameters      = self.params.alpha_beta,
        mask_parameters            = self.params.mask,
        target_weights             = self.target_weights,
        tan_b_iso_max              = self.tan_b_iso_max,
        macro_cycle                = self.macro_cycle,
        h_params                   = self.params.hydrogens,
        fmodels                    = self.fmodels,
        model                      = self.model,
        all_params                 = self.params,
        out                        = local_log,
        monitor                    = self.monitors.monitor_xray)
      if self.params.den.refine_adp:
        self.adp_refinement(log=local_log)
      cycle += 1
      r_free = self.fmodels.fmodel_xray().r_free()
      print >> local_log, "rfree at end of SA cycle: %f" % r_free
    r_free = self.fmodels.fmodel_xray().r_free()
    step_xray_structure = self.fmodels.fmodel_xray().\
      xray_structure.deep_copy_scatterers().scatterers()
    return (gamma_local,
            weight_local,
            r_free,
            step_xray_structure)

  def show_den_opt_summary_torsion(self, stdout_and_results):
    print_statistics.make_header(
      "DEN torsion weight optimization results", out=self.log)
    print >>self.log,"|---------------------------------------"+\
                "--------------------------------------|"
    print >>self.log,"|  Gamma    Weight    R-free            "+\
                "                                      |"
    for out, result in stdout_and_results:
      if result == None:
        raise RuntimeError("Parallel DEN job failed: %s" % str(out))
      cur_gamma = result[0]
      cur_weight = result[1]
      cur_r_free = result[2]
      print >> self.log, "| %6.1f    %6.1f    %6.4f              " %\
        (cur_gamma,
         cur_weight,
         cur_r_free)+\
                    "                                    |"
    print >>self.log,"|---------------------------------------"+\
                "--------------------------------------|"

  def show_den_opt_summary_cartesian(self, stdout_and_results):
    print_statistics.make_header(
      "DEN Cartesian weight optimization results", out=self.log)
    print >>self.log,"|---------------------------------------"+\
                "--------------------------------------|"
    print >>self.log,"|  Gamma    Weight    R-free            "+\
                "                                      |"
    for out, result in stdout_and_results:
      if result == None:
        raise RuntimeError("Parallel DEN job failed: %s" % str(out))
      cur_gamma = result[0]
      cur_weight = result[1]
      cur_r_free = result[2]
      print >> self.log, "| %6.1f    %6.1f    %6.4f              " %\
        (cur_gamma,
         cur_weight,
         cur_r_free)+\
                    "                                    |"
    print >>self.log,"|---------------------------------------"+\
                "--------------------------------------|"

  def adp_refinement(self, log):
    if log is None:
      log = sys.stdout
    save_xray_structure = self.fmodels.fmodel_xray().\
      xray_structure.deep_copy_scatterers().scatterers()
      ###> Make ADP of H/D sites equal
    self.model.reset_adp_of_hd_sites_to_be_equal()
    self.fmodels.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = True)
    adp_refinement_manager = adp_refinement.manager(
      fmodels                = self.fmodels,
      model                  = self.model,
      group_adp_selections   = self.model.refinement_flags.adp_group,
      group_adp_selections_h = self.model.refinement_flags.group_h,
      group_adp_params       = self.params.group_b_iso,
      tls_selections         = self.model.refinement_flags.adp_tls,
      all_params             = self.params,
      tls_params             = self.params.tls,
      individual_adp_params  = self.params.adp,
      adp_restraints_params  = self.params.adp_restraints,
      refine_adp_individual  = self.model.refinement_flags.individual_adp,
      refine_adp_group       = self.model.refinement_flags.group_adp,
      refine_tls             = self.model.refinement_flags.tls,
      tan_b_iso_max          = self.tan_b_iso_max,
      restraints_manager     = self.model.restraints_manager,
      macro_cycle            = self.macro_cycle,
      target_weights         = self.target_weights,
      log                    = log,
      h_params               = self.params.hydrogens)
    #self.monitors.collect(step    = str(self.macro_cycle)+"_adp:",
    #                 model   = self.model,
    #                 fmodels = self.fmodels)
