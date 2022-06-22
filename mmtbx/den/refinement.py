from __future__ import absolute_import, division, print_function
from mmtbx import utils
from mmtbx.dynamics import simulated_annealing
from mmtbx.refinement import tardy
from libtbx import easy_mp, Auto
from mmtbx.refinement import print_statistics
from mmtbx.refinement import adp_refinement
from cctbx.array_family import flex
from six.moves import cStringIO as StringIO
import sys, random
from mmtbx.refinement import geometry_minimization

class manager(object):
  def __init__(
            self,
            fmodels,
            model,
            params,
            target_weights,
            macro_cycle,
            ncs_manager=None,
            log=None):
    if log is None:
      log = sys.stdout
    # self.ncs_manager = ncs_manager
    self.nproc = params.main.nproc
    if self.nproc is Auto:
      self.nproc = 1
    self.verbose = params.den.verbose
    self.log = log
    self.fmodels = fmodels
    self.model = model
    self.params = params
    self.target_weights = target_weights
    self.adp_refinement_manager = None
    self.macro_cycle = macro_cycle
    self.tan_b_iso_max = 0
    self.random_seed = params.main.random_seed
    den_manager = model.restraints_manager. \
      geometry.den_manager
    print_statistics.make_header("DEN refinement", out=self.log)
    pdb_hierarchy = self.model.get_hierarchy()
    if den_manager.get_n_proxies() == 0:
      print_statistics.make_sub_header(
        "DEN restraint nework", out = self.log)
      den_manager.build_den_proxies(pdb_hierarchy=pdb_hierarchy)
      den_manager.build_den_restraints()
      den_manager.show_den_summary(
        sites_cart=self.model.get_sites_cart())
    if den_manager.params.output_kinemage:
      den_manager.output_kinemage(
        self.model.get_sites_cart())
    print_statistics.make_sub_header(
      "coordinate minimization before annealing", out=self.log)
    self.minimize(ca_only=self.params.den.minimize_c_alpha_only)
    self.save_scatterers_local = fmodels.fmodel_xray().\
      xray_structure.deep_copy_scatterers().scatterers()
    #DEN refinement start, turn on
    if params.den.optimize:
      grid = den_manager.get_optimization_grid()
      print("Running DEN torsion optimization on %d processors..." % \
        self.nproc, file=log)
    else:
      grid = [(params.den.gamma, params.den.weight)]
    grid_results = []
    grid_so = []

    if "torsion" in params.den.annealing_type:
      print("Running torsion simulated annealing", file=self.log)
      if ( (params.den.optimize) and
           ( (self.nproc is Auto) or (self.nproc > 1) )):
        stdout_and_results = easy_mp.pool_map(
          processes=params.main.nproc,
          fixed_func=self.try_den_weight_torsion,
          args=grid,
          func_wrapper="buffer_stdout_stderr")
        for so, r in stdout_and_results:
          if (r is None):
            raise RuntimeError(("DEN weight optimization failed:"+
              "\n%s\nThis is a "+
              "serious error; please contact help@phenix-online.org.") % so)
          grid_so.append(so)
          grid_results.append(r)
      else:
        for grid_pair in grid:
          result = self.try_den_weight_torsion(
                     grid_pair=grid_pair)
          grid_results.append(result)
      self.show_den_opt_summary_torsion(grid_results)
    elif "cartesian" in params.den.annealing_type:
      print("Running Cartesian simulated annealing", file=self.log)
      if ( (params.den.optimize) and
           ( (self.nproc is Auto) or (self.nproc > 1) )):
        stdout_and_results = easy_mp.pool_map(
          processes=params.main.nproc,
          fixed_func=self.try_den_weight_cartesian,
          args=grid,
          func_wrapper="buffer_stdout_stderr")
        for so, r in stdout_and_results:
          if (r is None):
            raise RuntimeError(("DEN weight optimization failed:"+
              "\n%s\nThis is a "+
              "serious error; please contact help@phenix-online.org.") % so)
          grid_so.append(so)
          grid_results.append(r)
      else:
        for grid_pair in grid:
          result = self.try_den_weight_cartesian(
                     grid_pair=grid_pair)
          grid_results.append(result)
      self.show_den_opt_summary_cartesian(grid_results)
    else:
      raise "error in DEN annealing type"
    low_r_free = 1.0
    best_xray_structure = None
    best_eq_distances = None
    best_gamma =  None
    best_weight = None
    best_so_i = None
    for i, result in enumerate(grid_results):
      cur_r_free = result[2]
      if cur_r_free < low_r_free:
        low_r_free = cur_r_free
        best_gamma = result[0]
        best_weight = result[1]
        best_xray_structure = result[3]
        best_eq_distances = result[4]
        best_so_i = i
    assert best_xray_structure is not None
    if params.den.optimize:
      print("\nbest gamma: %.1f" % best_gamma, file=self.log)
      print("best weight: %.1f\n" % best_weight, file=self.log)
      if params.den.verbose:
        if len(grid_so) >= (best_so_i+1):
          print("\nBest annealing results:\n", file=self.log)
          print(grid_so[best_so_i], file=self.log)
    fmodels.fmodel_xray().xray_structure.replace_scatterers(
      best_xray_structure.deep_copy())
    fmodels.update_xray_structure(
      xray_structure = fmodels.fmodel_xray().xray_structure,
      update_f_calc  = True)
    utils.assert_xray_structures_equal(
      x1 = fmodels.fmodel_xray().xray_structure,
      x2 = model.get_xray_structure())
    model.restraints_manager.geometry.\
      den_manager.import_eq_distances(eq_distances=best_eq_distances)
    self.model.torsion_ncs_restraints_update(log=self.log)
    #DEN refinement done, turn off

  def try_den_weight_torsion(self, grid_pair):
    #backup_k_rep = self.params.tardy.\
    #  prolsq_repulsion_function_changes.k_rep
    local_seed = int(self.random_seed+grid_pair[1])
    flex.set_random_seed(value=local_seed)
    random.seed(local_seed)
    self.fmodels.fmodel_xray().xray_structure.replace_scatterers(
      self.save_scatterers_local.deep_copy())
    self.fmodels.update_xray_structure(
      xray_structure = self.fmodels.fmodel_xray().xray_structure,
      update_f_calc  = True)
    utils.assert_xray_structures_equal(
      x1 = self.fmodels.fmodel_xray().xray_structure,
      x2 = self.model.get_xray_structure())
    gamma_local = grid_pair[0]
    weight_local = grid_pair[1]
    self.model.restraints_manager.geometry.\
        den_manager.gamma = gamma_local
    self.model.restraints_manager.geometry.\
        den_manager.weight = weight_local
    cycle = 0
    self.model.restraints_manager.geometry.\
        den_manager.current_cycle = cycle+1
    num_den_cycles = self.model.restraints_manager.geometry.\
        den_manager.num_cycles
    if self.params.den.optimize and \
       self.nproc != Auto and \
       self.nproc > 1:
      local_log = sys.stdout
    elif self.params.den.optimize and \
         self.nproc == 1:
      if self.verbose:
        local_log = self.log
      else:
        local_log = StringIO()
    else:
      local_log = self.log
    print("  ...trying gamma %.1f, weight %.1f" % (
      gamma_local, weight_local), file=self.log)
    while cycle < num_den_cycles:
      #if self.model.restraints_manager.geometry.\
      #     generic_restraints_manager.den_manager.current_cycle == \
      #     self.model.restraints_manager.geometry.\
      #     generic_restraints_manager.den_manager.torsion_mid_point+1:
      #  self.params.tardy.\
      #    prolsq_repulsion_function_changes.k_rep = 1.0
      print("DEN cycle %d" % (cycle+1), file=local_log)
      #print >> local_log, "Random seed: %d" % flex.get_random_seed()
      r_free = self.fmodels.fmodel_xray().r_free()
      print("rfree at start of SA cycle: %.4f" % r_free, file=local_log)
      print("k_rep = %.2f" % \
        self.params.tardy.\
          prolsq_repulsion_function_changes.k_rep, file=local_log)
      tardy.run(
        fmodels=self.fmodels,
        model=self.model,
        target_weights=self.target_weights,
        params=self.params.tardy,
        log=local_log,
        format_for_phenix_refine=True,
        call_back_after_step=False)
      if self.params.den.bulk_solvent_and_scale:
        self.bulk_solvent_and_scale(log=local_log)
        self.fmodels.fmodel_xray().xray_structure = self.model.get_xray_structure()
      if self.params.den.refine_adp:
        self.adp_refinement(log=local_log)
      self.model.torsion_ncs_restraints_update(log=local_log)
      cycle += 1
      self.model.restraints_manager.geometry.\
          den_manager.current_cycle += 1
      r_free = self.fmodels.fmodel_xray().r_free()
      print("rfree at end of SA cycle: %f" % r_free, file=local_log)
    r_free = self.fmodels.fmodel_xray().r_free()
    step_xray_structure = self.fmodels.fmodel_xray().\
      xray_structure.deep_copy_scatterers().scatterers()
    step_eq_distances = self.model.restraints_manager.geometry.\
        den_manager.get_current_eq_distances()
    return (gamma_local,
            weight_local,
            r_free,
            step_xray_structure,
            step_eq_distances)

  def try_den_weight_cartesian(self, grid_pair):
    local_seed = int(self.random_seed+grid_pair[1])
    flex.set_random_seed(value=local_seed)
    random.seed(local_seed)
    self.fmodels.fmodel_xray().xray_structure.replace_scatterers(
      self.save_scatterers_local.deep_copy())
    self.fmodels.update_xray_structure(
      xray_structure = self.fmodels.fmodel_xray().xray_structure,
      update_f_calc  = True)
    utils.assert_xray_structures_equal(
      x1 = self.fmodels.fmodel_xray().xray_structure,
      x2 = self.model.get_xray_structure())
    gamma_local = grid_pair[0]
    weight_local = grid_pair[1]
    self.model.restraints_manager.geometry.\
        den_manager.gamma = gamma_local
    self.model.restraints_manager.geometry.\
        den_manager.weight = weight_local
    cycle = 0
    self.model.restraints_manager.geometry.\
        den_manager.current_cycle = cycle+1
    num_den_cycles = self.model.restraints_manager.geometry.\
        den_manager.num_cycles
    if self.params.den.optimize and \
       self.nproc != Auto and \
       self.nproc > 1:
      local_log = sys.stdout
    elif self.params.den.optimize and \
         self.nproc == 1:
      if self.verbose:
        local_log = self.log
      else:
        local_log = StringIO()
    else:
      local_log = self.log
    print("  ...trying gamma %f, weight %f" % (
      gamma_local, weight_local), file=self.log)
    while cycle < num_den_cycles:
      print("DEN cycle %s" % (cycle+1), file=local_log)
      r_free = self.fmodels.fmodel_xray().r_free()
      print("rfree at start of SA cycle: %f" % r_free, file=local_log)
      simulated_annealing.manager(
        params         = self.params.simulated_annealing,
        target_weights = self.target_weights,
        macro_cycle    = self.macro_cycle,
        h_params       = self.params.hydrogens,
        fmodels        = self.fmodels,
        model          = self.model,
        all_params     = self.params,
        out            = local_log)
      if self.params.den.bulk_solvent_and_scale:
        self.bulk_solvent_and_scale(log=local_log)
      if self.params.den.refine_adp:
        self.adp_refinement(log=local_log)
      self.model.torsion_ncs_restraints_update(log=local_log)
      cycle += 1
      self.model.restraints_manager.geometry.\
          den_manager.current_cycle += 1
      r_free = self.fmodels.fmodel_xray().r_free()
      print("rfree at end of SA cycle: %f" % r_free, file=local_log)
    r_free = self.fmodels.fmodel_xray().r_free()
    step_xray_structure = self.fmodels.fmodel_xray().\
      xray_structure.deep_copy_scatterers().scatterers()
    step_eq_distances = self.model.restraints_manager.geometry.\
        den_manager.get_current_eq_distances()
    return (gamma_local,
            weight_local,
            r_free,
            step_xray_structure,
            step_eq_distances)

  def show_den_opt_summary_torsion(self, grid_results):
    print_statistics.make_header(
      "DEN torsion weight optimization results", out=self.log)
    print("|---------------------------------------"+\
                "--------------------------------------|", file=self.log)
    print("|  Gamma    Weight    R-free            "+\
                "                                      |", file=self.log)
    for result in grid_results:
      if result == None:
        raise RuntimeError("Parallel DEN job failed: %s" % str(out))
      cur_gamma = result[0]
      cur_weight = result[1]
      cur_r_free = result[2]
      print("| %6.1f    %6.1f    %6.4f              " %\
        (cur_gamma,
         cur_weight,
         cur_r_free)+\
                    "                                    |", file=self.log)
    print("|---------------------------------------"+\
                "--------------------------------------|", file=self.log)

  def show_den_opt_summary_cartesian(self, grid_results):
    print_statistics.make_header(
      "DEN Cartesian weight optimization results", out=self.log)
    print("|---------------------------------------"+\
                "--------------------------------------|", file=self.log)
    print("|  Gamma    Weight    R-free            "+\
                "                                      |", file=self.log)
    for result in grid_results:
      if result == None:
        raise RuntimeError("Parallel DEN job failed: %s" % str(out))
      cur_gamma = result[0]
      cur_weight = result[1]
      cur_r_free = result[2]
      print("| %6.1f    %6.1f    %6.4f              " %\
        (cur_gamma,
         cur_weight,
         cur_r_free)+\
                    "                                    |", file=self.log)
    print("|---------------------------------------"+\
                "--------------------------------------|", file=self.log)

  def bulk_solvent_and_scale(self, log):
    self.fmodels.update_all_scales(
      update_f_part1 = True,
      params = self.params.bulk_solvent_and_scale,
      optimize_mask = self.params.main.optimize_mask,
      force_update_f_mask = True,
      nproc=1,
      log=log)

  def adp_refinement(self, log):
    if log is None:
      log = sys.stdout
    save_xray_structure = self.fmodels.fmodel_xray().\
      xray_structure.deep_copy_scatterers().scatterers()
      ###> Make ADP of H/D sites equal
    self.model.reset_adp_of_hd_sites_to_be_equal()
    self.fmodels.update_xray_structure(
      xray_structure = self.model.get_xray_structure(),
      update_f_calc  = True)
    self.adp_refinement_manager = adp_refinement.manager(
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
      h_params               = self.params.hydrogens,
      nproc                  = 1)

  def minimize(self, ca_only=False):
    sel = "all"
    restraint_sites_cart = self.model.get_sites_cart().deep_copy()
    if ca_only:
      restraint_sites_cart = self.model.get_sites_cart().\
          deep_copy().select(self.model.get_atom_selection_cache().selection("name CA"))
      sel = "name CA"
    self.model.set_reference_coordinate_restraints(
        ref_model=self.model,
        selection=sel,
        exclude_outliers=False)
    ##### sanity check #####
    # How this check was supposed to hold up when ca_only=True???
    # When ca_only=True, number of restraints will be less than number of atoms.
    # Good that users never change defaults (minimize_c_alpha_only = False).
    #
    # XXX finally, removing this assertion. On top of c_alpha_only case,
    # it also fails (rightfully) when there are water molecules in the model,
    # they are excluded. The same goes to hydrogens, but it seems no one yet
    # tried den_refine for model with hydrogens on.
    # n_rcp = self.model.get_restraints_manager().geometry.get_n_reference_coordinate_proxies()
    # n_sc = len(restraint_sites_cart)
    # assert n_rcp >= n_sc, "%d >= %d fail" % (n_rcp, n_sc)
    ########################

    # selection = self.model.selection_moving
    geometry_minimization.run2(
        restraints_manager = self.model.get_restraints_manager(),
        pdb_hierarchy = self.model.get_hierarchy(),
        max_number_of_iterations       = 500,
        correct_special_position_tolerance=1.0,
        number_of_macro_cycles         = 1,
        bond                           = True,
        nonbonded                      = True,
        angle                          = True,
        dihedral                       = True,
        chirality                      = True,
        planarity                      = True)
    self.model.set_sites_cart_from_hierarchy()

    utils.assert_xray_structures_equal(
      x1 = self.fmodels.fmodel_xray().xray_structure,
      x2 = self.model.get_xray_structure())
    self.model.restraints_manager.geometry.\
        remove_reference_coordinate_restraints_in_place()
    # self.model.restraints_manager.geometry.generic_restraints_manager.\
    #      flags.reference_coordinate = False
