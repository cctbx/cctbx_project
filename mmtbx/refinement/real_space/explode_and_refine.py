from __future__ import division
from __future__ import print_function
from  mmtbx.refinement import geometry_minimization
import mmtbx.refinement.real_space.individual_sites
from cctbx import maptbx
from libtbx import adopt_init_args
import random, time
from scitbx.array_family import flex
from mmtbx.dynamics import simulated_annealing as sa
from cctbx import geometry_restraints
import scitbx.lbfgs
from mmtbx.building.merge_models import run as merge_models
import sys

if (0): # fixed random seed to avoid rare failures
  random.seed(1)
  flex.set_random_seed(1)

class run_sa(object):
  def __init__(
        self,
        xray_structure, # XXX redundant
        pdb_hierarchy,
        restraints_manager,
        map_data,
        number_of_trials,
        nproc,
        weight):
    adopt_init_args(self, locals())
    # Initialize states collector
    self.states = mmtbx.utils.states(
      xray_structure = self.xray_structure.deep_copy_scatterers(),
      pdb_hierarchy  = self.pdb_hierarchy.deep_copy())
    # SA params
    self.params = sa.master_params().extract()
    self.params.start_temperature=50000
    self.params.final_temperature=0
    self.params.cool_rate = 25000
    self.params.number_of_steps = 50
    # minimizer params
    self.grf = geometry_restraints.flags.flags(default=True)
    self.lbfgs_exception_handling_params = \
      scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = True,
        ignore_line_search_failed_step_at_upper_bound = True,
        ignore_line_search_failed_maxfev              = True)
    # pre-compute random seeds
    random_seeds = []
    for it in xrange(self.number_of_trials):
      random_seeds.append(random.randint(0,10000000))
    # run SA
    self.results = []
    if(self.nproc>1):
      from libtbx import easy_mp
      stdout_and_results = easy_mp.pool_map(
        processes    = self.nproc,
        fixed_func   = self.run,
        args         = random_seeds,
        func_wrapper = "buffer_stdout_stderr")
      for so, xrs in stdout_and_results :
        self.results.append(xrs)
        self.states.add(sites_cart = xrs.sites_cart())
    else:
      for random_seed in random_seeds:
        xrs = self.run(random_seed=random_seed).deep_copy_scatterers()
        self.results.append(xrs)
        self.states.add(sites_cart = xrs.sites_cart())
    assert len(self.results) == self.number_of_trials

  def ensemble_xrs(self):
    return self.results

  def show_states(self, file_name="SA_ensemble.pdb"):
    self.states.write(file_name=file_name)

  def run(self, random_seed):
    tmp = self.xray_structure.deep_copy_scatterers()
    # Shake and minimize to get variable starting points
    tmp.shake_sites_in_place(
      rms_difference = None,
      mean_distance  = 5)
    sites_cart_ = tmp.sites_cart()
    minimized = geometry_minimization.lbfgs(
      sites_cart                        = sites_cart_,
      correct_special_position_tolerance= 1.0,
      geometry_restraints_manager       = self.restraints_manager.geometry,
      geometry_restraints_flags         = self.grf,
      lbfgs_exception_handling_params   = self.lbfgs_exception_handling_params,
      lbfgs_termination_params          = scitbx.lbfgs.termination_parameters(
        max_iterations=50))
    tmp = tmp.replace_sites_cart(new_sites = sites_cart_)
    #
    random.seed(random_seed)
    flex.set_random_seed(random_seed)
    sa.run(
      params             = self.params,
      xray_structure     = tmp,
      real_space         = True,
      target_map         = self.map_data,
      restraints_manager = self.restraints_manager,
      wx                 = self.weight,
      wc                 = 1.,
      verbose            = False)
    return tmp

class ensemble_refine(object):
  def __init__(
        self,
        pdb_hierarchy,
        ensemble_xrs,
        restraints_manager,
        target_bond_rmsd,
        target_angle_rmsd,
        map_data,
        weight,
        nproc):
    adopt_init_args(self, locals())
    self.crystal_symmetry = self.ensemble_xrs[0].crystal_symmetry()
    # initialize states collector
    self.states = mmtbx.utils.states(
      pdb_hierarchy  = self.pdb_hierarchy.deep_copy())
    # run minimization
    if(self.nproc>1):
      from libtbx import easy_mp
      stdout_and_results = easy_mp.pool_map(
        processes    = self.nproc,
        fixed_func   = self.run,
        args         = self.ensemble_xrs,
        func_wrapper = "buffer_stdout_stderr")
      for so, sites_cart in stdout_and_results :
        self.states.add(sites_cart = sites_cart)
    else:
      for xrs in self.ensemble_xrs:
        sites_cart = self.run(xray_structure=xrs)
        self.states.add(sites_cart = sites_cart)

  def show_states(self, file_name="SA_ensemble_refined.pdb"):
    self.states.write(file_name=file_name,
      crystal_symmetry=self.crystal_symmetry)

  def ensemble_pdb_hierarchy_refined(self):
    return self.states.root

  def run(self, xray_structure):
    self.pdb_hierarchy.adopt_xray_structure(xray_structure)
    ro = mmtbx.refinement.real_space.individual_sites.easy(
      map_data                    = self.map_data,
      xray_structure              = xray_structure,
      pdb_hierarchy               = self.pdb_hierarchy,
      geometry_restraints_manager = self.restraints_manager,
      rms_bonds_limit             = self.target_bond_rmsd,
      rms_angles_limit            = self.target_angle_rmsd,
      max_iterations              = 250,
      selection                   = None,
      w                           = self.weight,
      log                         = None)
    return ro.xray_structure.sites_cart()

class run(object):
  def __init__(
        self,
        xray_structure,
        pdb_hierarchy,
        map_data,
        restraints_manager,
        resolution=None,
        mode="quick",
        target_bond_rmsd=0.02,
        target_angle_rmsd=2.0,
        number_of_trials=20,
        xyz_shake=5.0,
        nproc=1,
        states=None,
        show=True,
        log=None):
    adopt_init_args(self, locals())
    assert self.mode in ["quick", "thorough"]
    if(self.log is None): self.log = sys.stdout
    self.show_target(prefix="Target(start):")
    # Get refined starting model
    if(mode=="thorough"):
      weights = self.run_refine_flexible_rmsd_targets()
      er_weight = flex.mean(weights)
      sa_weight = weights[0]
    else:
      sa_weight = self.minimization(
        target_bond_rmsd  = 0.03,
        target_angle_rmsd = 3.5,
        weight            = None)
      er_weight = sa_weight
    self.show_target(prefix="Target(minimization):")
    # Generate SA ensemble of structures
    sao = run_sa(
      xray_structure     = self.xray_structure,
      pdb_hierarchy      = self.pdb_hierarchy,
      restraints_manager = self.restraints_manager,
      map_data           = self.map_data,
      number_of_trials   = self.number_of_trials,
      nproc              = self.nproc,
      weight             = sa_weight)
    self.ensemble_xrs = sao.ensemble_xrs()
    if(show): sao.show_states()
    # Refine each SA model
    self.ero = ensemble_refine(
      pdb_hierarchy       = self.pdb_hierarchy,
      ensemble_xrs        = self.ensemble_xrs,
      restraints_manager  = self.restraints_manager,
      target_bond_rmsd    = self.target_bond_rmsd,
      target_angle_rmsd   = self.target_angle_rmsd,
      map_data            = self.map_data,
      weight              = er_weight,
      nproc               = self.nproc)
    if(show): self.ero.show_states()
    # Merge SA refined models
    self.merge_models(pdb_hierarchy = self.ero.ensemble_pdb_hierarchy_refined())
    self.show_target(prefix="Target(merge):")
    # Finalize
    if(mode=="thorough"):
      self.run_refine_flexible_rmsd_targets(weights = weights)
    else:
      _ = self.minimization(
        target_bond_rmsd  = self.target_bond_rmsd,
        target_angle_rmsd = self.target_angle_rmsd,
        weight            = None)
    self.show_target(prefix="Target(final minimization):")

  def ensemble_pdb_hierarchy_refined(self):
    return self.ero.ensemble_pdb_hierarchy_refined()

  def pdb_hierarchy_overall_best(self):
    return self.pdb_hierarchy

  def merge_models(self, pdb_hierarchy):
    t0=time.time()
    assert pdb_hierarchy.models_size() == self.number_of_trials, \
      pdb_hierarchy.models_size()
    from mmtbx.building.merge_models import run as merge_models
    pdb_hierarchy_merged = merge_models(
      map_data         = self.map_data,
      resolution       = self.resolution,
      pdb_hierarchy    = pdb_hierarchy,
      crystal_symmetry = self.xray_structure.crystal_symmetry(),
      out              = self.log)
    if(self.show):
      pdb_hierarchy_merged.write_pdb_file(
        file_name="SA_ensemble_refined_merged.pdb")
    self.xray_structure = self.xray_structure.replace_sites_cart(
      new_sites=pdb_hierarchy_merged.atoms().extract_xyz())
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
    print("Time (merge): %10.3f"%(time.time()-t0), file=self.log)

  def minimization(self, target_bond_rmsd, target_angle_rmsd, weight):
    ro = mmtbx.refinement.real_space.individual_sites.easy(
      map_data                    = self.map_data,
      xray_structure              = self.xray_structure,
      pdb_hierarchy               = self.pdb_hierarchy,
      geometry_restraints_manager = self.restraints_manager,
      rms_bonds_limit             = target_bond_rmsd,
      rms_angles_limit            = target_angle_rmsd,
      max_iterations              = 50,
      selection                   = None,
      w                           = weight,
      states_accumulator          = self.states,
      log                         = None)
    self.xray_structure.replace_scatterers(
      scatterers=ro.xray_structure.scatterers())
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
    return ro.w

  def run_refine_flexible_rmsd_targets(self, weights=None):
    t0=time.time()
    ro=None
    result = flex.double()
    pairs = [(0.10,10.,), (0.05,5.0,), (0.03,3.5,), (0.025,2.5), (0.02,2.0)]
    for i, pair in enumerate(pairs):
      if(weights is None): w = None
      else:                w = weights[i-1]
      w_opt = self.minimization(
        target_bond_rmsd  = pair[0],
        target_angle_rmsd = pair[1],
        weight            = w)
      result.append(w_opt)
    print("Time (run_refine_flexible_rmsd_targets): %10.3f"%(
      time.time()-t0), file=self.log)
    return result

  def show_target(self, prefix):
    if(self.show):
      print(prefix, maptbx.real_space_target_simple(
          unit_cell   = self.xray_structure.unit_cell(),
          density_map = self.map_data,
          sites_cart  = self.xray_structure.sites_cart()), file=self.log)
