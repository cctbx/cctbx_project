from __future__ import division
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

if (0): # fixed random seed to avoid rare failures
  random.seed(1)
  flex.set_random_seed(1)


class scorer_whole_chain(object):
  def __init__(self, pdb_hierarchy, unit_cell, map_data):
    adopt_init_args(self, locals())
    self.sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    self.target = maptbx.real_space_target_simple(
      unit_cell   = self.unit_cell,
      density_map = self.map_data,
      sites_cart  = self.sites_cart)

  def update(self, sites_cart):
    target = maptbx.real_space_target_simple(
      unit_cell   = self.unit_cell,
      density_map = self.map_data,
      sites_cart  = sites_cart)
    if(target > self.target):
      self.target = target
      self.sites_cart = sites_cart
    #print self.target, target # XXX for debugging

class run(object):
  def __init__(
        self,
        xray_structure,
        pdb_hierarchy,
        map_data,
        restraints_manager,
        number_of_macro_cycles=5,
        target_bond_rmsd=0.02,
        target_angle_rmsd=2.0,
        number_of_trials=20,
        xyz_shake=10.0,
        number_of_sa_models=20,
        states=None,
        show=True):
    adopt_init_args(self, locals())
    self.ear_states=None
    self.show_target(prefix="Target(start):")
    # Get refined starting model
    weights = self.run_refine_flexible_rmsd_targets()
    self.show_target(prefix="Target(minimization):")
    # Weight
    self.weight = weights[len(weights)-1]
    # Generate SA ensemble of structures
    self.ensemble_xrs = self.simulated_annealing(weight = self.weight)
    # Explode and refine
    self.pdb_hierarchy.adopt_xray_structure(self.ensemble_xrs[0])
    sc = scorer_whole_chain(
      pdb_hierarchy = self.pdb_hierarchy,
      unit_cell     = self.xray_structure.unit_cell(),
      map_data      = map_data)
    self.ensemble_refine(scorer=sc, weight=flex.mean(weights))
    self.show_target(prefix="Target(ensemble_refine):")
    # Finalize
    self.run_refine_flexible_rmsd_targets(weights = weights)
    self.show_target(prefix="Target(minimization):")
    self.score=sc.target # save the score

  def simulated_annealing(self, weight, collect_states=True):
    t0 = time.time()
    params = sa.master_params().extract()
    params.start_temperature=50000
    params.final_temperature=0
    params.cool_rate = 25000
    params.number_of_steps = 50
    if(collect_states):
      states = mmtbx.utils.states(
        xray_structure = self.xray_structure,
        pdb_hierarchy  = self.pdb_hierarchy)
    grf = geometry_restraints.flags.flags(default=True)
    lbfgs_exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = True,
        ignore_line_search_failed_step_at_upper_bound = True,
        ignore_line_search_failed_maxfev              = True)
    result = []
    for it in xrange(self.number_of_sa_models):
      tmp = self.xray_structure.deep_copy_scatterers()
      # Shake and minimize to get variable starting points
      tmp.shake_sites_in_place(
        rms_difference = None,
        mean_distance  = self.xyz_shake)
      sites_cart_ = tmp.sites_cart()
      minimized = geometry_minimization.lbfgs(
        sites_cart                        = sites_cart_,
        correct_special_position_tolerance= 1.0,
        geometry_restraints_manager       = self.restraints_manager.geometry,
        geometry_restraints_flags         = grf,
        lbfgs_exception_handling_params   = lbfgs_exception_handling_params,
        lbfgs_termination_params          = scitbx.lbfgs.termination_parameters(
          max_iterations=50))
      tmp = tmp.replace_sites_cart(new_sites = sites_cart_)
      # Run SA
      ph = self.pdb_hierarchy.deep_copy() # need this?
      sa.run(
        params             = params,
        xray_structure     = tmp,
        real_space         = True,
        target_map         = self.map_data,
        restraints_manager = self.restraints_manager,
        wx                 = self.weight,
        wc                 = 1.,
        verbose            = False)
      ph.adopt_xray_structure(tmp)
      if(collect_states): states.add(sites_cart = tmp.sites_cart())
      result.append(tmp.deep_copy_scatterers())
    if(collect_states): states.write(file_name="SA_ensemble.pdb")
    print "Time (SA): %10.3f"%(time.time()-t0)
    return result

  def run_refine_flexible_rmsd_targets(self, weights=None):
    t0=time.time()
    ro=None
    result = flex.double()
    pairs = [(0.10,10.,), (0.05,5.0,), (0.03,3.5,), (0.025,2.5), (0.02,2.0)]
    for i, pair in enumerate(pairs):
      if(weights is None): w = None
      else:                w = weights[i-1]
      ro = mmtbx.refinement.real_space.individual_sites.easy(
        map_data                    = self.map_data,
        xray_structure              = self.xray_structure,
        pdb_hierarchy               = self.pdb_hierarchy,
        geometry_restraints_manager = self.restraints_manager,
        rms_bonds_limit             = pair[0],
        rms_angles_limit            = pair[1],
        max_iterations              = 50,
        selection                   = None,
        w                           = w,
        states_accumulator          = self.states,
        log                         = None)
      self.xray_structure.replace_scatterers(
       scatterers=ro.xray_structure.scatterers())
      result.append(ro.w)
      self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
    print "Time (run_refine_flexible_rmsd_targets): %10.3f"%(time.time()-t0)
    return result

  def ensemble_refine(self, scorer, weight, collect_states=True):
    if(collect_states):
      states = mmtbx.utils.states(
        xray_structure = self.xray_structure.deep_copy_scatterers(),
        pdb_hierarchy  = self.pdb_hierarchy.deep_copy())
    for xrs in self.ensemble_xrs:
      self.xray_structure.replace_scatterers(
         scatterers=xrs.scatterers())
      self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
      ro = mmtbx.refinement.real_space.individual_sites.easy(
        map_data                    = self.map_data,
        xray_structure              = self.xray_structure,
        pdb_hierarchy               = self.pdb_hierarchy,
        geometry_restraints_manager = self.restraints_manager,
        rms_bonds_limit             = self.target_bond_rmsd,
        rms_angles_limit            = self.target_angle_rmsd,
        max_iterations              = 250,
        selection                   = None,
        w                           = weight,
        log                         = None)
      if(collect_states):
        states.add(sites_cart = ro.xray_structure.sites_cart())
      self.xray_structure.replace_scatterers(
        scatterers=ro.xray_structure.scatterers())
      scorer.update(sites_cart = self.xray_structure.sites_cart())
      self.xray_structure = self.xray_structure.replace_sites_cart(
        new_sites=scorer.sites_cart)
      self.states.add(sites_cart = scorer.sites_cart)
      self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
    if(collect_states):
      states.add(sites_cart = self.xray_structure.sites_cart())
      states.write(file_name="SA_ensemble_refined.pdb")
    self.ear_states=states
    # Merge refined models
    from mmtbx.building.merge_models import run as merge_models
    pdb_hierarchy_merged = merge_models(
      map_data = self.map_data,
      pdb_hierarchy    = states.root,
      crystal_symmetry = self.xray_structure.crystal_symmetry())
    pdb_hierarchy_merged.write_pdb_file(file_name="SA_ensemble_refined_merged.pdb")
    self.xray_structure = self.xray_structure.replace_sites_cart(
      new_sites=pdb_hierarchy_merged.atoms().extract_xyz())
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)

  def show_target(self, prefix):
    if(self.show):
      print prefix, maptbx.real_space_target_simple(
          unit_cell   = self.xray_structure.unit_cell(),
          density_map = self.map_data,
          sites_cart  = self.xray_structure.sites_cart())
