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

def get_restraints_manager(raw_records):
      # Build geometry restraints
      # This is just because restraints_manager is not pickle-able.

      from mmtbx import monomer_library
      geo_params = monomer_library.pdb_interpretation.master_params.extract()
      processed_pdb_file = monomer_library.pdb_interpretation.process(
        mon_lib_srv              = monomer_library.server.server(),
        ener_lib                 = monomer_library.server.ener_lib(),
        raw_records              = raw_records,
        params                   = geo_params,
        strict_conflict_handling = True,
        force_symmetry           = True,
        log                      = None)
      geometry = processed_pdb_file.geometry_restraints_manager(
        show_energies                = False,
        plain_pairs_radius           = 5,
        assume_hydrogens_all_missing = True)
      restraints_manager = mmtbx.restraints.manager(
        geometry      = geometry,
        normalization = True)
      return restraints_manager

def simulated_annealing(
   xray_structure=None,
   pdb_hierarchy=None,
   restraints_manager=None,
   raw_records=None,
   map_data=None,
   weight=None,
   number_of_sa_models=None,
   xyz_shake=None,
   set_random_seed=None,
   collect_states=True,
   nproc=1
   ):

    t0 = time.time()
    params = sa.master_params().extract()
    params.start_temperature=50000
    params.final_temperature=0
    params.cool_rate = 25000
    params.number_of_steps = 50
    if(collect_states):
      states = mmtbx.utils.states(
        xray_structure = xray_structure,
        pdb_hierarchy  = pdb_hierarchy)
    grf = geometry_restraints.flags.flags(default=True)
    lbfgs_exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = True,
        ignore_line_search_failed_step_at_upper_bound = True,
        ignore_line_search_failed_maxfev              = True)
    kw_list = []

    if nproc==1:
      restraints_manager=restraints_manager
    else:
      restraints_manager=None # have to re-create it as not pickle-able

    number_of_models_list=[]
    remainder=number_of_sa_models-nproc*(number_of_sa_models//nproc)
    for it in xrange(nproc):
      nn=number_of_sa_models//nproc
      if it < remainder: nn+=1
      number_of_models_list.append(nn)

    random_seed_list=number_of_sa_models*[None]
    if set_random_seed:
      random_seed_list=[]
      for it in xrange(number_of_sa_models):
        local_random_seed=random.randint(0,10000000)
        random_seed_list.append(local_random_seed)

    for it,local_number_of_sa_models in zip(
        xrange(nproc),number_of_models_list):
      local_random_seed_list=random_seed_list[:local_number_of_sa_models]
      random_seed_list=random_seed_list[local_number_of_sa_models:]
      kw_list.append({
        'params':params,
        'xray_structure':xray_structure,
        'pdb_hierarchy':pdb_hierarchy,
        'raw_records':raw_records,
        'restraints_manager':restraints_manager,
        'grf':grf,
        'lbfgs_exception_handling_params':lbfgs_exception_handling_params,
        'map_data':map_data,
        'weight':weight,
        'xyz_shake':xyz_shake,
        'random_seed_list':local_random_seed_list,
        'number_of_sa_models':local_number_of_sa_models,
       })

    from libtbx.easy_mp import run_parallel
    sa_results=run_parallel(
     method='multiprocessing',
     qsub_command=None,
     nproc=nproc,
     target_function=run_sa,
     kw_list=kw_list)

    result = []
    for tmp_group in sa_results:
      for tmp in tmp_group:
        if(collect_states): states.add(sites_cart = tmp.sites_cart())
        result.append(tmp.deep_copy_scatterers())

    if(collect_states): states.write(file_name="SA_ensemble.pdb")
    print "Time (SA): %10.3f"%(time.time()-t0)
    return result

def run_sa(
   params=None,
   xray_structure=None,
   pdb_hierarchy=None,
   restraints_manager=None,
   raw_records=None,
   grf=None,
   lbfgs_exception_handling_params=None,
   map_data=None,
   weight=None,
   xyz_shake=None,
   random_seed_list=None,
   number_of_sa_models=None):

    if random_seed_list is None: random_seed_list=number_of_sa_models*[None]
    assert len(random_seed_list)==number_of_sa_models


    if not restraints_manager:
      restraints_manager=get_restraints_manager(raw_records)
    results=[]

    for it,random_seed in zip(
       xrange(number_of_sa_models),
       random_seed_list):

      if random_seed is not None: # XXX note: set out random seeds for each one
        random.seed(random_seed)
        flex.set_random_seed(random_seed)

      tmp = xray_structure.deep_copy_scatterers()

      # Shake and minimize to get variable starting points
      tmp.shake_sites_in_place(
        rms_difference = None,
        mean_distance  = xyz_shake)
      sites_cart_ = tmp.sites_cart()
      minimized = geometry_minimization.lbfgs(
        sites_cart                        = sites_cart_,
        correct_special_position_tolerance= 1.0,
        geometry_restraints_manager       = restraints_manager.geometry,
        geometry_restraints_flags         = grf,
        lbfgs_exception_handling_params   = lbfgs_exception_handling_params,
        lbfgs_termination_params          = scitbx.lbfgs.termination_parameters(
          max_iterations=50))
      tmp = tmp.replace_sites_cart(new_sites = sites_cart_)
      # Run SA
      sa.run(
        params             = params,
        xray_structure     = tmp,
        real_space         = True,
        target_map         = map_data,
        restraints_manager = restraints_manager,
        wx                 = weight,
        wc                 = 1.,
        verbose            = False)
      results.append(tmp)
    return results

def ensemble_refine(
    ensemble_xrs=None,
    xray_structure=None,
    pdb_hierarchy=None,
    restraints_manager=None,
    map_data=None,
    raw_records=None,
    scorer=None,
    weight=None,
    target_bond_rmsd=None,
    target_angle_rmsd=None,
    states=None,
    collect_states=True,
    nproc=1):


    t0 = time.time()
    if(collect_states):
      local_states = mmtbx.utils.states(
        xray_structure = xray_structure.deep_copy_scatterers(),
        pdb_hierarchy  = pdb_hierarchy.deep_copy())

    # Set up call to multiprocessing

    kw_list=[]

    if nproc==1:
      restraints_manager=restraints_manager
    else:
      restraints_manager=None # have to re-create it as not pickle-able

    # number of jobs per process
    nn=len(ensemble_xrs)//nproc
    # number of processes that get an extra ensemble
    remainder=len(ensemble_xrs)-nproc*(len(ensemble_xrs)//nproc)

    for it in xrange(nproc):
      if it < remainder:
        n_models=nn+1
      else:
        n_models=nn
      local_ensemble_xrs=ensemble_xrs[:n_models]
      ensemble_xrs=ensemble_xrs[n_models:]

      kw_list.append({
        'ensemble_xrs':local_ensemble_xrs,
        'xray_structure':xray_structure.deep_copy_scatterers(),
        'pdb_hierarchy':pdb_hierarchy,
        'raw_records':raw_records,
        'restraints_manager':restraints_manager,
        'map_data':map_data,
        'weight':weight,
        'target_bond_rmsd':target_bond_rmsd,
        'target_angle_rmsd':target_angle_rmsd,
       })


    from libtbx.easy_mp import run_parallel
    ensemble_results=run_parallel(
     method='multiprocessing',
     qsub_command=None,
     nproc=nproc,
     target_function=run_ensemble_refine,
     kw_list=kw_list)

    result = []
    for tmp_group in ensemble_results:
      for sc in tmp_group: # tmp is a set of scatterers
        if(collect_states):
          local_states.add(sites_cart = sc.sites_cart())
        states.add(sites_cart = sc.sites_cart())
        scorer.update(sites_cart = sc.sites_cart())
        xray_structure = xray_structure.replace_sites_cart(
          new_sites=scorer.sites_cart)

    if(collect_states):
      local_states.write(file_name="SA_ensemble_refined.pdb")
    print "Time (ensemble_refine): %10.3f"%(time.time()-t0)
    # Merge refined models
    from mmtbx.building.merge_models import run as merge_models
    pdb_hierarchy_merged = merge_models(
      map_data = map_data,
      pdb_hierarchy    = local_states.root,
      crystal_symmetry = xray_structure.crystal_symmetry())
    pdb_hierarchy_merged.write_pdb_file(file_name="SA_ensemble_refined_merged.pdb")
    xray_structure = xray_structure.replace_sites_cart(
      new_sites=pdb_hierarchy_merged.atoms().extract_xyz())
    print "Time (ensemble_refine+merge): %10.3f"%(time.time()-t0)
    return xray_structure,local_states

def run_ensemble_refine(
    ensemble_xrs=None,
    xray_structure=None,
    pdb_hierarchy=None,
    restraints_manager=None,
    map_data=None,
    raw_records=None,
    weight=None,
    target_bond_rmsd=None,
    target_angle_rmsd=None):

    results=[]
    if not restraints_manager:
      restraints_manager=get_restraints_manager(raw_records)

    for xrs in ensemble_xrs:
      xray_structure.replace_scatterers(
         scatterers=xrs.scatterers())
      pdb_hierarchy.adopt_xray_structure(xray_structure)
      ro = mmtbx.refinement.real_space.individual_sites.easy(
        map_data                    = map_data,
        xray_structure              = xray_structure,
        pdb_hierarchy               = pdb_hierarchy,
        geometry_restraints_manager = restraints_manager,
        rms_bonds_limit             = target_bond_rmsd,
        rms_angles_limit            = target_angle_rmsd,
        max_iterations              = 250,
        selection                   = None,
        w                           = weight,
        log                         = None)

      results.append(ro.xray_structure.deep_copy_scatterers())
    return results

class run(object):
  def __init__(
        self,
        xray_structure,
        pdb_hierarchy,
        map_data,
        raw_records, # XXX as restraints_manager is not pickle-able
        number_of_macro_cycles=5,
        target_bond_rmsd=0.02,
        target_angle_rmsd=2.0,
        number_of_trials=20,
        xyz_shake=10.0,
        number_of_sa_models=20,
        states=None,
        random_seed=None,
        show=True,
        nproc=1):
    adopt_init_args(self, locals())
    self.ear_states=None
    self.show_target(prefix="Target(start):")

    # Set random seed if supplied
    if random_seed is not None:
      random.seed(random_seed)
      flex.set_random_seed(random_seed)
      set_random_seed=True
    else:
      set_random_seed=False

    # Get restraints_manager
    restraints_manager=get_restraints_manager(raw_records)
    self.restraints_manager=restraints_manager

    # Get refined starting model
    weights = self.run_refine_flexible_rmsd_targets()
    self.show_target(prefix="Target(minimization):")
    # Weight
    self.weight = weights[len(weights)-1]
    # Generate SA ensemble of structures
    ensemble_xrs = simulated_annealing(
      xray_structure=xray_structure,
      pdb_hierarchy=pdb_hierarchy,
      restraints_manager=restraints_manager,
      raw_records=raw_records,
      map_data=map_data,
      weight=self.weight,
      number_of_sa_models=number_of_sa_models,
      xyz_shake=xyz_shake,
      set_random_seed=set_random_seed,
      nproc=nproc)

    # Explode and refine
    self.pdb_hierarchy.adopt_xray_structure(ensemble_xrs[0])
    sc = scorer_whole_chain(
      pdb_hierarchy = self.pdb_hierarchy,
      unit_cell     = self.xray_structure.unit_cell(),
      map_data      = map_data)
    self.xray_structure,ear_states=ensemble_refine(
      xray_structure=xray_structure,
      pdb_hierarchy=pdb_hierarchy,
      ensemble_xrs=ensemble_xrs,
      map_data=map_data,
      restraints_manager=restraints_manager,
      raw_records=raw_records,
      scorer=sc,
      weight=flex.mean(weights),
      target_bond_rmsd=target_angle_rmsd,
      target_angle_rmsd=target_angle_rmsd,
      states=self.states,
      nproc=nproc
      )
    self.ear_states=ear_states
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
    self.show_target(prefix="Target(ensemble_refine):")
    # Finalize
    self.run_refine_flexible_rmsd_targets(weights = weights)
    self.show_target(prefix="Target(minimization):")
    self.score=sc.target # save the score


  def run_refine_flexible_rmsd_targets(self, weights=None):
    t0=time.time()
    ro=None
    result = flex.double()
    pairs = [(0.10,10.,), (0.05,5.0,), (0.03,3.5,), (0.025,2.5), (0.02,2.0)]
    for i, pair in enumerate(pairs):
      if(weights is None): w = None
      else:                w = weights[i-1]
      if not self.restraints_manager:
        self.restraints_manager=get_restraints_manager(self.raw_records)
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


  def show_target(self, prefix):
    if(self.show):
      print prefix, maptbx.real_space_target_simple(
          unit_cell   = self.xray_structure.unit_cell(),
          density_map = self.map_data,
          sites_cart  = self.xray_structure.sites_cart())
