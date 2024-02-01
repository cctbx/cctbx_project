from __future__ import absolute_import, division, print_function
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
from six.moves import range
from libtbx.test_utils import approx_equal

if (0): # fixed random seed to avoid rare failures
  random.seed(1)
  flex.set_random_seed(1)

def sa_simple(rm, xrs, ph, map_data, log):
  tmp_xrs = xrs.deep_copy_scatterers()
  ro = mmtbx.refinement.real_space.individual_sites.easy(
    map_data                    = map_data,
    xray_structure              = tmp_xrs,
    pdb_hierarchy               = ph.deep_copy(),
    geometry_restraints_manager = rm,
    rms_bonds_limit             = 0.01,
    rms_angles_limit            = 1.0,
    selection                   = None, #TODO
    log                         = log)
  weight = ro.w
  #
  from mmtbx.dynamics import simulated_annealing as sa
  tmp = xrs.deep_copy_scatterers()
  params = sa.master_params().extract()
  params.start_temperature=5000
  params.cool_rate=500
  sa.run(
    params             = params,
    xray_structure     = tmp,
    real_space         = True,
    target_map         = map_data,
    restraints_manager = rm,
    wx                 = weight,
    wc                 = 1.,
    verbose            = False,
    log                = log)
  return tmp.sites_cart()

class scorer(object):
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
      pdb_hierarchy = self.pdb_hierarchy.deep_copy())
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
    for it in range(self.number_of_trials):
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
        score_method=["merge_models",],
        resolution=None,
        mode="quick",
        target_bond_rmsd=0.02,
        target_angle_rmsd=2.0,
        number_of_trials=20,
        nproc=1,
        states=None,
        map_data_ref=None,
        fragments=None,
        show=True,
        log=None):
    if not map_data_ref:
      map_data_ref = map_data # XXX in case map_data_ref is not supplied
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
    #
    assert approx_equal(self.xray_structure.sites_cart(),
                        self.pdb_hierarchy.atoms().extract_xyz())
    self.final_cc = None
    if(self.resolution is not None):
      self.final_cc = self.get_cc(xrs=self.xray_structure)
    self.show_target(prefix="Target(final minimization):")

  def ensemble_pdb_hierarchy_refined(self):
    return self.ero.ensemble_pdb_hierarchy_refined()

  def pdb_hierarchy_overall_best(self):
    return self.pdb_hierarchy

  def get_cc(self, xrs=None, sites_cart=None):
    if(sites_cart is not None):
      xrs = self.xray_structure.deep_copy_scatterers()
      xrs.set_sites_cart(sites_cart)
    fc = xrs.structure_factors(d_min=self.resolution).f_calc()
    fo = fc.structure_factors_from_map(map=self.map_data_ref, use_sg=False)
    return fc.map_correlation(other=fo)

  def _score_by_geometry(self, pdb_hierarchy):
    uc = self.xray_structure.unit_cell()
    sites_cart = flex.vec3_double(self.xray_structure.scatterers().size())
    # Planes
    sel = flex.size_t(self.fragments.planes_all)
    SC = scorer(
      pdb_hierarchy = self.pdb_hierarchy.select(sel),
      unit_cell     = uc,
      map_data      = self.map_data_ref)
    for model in pdb_hierarchy.models():
      SC.update(sites_cart = model.atoms().extract_xyz().select(sel))
    sites_cart = sites_cart.set_selected(sel, SC.sites_cart)
    # Chirals
    sel = flex.size_t(self.fragments.chirals_all)
    SC = scorer(
      pdb_hierarchy = self.pdb_hierarchy.select(sel),
      unit_cell     = uc,
      map_data      = self.map_data)
    for model in pdb_hierarchy.models():
      SC.update(sites_cart = model.atoms().extract_xyz().select(sel))
    sel = flex.size_t(self.fragments.chirals_unique)
    tmp = SC.sites_cart.select(self.fragments.chirals_mapping)
    sites_cart = sites_cart.set_selected(sel, tmp)
    # Dihedrals
    sel = flex.size_t(self.fragments.dihedrals_all)
    SC = scorer(
      pdb_hierarchy = self.pdb_hierarchy.select(sel),
      unit_cell     = uc,
      map_data      = self.map_data)
    for model in pdb_hierarchy.models():
      SC.update(sites_cart = model.atoms().extract_xyz().select(sel))
    sel = flex.size_t(self.fragments.dihedrals_unique)
    tmp = SC.sites_cart.select(self.fragments.dihedrals_mapping)
    sites_cart = sites_cart.set_selected(sel, tmp)
    #
    self.xray_structure.set_sites_cart(sites_cart = sites_cart)
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
    #
    sites_cart = sa_simple(
      rm       = self.restraints_manager, xrs=self.xray_structure,
      ph       = self.pdb_hierarchy,
      map_data = self.map_data,
      log      = None)
    self.xray_structure.set_sites_cart(sites_cart = sites_cart)
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
    sites_cart = sa_simple(
        rm=self.restraints_manager, xrs=self.xray_structure,
        ph=self.pdb_hierarchy,
        map_data=self.map_data,
        log=None)
    return sites_cart

  def _score_by_cc(self, pdb_hierarchy):
    t0 = time.time()
    xrs = self.xray_structure.deep_copy_scatterers()
    cc=-1
    sites_cart = None
    for model in pdb_hierarchy.models():
      sites_cart_ = model.atoms().extract_xyz()
      xrs.set_sites_cart(sites_cart = sites_cart_)
      cc_ = self.get_cc(xrs=xrs)
      if(cc_>cc):
        cc = cc_
        sites_cart = sites_cart_.deep_copy()
    return sites_cart

  def merge_models(self, pdb_hierarchy):
    #
    t0=time.time()
    sites_cart = []
    if("geometry" in self.score_method or "cc" in self.score_method):
      if("geometry" in self.score_method):
        sites_cart.append(self._score_by_geometry(pdb_hierarchy = pdb_hierarchy))
      if("cc" in self.score_method):
        sites_cart.append(self._score_by_cc(pdb_hierarchy = pdb_hierarchy))
      sites_cart_best = None
      cc_best = -1
      for sc in sites_cart:
        cc = self.get_cc(sites_cart=sc)
        if(cc>cc_best):
          cc_best = cc
          sites_cart_best = sc.deep_copy()
      self.xray_structure.set_sites_cart(sites_cart = sites_cart_best)
      self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
      self.pdb_hierarchy.write_pdb_file(file_name="merged.pdb")
    #
    if("merge_models" in self.score_method):
      from mmtbx.building.merge_models import run as merge_models
      pdb_hierarchy_merged, pdb_out = merge_models(
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
