from __future__ import division
import iotbx.phil
from mmtbx.refinement import print_statistics
from cctbx.maptbx import real_space_target_and_gradients
from libtbx import adopt_init_args
from libtbx.str_utils import format_value
import scitbx.lbfgs
from cctbx import maptbx
from cctbx.array_family import flex
from mmtbx import utils
from mmtbx.validation import clashscore
from cctbx import miller
from libtbx.test_utils import approx_equal

master_params_str = """\
real_space_refinement
  .short_caption = Real-space refinement
  .style = menu_item auto_align
{
    mode = simple diff_map *lockit
      .type = choice(multi=False)
      .help = Real space refinement method (diff_map is much slower but might \
              have larger convergence radius especially at low resolution)
    target_map_name = 2mFo-DFc
      .type = str
      .help = Target map type to refine against
    model_map_name = Fc
      .type = str
      .help = Current model map type (used only if diff_map is selected)
    grid_resolution_factor = 1./4
      .type = float
      .help = Defines coarseness of the map
    target_weights = value grid_search *gradients_ratio
      .type = choice(multi=False)
      .help = Method for relative target weight determination
    real_space_target_weight = None
      .type = float
      .help = Weight for the real-space target (used only if diff_map is selected)
    restraints_target_weight = None
      .type = float
      .help = Weight for the restraints target (used only if diff_map is selected)
    number_of_cycles = 1
      .type = int
      .help = Number of minimization cycles
    rmsd_max_bonds = 0.07
      .type = float
      .short_caption = Max. RMSD(bonds)
      .help = Refinement result is ignored if this max allowable deviation exceeded
    rmsd_max_angles = 5.0
      .type = float
      .short_caption = Max. RMSD(angles)
      .help = Refinement result is ignored if this max allowable deviation exceeded
    grid_search_scales = 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.5 3.0
      .type = floats
      .help = Defines the range of values for the weight grid search
    verbose = 1
      .type = int
      .help = All output is supressed if it is negative
    lockit_parameters
      .expert_level=3
    {
      include scope \
        mmtbx.command_line.lockit.coordinate_refinement_export_phil_str
    }
  }
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

class run(object):

  def __init__(self, fmodel, model, params=None, log=None):
    adopt_init_args(self, locals())
    if(params is None):
      params = master_params().extract()
      params = params.real_space_refinement
      self.params = params
    if(log is None): log = sys.stdout
    if(params.verbose):
      print_statistics.make_header("Real-space coordinate refinement", out=log)
    fft_map_target = self.compute_map(map_type = params.target_map_name)
    geom = self.show()
    step = self.fmodel.f_obs_work().d_min()*\
      params.grid_resolution_factor
    if(params.real_space_target_weight is not None):
      real_space_target_weight = params.real_space_target_weight
    else:
      real_space_target_weight = 1.
    if(params.target_weights == "grid_search"):
      restraints_target_weights = params.grid_search_scales
    elif(params.target_weights == "value"):
      assert [params.real_space_target_weight,
              params.restraints_target_weight].count(None)==0
      if(params.real_space_target_weight is not None):
        real_space_target_weight = params.real_space_target_weight
      if(params.restraints_target_weight is not None):
        restraints_target_weights = [params.restraints_target_weight]
    elif(params.target_weights == "gradients_ratio"):
      real_space_target_weight = 1.
      restraints_target_weights = [1.,]
    if(params.target_weights == "gradients_ratio" or
       params.target_weights == "grid_search"):
      map_current = None
      if(params.mode == "diff_map"):
          map_current = self.compute_map(map_type =
            params.model_map_name).real_map_unpadded()
      restraints_target_weight_ = \
        real_space_target_and_gradients.target_and_gradients(
          unit_cell                   = self.fmodel.xray_structure.unit_cell(),
          map_target                  = fft_map_target.real_map_unpadded(),
          map_current                 = map_current,
          real_space_gradients_delta  = step,
          sites_frac                  = self.fmodel.xray_structure.sites_frac(),
          geometry_restraints_manager = model.restraints_manager,
          real_space_target_weight    = 1,
          restraints_target_weight    = 1,
          sites_cart                  = self.fmodel.xray_structure.sites_cart(),
          target_type                 = params.mode).weight()
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True,
      update_f_mask  = True)
    xrs_start = self.fmodel.xray_structure.deep_copy_scatterers()
    geom = model.geometry_statistics(ignore_hd = True) # XXX
    b_start, a_start = geom.b_mean, geom.a_mean
    best_weight_scale = None
    for i_w, restraints_target_weight in enumerate(restraints_target_weights):
      self.fmodel.update_xray_structure(
        xray_structure = xrs_start,
        update_f_calc  = True,
        update_f_mask  = True)
      model.xray_structure = self.fmodel.xray_structure
      for cycle in range(1,params.number_of_cycles+1):
        map_current = None
        if(params.mode == "diff_map"):
            map_current = self.compute_map(map_type =
              params.model_map_name).real_map_unpadded()
        if(params.target_weights == "gradients_ratio" or
           params.target_weights == "grid_search"):
          if(params.target_weights == "grid_search"):
            restraints_target_weight = restraints_target_weight_*\
              restraints_target_weight
          else: restraints_target_weight = restraints_target_weight_
        minimized = real_space_target_and_gradients.minimization(
          xray_structure              = xrs_start.deep_copy_scatterers(),
          miller_array                = self.fmodel.f_obs_work(),
          crystal_gridding            = fft_map_target,
          map_target                  = fft_map_target.real_map_unpadded(),
          step                        = step,
          restraints_target_weight    = restraints_target_weight,
          target_type                 = params.mode,
          geometry_restraints_manager = model.restraints_manager)
        self.fmodel.update_xray_structure(
          xray_structure = minimized.xray_structure,
          update_f_calc  = True,
          update_f_mask  = True)
        model.xray_structure = self.fmodel.xray_structure
        self.show(weight=restraints_target_weight)
        if(i_w == 0):
          xrs_best = self.fmodel.xray_structure.deep_copy_scatterers()
          r_free_best = self.fmodel.r_free()
        else:
          r_free = self.fmodel.r_free()
          geom = model.geometry_statistics(ignore_hd = True) # XXX
          if(r_free < r_free_best and ((geom.b_mean <= params.rmsd_max_bonds and
             geom.a_mean <= params.rmsd_max_angles) or
             ((geom.b_mean > params.rmsd_max_bonds or geom.a_mean > params.rmsd_max_angles)
             and (geom.b_mean <= b_start or geom.a_mean <= a_start)) )):
            r_free_best = r_free
            xrs_best = minimized.xray_structure.deep_copy_scatterers()
            best_weight_scale = restraints_target_weight/restraints_target_weight_
    self.fmodel.update_xray_structure(
      xray_structure = xrs_best,
      update_f_calc  = True,
      update_f_mask  = True)
    model.xray_structure = self.fmodel.xray_structure
    self.show()
    print >> self.log, "Best weight scale: %s"%format_value("%8.4f",best_weight_scale)

  def compute_map(self, map_type=None, use_all_data=False):
    e_map_manager = self.fmodel.electron_density_map() # XXX pass in map filling options
    fft_map = e_map_manager.fft_map(
      resolution_factor = self.params.grid_resolution_factor,
      map_type          = map_type,
      use_all_data      = use_all_data)
    fft_map.apply_sigma_scaling()
    return fft_map

  def show(self, weight=None):
    geom = self.model.geometry_statistics(ignore_hd = True) # XXX
    if(self.params.verbose):
      p1 = "|-real space refinement (target=%s)"%self.params.mode
      p1 +="-"*(78-len(p1))+"|"
      print >> self.log, p1
      if(weight is not None):
        p2 = "| r_work=%s r_free=%s  rmsd bonds=%s angles=%s  weight=%s  |"%(
          format_value("%6.4f", self.fmodel.r_work()),
          format_value("%6.4f", self.fmodel.r_free()),
          format_value("%5.3f", geom.b_mean),
          format_value("%5.2f", geom.a_mean),
          format_value("%7.2f", weight))
      else:
        p2 = "| r_work=%s r_free=%s rmsd bonds=%s angles=%s%s|"%(
          format_value("%6.4f", self.fmodel.r_work()),
          format_value("%6.4f", self.fmodel.r_free()),
          format_value("%5.3f", geom.b_mean),
          format_value("%5.2f", geom.a_mean), " "*19)
      print >> self.log, p2
      print >> self.log, "|"+"-"*77+"|"
    return geom

class simple(object):
  def __init__(self,
               target_map,
               selection,
               real_space_gradients_delta,
               geometry_restraints_manager=None,
               max_iterations=150):
    adopt_init_args(self, locals())
    self.lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations = max_iterations)
    self.lbfgs_exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True,
      ignore_line_search_failed_step_at_upper_bound = True,
      ignore_line_search_failed_maxfev              = True)
    self.refined = None

  def refine(self, weight, xray_structure):
    assert self.selection.size() == xray_structure.scatterers().size()
    self.refined = maptbx.real_space_refinement_simple.lbfgs(
      selection_variable              = self.selection,
      sites_cart                      = xray_structure.sites_cart(),
      density_map                     = self.target_map,
      geometry_restraints_manager     = self.geometry_restraints_manager,
      real_space_target_weight        = weight,
      real_space_gradients_delta      = self.real_space_gradients_delta,
      lbfgs_termination_params        = self.lbfgs_termination_params,
      lbfgs_exception_handling_params = self.lbfgs_exception_handling_params)

  def sites_cart(self):
    assert self.refined is not None
    sites_cart = self.refined.sites_cart
    sites_cart.set_selected(self.selection, self.refined.sites_cart_variable)
    return sites_cart

class diff_map(object):
  def __init__(self,
               miller_array,
               crystal_gridding,
               map_target,
               geometry_restraints_manager,
               restraints_target_weight = 1,
               max_iterations = 500,
               min_iterations = 500):
    adopt_init_args(self, locals())
    self.step = miller_array.d_min()/4.
    self.refined = None

  def refine(self, weight, sites_cart=None, xray_structure=None):
    assert xray_structure is not None and [sites_cart,xray_structure].count(None)==1
    self.refined = real_space_target_and_gradients.minimization(
      xray_structure              = xray_structure,
      miller_array                = self.miller_array,
      crystal_gridding            = self.crystal_gridding,
      map_target                  = self.map_target,
      max_iterations              = self.max_iterations,
      min_iterations              = self.min_iterations,
      step                        = self.step,
      real_space_target_weight    = weight,
      restraints_target_weight    = self.restraints_target_weight,
      geometry_restraints_manager = self.geometry_restraints_manager,
      target_type                 = "diff_map")

  def sites_cart(self):
    assert self.refined is not None
    return self.refined.xray_structure.sites_cart()

class refinery(object):
  def __init__(self,
               refiner,
               xray_structure,
               start_trial_weight_value = 50.,
               weight_sample_rate = 10,
               rms_bonds_limit = 0.03,
               rms_angles_limit = 3.0,
               optimize_weight = True):
    self.rms_angles_start = None
    self.rms_bonds_start = None
    self.refiner = refiner
    self.weight_start=start_trial_weight_value
    self.rms_bonds_start, self.rms_angles_start  = \
      self.rmsds(sites_cart=xray_structure.sites_cart())
    self.weight_sample_rate = weight_sample_rate
    # results
    self.weight_final = None
    self.sites_cart_result = None
    self.rms_bonds_final,self.rms_angles_final = None,None
    #
    pool = {}
    bonds = flex.double()
    angles = flex.double()
    weights = flex.double()
    #
    weight = start_trial_weight_value
    weight_last = weight
    self.adjust_weight_sample_rate(weight=weight)
    if(optimize_weight):
      while True:
        self.adjust_weight_sample_rate(weight=weight_last)
        refiner.refine(
          xray_structure = xray_structure.deep_copy_scatterers(), # XXX
          weight     = weight)
        sites_cart_result = refiner.sites_cart()
        bd, ad = self.rmsds(sites_cart=sites_cart_result)
        bonds.append(bd)
        angles.append(ad)
        weights.append(weight)
        pool.setdefault(weight,[]).append([sites_cart_result.deep_copy(),bd,ad])
        if(refiner.geometry_restraints_manager is None): break
        weight_last = weight
        if(ad>rms_angles_limit or bd > rms_bonds_limit):
          weight -= self.weight_sample_rate
        else:
          weight += self.weight_sample_rate
        #print ">>> ", "%6.2f %6.2f"%(weight, weight_last), "%6.4f %5.2f"%(bd, ad)
        if((weight<0 or weight>1000) or weight in weights): break
    else:
      refiner.refine(
        xray_structure = xray_structure.deep_copy_scatterers(), # XXX
        weight     = weight)
      sites_cart_result = refiner.sites_cart()
    # select results
    sel  = bonds <= rms_bonds_limit
    sel &= angles <= rms_angles_limit
    bonds   = bonds  .select(sel)
    angles  = angles .select(sel)
    weights = weights.select(sel)
    if(sel.count(True)>0):
      bond_max = flex.max(bonds)
      ind = None
      for i, b in enumerate(bonds):
        if(b==bond_max):
          ind = i
          break
      assert ind is not None
      self.weight_final = weights[ind]
      self.sites_cart_result = pool[self.weight_final][0][0]
      self.rms_bonds_final,self.rms_angles_final = \
        self.rmsds(sites_cart=self.sites_cart_result)
      assert approx_equal(pool[self.weight_final][0][2], angles[ind])
      assert approx_equal(pool[self.weight_final][0][1], bonds[ind])
      assert approx_equal(self.rms_angles_final, angles[ind])
      assert approx_equal(self.rms_bonds_final, bonds[ind])

  def rmsds(self, sites_cart):
    b,a = None,None
    if(self.refiner.geometry_restraints_manager is not None):
      es = self.refiner.geometry_restraints_manager.energies_sites(
        sites_cart = sites_cart)
      a = es.angle_deviations()[2]
      b = es.bond_deviations()[2]
    return b,a

  def adjust_weight_sample_rate(self, weight):
    if(  weight <= 1000.): self.weight_sample_rate=100.
    if(  weight <= 100.):  self.weight_sample_rate=10.
    if(  weight <= 10.):   self.weight_sample_rate=1.
    elif(weight <= 1.0):   self.weight_sample_rate=0.1
    elif(weight <= 0.1):   self.weight_sample_rate=0.01


class monitor(object):
  def __init__(self,
               xray_structure,
               target_map,
               geometry_restraints_manager,
               xray_structure_reference = None):
    adopt_init_args(self, locals())
    self.xray_structure_start = xray_structure.deep_copy_scatterers()
    self.map_cc = None
    self.rmsd_b = None
    self.rmsd_a = None
    self.dist_from_start = None
    self.dist_from_ref = None
    self.update(xray_structure=xray_structure)

  def update(self, xray_structure, accept_any=False):
    def cc(xray_structure):
      mc = self.target_map.miller_array.structure_factors_from_scatterers(
        xray_structure = xray_structure).f_calc()
      fft_map = miller.fft_map(
        crystal_gridding     = self.target_map.crystal_gridding,
        fourier_coefficients = mc)
      fft_map.apply_sigma_scaling()
      current_map = fft_map.real_map_unpadded()
      return flex.linear_correlation(x=current_map.as_1d(),
        y=self.target_map.data.as_1d()).coefficient()
    map_cc = cc(xray_structure=xray_structure)
    if(accept_any or (self.map_cc is None or self.map_cc < map_cc)):
      self.map_cc = map_cc
      es = self.geometry_restraints_manager.energies_sites(
        sites_cart = xray_structure.sites_cart())
      self.rmsd_a = es.angle_deviations()[2]
      self.rmsd_b = es.bond_deviations()[2]
      self.dist_from_start = flex.mean(self.xray_structure_start.distances(
        other = xray_structure))
      if(self.xray_structure_reference is not None):
        self.dist_from_ref = flex.mean(self.xray_structure_reference.distances(
          other = xray_structure))
      self.xray_structure = xray_structure # must be last line in this function

  def show(self, suffix=""):
    if(self.dist_from_ref is None):
      f="cc: %6.4f rmsd_b: %6.4f rmsd_a: %5.2f d(start): %6.3f"
      print f%(self.map_cc,self.rmsd_b,self.rmsd_a,self.dist_from_start),suffix
    else:
      f="cc: %6.4f rmsd_b: %6.4f rmsd_a: %5.2f d(start): %6.3f d(ref): %6.3f"
      print f%(self.map_cc,self.rmsd_b,self.rmsd_a,self.dist_from_start,
        self.dist_from_ref),suffix

def run_tmp(target_map,
            xray_structure,
            geometry_restraints_manager,
            xray_structure_reference = None,
            max_iterations = 500,
            macro_cycles   = 10): # XXX Replace run() above with this run_tmp
  sel = flex.bool(xray_structure.scatterers().size(), True)
  d_min = target_map.miller_array.d_min()
  #
  rsr_simple_refiner = simple(
    target_map                  = target_map.data,
    selection                   = sel,
    real_space_gradients_delta  = d_min/4,
    max_iterations              = max_iterations,
    geometry_restraints_manager = geometry_restraints_manager.geometry)

  rsr_diff_map_refiner = diff_map(
    miller_array                = target_map.miller_array,
    crystal_gridding            = target_map.crystal_gridding,
    map_target                  = target_map.data,
    geometry_restraints_manager = geometry_restraints_manager.geometry,
    restraints_target_weight    = 1,
    max_iterations              = max_iterations,
    min_iterations              = max_iterations)
  #
  #restraints_manager.geometry.remove_dihedrals_in_place(sel)
  if(xray_structure_reference is not None):
    xray_structure_reference = xray_structure_reference.deep_copy_scatterers()
  monitor_object = monitor(
    xray_structure = xray_structure.deep_copy_scatterers(),
    target_map = target_map,
    geometry_restraints_manager = geometry_restraints_manager.geometry,
    xray_structure_reference = xray_structure_reference)
  monitor_object.show(suffix="start")
  #
  tmp = monitor_object.xray_structure.deep_copy_scatterers()
  weight_d = 50
  weight_s = 50
  #
  b_upper_limit = 0.03#0.2
  b_lower_limit = 0.03#0.02
  b_inc = (b_upper_limit-b_lower_limit)/macro_cycles
  a_upper_limit = 3.0#20.0
  a_lower_limit = 3.0#2.5
  a_inc = (a_upper_limit-a_lower_limit)/macro_cycles
  #
  rms_bonds_limit  = b_upper_limit
  rms_angles_limit = a_upper_limit
  for i in range(macro_cycles):
    if(i>3 and i%2==0):
      tmp.shake_sites_in_place(mean_distance=1)

    if(0):
      target_type = "diff_map"
      refined = refinery(
        refiner          = rsr_diff_map_refiner,
        xray_structure   = tmp,
        start_trial_weight_value = weight_d,
        rms_bonds_limit  = rms_bonds_limit,
        rms_angles_limit = rms_angles_limit)
      if(refined.sites_cart_result is not None):
        tmp = tmp.replace_sites_cart(refined.sites_cart_result)
        weight_d = refined.weight_final
        monitor_object.update(xray_structure=tmp)#, accept_any=True)
        monitor_object.show(suffix=" | target=%s weight: %s"%(target_type, str(weight_d)))
        tmp = monitor_object.xray_structure.deep_copy_scatterers()
        if(weight_d<0.1): weight_d=1
      else: print "Refinement failed."
    #
    if(1):
      target_type = "simple"
      refined = refinery(
        refiner          = rsr_simple_refiner,
        xray_structure   = tmp,
        start_trial_weight_value = weight_s,
        rms_bonds_limit  = rms_bonds_limit,
        rms_angles_limit = rms_angles_limit)
      if(refined.sites_cart_result is not None):
        tmp = tmp.replace_sites_cart(refined.sites_cart_result)
        weight_s = refined.weight_final
        monitor_object.update(xray_structure=tmp)#, accept_any=True)
        monitor_object.show(suffix=" | target=%s weight: %s"%(target_type, str(weight_s)))
        tmp = monitor_object.xray_structure.deep_copy_scatterers()
      else: print "Refinement failed."
    rms_bonds_limit -= b_inc
    rms_angles_limit-= a_inc
  return monitor_object.xray_structure


class box_refinement_manager(object):
  def __init__(self,
               xray_structure,
               pdb_hierarchy,
               target_map,
               geometry_restraints_manager,
               real_space_gradients_delta=1./4,
               max_iterations = 50):
    self.xray_structure = xray_structure
    self.sites_cart = xray_structure.sites_cart()
    self.pdb_hierarchy = pdb_hierarchy
    self.target_map = target_map
    self.geometry_restraints_manager = geometry_restraints_manager
    self.max_iterations=max_iterations
    self.real_space_gradients_delta = real_space_gradients_delta

  def update_xray_structure(self, new_xray_structure):
    self.xray_structure = new_xray_structure

  def update_target_map(self, new_target_map):
    self.target_map = new_target_map

  def refine(self,
             selection,
             selection_buffer_radius=5,
             box_cushion=2,
             monitor_clashscore=False):
    sites_cart_moving = self.sites_cart
    selection_within = self.xray_structure.selection_within(
      radius    = selection_buffer_radius,
      selection = selection)
    sel = selection.select(selection_within)
    iselection = flex.size_t()
    for i, state in enumerate(selection):
      if state:
        iselection.append(i)
    if monitor_clashscore:
      pdb_string = utils.write_pdb_file(
                     xray_structure=self.xray_structure,
                     pdb_hierarchy=self.pdb_hierarchy,
                     write_cryst1_record = False,
                     selection = selection_within,
                     return_pdb_string = True)
      csm = clashscore.probe_clashscore_manager(
              pdb_string=pdb_string)
      self.clashscore = csm.clashscore
    box = utils.extract_box_around_model_and_map(
            xray_structure   = self.xray_structure,
            pdb_hierarchy    = self.pdb_hierarchy,
            map_data         = self.target_map,
            selection        = selection_within,
            box_cushion      = box_cushion)
    new_unit_cell = box.xray_structure_box.unit_cell()
    geo_box = \
      self.geometry_restraints_manager.select(box.selection_within)
    geo_box.discard_symmetry(new_unit_cell=new_unit_cell)
    map_box = box.map_box
    sites_cart_box = box.xray_structure_box.sites_cart()
    selection = flex.bool(sites_cart_box.size(), True)
    rsr_simple_refiner = simple(
      target_map                  = map_box,
      selection                   = sel,
      real_space_gradients_delta  = self.real_space_gradients_delta,
      max_iterations              = self.max_iterations,
      geometry_restraints_manager = geo_box)
    real_space_result = refinery(
      refiner          = rsr_simple_refiner,
      xray_structure   = box.xray_structure_box)
    sites_cart_box_refined = real_space_result.sites_cart_result
    sites_cart_box_refined_shifted_back = \
      sites_cart_box_refined + box.shift_to_map_boxed_sites_back
    sites_cart_refined = sites_cart_box_refined_shifted_back.select(
                           sel)
    sites_cart_moving = sites_cart_moving.set_selected(
      iselection, sites_cart_refined)
    self.xray_structure.set_sites_cart(sites_cart_moving)
    self.sites_cart = self.xray_structure.sites_cart()
    if monitor_clashscore:
      pdb_string = utils.write_pdb_file(
                     xray_structure=self.xray_structure,
                     pdb_hierarchy=self.pdb_hierarchy,
                     write_cryst1_record = False,
                     selection = selection_within,
                     return_pdb_string = True)
      csm = clashscore.probe_clashscore_manager(
              pdb_string=pdb_string)
      self.clashscore_refined = csm.clashscore
