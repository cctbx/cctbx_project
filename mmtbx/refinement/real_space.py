import iotbx.phil
from mmtbx.refinement import print_statistics
from cctbx.maptbx import real_space_target_and_gradients
from libtbx import adopt_init_args
from libtbx.str_utils import format_value

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

    lockit_run_conditions
      .expert_level=3
    {
      always = False
        .type = bool
      acceptable_r_factor_increase = 2
        .type = float
      high_resolution_cutoff = 1.5
        .type = float
      low_resolution_cutoff = 3.5
        .type = float
      min_r_work = 0.25
        .type = float
      r_free_r_work_diff = 5
        .type = float
      debug = False
        .type = bool
    }

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
