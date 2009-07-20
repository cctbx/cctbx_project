from cctbx.array_family import flex
import iotbx.phil
from mmtbx.refinement import print_statistics
from libtbx.test_utils import approx_equal
from cctbx.maptbx import real_space_target_and_gradients
from libtbx import adopt_init_args

master_params_str = """\
real_space_refinement
  {
    mode = simple *diff_map
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
    number_of_cycles = 2
      .type = int
      .help = Number of minimization cycles
    rmsd_max_bonds = 0.07
      .type = float
      .help = Refinement result is ignored if this max allowable deviation exceeded
    rmsd_max_angles = 5.0
      .type = float
      .help = Refinement result is ignored if this max allowable deviation exceeded
    grid_search_weight_min_max_step = 0 105 5
      .type = ints
      .help = Defines the range of values for the weight grid search
    verbose = 1
      .type = int
      .help = All output is supressed if it is negative
  }
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

class run(object):

  def __init__(self, fmodel, model, params=None, log=None):
    adopt_init_args(self, locals())
    if(params is None):
      params = master_params().extract()
      self.params = params
    if(log is None): log = sys.stdout
    if(params.verbose):
      print_statistics.make_header("Real-space coordinate refinement", out=log)
    fft_map_target = self.compute_map(map_type = params.target_map_name)
    geom = self.show(prefix = "Start")
    step = fmodel.f_obs_w.d_min()*\
      params.grid_resolution_factor
    if(params.real_space_target_weight is not None):
      real_space_target_weight = params.real_space_target_weight
    else:
      real_space_target_weight = 1.
    if(params.target_weights == "grid_search"):
      gsw = params.grid_search_weight_min_max_step
      assert gsw[0] < gsw[1]
      assert gsw[0]*gsw[1]*gsw[2] >= 0
      restraints_target_weights = [float(i) for i in range(gsw[0],gsw[1],gsw[2])]
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
    xrs_start = fmodel.xray_structure.deep_copy_scatterers()
    xrs_best = fmodel.xray_structure.deep_copy_scatterers()
    r_work_best = fmodel.r_work()
    geom = model.geometry_statistics(ignore_hd = True) # XXX
    b_start, a_start = geom.b_mean, geom.a_mean
    for restraints_target_weight in restraints_target_weights:
      fmodel.update_xray_structure(
        xray_structure = xrs_start,
        update_f_calc  = True,
        update_f_mask  = True)
      model.xray_structure = fmodel.xray_structure
      for cycle in range(1,params.number_of_cycles+1):
        if(params.target_weights == "gradients_ratio"):
          map_current = None
          if(params.mode == "diff_map"):
            map_current = self.compute_map(map_type =
              params.model_map_name).real_map_unpadded()
          restraints_target_weight = \
            real_space_target_and_gradients.target_and_gradients(
              unit_cell                   = fmodel.xray_structure.unit_cell(),
              map_target                  = fft_map_target.real_map_unpadded(),
              map_current                 = map_current,
              real_space_gradients_delta  = step,
              sites_frac                  = fmodel.xray_structure.sites_frac(),
              geometry_restraints_manager = model.restraints_manager,
              real_space_target_weight    = 1,
              restraints_target_weight    = 1,
              sites_cart                  = fmodel.xray_structure.sites_cart(),
              target_type                 = params.mode).weight()
        # ??? if(cycle==1 and restraints_target_weight > 100.):
        # ???   restraints_target_weight=100.0
        minimized = real_space_target_and_gradients.minimization(
          xray_structure              = xrs_start.deep_copy_scatterers(),
          miller_array                = fmodel.f_obs_w,
          crystal_gridding            = fft_map_target,
          map_target                  = fft_map_target.real_map_unpadded(),
          step                        = step,
          restraints_target_weight    = restraints_target_weight,
          target_type                 = params.mode,
          geometry_restraints_manager = model.restraints_manager)
        fmodel.update_xray_structure(
          xray_structure = minimized.xray_structure,
          update_f_calc  = True,
          update_f_mask  = True)
        model.xray_structure = fmodel.xray_structure
        if(params.verbose):
          print >> log, "  weight=%10.3f"%restraints_target_weight
        self.show(prefix = "mc %s:"%str(cycle))
        r_work = fmodel.r_work()
        geom = model.geometry_statistics(ignore_hd = True) # XXX
        if(r_work < r_work_best and ((geom.b_mean <= params.rmsd_max_bonds and
           geom.a_mean <= params.rmsd_max_angles) or
           ((geom.b_mean > params.rmsd_max_bonds or geom.a_mean > params.rmsd_max_angles)
           and (geom.b_mean <= b_start or geom.a_mean <= a_start)) )):
          r_work_best = r_work
          xrs_best = minimized.xray_structure.deep_copy_scatterers()
    fmodel.update_xray_structure(
      xray_structure = xrs_best,
      update_f_calc  = True,
      update_f_mask  = True)
    model.xray_structure = fmodel.xray_structure
    self.show(prefix = "Final")

  def compute_map(self, map_type=None, use_all_data=False):
    e_map_manager = self.fmodel.electron_density_map() # XXX pass in map filling options
    fft_map = e_map_manager.fft_map(
      resolution_factor = self.params.grid_resolution_factor,
      map_type          = map_type,
      use_all_data      = use_all_data)
    fft_map.apply_sigma_scaling()
    return fft_map

  def show(self, prefix):
    geom = self.model.geometry_statistics(ignore_hd = True) # XXX
    if(self.params.verbose):
      print >> self.log, \
        "%s r_work=%6.4f r_free=%6.4f rmsd bonds=%6.3f angles=%6.2f"%(
        prefix, self.fmodel.r_work(), self.fmodel.r_free(), geom.b_mean,
        geom.a_mean)
    return geom
